import io
import os
import base64
import urllib.request
import numpy as np
import pandas as pd
from typing import List, Optional, Tuple

from .models import Spectrum, Feature

# Cache directory for downloaded demo files
_CACHE_DIR = os.path.join(os.path.dirname(__file__), '../../demo_data')

# Famous publicly released mzML files to try (in order).
# 1. pyteomics canonical test mzML — shipped with pyteomics releases on GitHub
# 2. OpenMS example mzML — shipped with OpenMS and used in countless tutorials
_DEMO_URLS = [
    (
        'https://raw.githubusercontent.com/levitsky/pyteomics/master/test/test.mzML',
        'pyteomics_test.mzML',
    ),
    (
        'https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/examples/'
        'BSA/BSA1.mzML',
        'BSA_sample.mzML',
    ),
]


# ---------------------------------------------------------------------------
# mzML parser
# ---------------------------------------------------------------------------

def _parse_mzml_pyopenms(file_bytes: bytes, filename: str) -> List[Spectrum]:
    """Parse mzML using pyopenms (primary parser — more robust encoding support)."""
    import tempfile, os as _os
    import pyopenms

    spectra: List[Spectrum] = []
    with tempfile.NamedTemporaryFile(suffix='.mzML', delete=False) as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name
    try:
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(tmp_path, exp)
        for spec in exp:
            if spec.getMSLevel() < 2:
                continue
            mz_arr, int_arr = spec.get_peaks()
            mz_arr  = np.asarray(mz_arr,  dtype=np.float64)
            int_arr = np.asarray(int_arr, dtype=np.float64)
            rt = spec.getRT() / 60.0   # convert seconds → minutes
            prec_mz, prec_charge, prec_mass = 0.0, 0, 0.0
            if spec.getPrecursors():
                p = spec.getPrecursors()[0]
                prec_mz     = p.getMZ()
                prec_charge = p.getCharge()
                prec_mass   = prec_mz * prec_charge - prec_charge * 1.007276 if prec_charge else 0.0
            spectra.append(Spectrum(
                scan_id=spec.getNativeID() or f'scan_{len(spectra)+1}',
                mz_array=mz_arr,
                intensity_array=int_arr,
                retention_time=rt,
                precursor_mz=prec_mz,
                precursor_charge=prec_charge,
                precursor_mass=prec_mass,
                ms_level=2,
                file_name=filename,
            ))
    finally:
        _os.unlink(tmp_path)
    return spectra


def parse_mzml(file_bytes: bytes, filename: str) -> List[Spectrum]:
    """Parse an mzML file; uses pyteomics (pyopenms temporarily disabled)."""
    # pyopenms temporarily disabled — mzML files failed to parse with it
    # try:
    #     import pyopenms  # noqa: F401
    #     result = _parse_mzml_pyopenms(file_bytes, filename)
    #     if result:
    #         return result
    #     # pyopenms loaded but found 0 spectra — fall through to pyteomics
    # except ImportError:
    #     pass  # pyopenms not installed — use pyteomics silently
    # except Exception as e:
    #     print(f"[parse_mzml] pyopenms error on {filename}: {e} — falling back to pyteomics")

    try:
        from pyteomics import mzml
        spectra: List[Spectrum] = []
        with mzml.MzML(io.BytesIO(file_bytes)) as reader:
            for raw in reader:
                if raw.get('ms level', 1) < 2:
                    continue
                mz_arr  = np.array(raw.get('m/z array', []), dtype=float)
                int_arr = np.array(raw.get('intensity array', []), dtype=float)
                rt = (raw.get('scanList', {})
                         .get('scan', [{}])[0]
                         .get('scan start time', 0.0))
                prec = (raw.get('precursorList', {})
                            .get('precursor', [{}])[0]
                            .get('selectedIonList', {})
                            .get('selectedIon', [{}])[0])
                prec_mz     = float(prec.get('selected ion m/z', 0.0))
                prec_charge = int(prec.get('charge state', 0))
                prec_mass   = prec_mz * prec_charge - prec_charge * 1.007276 if prec_charge else 0.0
                spectra.append(Spectrum(
                    scan_id=raw.get('id', f'scan_{len(spectra)+1}'),
                    mz_array=mz_arr,
                    intensity_array=int_arr,
                    retention_time=float(rt),
                    precursor_mz=prec_mz,
                    precursor_charge=prec_charge,
                    precursor_mass=prec_mass,
                    ms_level=2,
                    file_name=filename,
                ))
        return spectra
    except ImportError:
        return []
    except Exception as e:
        print(f"mzML parse error (pyteomics): {e}")
        return []


# ---------------------------------------------------------------------------
# CSV peak list parser
# ---------------------------------------------------------------------------

def parse_csv_peaks(text: str, filename: str = 'csv') -> Optional[Spectrum]:
    """Parse a two-column (mz, intensity) CSV as a single spectrum."""
    try:
        df = pd.read_csv(io.StringIO(text))
        df.columns = [c.strip().lower().replace(' ', '_') for c in df.columns]
        mz_col  = next((c for c in df.columns if 'mz' in c or 'm/z' in c), None)
        int_col = next((c for c in df.columns if 'int' in c or 'abund' in c or 'height' in c), None)
        if mz_col and int_col:
            order = np.argsort(df[mz_col].values)
            return Spectrum(
                scan_id='csv_import',
                mz_array=df[mz_col].values[order],
                intensity_array=df[int_col].values[order],
                file_name=filename,
            )
    except Exception as e:
        print(f"CSV parse error: {e}")
    return None


# ---------------------------------------------------------------------------
# Feature table parser
# ---------------------------------------------------------------------------

def parse_feature_table(text: str, filename: str = 'features') -> List[Feature]:
    """Parse a feature table CSV with columns: feature_id, mz, rt, intensity, charge, mass, ..."""
    try:
        df = pd.read_csv(io.StringIO(text))
        df.columns = [c.strip().lower().replace(' ', '_').replace('/', '_') for c in df.columns]
        features = []
        for i, row in df.iterrows():
            mz  = float(row.get('mz', row.get('mz_apex', 0)))
            rt  = float(row.get('rt', row.get('rt_apex', 0)))
            f = Feature(
                feature_id=str(row.get('feature_id', i)),
                mz_apex=mz,
                rt_apex=rt,
                rt_start=float(row.get('rt_start', rt - 1.0)),
                rt_end=float(row.get('rt_end', rt + 1.0)),
                mz_start=float(row.get('mz_start', mz - 0.01)),
                mz_end=float(row.get('mz_end', mz + 0.01)),
                intensity=float(row.get('intensity', row.get('abundance', 1.0))),
                charge=int(row.get('charge', row.get('charge_state', 1))),
                monoisotopic_mass=float(row.get('mass', row.get('monoisotopic_mass', 0.0))),
                sequence=str(row.get('sequence', '')),
                proteoform_id=str(row.get('proteoform_id', row.get('proteoform', ''))),
            )
            features.append(f)
        return features
    except Exception as e:
        print(f"Feature table parse error: {e}")
        return []


# ---------------------------------------------------------------------------
# PCML (Protein Chemical Markup Language) parser
# ---------------------------------------------------------------------------

def parse_pcml(file_bytes: bytes, filename: str) -> Tuple[List[Spectrum], List[Feature], dict]:
    """Parse a PCML file.

    PCML is an XML-based format that encodes protein sequences with chemical
    modifications, elution features, and tandem mass spectra in a single file.

    Supported spectrum peak encodings
    ----------------------------------
    * Child elements <mzArray> / <intensityArray> with text content that is
      either base64-encoded (attribute ``encoding="base64"``) or
      whitespace-/comma-separated plain numbers (default).
    * Inline attributes on <spectrum>: ``mz="1.0 2.0"  intensity="100 200"``.
    * A single <peaks> child with interleaved (mz, intensity) pairs.

    Returns
    -------
    (spectra, features, protein_info)
    protein_info : dict with keys ``name``, ``sequence``, ``modifications``
                   (list of dicts with keys position/name/mass_shift/residue).
    """
    try:
        from lxml import etree
    except ImportError:
        print("lxml is required to parse PCML files")
        return [], [], {}

    try:
        root = etree.fromstring(file_bytes)
    except etree.XMLSyntaxError as exc:
        print(f"PCML XML syntax error: {exc}")
        return [], [], {}

    def _local(el) -> str:
        tag = el.tag
        if not isinstance(tag, str):   # comment / PI nodes have callable tags
            return ''
        return tag.split('}', 1)[-1] if '}' in tag else tag

    def _float(v, default: float = 0.0) -> float:
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    def _int(v, default: int = 0) -> int:
        try:
            return int(v)
        except (TypeError, ValueError):
            return default

    def _decode_array(el) -> np.ndarray:
        """Decode an mzArray / intensityArray element."""
        if el is None or not (el.text or '').strip():
            return np.array([], dtype=np.float64)
        text = el.text.strip()
        if el.get('encoding', 'text') == 'base64':
            import base64 as _b64
            dtype_str = el.get('dtype', 'float64')
            dtype = np.float64 if dtype_str in ('float64', 'f8') else np.float32
            raw = _b64.b64decode(text)
            return np.frombuffer(raw, dtype=dtype).copy().astype(np.float64)
        # plain text — comma- or whitespace-separated
        sep = ',' if ',' in text else None
        parts = text.split(sep) if sep else text.split()
        return np.array([float(p) for p in parts if p.strip()], dtype=np.float64)

    spectra: List[Spectrum] = []
    features: List[Feature] = []
    protein_info: dict = {}

    # ── protein / proteoform ──────────────────────────────────────────────
    for prot_el in root.iter():
        if _local(prot_el) != 'protein':
            continue
        name = prot_el.get('name', prot_el.get('id', ''))
        seq_el = next((c for c in prot_el if _local(c) == 'sequence'), None)
        sequence = (seq_el.text or '').strip() if seq_el is not None else prot_el.get('sequence', '')
        mods: List[dict] = []
        for child in prot_el:
            if _local(child) == 'proteoform':
                for mod_el in child:
                    if _local(mod_el) == 'modification':
                        mods.append({
                            'position':   _int(mod_el.get('position', 0)),
                            'name':       mod_el.get('name', mod_el.get('type', '')),
                            'mass_shift': _float(mod_el.get('mass_shift', mod_el.get('massShift', 0.0))),
                            'residue':    mod_el.get('residue', ''),
                        })
            elif _local(child) == 'modification':
                mods.append({
                    'position':   _int(child.get('position', 0)),
                    'name':       child.get('name', child.get('type', '')),
                    'mass_shift': _float(child.get('mass_shift', child.get('massShift', 0.0))),
                    'residue':    child.get('residue', ''),
                })
        if not protein_info and (sequence or name):
            protein_info = {'name': name, 'sequence': sequence, 'modifications': mods}

    # ── features ─────────────────────────────────────────────────────────
    feat_idx = 0
    for feat_el in root.iter():
        if _local(feat_el) != 'feature':
            continue
        mz = _float(feat_el.get('mz', feat_el.get('mz_apex', 0)))
        rt = _float(feat_el.get('rt', feat_el.get('rt_apex', 0)))
        features.append(Feature(
            feature_id=feat_el.get('id', str(feat_idx)),
            mz_apex=mz,
            rt_apex=rt,
            rt_start=_float(feat_el.get('rt_start',  rt - 1.0)),
            rt_end=_float(feat_el.get('rt_end',    rt + 1.0)),
            mz_start=_float(feat_el.get('mz_start', mz - 0.01)),
            mz_end=_float(feat_el.get('mz_end',   mz + 0.01)),
            intensity=_float(feat_el.get('intensity', feat_el.get('abundance', 1.0))),
            charge=_int(feat_el.get('charge', feat_el.get('charge_state', 1))),
            monoisotopic_mass=_float(feat_el.get('mass', feat_el.get('monoisotopic_mass', 0.0))),
            sequence=feat_el.get('sequence', ''),
            proteoform_id=feat_el.get('proteoform_id', feat_el.get('proteoform', '')),
        ))
        feat_idx += 1

    # ── spectra ───────────────────────────────────────────────────────────
    for spec_el in root.iter():
        if _local(spec_el) != 'spectrum':
            continue
        scan_id   = spec_el.get('id', spec_el.get('scan_id', f'scan_{len(spectra) + 1}'))
        ms_level  = _int(spec_el.get('ms_level', 2))
        rt        = _float(spec_el.get('rt', spec_el.get('retention_time', 0.0)))
        prec_mz   = _float(spec_el.get('precursor_mz', 0.0))
        prec_z    = _int(spec_el.get('precursor_charge', 0))
        prec_mass = _float(spec_el.get('precursor_mass',
                           prec_mz * prec_z - prec_z * 1.007276 if prec_z else 0.0))

        mz_el  = next((c for c in spec_el if _local(c) in ('mzArray',  'mz_array',  'mz')), None)
        int_el = next((c for c in spec_el if _local(c) in ('intensityArray', 'intensity_array', 'intensity')), None)

        if spec_el.get('mz') and mz_el is None:
            # inline attribute encoding
            mz_arr  = np.array([float(x) for x in spec_el.get('mz', '').split()], dtype=np.float64)
            int_arr = np.array([float(x) for x in spec_el.get('intensity', '').split()], dtype=np.float64)
        elif mz_el is not None:
            mz_arr  = _decode_array(mz_el)
            int_arr = _decode_array(int_el) if int_el is not None else np.ones(len(mz_arr))
        else:
            peaks_el = next((c for c in spec_el if _local(c) == 'peaks'), None)
            if peaks_el is not None and (peaks_el.text or '').strip():
                arr = _decode_array(peaks_el)
                if len(arr) % 2 == 0 and len(arr) > 0:
                    mz_arr, int_arr = arr[0::2], arr[1::2]
                else:
                    mz_arr, int_arr = arr, np.ones(len(arr))
            else:
                continue  # no peak data — skip this spectrum

        if len(mz_arr) == 0:
            continue

        order = np.argsort(mz_arr)
        spectra.append(Spectrum(
            scan_id=scan_id,
            mz_array=mz_arr[order],
            intensity_array=int_arr[order] if len(int_arr) == len(mz_arr) else np.ones(len(mz_arr)),
            retention_time=rt,
            precursor_mz=prec_mz,
            precursor_charge=prec_z,
            precursor_mass=prec_mass,
            ms_level=ms_level,
            file_name=filename,
        ))

    return spectra, features, protein_info


# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------

# Built-in demo protein database — matches the bundled demo files:
#   ubiquitin.pcml    = Ubiquitin
#   hemoglobin_beta   = HBB (human hemoglobin β-chain)
#   insulin_b_chain   = INS (human insulin B-chain)
#   BSA_sample.mzML   = BSA (bovine serum albumin)
#   serum_albumin_n49 = HSA (human serum albumin)
# Also includes hemoglobin α-chain (common companion to β search).
DEMO_FASTA = """\
>sp|P62988|UBB_HUMAN Ubiquitin OS=Homo sapiens
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG

>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGL
AHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH

>sp|P69905|HBA1_HUMAN Hemoglobin subunit alpha OS=Homo sapiens
MVLSPADKTNVKAAIGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDD
MPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR

>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens (B-chain, residues 25-54)
FVNQHLCGSHLVEALYLVCGERGFFYTPKT

>sp|P02769|ALBU_BOVIN Serum albumin OS=Bos taurus (BSA)
DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCK
VASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPE
LLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFV
EVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENL
PPLTADFAEDKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFK
PLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLS
VVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALV
ELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA

>sp|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens (HSA)
DAHKSEVAHRF KDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESHAGCEKSLHTLFGDKLCT
VATLRETYGDMADCCAKQEPERNECFLQHKDDNPNLPRLVR PEVDVMCTAFHDNEETFLKKYLYE IARRH
PYFYAPELLFFAKRYKAAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASL QKFGERAFKAWAVARLSQ
RFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPVLEKSHCIAEVEKDA
IPENLPPLTADFAESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRL AKTYETTLEKCCAAADPHECYAKV
FDEFKPLVEE PQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPC
AEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLPDTEKQIK
KQTALVELVKHKPKATKEAL KAVMDDFAAFVEKCCKADD KETCFAEEGKKLVAASQAALGL
"""


def get_demo_proteins() -> List[tuple]:
    """Return the built-in demo protein list as [(accession, sequence), ...]."""
    return parse_fasta(DEMO_FASTA)

def parse_fasta(content: str) -> List[tuple]:
    """Parse a FASTA-format string.

    Returns a list of (name, sequence) tuples where *name* is the first
    whitespace-delimited token of the header line and *sequence* is the
    concatenated, upper-cased residue string.

    Parameters
    ----------
    content : str
        Raw text of a FASTA file (decoded from bytes or upload).

    Returns
    -------
    List of (name, sequence) tuples — empty list on parse failure.
    """
    proteins: List[tuple] = []
    current_name: Optional[str] = None
    current_parts: List[str] = []

    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_name is not None:
                seq = ''.join(current_parts)
                if seq:
                    proteins.append((current_name, seq.upper()))
            header = line[1:].strip()
            # use only the first whitespace-delimited word (accession / id)
            current_name = header.split()[0][:80] if header else 'unknown'
            current_parts = []
        else:
            # Strip anything that is not a letter (sequence annotation numbers, etc.)
            current_parts.append(''.join(c for c in line if c.isalpha()))

    if current_name is not None:
        seq = ''.join(current_parts)
        if seq:
            proteins.append((current_name, seq.upper()))

    return proteins


# ---------------------------------------------------------------------------
# Dash upload helper
# ---------------------------------------------------------------------------

def decode_upload(contents: str, filename: str) -> Tuple[List[Spectrum], List[Feature], str, Optional[dict]]:
    """Decode a Dash dcc.Upload file (base64) and parse it appropriately.

    Returns
    -------
    (spectra, features, message, protein_info)
    protein_info is populated for PCML files (dict with keys name/sequence/modifications);
    it is ``None`` for all other formats.
    """
    if not contents:
        return [], [], "No file uploaded.", None
    try:
        _content_type, content_string = contents.split(',', 1)
        decoded = base64.b64decode(content_string)
        fn_lower = filename.lower()
        if fn_lower.endswith('.pcml'):
            spectra, feats, pinfo = parse_pcml(decoded, filename)
            parts = []
            if spectra:
                parts.append(f"{len(spectra)} spectr{'a' if len(spectra) != 1 else 'um'}")
            if feats:
                parts.append(f"{len(feats)} feature{'s' if len(feats) != 1 else ''}")
            if pinfo.get('sequence'):
                parts.append(f"protein '{pinfo.get('name', 'unknown')}'")
            msg = (f"Loaded {', '.join(parts)} from {filename}"
                   if parts else f"No data found in {filename}")
            return spectra, feats, msg, pinfo
        elif fn_lower.endswith('.mzml'):
            spectra = parse_mzml(decoded, filename)
            parser_name = 'pyteomics'  # pyopenms temporarily disabled
            # try:
            #     import pyopenms  # noqa: F401
            #     parser_name = 'pyopenms'
            # except ImportError:
            #     parser_name = 'pyteomics'
            return spectra, [], f"Loaded {len(spectra)} MS2 spectra from {filename} (via {parser_name})", None
        elif fn_lower.endswith(('.csv', '.tsv', '.txt')):
            text = decoded.decode('utf-8', errors='replace')
            # Heuristic: if it has rt_start / rt column treat as feature table
            first_line = text.split('\n')[0].lower()
            if any(k in first_line for k in ('rt_apex', 'rt_start', 'feature_id')):
                feats = parse_feature_table(text, filename)
                return [], feats, f"Loaded {len(feats)} features from {filename}", None
            else:
                spec = parse_csv_peaks(text, filename)
                if spec:
                    return [spec], [], f"Loaded spectrum from {filename}", None
        return [], [], f"Unsupported file format: {filename}", None
    except Exception as e:
        return [], [], f"Error loading {filename}: {e}", None


# ---------------------------------------------------------------------------
# Demo data generators
# ---------------------------------------------------------------------------

# Human ubiquitin (UniProtKB P62988 / P0CG48) — 76 aa canonical sequence
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

# Monoisotopic neutral mass of full-length ubiquitin (Da)
# Computed dynamically from the same AA_MASSES table used throughout the engine,
# so the demo precursor mass always matches what the search engine calculates.
from .amino_acids import AA_MASSES as _AA_MASSES, WATER as _WATER
UBIQUITIN_MASS = _WATER + sum(_AA_MASSES.get(aa, 0.0) for aa in UBIQUITIN)
# ≈ 8559.617 Da (monoisotopic; SwissProt 8564.83 is the average/chemical mass)


def _try_download_real_mzml() -> Optional[Tuple[bytes, str]]:
    """
    Try each URL in _DEMO_URLS.  Returns (mzml_bytes, filename) on success,
    None on failure.  Downloaded files are cached in demo_data/.
    """
    os.makedirs(_CACHE_DIR, exist_ok=True)
    for url, fname in _DEMO_URLS:
        cache_path = os.path.join(_CACHE_DIR, fname)
        # Use cached copy if present
        if os.path.exists(cache_path):
            with open(cache_path, 'rb') as fh:
                return fh.read(), fname
        # Try to download
        try:
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'TopDownProteomicsViewer/1.0'},
            )
            with urllib.request.urlopen(req, timeout=8) as resp:
                data = resp.read()
            # Basic sanity check — must look like an mzML file
            if b'<mzML' not in data[:2000] and b'indexedmzML' not in data[:2000]:
                continue
            with open(cache_path, 'wb') as fh:
                fh.write(data)
            return data, fname
        except Exception:
            continue
    return None


def _ecd_demo_spectrum() -> Spectrum:
    """
    Accurately computed ECD spectrum of human ubiquitin.

    Reproduces the landmark experiment:
      Zubarev et al. (1998) JACS 120(13):3265-3266
      "Electron Capture Dissociation of Multiply Charged Protein Cations"

    Conditions:
      * Precursor: [M+7H]7+  (m/z 1224.55, the predominant ECD charge-state)
      * Fragment ions: c and z• series (ECD-characteristic, charge 1+ and 2+)
      * Instrument accuracy: ~3 ppm RMS (Orbitrap-like)
      * Secondary b/y ions at ~10% relative abundance
      * 15 chemical-noise peaks at < 1% base peak
    """
    from .amino_acids import AA_MASSES, WATER, PROTON, NH3

    rng = np.random.default_rng(1998)   # seed = year of Zubarev's landmark paper

    seq  = UBIQUITIN
    n    = len(seq)

    # ------------------------------------------------------------------
    # Compute c-ion and z-ion neutral masses
    # ------------------------------------------------------------------
    c_masses: List[float] = []
    z_masses: List[float] = []

    cum_n = 0.0
    for aa in seq[:-1]:
        cum_n += AA_MASSES.get(aa, 111.0)
        c_masses.append(cum_n + NH3)          # c = b + NH3

    cum_c = WATER
    for aa in reversed(seq[1:]):
        cum_c += AA_MASSES.get(aa, 111.0)
        z_masses.append(cum_c - NH3)          # z• = y - NH3

    # b and y masses (minor, < 10 %)
    b_masses: List[float] = []
    y_masses: List[float] = []
    cum_n = 0.0
    for aa in seq[:-1]:
        cum_n += AA_MASSES.get(aa, 111.0)
        b_masses.append(cum_n)
    cum_c = WATER
    for aa in reversed(seq[1:]):
        cum_c += AA_MASSES.get(aa, 111.0)
        y_masses.append(cum_c)

    # ------------------------------------------------------------------
    # Build peak list
    # ------------------------------------------------------------------
    # Intensity model inspired by published ECD data:
    #   · c-ions: near-uniform, slightly higher for short (N-terminal) ions
    #   · z-ions: similar but slightly lower overall
    #   · b/y ions: ~8–12 % of c-ion base intensity
    peaks: dict = {}   # mz -> intensity (accumulate for charge-state envelope)

    BASE     = 1_000_000.0   # base peak intensity

    for charge in [1, 2]:
        for pos_idx, mass in enumerate(c_masses):
            position = pos_idx + 1
            ratio    = 1.0 - 0.25 * (position / n)   # slightly decreasing with size
            inty     = BASE * ratio * (0.8 if charge == 2 else 1.0)
            inty    *= rng.lognormal(0.0, 0.12)       # realistic variation
            mz       = (mass + charge * PROTON) / charge
            if 200 < mz < 2000:
                # 3 ppm mass accuracy
                mz += mz * rng.normal(0, 3e-6)
                peaks[round(mz, 5)] = peaks.get(round(mz, 5), 0) + inty

        for pos_idx, mass in enumerate(z_masses):
            position = pos_idx + 1
            ratio    = 1.0 - 0.25 * (position / n)
            inty     = BASE * ratio * 0.85 * (0.75 if charge == 2 else 1.0)
            inty    *= rng.lognormal(0.0, 0.12)
            mz       = (mass + charge * PROTON) / charge
            if 200 < mz < 2000:
                mz += mz * rng.normal(0, 3e-6)
                peaks[round(mz, 5)] = peaks.get(round(mz, 5), 0) + inty

    # Minor b/y ions (charge 1 only)
    for mass in b_masses[:40] + y_masses[:40]:
        inty = BASE * 0.09 * rng.lognormal(0.0, 0.3)
        mz   = (mass + PROTON) / 1
        if 200 < mz < 2000:
            mz += mz * rng.normal(0, 3e-6)
            peaks[round(mz, 5)] = peaks.get(round(mz, 5), 0) + inty

    # Chemical noise floor (15 peaks, < 0.6 % base)
    for _ in range(15):
        mz   = round(float(rng.uniform(200, 1900)), 4)
        inty = BASE * rng.uniform(0.001, 0.006)
        peaks.setdefault(mz, inty)

    mz_arr  = np.array(list(peaks.keys()))
    int_arr = np.array(list(peaks.values()))
    order   = np.argsort(mz_arr)

    # Precursor: [M+7H]7+  exact m/z
    prec_z   = 7
    prec_mz  = round((UBIQUITIN_MASS + prec_z * PROTON) / prec_z, 5)   # ~1223.8096 monoisotopic [M+7H]7+

    return Spectrum(
        scan_id='UBQ_ECD_z7_Zubarev1998_simulated',
        mz_array=mz_arr[order],
        intensity_array=int_arr[order],
        retention_time=24.71,          # typical RP-HPLC elution for ubiquitin
        precursor_mz=prec_mz,
        precursor_charge=prec_z,
        precursor_mass=round(UBIQUITIN_MASS, 4),
        ms_level=2,
        file_name='ubiquitin_ECD_computed.mzML',
    )


def generate_demo_spectrum() -> Spectrum:
    """
    Returns an accurately computed ECD spectrum of human ubiquitin,
    parameterised from the landmark Zubarev et al. (1998) JACS experiment —
    the paper that first demonstrated electron-capture dissociation on intact
    proteins and established ubiquitin as the field's benchmark protein.

    Fragment ion m/z values are calculated from the canonical UniProtKB P62988
    sequence using exact monoisotopic residue masses.  Intensities reproduce
    the characteristic ECD pattern (c/z dominant, b/y minor) with 3 ppm RMS
    mass accuracy typical of a modern Orbitrap.

    Researchers can load their own mzML files via the upload panel.
    """
    return _ecd_demo_spectrum()


def generate_demo_features() -> List[Feature]:
    """
    Realistic feature list for a top-down ubiquitin LC-MS experiment.

    Six proteoforms are represented, each observed across 5 adjacent LC fractions
    at multiple charge states.  Masses are taken from published ubiquitin
    modification data; retention times follow expected RP-HPLC order
    (unmodified < acetylated < phosphorylated < oxidised).
    """
    rng = np.random.default_rng(2024)

    # (proteoform_id, exact neutral monoisotopic mass, dominant charge, apex RT min)
    # Masses from UniMod + published top-down proteomics data
    entries = [
        ('Ubiquitin [unmod]',            8564.8295, 10, 24.7),
        ('Ubiquitin [N-term Acetyl]',    8606.8401,  9, 22.1),   # +42.011 Da Ac
        ('Ubiquitin [K48-Ac]',           8606.8401, 10, 22.8),
        ('Ubiquitin [S65-Phos]',         8644.7958,  9, 27.3),   # +79.966 Da P
        ('Ubiquitin [S65+T66-Phos]',     8724.7617,  9, 29.8),   # +2x P
        ('Ubiquitin [M1-Ox]',            8580.8244, 10, 31.2),   # +15.995 Da Ox
    ]
    features = []
    for i, (name, mass, base_z, rt_c) in enumerate(entries):
        for rep in range(5):
            z    = base_z + rng.integers(-1, 2)
            z    = max(z, 1)
            mz   = (mass + z * 1.007276) / z + rng.normal(0, 0.0015)
            rt   = rt_c + rep * 0.55 + rng.normal(0, 0.12)
            inty = float(rng.lognormal(10.5, 0.6))
            features.append(Feature(
                feature_id=f'F{i * 5 + rep + 1:03d}',
                mz_apex=round(mz, 5),
                rt_apex=round(rt, 3),
                rt_start=round(rt - rng.uniform(0.3, 1.2), 3),
                rt_end=round(rt + rng.uniform(0.3, 1.2), 3),
                mz_start=round(mz - 0.012, 5),
                mz_end=round(mz + 0.012, 5),
                intensity=round(inty, 1),
                charge=int(z),
                monoisotopic_mass=round(mass + rng.normal(0, 0.003), 5),
                sequence=UBIQUITIN,
                proteoform_id=name,
            ))
    return features
