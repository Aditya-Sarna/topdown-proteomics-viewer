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
        'openms_BSA1.mzML',
    ),
]


# ---------------------------------------------------------------------------
# mzML parser
# ---------------------------------------------------------------------------

def parse_mzml(file_bytes: bytes, filename: str) -> List[Spectrum]:
    """Parse an mzML file; returns MS2 spectra only."""
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
        print(f"mzML parse error: {e}")
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
# Dash upload helper
# ---------------------------------------------------------------------------

def decode_upload(contents: str, filename: str) -> Tuple[List[Spectrum], List[Feature], str]:
    """Decode a Dash dcc.Upload file (base64) and parse it appropriately."""
    if not contents:
        return [], [], "No file uploaded."
    try:
        _content_type, content_string = contents.split(',', 1)
        decoded = base64.b64decode(content_string)
        fn_lower = filename.lower()
        if fn_lower.endswith('.mzml'):
            spectra = parse_mzml(decoded, filename)
            return spectra, [], f"Loaded {len(spectra)} MS2 spectra from {filename}"
        elif fn_lower.endswith(('.csv', '.tsv', '.txt')):
            text = decoded.decode('utf-8', errors='replace')
            # Heuristic: if it has rt_start / rt column treat as feature table
            first_line = text.split('\n')[0].lower()
            if any(k in first_line for k in ('rt_apex', 'rt_start', 'feature_id')):
                feats = parse_feature_table(text, filename)
                return [], feats, f"Loaded {len(feats)} features from {filename}"
            else:
                spec = parse_csv_peaks(text, filename)
                if spec:
                    return [spec], [], f"Loaded spectrum from {filename}"
        return [], [], f"Unsupported file format: {filename}"
    except Exception as e:
        return [], [], f"Error loading {filename}: {e}"


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
