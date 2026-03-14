"""Quick smoke test — run with: python test_pipeline.py"""
import pathlib
from src.data.parsers import (
    generate_demo_spectrum, generate_demo_features, UBIQUITIN,
    parse_pcml, parse_mzml,
)
from src.analysis.fragment_ions import calc_ions
from src.analysis.peak_matching import match_peaks, sequence_coverage_pct
from src.analysis.proteoform_search import run_targeted_search

DEMO_DATA = pathlib.Path(__file__).parent / "demo_data"

# ── 1. Basic demo spectrum / features ─────────────────────────────────────
spec  = generate_demo_spectrum()
feats = generate_demo_features()
print(f"Demo spectrum : {len(spec.mz_array)} peaks | RT={spec.retention_time} min")
print(f"Demo features : {len(feats)} features")

# ── 2. Fragment ions & peak matching ──────────────────────────────────────
ions    = calc_ions(UBIQUITIN[:30], ['b','y','c','z'], {}, max_charge=3)
matched = match_peaks(ions, spec, tolerance_ppm=10)
n_match = sum(1 for i in matched if i.matched)
cov     = sequence_coverage_pct(UBIQUITIN[:30], matched)
print(f"Fragment ions : {len(ions)} theoretical, {n_match} matched, coverage {cov:.1f}%")

# ── 3. Proteoform search ──────────────────────────────────────────────────
results = run_targeted_search(
    spec, UBIQUITIN, 'Ubiquitin',
    ion_types=['b','y','c','z'], max_charge=3, tolerance_ppm=10,
    search_truncations=True, search_modifications=True,
    variable_mods=['Phosphorylation','Acetylation'],
)
print(f"Search results: {len(results)} hits")
if results:
    pf = results[0].proteoform
    print(f"  Top hit      : score={pf.score}, "
          f"coverage={results[0].sequence_coverage}%, "
          f"Δppm={pf.mass_error_ppm}, PTMs={[m.name for m in pf.modifications]}")


# ── 5. PCML parser — hemoglobin_beta.pcml ─────────────────────────────────
hbb_bytes = (DEMO_DATA / "hemoglobin_beta.pcml").read_bytes()
hbb_spectra, hbb_feats, hbb_info = parse_pcml(hbb_bytes, "hemoglobin_beta.pcml")
assert hbb_info.get("sequence"), "hemoglobin_beta.pcml: protein sequence missing"
print(f"PCML HBB      : protein='{hbb_info['name']}' seq_len={len(hbb_info['sequence'])} "
      f"| {len(hbb_spectra)} spectrum/a | {len(hbb_feats)} features "
      f"| {len(hbb_info['modifications'])} mods")

# ── 6. PCML parser — insulin_b_chain.pcml ────────────────────────────────
ins_bytes = (DEMO_DATA / "insulin_b_chain.pcml").read_bytes()
ins_spectra, ins_feats, ins_info = parse_pcml(ins_bytes, "insulin_b_chain.pcml")
assert ins_info.get("sequence"), "insulin_b_chain.pcml: protein sequence missing"
print(f"PCML Insulin  : protein='{ins_info['name']}' seq_len={len(ins_info['sequence'])} "
      f"| {len(ins_spectra)} spectrum/a | {len(ins_feats)} features "
      f"| {len(ins_info['modifications'])} mods")

# ── 7. PCML parser — serum_albumin_n49.pcml ──────────────────────────────
hsa_bytes = (DEMO_DATA / "serum_albumin_n49.pcml").read_bytes()
hsa_spectra, hsa_feats, hsa_info = parse_pcml(hsa_bytes, "serum_albumin_n49.pcml")
assert hsa_info.get("sequence"), "serum_albumin_n49.pcml: protein sequence missing"
print(f"PCML HSA-N49  : protein='{hsa_info['name']}' seq_len={len(hsa_info['sequence'])} "
      f"| {len(hsa_spectra)} spectrum/a | {len(hsa_feats)} features "
      f"| {len(hsa_info['modifications'])} mods")

# ── 8. mzML parser — BSA_sample.mzML ────────────────────────────────────
mzml_bytes = (DEMO_DATA / "BSA_sample.mzML").read_bytes()
mzml_spectra = parse_mzml(mzml_bytes, "BSA_sample.mzML")
assert len(mzml_spectra) >= 1, "BSA_sample.mzML: no spectra parsed"
print(f"mzML BSA1     : {len(mzml_spectra)} spectra")

# ── 9. PCML robustness — XML comment nodes must not crash parser ──────────
_xml_with_comments = b"""<?xml version="1.0"?>
<!-- top-level comment -->
<pcml version="1.0">
  <!-- section comment -->
  <protein id="P1" name="TestProt">
    <sequence>ACDEFGHIKLM</sequence>
    <!-- mod comment -->
  </protein>
  <spectrumList>
    <spectrum id="s1" ms_level="2" rt="1.0" precursor_mz="500.0" precursor_charge="2">
      <mzArray>100.0 200.0 300.0</mzArray>
      <intensityArray>1000 2000 3000</intensityArray>
    </spectrum>
  </spectrumList>
</pcml>"""
_c_spectra, _c_feats, _c_info = parse_pcml(_xml_with_comments, "comment_test.pcml")
assert _c_info.get("sequence") == "ACDEFGHIKLM", "comment-node robustness: wrong sequence"
assert len(_c_spectra) == 1,                      "comment-node robustness: wrong spectrum count"
print("PCML comments : XML comment-node robustness OK")

print("\nALL TESTS PASS")
