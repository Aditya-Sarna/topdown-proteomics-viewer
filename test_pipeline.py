"""Quick smoke test — run with: python test_pipeline.py"""
from src.data.parsers import generate_demo_spectrum, generate_demo_features, UBIQUITIN
from src.analysis.fragment_ions import calc_ions
from src.analysis.peak_matching import match_peaks, sequence_coverage_pct
from src.analysis.proteoform_search import run_targeted_search

spec  = generate_demo_spectrum()
feats = generate_demo_features()
print(f"Demo spectrum : {len(spec.mz_array)} peaks | RT={spec.retention_time} min")
print(f"Demo features : {len(feats)} features")

ions    = calc_ions(UBIQUITIN[:30], ['b','y','c','z'], {}, max_charge=3)
matched = match_peaks(ions, spec, tolerance_ppm=10)
n_match = sum(1 for i in matched if i.matched)
cov     = sequence_coverage_pct(UBIQUITIN[:30], matched)
print(f"Fragment ions : {len(ions)} theoretical, {n_match} matched, coverage {cov:.1f}%")

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
print("PASS")
