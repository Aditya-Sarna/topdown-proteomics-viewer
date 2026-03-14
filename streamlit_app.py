
import io
import os
import csv
import json
import base64

import numpy as np
import pandas as pd
import streamlit as st

from src.data.parsers import (
    parse_mzml, parse_csv_peaks, parse_feature_table, parse_pcml,
)
from src.data.models import Spectrum, Feature, Proteoform, FragmentIon, SearchResult
from src.analysis.mass_utils import calc_sequence_mass, suggest_modifications, ppm_error
from src.analysis.proteoform_search import run_targeted_search
from src.analysis.peak_matching import coverage_map
from src.viz.spectrum_plots import create_spectrum_plot
from src.viz.sequence_plots import create_sequence_plot
from src.viz.mirror_plots import create_mirror_plot
from src.viz.feature_plots import create_feature_map, create_intensity_trace
from src.data.amino_acids import AA_MASSES, PTM_DATABASE

st.set_page_config(
    page_title="Top-Down Proteomics Viewer",
    layout="wide",
    initial_sidebar_state="expanded",
)

_DEFAULTS = {
    "spectra":         [],
    "features":        [],
    "scan_idx":        0,
    "protein":         {},
    "search_results":  [],
    "matched_ions":    [],
    "selected_result": None,
    "_prot_name":      "",
    "_prot_seq":       "",
}
for _k, _v in _DEFAULTS.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v


def _clean_seq(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in AA_MASSES)


_DEMO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_data")


def _load_demo_file(filename: str) -> tuple[list, list, dict, str]:
    """Read a file from demo_data/ by name and return (spectra, feats, pinfo, msg)."""
    path = os.path.join(_DEMO_DIR, filename)
    with open(path, "rb") as _fh:
        raw = _fh.read()
    fname = filename.lower()
    if fname.endswith(".pcml"):
        spectra, feats, pinfo = parse_pcml(raw, filename)
        parts = []
        if spectra:
            parts.append(f"{len(spectra)} spectr{'a' if len(spectra) != 1 else 'um'}")
        if feats:
            parts.append(f"{len(feats)} feature{'s' if len(feats) != 1 else ''}")
        if pinfo.get('sequence'):
            parts.append(f"protein '{pinfo.get('name', 'unknown')}'")
        msg = (f"Loaded {', '.join(parts)} from {filename}"
               if parts else f"No data found in {filename}")
        return spectra, feats, pinfo, msg
    if fname.endswith(".mzml"):
        spectra = parse_mzml(raw, filename)
        # Auto-load paired feature CSV if present (e.g. hemoglobin_beta_features.csv)
        feats: list = []
        base = filename.rsplit(".", 1)[0]
        csv_path = os.path.join(_DEMO_DIR, base + "_features.csv")
        if os.path.exists(csv_path):
            with open(csv_path, "r", encoding="utf-8") as _cf:
                feats = parse_feature_table(_cf.read(), base + "_features.csv")
        n_feats = f" + {len(feats)} features" if feats else ""
        return spectra, feats, {}, f"Loaded {len(spectra)} MS2 scan(s){n_feats} from {filename}"
    if fname.endswith((".csv", ".tsv", ".txt")):
        text = raw.decode("utf-8", errors="replace")
        spec = parse_csv_peaks(text, filename)
        if spec:
            return [spec], [], {}, f"Loaded spectrum from {filename}"
        return [], [], {}, f"Could not parse peak list from {filename}"
    return [], [], {}, f"Unsupported format: {filename}"


def _load_spectrum_file(uploaded) -> tuple[list, list, dict, str]:
    """Return (spectra, features, protein_info, message)."""
    raw = uploaded.read()
    fname = uploaded.name.lower()
    if fname.endswith(".pcml"):
        spectra, feats, pinfo = parse_pcml(raw, uploaded.name)
        parts = []
        if spectra:
            parts.append(f"{len(spectra)} spectr{'a' if len(spectra) != 1 else 'um'}")
        if feats:
            parts.append(f"{len(feats)} feature{'s' if len(feats) != 1 else ''}")
        if pinfo.get('sequence'):
            parts.append(f"protein '{pinfo.get('name', 'unknown')}'")
        msg = (f"Loaded {', '.join(parts)} from {uploaded.name}"
               if parts else f"No data found in {uploaded.name}")
        return spectra, feats, pinfo, msg
    if fname.endswith(".mzml"):
        spectra = parse_mzml(raw, uploaded.name)
        return spectra, [], {}, f"Loaded {len(spectra)} MS2 scan(s) from {uploaded.name}"
    elif fname.endswith((".csv", ".tsv", ".txt")):
        text = raw.decode("utf-8", errors="replace")
        spec = parse_csv_peaks(text, uploaded.name)
        if spec:
            return [spec], [], {}, f"Loaded spectrum from {uploaded.name}"
        return [], [], {}, f"Could not parse peak list from {uploaded.name}"
    return [], [], {}, f"Unsupported format: {uploaded.name}"


def _load_features_file(uploaded) -> tuple[list, str]:
    text = uploaded.read().decode("utf-8", errors="replace")
    feats = parse_feature_table(text, uploaded.name)
    return feats, f"Loaded {len(feats)} feature(s) from {uploaded.name}"


with st.sidebar:
    st.markdown("## Top-Down Proteomics Viewer")
    st.caption("v1.0 — Single-spectrum proteoform analysis")
    st.divider()

    st.markdown("##### TEST DATASETS / FILES")

    # Hardcoded demo entries — label → filename in demo_data/
    _DEMO_OPTIONS: dict[str, str | None] = {
        "— select a dataset —": None,
        # ── PCML temporarily removed ──────────────────────────────────────
        # "Hemoglobin Beta · PCML (spectra + features + sequence)": "hemoglobin_beta.pcml",
        # "Insulin B Chain · PCML (spectra + features + sequence)": "insulin_b_chain.pcml",
        # "Serum Albumin N49 · PCML (spectra + features + sequence)": "serum_albumin_n49.pcml",
        # ── mzML + auto-loaded features CSV ───────────────────────────────
        "Hemoglobin Beta · mzML (3 scans, 3 features)": "hemoglobin_beta.mzML",
        "Hemoglobin Alpha · mzML (3 scans, 3 features)": "hemoglobin_alpha.mzML",
        "Insulin B Chain · mzML (3 scans, 3 features)": "insulin_b_chain.mzML",
        "Serum Albumin N49 · mzML (3 scans, 3 features)": "serum_albumin_n49.mzML",
        "Ubiquitin · mzML (3 scans, 3 features)": "ubiquitin.mzML",
        "Cytochrome C · mzML (3 scans, 3 features)": "cytochrome_c.mzML",
        "Thioredoxin · mzML (3 scans, 3 features)": "thioredoxin.mzML",
    }
    _selected_demo = st.selectbox(
        "Select dataset",
        options=list(_DEMO_OPTIONS.keys()),
        key="demo_file_select",
        label_visibility="collapsed",
    )
    if st.button("Load", use_container_width=True, type="secondary"):
        _chosen = st.session_state.get("demo_file_select", "— select a dataset —")
        if _DEMO_OPTIONS.get(_chosen) is None:
            st.warning("Please select a dataset first.")
        else:
            _d_spectra, _d_feats, _d_pinfo, _d_msg = _load_demo_file(_DEMO_OPTIONS[_chosen])
            if _d_spectra or _d_feats or _d_pinfo.get('sequence'):
                if _d_spectra:
                    st.session_state.spectra = _d_spectra
                    st.session_state.scan_idx = 0
                if _d_feats:
                    st.session_state.features = _d_feats
                if _d_pinfo.get('sequence'):
                    _clean = _clean_seq(_d_pinfo['sequence'])
                    st.session_state["_prot_name"] = _d_pinfo.get('name', '')
                    st.session_state["_prot_seq"]  = _clean
                    st.session_state.protein = {
                        'name': _d_pinfo.get('name', ''),
                        'sequence': _clean,
                        'mass': calc_sequence_mass(_clean),
                    }
                else:
                    st.session_state["_prot_name"] = ""
                    st.session_state["_prot_seq"]  = ""
                    st.session_state.protein = {}
                st.session_state.search_results = []
                st.session_state.matched_ions = []
                st.session_state.selected_result = None
                st.session_state["_demo_msg"] = ("success", _d_msg)
            else:
                st.session_state["_demo_msg"] = ("error", _d_msg)
            st.rerun()
    _dm = st.session_state.pop("_demo_msg", None)
    if _dm:
        if _dm[0] == "success":
            st.success(_dm[1])
        else:
            st.error(_dm[1])


    # File upload disabled — use the demo datasets above
    # uploaded_spec = st.file_uploader(
    #     "Upload Spectrum (.mzML / .pcml / peak CSV)",
    #     type=["pcml", "mzml", "csv", "tsv", "txt"],
    #     key="upload_spectrum",
    # )
    # if uploaded_spec:
    #     spectra, feats, pinfo, msg = _load_spectrum_file(uploaded_spec)
    #     ...

    # uploaded_feats = st.file_uploader(
    #     "Upload Feature Table (.csv)",
    #     type=["csv", "tsv", "txt"],
    #     key="upload_features",
    # )
    # if uploaded_feats:
    #     feats, msg = _load_features_file(uploaded_feats)
    #     ...

    st.divider()

    st.markdown("##### SCAN")
    spectra: list[Spectrum] = st.session_state.spectra
    if spectra:
        scan_labels = [
            f"{s.scan_id}  RT={s.retention_time:.2f} min  "
            f"precursor={s.precursor_mz:.4f} m/z  z={s.precursor_charge}"
            for s in spectra
        ]
        chosen = st.selectbox(
            "Select Scan", options=range(len(spectra)),
            format_func=lambda i: scan_labels[i],
            index=st.session_state.scan_idx,
            label_visibility="collapsed",
        )
        st.session_state.scan_idx = chosen
        sel = spectra[chosen]
        st.caption(
            f"Scan {sel.scan_id} | {len(sel.mz_array)} peaks | "
            f"RT {sel.retention_time:.2f} min | "
            f"Prec {sel.precursor_mz:.3f} m/z z={sel.precursor_charge}"
        )
    else:
        st.caption("_Load a spectrum file first…_")

    st.divider()

    st.markdown("##### PROTEIN")
    prot_name = st.text_input(
        "Protein Name",
        key="_prot_name",
        placeholder="e.g. Ubiquitin",
    )
    prot_seq_raw = st.text_area(
        "Sequence (single-letter)",
        key="_prot_seq",
        placeholder="MQIFVKTLTGK…",
        height=100,
    )
    prot_seq = _clean_seq(prot_seq_raw)
    if prot_seq:
        prot_mass = calc_sequence_mass(prot_seq)
        st.session_state.protein = {"name": prot_name or "Protein", "sequence": prot_seq, "mass": prot_mass}
        st.caption(f"{len(prot_seq)} aa | Monoisotopic mass: {prot_mass:.4f} Da")
    else:
        st.session_state.protein = {}

    st.divider()

    st.markdown("##### SEARCH PARAMETERS")
    col1, col2 = st.columns(2)
    with col1:
        tol_ppm = st.number_input("Tolerance (ppm)", min_value=1, max_value=100, value=10, step=1)
    with col2:
        max_z = st.number_input("Max z", min_value=1, max_value=20, value=4, step=1)

    ion_types = st.multiselect(
        "Ion Types",
        options=["b", "y", "c", "z", "a"],
        default=["b", "y", "c", "z"],
    )

    col3, col4 = st.columns(2)
    with col3:
        search_trunc = st.checkbox("Truncations", value=True)
    with col4:
        search_mods = st.checkbox("Modifications", value=True)

    ptm_options = sorted(PTM_DATABASE.keys())
    variable_mods = st.multiselect(
        "Variable PTMs",
        options=ptm_options,
        default=["Phosphorylation", "Acetylation", "Oxidation"],
    )

    st.divider()

    st.markdown("##### MANUAL MASS INPUT")
    manual_mass = st.number_input(
        "Observed Neutral Mass (Da)",
        min_value=0.0, value=0.0, step=0.001,
        format="%.3f",
        label_visibility="collapsed",
        placeholder="e.g. 8564.83",
        key="manual_mass_input",
    )
    if st.button("Calculate Δ mass", use_container_width=True):
        prot_data = st.session_state.protein
        if manual_mass and prot_data:
            obs = float(manual_mass)
            th = prot_data.get("mass", 0)
            diff = obs - th
            ppm = ppm_error(obs, th) if th else 0.0
            st.info(f"Δmass = {diff:+.4f} Da  ({ppm:+.1f} ppm)")
            sugg = suggest_modifications(diff)
            if sugg:
                for s in sugg[:6]:
                    st.caption(f"• {s['name']}  ({s['mass_shift']:+.4f} Da,  residual {s['residual_da']:+.5f} Da)")
            else:
                st.caption("No matching PTMs found.")
        else:
            st.warning("Enter a protein sequence and observed mass first.")

    st.divider()

    run_clicked = st.button("Run Search", use_container_width=True, type="primary")
    if run_clicked:
        if not st.session_state.spectra:
            st.error("Load a spectrum first.")
        elif not st.session_state.protein.get("sequence"):
            st.error("Enter a protein sequence.")
        else:
            with st.spinner("Running search…"):
                spectrum = st.session_state.spectra[st.session_state.scan_idx]
                seq = st.session_state.protein["sequence"]
                name = st.session_state.protein.get("name", "Protein")
                results = run_targeted_search(
                    spectrum=spectrum,
                    protein_sequence=seq,
                    protein_name=name,
                    ion_types=ion_types or ["b", "y", "c", "z"],
                    max_charge=int(max_z),
                    tolerance_ppm=float(tol_ppm),
                    search_truncations=search_trunc,
                    search_modifications=search_mods,
                    variable_mods=variable_mods,
                )
            if results:
                st.session_state.search_results = results
                st.session_state.matched_ions = results[0].fragment_ions
                st.session_state.selected_result = results[0].proteoform
                st.success(f"{len(results)} hit(s) found")
            else:
                st.session_state.search_results = []
                st.session_state.matched_ions = []
                st.session_state.selected_result = None
                st.warning("No results found.")

    if st.session_state.search_results:
        st.divider()
        st.markdown("##### TOP HIT SUMMARY")
        top = st.session_state.search_results[0]
        pf0 = top.proteoform
        st.caption(
            f"Score: **{pf0.score:.2f}**  \n"
            f"Proteoform: {pf0.start_pos}–{pf0.end_pos}  \n"
            f"Ions matched: {pf0.matched_ions}/{pf0.total_ions}  \n"
            f"Coverage: {top.sequence_coverage:.1f}%  \n"
            f"Th. mass: {pf0.theoretical_mass:.4f} Da  \n"
            f"Delta mass: {pf0.mass_error_ppm:+.2f} ppm"
            + (f"  \nPTMs: {'; '.join(m.name for m in pf0.modifications)}"
               if pf0.modifications else "")
        )


st.title("Top-Down Proteomics Viewer")

tab_spec, tab_seq, tab_mirror, tab_results, tab_features = st.tabs(
    ["Spectrum", "Sequence", "Mirror Plot", "Search Results", "Feature Map"]
)

with tab_spec:
    spectra = st.session_state.spectra
    if spectra:
        spectrum = spectra[st.session_state.scan_idx]
        ions = st.session_state.matched_ions or None
        fig = create_spectrum_plot(spectrum, ions if ions else None)
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Upload a spectrum file or select a dataset from **Test Datasets / Files** to begin.")

with tab_seq:
    show_cleavage = st.checkbox("Show cleavage ticks", value=True, key="show_cleavage")

    prot_data = st.session_state.protein
    selected_result = st.session_state.selected_result
    ions = st.session_state.matched_ions

    if selected_result:
        pf = selected_result
    elif prot_data and prot_data.get("sequence"):
        pf = Proteoform(sequence=prot_data["sequence"], protein_name=prot_data.get("name", "Protein"))
    else:
        pf = Proteoform(sequence="", protein_name="")

    if pf.sequence:
        cov = coverage_map(pf.sequence, ions)
        n_cov = sum(1 for v in cov.values() if v)
        pct = n_cov / len(pf.sequence) * 100
        st.caption(f"Sequence coverage: {n_cov}/{len(pf.sequence)} aa ({pct:.1f}%)")

    fig_seq = create_sequence_plot(pf, ions, show_cleavage=show_cleavage)
    st.plotly_chart(fig_seq, use_container_width=True)

with tab_mirror:
    col_a, col_b = st.columns(2)
    with col_a:
        mz_min = st.number_input("m/z Min", value=200.0, step=10.0)
    with col_b:
        mz_max = st.number_input("m/z Max", value=2000.0, step=10.0)

    if st.session_state.spectra:
        spectrum = st.session_state.spectra[st.session_state.scan_idx]
        ions = st.session_state.matched_ions or None
        fig_mirror = create_mirror_plot(spectrum, ions, mz_range=(mz_min, mz_max))
        st.plotly_chart(fig_mirror, use_container_width=True)
    else:
        st.info("Load a spectrum first.")

with tab_results:
    results: list[SearchResult] = st.session_state.search_results
    if not results:
        st.info("Run a search to see results here.")
    else:
        spectrum = st.session_state.spectra[st.session_state.scan_idx]
        st.caption(
            f"{len(results)} candidate(s) | Scan {spectrum.scan_id} | "
            f"Tolerance {tol_ppm} ppm"
        )

        rows = []
        for rank, r in enumerate(results, 1):
            pf = r.proteoform
            mods_str = "; ".join(f"{m.name}@{m.position}" for m in pf.modifications) or "—"
            rows.append({
                "Rank": rank,
                "Sequence": pf.sequence[:30] + ("…" if len(pf.sequence) > 30 else ""),
                "Range": f"{pf.start_pos}–{pf.end_pos}",
                "Score": round(pf.score, 2),
                "Matched": f"{pf.matched_ions}/{pf.total_ions}",
                "Coverage %": round(r.sequence_coverage, 1),
                "Th. mass (Da)": round(pf.theoretical_mass, 4),
                "Obs. mass (Da)": round(pf.observed_mass, 4),
                "Δ ppm": round(pf.mass_error_ppm, 2),
                "PTMs": mods_str,
            })
        df_results = pd.DataFrame(rows)
        selected_rows = st.dataframe(
            df_results,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row",
            key="results_table",
        )

        sel_indices = selected_rows.selection.rows if selected_rows else []
        if sel_indices:
            sel_idx = sel_indices[0]
            chosen_result = results[sel_idx]
            st.session_state.matched_ions = chosen_result.fragment_ions
            st.session_state.selected_result = chosen_result.proteoform

        display_result = (
            results[sel_indices[0]] if sel_indices else results[0]
        )
        st.markdown("**Fragment Ion Matches**")
        ion_rows = []
        for ion in sorted(display_result.fragment_ions, key=lambda x: (x.ion_type, x.position)):
            ion_rows.append({
                "Ion": ion.label(),
                "Th. m/z": round(ion.mz, 4),
                "Obs. m/z": round(ion.observed_mz, 4) if ion.matched else "—",
                "Δ ppm": round(ion.mass_error_ppm, 2) if ion.matched else "—",
                "z": ion.charge,
                "Matched": "yes" if ion.matched else "no",
                "Sequence": ion.sequence[:15] if ion.sequence else "",
            })
        st.dataframe(pd.DataFrame(ion_rows), use_container_width=True, hide_index=True)

        st.markdown("**Export**")
        exp_col1, exp_col2, exp_col3, exp_col4 = st.columns(4)

        with exp_col1:
            buf = io.StringIO()
            fields = ["rank","protein_name","sequence","start_pos","end_pos",
                      "score","matched_ions","total_ions","sequence_coverage",
                      "theoretical_mass","observed_mass","mass_error_da","mass_error_ppm","modifications"]
            w = csv.DictWriter(buf, fieldnames=fields)
            w.writeheader()
            for rank, r in enumerate(results, 1):
                pf = r.proteoform
                w.writerow({
                    "rank": rank, "protein_name": pf.protein_name,
                    "sequence": pf.sequence, "start_pos": pf.start_pos, "end_pos": pf.end_pos,
                    "score": pf.score, "matched_ions": pf.matched_ions, "total_ions": pf.total_ions,
                    "sequence_coverage": r.sequence_coverage,
                    "theoretical_mass": pf.theoretical_mass, "observed_mass": pf.observed_mass,
                    "mass_error_da": pf.mass_error_da, "mass_error_ppm": pf.mass_error_ppm,
                    "modifications": "; ".join(f"{m.name}@{m.position}" for m in pf.modifications),
                })
            st.download_button("Results CSV", buf.getvalue(), "proform_results.csv", "text/csv")

        with exp_col2:
            buf2 = io.StringIO()
            ions_all = st.session_state.matched_ions
            fields2 = ["ion_type","position","charge","th_mz","obs_mz","mass_error_da","mass_error_ppm","matched","sequence"]
            w2 = csv.DictWriter(buf2, fieldnames=fields2)
            w2.writeheader()
            for ion in ions_all:
                w2.writerow({
                    "ion_type": ion.ion_type, "position": ion.position, "charge": ion.charge,
                    "th_mz": ion.mz, "obs_mz": ion.observed_mz if ion.matched else "",
                    "mass_error_da": ion.mass_error_da if ion.matched else "",
                    "mass_error_ppm": ion.mass_error_ppm if ion.matched else "",
                    "matched": ion.matched, "sequence": ion.sequence or "",
                })
            st.download_button("Ions CSV", buf2.getvalue(), "proform_ions.csv", "text/csv")

        with exp_col3:
            feats = st.session_state.features
            if feats:
                buf3 = io.StringIO()
                rows_f = [f.__dict__ for f in feats]
                w3 = csv.DictWriter(buf3, fieldnames=list(rows_f[0].keys()))
                w3.writeheader()
                w3.writerows(rows_f)
                st.download_button("Features CSV", buf3.getvalue(), "proform_features.csv", "text/csv")
            else:
                st.button("Features CSV", disabled=True)

        with exp_col4:
            session_json = json.dumps({
                "search_results": [r.to_dict() for r in results],
                "matched_ions":   [i.to_dict() for i in st.session_state.matched_ions],
                "features":       [f.__dict__ for f in st.session_state.features],
            }, indent=2, default=str)
            st.download_button("JSON Session", session_json, "proform_session.json", "application/json")

with tab_features:
    features: list[Feature] = st.session_state.features
    if features:
        col_f1, col_f2 = st.columns([2, 3])
        with col_f1:
            color_by = st.radio("Color by", ["charge", "proteoform"], horizontal=True)
            fid_filter = st.text_input("Feature ID filter", placeholder="e.g. F001")

        with col_f2:
            col_m1, col_m2 = st.columns(2)
            with col_m1:
                th_mass_input = st.number_input("Theoretical mass overlay (Da)", value=0.0, step=0.001, format="%.3f")
            with col_m2:
                th_charge_input = st.number_input("Charge for trace", value=10, min_value=1, step=1)

        proteoforms = []
        if th_mass_input:
            proteoforms = [Proteoform(sequence="", protein_name="Target", theoretical_mass=float(th_mass_input))]
        elif st.session_state.search_results:
            for r in st.session_state.search_results[:3]:
                proteoforms.append(r.proteoform)

        selected_id = fid_filter.strip() if fid_filter else None

        fig_feat = create_feature_map(features, proteoforms or None, selected_id, color_by)
        st.plotly_chart(fig_feat, use_container_width=True, key="feature_map_chart")

        fig_trace = create_intensity_trace(
            features,
            selected_id=selected_id,
            theoretical_mass=float(th_mass_input) if th_mass_input else 0.0,
            theoretical_charge=int(th_charge_input) if th_charge_input else 0,
        )
        st.plotly_chart(fig_trace, use_container_width=True, key="intensity_trace_chart")
    else:
        st.info("Upload a feature table CSV or load one of the **Test Datasets / Files** to see the feature map.")
