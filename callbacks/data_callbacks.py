"""
Data-loading callbacks:
  - Upload spectrum file  → store-spectra, scan-selector options
  - Upload feature table  → store-features
  - Load demo button      → both stores
  - Scan selector change  → store-selected-scan-idx, scan-info
  - Protein text input    → store-protein, protein-mass-display
  - Manual mass input     → mass-diff-display, mod-suggestions
"""
import json
import os
import functools
from dash import Input, Output, State, no_update

from src.data.parsers import (decode_upload, generate_demo_spectrum,
                               generate_demo_features, UBIQUITIN, UBIQUITIN_MASS,
                               parse_pcml, parse_mzml, parse_fasta,
                               get_demo_proteins)

_DEMO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo_data')


@functools.lru_cache(maxsize=8)
def _load_demo_file(filename: str):
    """Return (spectra, feats, pinfo) for a file in demo_data/."""
    path = os.path.join(_DEMO_DIR, filename)
    with open(path, 'rb') as fh:
        raw = fh.read()
    if filename.lower().endswith('.pcml'):
        return parse_pcml(raw, filename)
    if filename.lower().endswith('.mzml'):
        spectra = parse_mzml(raw, filename)
        # Load matching _features.csv if it exists
        base = os.path.splitext(filename)[0]
        feat_path = os.path.join(_DEMO_DIR, f'{base}_features.csv')
        feats = []
        pinfo = {}
        if os.path.exists(feat_path):
            from src.data.parsers import parse_feature_table
            with open(feat_path, 'r') as fh:
                feats = parse_feature_table(fh.read(), feat_path)
            # Extract protein name / sequence from first feature that has one
            for f in feats:
                if f.sequence:
                    pinfo = {
                        'name': base.replace('_', ' ').title(),
                        'sequence': f.sequence,
                    }
                    break
        if not pinfo:
            pinfo = {'name': base.replace('_', ' ').title(), 'sequence': ''}
        return spectra, feats, pinfo
    return [], [], {}
from src.data.models import Spectrum
from src.analysis.mass_utils import (calc_sequence_mass, suggest_modifications,
                                      ppm_error)
from src.data.amino_acids import AA_MASSES


def register_callbacks(app):

    # ── Load dataset — fills spectrum, features, protein name+sequence ─────
    @app.callback(
        Output('store-spectra',         'data'),
        Output('scan-selector',         'options'),
        Output('scan-selector',         'value'),
        Output('upload-spectrum-status','children'),
        Output('protein-name',          'value'),
        Output('protein-sequence',      'value'),
        Input('demo-file-select', 'value'),
        Input('load-demo-btn', 'n_clicks'),
        prevent_initial_call=True,
    )
    def load_demo(selected_file, _):
        if selected_file:
            spectra, feats, pinfo = _load_demo_file(selected_file)
            if not spectra and not feats and not (pinfo or {}).get('sequence'):
                return [], [], None, f'No data found in {selected_file}', no_update, no_update
            opts = [{'label': f"{s.scan_id}  RT={s.retention_time:.2f} min  "
                              f"precursor={s.precursor_mz:.4f} m/z  z={s.precursor_charge}",
                     'value': i}
                    for i, s in enumerate(spectra)]
            scan_val = 0 if spectra else no_update
            prot_name = pinfo.get('name', '') if pinfo else no_update
            prot_seq  = pinfo.get('sequence', '') if pinfo else no_update
            msg = f'Loaded {len(spectra)} spectrum/a from {selected_file}'
            return [s.to_dict() for s in spectra], opts, scan_val, msg, prot_name, prot_seq
        # fallback: original ubiquitin demo
        demo = generate_demo_spectrum()
        spectra_data = [demo.to_dict()]
        opts = [{'label': f"{demo.scan_id}  RT={demo.retention_time:.2f} min  "
                          f"precursor={demo.precursor_mz:.4f} m/z  z={demo.precursor_charge}",
                 'value': 0}]
        src_note = (
            'Loaded: real instrument file (online)'
            if 'computed' not in demo.file_name
            else 'Loaded: computed ECD spectrum of human ubiquitin (Zubarev 1998)'
        )
        return spectra_data, opts, 0, src_note, 'Ubiquitin', UBIQUITIN

    # ── Upload spectrum (also handles PCML files) ────────────────────────
    @app.callback(
        Output('store-spectra',          'data', allow_duplicate=True),
        Output('scan-selector',          'options', allow_duplicate=True),
        Output('scan-selector',          'value', allow_duplicate=True),
        Output('upload-spectrum-status', 'children', allow_duplicate=True),
        Output('store-features',         'data', allow_duplicate=True),
        Output('upload-features-status', 'children', allow_duplicate=True),
        Output('protein-name',           'value', allow_duplicate=True),
        Output('protein-sequence',       'value', allow_duplicate=True),
        Input('upload-spectrum', 'contents'),
        State('upload-spectrum', 'filename'),
        prevent_initial_call=True,
    )
    def load_spectrum(contents, filename):
        if not contents:
            return (no_update,) * 8

        spectra, feats, msg, pinfo = decode_upload(contents, filename or 'file')
        if not spectra and not feats and not (pinfo or {}).get('sequence'):
            return [], [], None, f'Error: {msg}', no_update, no_update, no_update, no_update

        # spectra
        data = [s.to_dict() for s in spectra]
        opts = [{'label': f"{s.scan_id}  RT={s.retention_time:.2f}  "
                          f"precursor={s.precursor_mz:.3f}",
                 'value': i}
                for i, s in enumerate(spectra)]
        scan_val = 0 if spectra else no_update

        # features (only update store when PCML contained features)
        feats_data   = [f.to_dict() for f in feats] if feats else no_update
        feats_status = f'{len(feats)} features loaded (PCML)' if feats else no_update

        # protein (only populate when PCML contained protein info)
        prot_name = pinfo.get('name', '') if pinfo else no_update
        prot_seq  = pinfo.get('sequence', '') if pinfo else no_update

        return data, opts, scan_val, f'Loaded: {msg}', feats_data, feats_status, prot_name, prot_seq

    # ── Load dataset — also fills features ──────────────────────────────────
    @app.callback(
        Output('store-features', 'data'),
        Output('upload-features-status', 'children'),
        Input('demo-file-select', 'value'),
        Input('load-demo-btn', 'n_clicks'),
        prevent_initial_call=True,
    )
    def load_demo_features(selected_file, _):
        if selected_file:
            _, feats, _ = _load_demo_file(selected_file)
            if feats:
                return [f.to_dict() for f in feats], f'{len(feats)} features loaded from {selected_file}'
            return [], f'No features in {selected_file}'
        feats = generate_demo_features()
        return [f.to_dict() for f in feats], f'{len(feats)} features loaded (demo)'

    # ── Upload features ────────────────────────────────────────────────────
    @app.callback(
        Output('store-features',          'data', allow_duplicate=True),
        Output('upload-features-status',  'children', allow_duplicate=True),
        Input('upload-features', 'contents'),
        State('upload-features', 'filename'),
        prevent_initial_call=True,
    )
    def load_features(contents, filename):
        if not contents:
            return no_update, no_update
        _, features, msg, _ = decode_upload(contents, filename or 'file')
        if not features:
            return [], f'Error: {msg}'
        return [f.to_dict() for f in features], f'Loaded: {msg}'

    # ── Scan info on selection ─────────────────────────────────────────────
    @app.callback(
        Output('store-selected-scan-idx', 'data'),
        Output('scan-info', 'children'),
        Input('scan-selector', 'value'),
        State('store-spectra', 'data'),
        prevent_initial_call=True,
    )
    def on_scan_select(idx, spectra_data):
        if idx is None or not spectra_data:
            return 0, ''
        d = spectra_data[idx]
        n_peaks = len(d['mz'])
        info = (f"Scan {d['scan_id']} | {n_peaks} peaks | "
                f"RT {d['retention_time']:.2f} min | "
                f"Prec {d['precursor_mz']:.3f} m/z z={d['precursor_charge']}")
        return idx, info

    # ── Protein input → store-protein, mass display ────────────────────────
    @app.callback(
        Output('store-protein', 'data'),
        Output('protein-mass-display', 'children'),
        Input('protein-sequence', 'value'),
        Input('protein-name', 'value'),
    )
    def on_protein_input(seq, name):
        seq  = (seq or '').strip().upper()
        name = (name or '').strip()
        if not seq:
            return {}, ''
        # Clean sequence
        seq = ''.join(c for c in seq if c in AA_MASSES)
        mass = calc_sequence_mass(seq)
        prot_data = {'name': name or 'Protein', 'sequence': seq, 'mass': mass}
        info = (f"{len(seq)} aa | Monoisotopic mass: {mass:.4f} Da")
        return prot_data, info

    # ── N-terminal modification quick-picker ──────────────────────────────
    _NTERM_PRESETS = {
        'nterm-mod-none':  ('No Modification',  0.0),
        'nterm-mod-ac':    ('Acetylation',      42.0106),
        'nterm-mod-fo':    ('Formylation',      27.9949),
        'nterm-mod-tri':   ('Trimethylation',   42.0470),
        'nterm-mod-palm':  ('Palmitate',       238.2297),
    }

    @app.callback(
        Output('store-nterm-mod',   'data'),
        Output('nterm-mod-display', 'children'),
        Input('nterm-mod-none',        'n_clicks'),
        Input('nterm-mod-ac',          'n_clicks'),
        Input('nterm-mod-fo',          'n_clicks'),
        Input('nterm-mod-tri',         'n_clicks'),
        Input('nterm-mod-palm',        'n_clicks'),
        Input('nterm-mod-custom-btn',  'n_clicks'),
        State('nterm-mod-custom-mass', 'value'),
        prevent_initial_call=True,
    )
    def set_nterm_mod(*args):
        from dash import ctx
        custom_mass = args[-1]
        btn_id = ctx.triggered_id
        if btn_id == 'nterm-mod-custom-btn':
            shift = float(custom_mass) if custom_mass is not None else 0.0
            name  = f'Custom ({shift:+.4f} Da)'
        else:
            name, shift = _NTERM_PRESETS.get(btn_id, ('No Modification', 0.0))
        store = {'name': name, 'mass_shift': shift}
        if shift == 0.0:
            display = 'N-term mod: None'
        else:
            display = f'N-term mod: {name}  ({shift:+.4f} Da)'
        return store, display

    # ── Manual mass → mass-diff-display, mod suggestions ─────────────────
    @app.callback(
        Output('mass-diff-display', 'children'),
        Output('mod-suggestions', 'children'),
        Input('calc-mass-diff-btn', 'n_clicks'),
        State('manual-mass', 'value'),
        State('store-protein', 'data'),
        prevent_initial_call=True,
    )
    def calc_mass_diff(n_clicks, obs_mass_val, prot_data):
        if obs_mass_val is None or not prot_data:
            return 'Enter a protein sequence and observed mass first.', ''
        obs  = float(obs_mass_val)
        th   = prot_data.get('mass', 0)
        diff = obs - th
        ppm  = ppm_error(obs, th) if th else 0.0
        diff_text = (f"Δmass = {diff:+.4f} Da  ({ppm:+.1f} ppm)")
        sugg = suggest_modifications(diff)
        if not sugg:
            sugg_el = 'No matching PTMs found.'
        else:
            items = [
                f"• {s['name']}  ({s['mass_shift']:+.4f} Da,  "
                f"residual {s['residual_da']:+.5f} Da)"
                for s in sugg[:6]
            ]
            sugg_el = html_list(items)
        return diff_text, sugg_el

    # ── Load Demo DB button → store-fasta-proteins ───────────────────────
    @app.callback(
        Output('store-fasta-proteins', 'data',   allow_duplicate=True),
        Output('upload-fasta-status',  'children', allow_duplicate=True),
        Input('load-demo-db-btn', 'n_clicks'),
        prevent_initial_call=True,
    )
    def load_demo_db(_):
        proteins = get_demo_proteins()
        data = [[name, seq] for name, seq in proteins]
        names = ', '.join(p[0].split('|')[2].split('_')[0] if '|' in p[0] else p[0]
                          for p in proteins)
        return data, f'✓ {len(proteins)} proteins loaded: {names}'

    # ── FASTA upload → store-fasta-proteins ───────────────────────────────
    @app.callback(
        Output('store-fasta-proteins', 'data'),
        Output('upload-fasta-status',  'children'),
        Input('upload-fasta', 'contents'),
        State('upload-fasta', 'filename'),
        prevent_initial_call=True,
    )
    def load_fasta_db(contents, filename):
        if not contents:
            return no_update, no_update
        import base64 as _b64
        try:
            _ctype, content_str = contents.split(',', 1)
            decoded = _b64.b64decode(content_str).decode('utf-8', errors='replace')
        except Exception as e:
            return [], f'Error decoding {filename}: {e}'
        proteins = parse_fasta(decoded)
        if not proteins:
            return [], f'No protein entries found in {filename}'
        # Hard-cap at 500 proteins for browser performance
        capped = proteins[:500]
        cap_note = f' (capped at 500)' if len(proteins) > 500 else ''
        data = [[name, seq] for name, seq in capped]
        return data, f'{len(capped)} protein(s) loaded{cap_note} from {filename}'

    # ── Mode toggle → show/hide targeted vs database controls ─────────────
    @app.callback(
        Output('targeted-protein-controls', 'style'),
        Output('database-protein-controls', 'style'),
        Input('search-mode', 'value'),
    )
    def toggle_search_mode(mode):
        if mode == 'database':
            return {'display': 'none'}, {}
        return {}, {'display': 'none'}


def html_list(items):
    from dash import html
    return html.Ul([html.Li(i, style={'marginBottom': '2px'}) for i in items],
                   style={'paddingLeft': '12px', 'margin': '0'})
