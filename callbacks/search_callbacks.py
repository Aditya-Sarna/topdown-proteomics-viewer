"""  
Search callbacks:
  - Run Search button → store-search-results, store-matched-ions,
                        results-table, search-result-summary, top-hit-summary
  - Results table row click → store-selected-result, ion-table
"""
from dash import Input, Output, State, no_update, html
import dash

from src.data.models import Spectrum, SearchResult, FragmentIon
from src.analysis.proteoform_search import run_targeted_search, run_database_search
from src.analysis.deconvolution import deconvolute_spectrum
from src.data.amino_acids import AA_MASSES

_PROTON = 1.007276


def _build_ion_row(ion):
    """Convert a FragmentIon to a table row matching the reference column layout."""
    obs_mass = (
        ion.observed_mz * ion.charge - ion.charge * _PROTON
        if (ion.matched and ion.charge > 0) else 0.0
    )
    return {
        'name':     ion.label(),
        'ion_type': ion.ion_type,
        'ion_num':  ion.position,
        'th_mass':  f"{ion.mass:.4f}",
        'obs_mass': f"{obs_mass:.4f}" if ion.matched else '—',
        'da_err':   f"{ion.mass_error_da:+.4f}" if ion.matched else '—',
        'ppm_err':  f"{ion.mass_error_ppm:+.2f}" if ion.matched else '—',
        'matched':  '1' if ion.matched else '0',  # hidden; used for conditional styling
    }


def register_callbacks(app):

    # ── Run targeted search ────────────────────────────────────────────────
    @app.callback(
        Output('store-search-results',  'data'),
        Output('store-matched-ions',    'data'),
        Output('results-table',         'data'),
        Output('ion-table',             'data', allow_duplicate=True),
        Output('search-result-summary', 'children'),
        Output('search-status',         'children'),
        Output('top-hit-summary',       'children'),
        Output('store-selected-result', 'data', allow_duplicate=True),
        Input('run-search-btn', 'n_clicks'),
        State('store-spectra',           'data'),
        State('store-selected-scan-idx', 'data'),
        State('store-protein',           'data'),
        State('tolerance-ppm',           'value'),
        State('max-charge',              'value'),
        State('ion-types',               'value'),
        State('variable-mods',           'value'),
        State('search-truncations',      'value'),
        State('search-mods',             'value'),
        State('search-mode',             'value'),
        State('store-fasta-proteins',    'data'),
        State('deconvolute-spectrum',    'value'),
        State('manual-mass',             'value'),
        State('store-nterm-mod',         'data'),
        prevent_initial_call=True,
    )
    def run_search(n_clicks, spectra_data, scan_idx, prot_data,
                   tol, max_z, ion_types, vmods, do_trunc, do_mods,
                   search_mode, fasta_proteins, do_deconv, manual_mass, nterm_mod):
        if not spectra_data:
            return (no_update,) * 7 + ('⚠ Load a spectrum first.',)

        scan_idx = scan_idx or 0
        spectrum = Spectrum.from_dict(spectra_data[scan_idx])

        # ── Optional charge deconvolution via pyopenms.Deisotoper ──────────
        deconv_note = ''
        if do_deconv:
            dr = deconvolute_spectrum(
                spectrum,
                fragment_tolerance_ppm=float(tol or 10),
                max_charge=int(max_z or 20),
            )
            spectrum = dr.spectrum
            if dr.used_openms:
                deconv_note = (f' [Deconvoluted: {dr.n_original_peaks}→'
                               f'{dr.n_deconvoluted_peaks} peaks via OpenMS]')

        # ── Database search mode ────────────────────────────────────────────
        if search_mode == 'database':
            if not fasta_proteins:
                return (no_update,) * 7 + ('⚠ Upload a FASTA file first.',)
            proteins = [(p[0], p[1]) for p in fasta_proteins]
            results = run_database_search(
                spectrum       = spectrum,
                proteins       = proteins,
                ion_types      = ion_types or ['b', 'y', 'c', 'z'],
                max_charge     = int(max_z or 4),
                tolerance_ppm  = float(tol or 10),
                search_truncations = bool(do_trunc),
                top_n          = 25,
            )
            search_label = f'{len(proteins)} proteins (database)'
        else:
            # ── Targeted search mode ────────────────────────────────────────
            if not prot_data or not prot_data.get('sequence'):
                return (no_update,) * 7 + ('⚠ Enter a protein sequence.',)
            seq  = ''.join(c for c in prot_data['sequence'] if c in AA_MASSES)
            name = prot_data.get('name', 'Protein')
            nterm_shift = float((nterm_mod or {}).get('mass_shift', 0.0))
            nterm_name  = (nterm_mod or {}).get('name', '')
            results = run_targeted_search(
                spectrum        = spectrum,
                protein_sequence= seq,
                protein_name    = name,
                ion_types       = ion_types or ['b', 'y', 'c', 'z'],
                max_charge      = int(max_z or 4),
                tolerance_ppm   = float(tol or 10),
                search_truncations   = bool(do_trunc),
                search_modifications = bool(do_mods),
                variable_mods   = vmods or [],
                obs_mass_override= float(manual_mass) if manual_mass else 0.0,
            )
            search_label = name

        if not results:
            return [], [], [], [], 'No results found.', 'No matches.', 'No results.', {}

        # ── Apply N-term mod as annotation overlay on the top result ───────
        if nterm_shift != 0.0 and nterm_name and results:
            from src.data.models import Modification
            from src.analysis.mass_utils import calc_sequence_mass, ppm_error as _ppm_err
            from src.analysis.fragment_ions import calc_ions
            from src.analysis.peak_matching import match_peaks as _match
            top_pf = results[0].proteoform
            if top_pf.sequence:
                _mods = [m for m in top_pf.modifications if m.position != 1]
                _mods.append(Modification(1, nterm_name, nterm_shift, top_pf.sequence[0]))
                top_pf.modifications = _mods
                top_pf.theoretical_mass = round(calc_sequence_mass(top_pf.sequence, _mods), 4)
                top_pf.mass_error_da = round(top_pf.observed_mass - top_pf.theoretical_mass, 4)
                if top_pf.theoretical_mass:
                    top_pf.mass_error_ppm = round(_ppm_err(top_pf.observed_mass, top_pf.theoretical_mass), 2)
                _mod_map = {m.position: m.mass_shift for m in _mods}
                _itypes = ion_types or ['b', 'y', 'c', 'z']
                _new_ions = calc_ions(top_pf.sequence, _itypes, _mod_map, int(max_z or 4))
                results[0].fragment_ions = list(_match(_new_ions, spectrum, float(tol or 10)))

        # Store results
        results_store = [r.to_dict() for r in results]

        # Matched ions of top hit
        top       = results[0]
        ions_data = [ion.to_dict() for ion in top.fragment_ions]

        # Table rows
        table_rows = []
        for rank, r in enumerate(results, 1):
            pf = r.proteoform
            mods_str = '; '.join(
                f"{m.name}@{m.position}" for m in pf.modifications
            ) or '—'
            # Format e-value in scientific notation
            ev = r.e_value
            if ev < 0.001:
                ev_str = f'{ev:.2e}'
            else:
                ev_str = f'{ev:.4f}'
            # Trim [DECOY] suffix from protein name display
            prot_display = pf.protein_name.replace(' [DECOY]', '')
            table_rows.append({
                'rank':       rank,
                'protein':    prot_display[:25] + ('…' if len(prot_display) > 25 else ''),
                'sequence':   pf.sequence[:30] + ('…' if len(pf.sequence) > 30 else ''),
                'range':      f"{pf.start_pos}–{pf.end_pos}",
                'matched_aa': r.matched_aa,
                'n_tags':     r.n_tags,
                'score':      f"{pf.score:.2f}",
                'e_value':    ev_str,
                'q_value':    f"{r.q_value:.4f}",
                'matched':    f"{pf.matched_ions}/{pf.total_ions}",
                'coverage':   f"{r.sequence_coverage:.1f}",
                'th_mass':    f"{pf.theoretical_mass:.4f}",
                'obs_mass':   f"{pf.observed_mass:.4f}",
                'ppm_error':  f"{pf.mass_error_ppm:+.2f}",
                'ptms':       mods_str,
            })

        # Summary
        n = len(results)
        summary = (f"{n} candidate(s) found for {search_label} | "
                   f"Scan {spectrum.scan_id} | "
                   f"Tolerance {tol} ppm"
                   f"{deconv_note}")

        # Top hit summary (sidebar)
        pf0 = top.proteoform
        ev0 = top.e_value
        ev0_str = f'{ev0:.2e}' if ev0 < 0.001 else f'{ev0:.4f}'
        qv_color = '#66BB6A' if top.q_value < 0.01 else ('#FFD54F' if top.q_value < 0.05 else '#EF5350')
        top_text = html.Div([
            html.Div(f"Score: {pf0.score:.2f}",           className='mb-1'),
            html.Div(f"E-value: {ev0_str}",
                     className='mb-1',
                     style={'color': '#66BB6A' if ev0 < 0.01 else '#FFD54F' if ev0 < 0.1 else '#EF5350'}),
            html.Div(f"q-value (FDR): {top.q_value:.4f}",
                     className='mb-1',
                     style={'color': qv_color}),
            html.Div(f"Proteoform: {pf0.start_pos}–{pf0.end_pos}",  className='mb-1'),
            html.Div(f"Ions matched: {pf0.matched_ions}/{pf0.total_ions}", className='mb-1'),
            html.Div(f"Coverage: {top.sequence_coverage:.1f}%",      className='mb-1'),
            html.Div(f"Th. mass: {pf0.theoretical_mass:.4f} Da",     className='mb-1'),
            html.Div(f"Δmass: {pf0.mass_error_ppm:+.2f} ppm",
                     style={'color': '#EF5350' if abs(pf0.mass_error_ppm) > 5 else '#66BB6A'}),
            html.Div(
                f"PTMs: {'; '.join(m.name for m in pf0.modifications) or 'none'}",
                className='mt-1', style={'color': '#FFD54F'}
            ) if pf0.modifications else html.Div(),
        ])

        # Ion rows for the top hit (auto-populated without requiring a row click)
        top_ion_rows = [
            _build_ion_row(ion)
            for ion in sorted(top.fragment_ions, key=lambda x: (x.ion_type, x.position))
        ]

        return results_store, ions_data, table_rows, top_ion_rows, summary, f'✓ {n} hits', top_text, top.proteoform.to_dict()

    # ── Select row in results table → update selected proteoform ──────────
    @app.callback(
        Output('store-selected-result', 'data'),
        Output('store-matched-ions',    'data', allow_duplicate=True),
        Output('ion-table',             'data'),
        Input('results-table', 'selected_rows'),
        State('store-search-results',  'data'),
        State('store-nterm-mod',       'data'),
        prevent_initial_call=True,
    )
    def on_result_select(selected_rows, results_data, nterm_mod):
        if not selected_rows or not results_data:
            return no_update, no_update, no_update

        idx = selected_rows[0]
        if idx >= len(results_data):
            return no_update, no_update, no_update

        sr   = SearchResult.from_dict(results_data[idx])
        pf   = sr.proteoform

        # Annotate proteoform with the currently-selected N-term mod (sequence view)
        nterm_shift = float((nterm_mod or {}).get('mass_shift', 0.0))
        nterm_name  = (nterm_mod or {}).get('name', '')
        if nterm_shift != 0.0 and nterm_name and pf.sequence:
            from src.data.models import Modification
            from src.analysis.mass_utils import calc_sequence_mass
            _mods = [m for m in pf.modifications if m.position != 1]
            _mods.append(Modification(1, nterm_name, nterm_shift, pf.sequence[0]))
            pf.modifications = _mods
            pf.theoretical_mass = round(calc_sequence_mass(pf.sequence, _mods), 4)

        ions_data = [ion.to_dict() for ion in sr.fragment_ions]

        # Ion table
        ion_rows = [
            _build_ion_row(ion)
            for ion in sorted(sr.fragment_ions, key=lambda x: (x.ion_type, x.position))
        ]

        return pf.to_dict(), ions_data, ion_rows

    # ── N-term mod overlay: immediately re-annotate current result ─────────
    @app.callback(
        Output('store-selected-result', 'data', allow_duplicate=True),
        Output('store-matched-ions',    'data', allow_duplicate=True),
        Input('store-nterm-mod', 'data'),
        State('store-selected-result',   'data'),
        State('store-spectra',           'data'),
        State('store-selected-scan-idx', 'data'),
        State('tolerance-ppm',           'value'),
        State('max-charge',              'value'),
        State('ion-types',               'value'),
        prevent_initial_call=True,
    )
    def apply_nterm_mod_overlay(nterm_mod, selected_result, spectra_data, scan_idx,
                                tol, max_z, ion_types):
        """Apply (or remove) the N-term mod on the currently displayed proteoform
        without re-running a full database/targeted search.  Gives instant visual
        feedback in the Sequence and Spectrum tabs."""
        if not selected_result or not spectra_data:
            return no_update, no_update

        from src.data.models import Proteoform, Modification
        from src.analysis.fragment_ions import calc_ions
        from src.analysis.peak_matching import match_peaks
        from src.analysis.mass_utils import calc_sequence_mass

        pf = Proteoform.from_dict(selected_result)
        if not pf.sequence:
            return no_update, no_update

        # Rebuild mod list: keep any mods that are NOT at position 1, then
        # optionally add the new N-term mod.
        mods = [m for m in pf.modifications if m.position != 1]
        nterm_shift = float((nterm_mod or {}).get('mass_shift', 0.0))
        nterm_name  = (nterm_mod or {}).get('name', '')
        if nterm_shift != 0.0 and nterm_name:
            mods.append(Modification(
                position=1, name=nterm_name,
                mass_shift=nterm_shift, residue=pf.sequence[0],
            ))

        # Update proteoform in-place (dataclass is mutable)
        pf.modifications = mods
        pf.theoretical_mass = round(calc_sequence_mass(pf.sequence, mods), 4)
        pf.mass_error_da  = round(pf.observed_mass - pf.theoretical_mass, 4)
        if pf.theoretical_mass:
            from src.analysis.mass_utils import ppm_error
            pf.mass_error_ppm = round(ppm_error(pf.observed_mass, pf.theoretical_mass), 2)

        # Recompute fragment ions with the updated mod map
        mod_map   = {m.position: m.mass_shift for m in mods}
        _ion_types = ion_types or ['b', 'y', 'c', 'z']
        ions = calc_ions(pf.sequence, _ion_types, mod_map, int(max_z or 4))

        spectrum = Spectrum.from_dict(spectra_data[scan_idx or 0])
        matched  = match_peaks(ions, spectrum, float(tol or 10))

        return pf.to_dict(), [ion.to_dict() for ion in matched]

    # ── Score distribution (updates whenever search results change) ─────────
    @app.callback(
        Output('score-dist-graph', 'figure'),
        Input('store-search-results', 'data'),
    )
    def update_score_dist(results_data):
        import plotly.graph_objects as go
        fig = go.Figure()
        if not results_data:
            fig.update_layout(title='Run search to see score distribution',
                              template='plotly_white', paper_bgcolor='#ffffff',
                              font=dict(color='#aaaaaa'), height=220,
                              margin=dict(t=35, b=35, l=40, r=10))
            return fig

        from src.data.models import SearchResult
        scores = [SearchResult.from_dict(r).proteoform.score for r in results_data]
        fig.add_trace(go.Histogram(
            x=scores, nbinsx=max(10, len(scores) // 2),
            marker_color='#1a73e8', opacity=0.8, name='Target',
        ))
        fig.update_layout(
            title='Score Distribution',
            xaxis_title='Score', yaxis_title='Count',
            template='plotly_white', paper_bgcolor='#ffffff',
            font=dict(color='#111111'), height=220,
            margin=dict(t=35, b=35, l=40, r=10),
            bargap=0.05,
        )
        return fig

    # ── Internal fragment map (updates on row selection or new search) ──────
    @app.callback(
        Output('internal-frag-graph', 'figure'),
        Input('store-selected-result', 'data'),
        Input('store-search-results',  'data'),
        State('store-spectra',           'data'),
        State('store-selected-scan-idx', 'data'),
    )
    def update_internal_frag_map(selected_result, results_data, spectra_data, scan_idx):
        import plotly.graph_objects as go
        from src.viz.sequence_plots import create_internal_fragment_map
        from src.analysis.fragment_ions import calc_internal_ions
        import numpy as np

        # Determine which proteoform/spectrum to use
        pf_data = selected_result
        if not pf_data and results_data:
            from src.data.models import SearchResult
            pf_data = SearchResult.from_dict(results_data[0]).proteoform.to_dict()

        if not pf_data or not pf_data.get('sequence') or not spectra_data:
            return create_internal_fragment_map('', [])

        seq      = pf_data['sequence']
        idx      = scan_idx or 0
        spectrum = Spectrum.from_dict(spectra_data[idx])
        obs_mz   = np.sort(spectrum.mz_array)

        internal = calc_internal_ions(seq, obs_mz, tolerance_ppm=20.0)
        return create_internal_fragment_map(seq, internal)

    # ── Fragment stats header (above ion table) ─────────────────────────────
    @app.callback(
        Output('fragment-stats-header', 'children'),
        Input('store-selected-result', 'data'),
        Input('store-search-results',  'data'),
    )
    def update_fragment_stats(selected_result, results_data):
        import dash_bootstrap_components as dbc

        matched  = 0
        coverage = 0.0

        if selected_result:
            from src.data.models import Proteoform
            pf = Proteoform.from_dict(selected_result)
            matched = pf.matched_ions
            # find matching SearchResult for coverage
            if results_data:
                for r in results_data:
                    sr = SearchResult.from_dict(r)
                    if sr.proteoform.sequence == pf.sequence:
                        coverage = sr.sequence_coverage
                        break
        elif results_data:
            sr       = SearchResult.from_dict(results_data[0])
            matched  = sr.proteoform.matched_ions
            coverage = sr.sequence_coverage

        if matched == 0 and not results_data:
            return html.Small('Run search to see fragment matches.',
                               className='text-muted')

        return dbc.Row([
            dbc.Col(
                html.Span(f"Matching fragments (# {matched})",
                          style={'fontWeight': '600', 'fontSize': '0.82rem',
                                 'color': '#111111'}),
                width='auto',
            ),
            dbc.Col(
                html.Span(f"% Residue cleavage: {coverage:.3f}%",
                          style={'fontSize': '0.82rem', 'color': '#555555'}),
                width='auto', className='ms-auto',
            ),
        ], className='align-items-center mb-1')

    # ── Inline row-detail expansion (Mass / Start / End / Description) ─────
    @app.callback(
        Output('result-row-detail', 'children'),
        Output('result-row-detail', 'style'),
        Input('results-table', 'selected_rows'),
        State('store-search-results', 'data'),
        prevent_initial_call=True,
    )
    def update_row_detail(selected_rows, results_data):
        import dash_bootstrap_components as dbc

        _hidden = {'display': 'none'}
        _visible = {
            'display': 'block',
            'backgroundColor': '#e8f0fe',
            'borderLeft': '3px solid #1a73e8',
            'borderBottom': '1px solid #b0c4de',
            'padding': '7px 14px',
            'marginBottom': '6px',
            'fontSize': '0.78rem',
            'borderRadius': '0 0 4px 0',
        }

        if not selected_rows or not results_data:
            return '', _hidden

        idx = selected_rows[0]
        if idx >= len(results_data):
            return '', _hidden

        sr = SearchResult.from_dict(results_data[idx])
        pf = sr.proteoform

        # Parse the protein name for a clean accession + description split
        full_name = pf.protein_name.replace(' [DECOY]', '')
        parts = full_name.split(' ', 1)
        accession   = parts[0]
        description = parts[1] if len(parts) > 1 else full_name

        def _field(label, value):
            return dbc.Col([
                html.Span(label, style={'fontWeight': '600', 'marginRight': '4px',
                                        'color': '#1a73e8'}),
                html.Span(str(value), style={'color': '#111111'}),
            ], width='auto', className='me-4 mb-1')

        return dbc.Row([
            _field('Mass',           f'{pf.theoretical_mass:.7f}'),
            _field('Start Position', pf.start_pos),
            _field('End Position',   pf.end_pos),
            dbc.Col([
                html.Span('Description', style={'fontWeight': '600', 'marginRight': '4px',
                                                'color': '#1a73e8'}),
                html.Span(description,
                          style={'color': '#111111', 'wordBreak': 'break-word'}),
            ], width=True),
        ], className='gx-2 align-items-start flex-wrap'), _visible
