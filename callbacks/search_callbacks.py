"""
Search callbacks:
  - Run Search button → store-search-results, store-matched-ions,
                        results-table, search-result-summary, top-hit-summary
  - Results table row click → store-selected-result, ion-table
"""
from dash import Input, Output, State, no_update, html
import dash

from src.data.models import Spectrum, SearchResult, FragmentIon
from src.analysis.proteoform_search import run_targeted_search
from src.data.amino_acids import AA_MASSES


def register_callbacks(app):

    # ── Run targeted search ────────────────────────────────────────────────
    @app.long_callback(
        Output('store-search-results',  'data'),
        Output('store-matched-ions',    'data'),
        Output('results-table',         'data'),
        Output('search-result-summary', 'children'),
        Output('search-status',         'children'),
        Output('top-hit-summary',       'children'),
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
        progress=[
            Output('search-progress-bar', 'value'),
            Output('search-status',       'children'),
        ],
        running=[
            (Output('run-search-btn',      'disabled'), True,  False),
            (Output('search-progress-bar', 'value'),    0,     0),
            (Output('search-progress-bar', 'style'),
             {'height': '10px', 'display': 'block', 'marginBottom': '4px'},
             {'height': '10px', 'display': 'none'}),
        ],
        prevent_initial_call=True,
    )
    def run_search(set_progress, n_clicks, spectra_data, scan_idx, prot_data,
                   tol, max_z, ion_types, vmods, do_trunc, do_mods):

        def _prog(pct, label=''):
            set_progress((pct, f'Searching… {pct}%'))

        if not spectra_data:
            return (no_update,) * 5 + ('⚠ Load a spectrum first.',)
        if not prot_data or not prot_data.get('sequence'):
            return (no_update,) * 5 + ('⚠ Enter a protein sequence.',)

        _prog(5, '5%')
        scan_idx = scan_idx or 0
        spectrum = Spectrum.from_dict(spectra_data[scan_idx])
        seq      = ''.join(c for c in prot_data['sequence'] if c in AA_MASSES)
        name     = prot_data.get('name', 'Protein')

        _prog(15, '15%')
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
            progress_cb     = _prog,
        )
        _prog(90, '90%')

        if not results:
            return ([], [], [], 'No results found.', 'No matches.', 'No results.')

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
            table_rows.append({
                'rank':     rank,
                'sequence': pf.sequence[:30] + ('…' if len(pf.sequence) > 30 else ''),
                'range':    f"{pf.start_pos}–{pf.end_pos}",
                'score':    f"{pf.score:.2f}",
                'matched':  f"{pf.matched_ions}/{pf.total_ions}",
                'coverage': f"{r.sequence_coverage:.1f}",
                'th_mass':  f"{pf.theoretical_mass:.4f}",
                'obs_mass': f"{pf.observed_mass:.4f}",
                'ppm_error': f"{pf.mass_error_ppm:+.2f}",
                'ptms':     mods_str,
            })

        # Summary
        n = len(results)
        summary = (f"{n} candidate(s) found for {name} | "
                   f"Scan {spectrum.scan_id} | "
                   f"Tolerance {tol} ppm")

        # Top hit summary (sidebar)
        pf0 = top.proteoform
        top_text = html.Div([
            html.Div(f"Score: {pf0.score:.2f}",           className='mb-1'),
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

        return results_store, ions_data, table_rows, summary, f'✓ {n} hits', top_text

    # ── Select row in results table → update selected proteoform ──────────
    @app.callback(
        Output('store-selected-result', 'data'),
        Output('store-matched-ions',    'data', allow_duplicate=True),
        Output('ion-table',             'data'),
        Input('results-table', 'selected_rows'),
        State('store-search-results',  'data'),
        prevent_initial_call=True,
    )
    def on_result_select(selected_rows, results_data):
        if not selected_rows or not results_data:
            return no_update, no_update, no_update

        idx = selected_rows[0]
        if idx >= len(results_data):
            return no_update, no_update, no_update

        sr   = SearchResult.from_dict(results_data[idx])
        pf   = sr.proteoform

        ions_data = [ion.to_dict() for ion in sr.fragment_ions]

        # Ion table
        ion_rows = []
        for ion in sorted(sr.fragment_ions, key=lambda x: (x.ion_type, x.position)):
            ion_rows.append({
                'ion':     ion.label(),
                'th_mz':  f"{ion.mz:.4f}",
                'obs_mz': f"{ion.observed_mz:.4f}" if ion.matched else '—',
                'ppm':    f"{ion.mass_error_ppm:+.2f}" if ion.matched else '—',
                'charge': str(ion.charge),
                'matched': '✓' if ion.matched else '✗',
                'seq':     ion.sequence[:15] if ion.sequence else '',
            })

        return pf.to_dict(), ions_data, ion_rows
