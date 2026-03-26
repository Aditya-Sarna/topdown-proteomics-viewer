"""
Diagnostics tab callbacks:
  - TIC graph               (store-spectra + selected scan idx)
  - TIC click → scan select (updates store-selected-scan-idx)
  - Spectrum QC 4-panel     (current scan + matched ions)
  - Ion type breakdown      (matched ions)
  - Precursor envelope      (current scan precursor m/z + z)
  - FDR curve               (search results)
  - Sequence tag map        (selected proteoform + matched ions)
  - Spectrum stats bar      (metrics pill row above spectrum graph)
"""
from dash import Input, Output, State, no_update, html
import dash_bootstrap_components as dbc

from src.viz.diagnostics_plots import (
    create_tic_plot,
    create_spectrum_qc,
    create_ion_breakdown,
    create_fdr_curve,
    create_precursor_envelope,
    create_sequence_tag_stats,
)


def _metric_pill(label: str, value: str, color: str = '#1a73e8') -> dbc.Badge:
    return html.Span([
        html.Span(label + ': ', style={'color': '#5f6368', 'fontSize': '0.72rem'}),
        html.Span(value,        style={'color': color,    'fontSize': '0.72rem',
                                       'fontWeight': '600'}),
    ], className='me-3')


def register_callbacks(app):

    # ── TIC graph ────────────────────────────────────────────────────────────
    @app.callback(
        Output('tic-graph', 'figure'),
        Input('store-spectra',            'data'),
        Input('store-selected-scan-idx',  'data'),
    )
    def update_tic(spectra_data, scan_idx):
        if not spectra_data:
            return create_tic_plot([])
        scan_idx = scan_idx or 0
        rt_sel   = None
        if 0 <= scan_idx < len(spectra_data):
            rt_sel = spectra_data[scan_idx].get('retention_time')
        return create_tic_plot(spectra_data, selected_rt=rt_sel, selected_idx=scan_idx)

    # ── TIC click → select that scan ─────────────────────────────────────────
    @app.callback(
        Output('store-selected-scan-idx', 'data', allow_duplicate=True),
        Output('scan-selector',           'value',  allow_duplicate=True),
        Input('tic-graph', 'clickData'),
        State('store-spectra', 'data'),
        prevent_initial_call=True,
    )
    def tic_select_scan(click_data, spectra_data):
        if not click_data or not spectra_data:
            return no_update, no_update
        clicked_rt = click_data['points'][0].get('x', None)
        if clicked_rt is None:
            return no_update, no_update
        rts   = [s.get('retention_time', 0.0) for s in spectra_data]
        import numpy as np
        idx   = int(np.argmin(np.abs(np.array(rts) - float(clicked_rt))))
        return idx, idx

    # ── Spectrum QC panel ────────────────────────────────────────────────────
    @app.callback(
        Output('diag-qc-graph', 'figure'),
        Input('store-selected-scan-idx', 'data'),
        Input('store-spectra',           'data'),
        Input('store-matched-ions',      'data'),
    )
    def update_qc(scan_idx, spectra_data, ions_data):
        if not spectra_data:
            return create_spectrum_qc({})
        scan_idx = scan_idx or 0
        s_dict   = spectra_data[scan_idx] if scan_idx < len(spectra_data) else {}
        return create_spectrum_qc(s_dict, ions_data or [])

    # ── Ion type breakdown ───────────────────────────────────────────────────
    @app.callback(
        Output('diag-ion-breakdown', 'figure'),
        Input('store-matched-ions', 'data'),
    )
    def update_ion_breakdown(ions_data):
        return create_ion_breakdown(ions_data or [])

    # ── Precursor envelope ───────────────────────────────────────────────────
    @app.callback(
        Output('diag-precursor-envelope', 'figure'),
        Input('store-spectra',           'data'),
        Input('store-selected-scan-idx', 'data'),
    )
    def update_precursor_envelope(spectra_data, scan_idx):
        if not spectra_data:
            return create_precursor_envelope({})
        scan_idx = scan_idx or 0
        s_dict   = spectra_data[scan_idx] if scan_idx < len(spectra_data) else {}
        prec_mz  = float(s_dict.get('precursor_mz',     0.0))
        prec_z   = int(s_dict.get('precursor_charge', 0))
        return create_precursor_envelope(s_dict, prec_mz, prec_z)

    # ── FDR curve ────────────────────────────────────────────────────────────
    @app.callback(
        Output('fdr-curve-graph', 'figure'),
        Input('store-search-results', 'data'),
    )
    def update_fdr_curve(results_data):
        return create_fdr_curve(results_data or [])

    # ── Sequence tag map (diagnostics tab) ───────────────────────────────────
    @app.callback(
        Output('diag-seq-tag-map', 'figure'),
        Input('store-selected-result', 'data'),
        Input('store-matched-ions',    'data'),
    )
    def update_seq_tag_map(selected_result, ions_data):
        seq = ''
        if selected_result:
            seq = selected_result.get('sequence', '')
        return create_sequence_tag_stats(seq, ions_data or [])

    # ── Spectrum stats bar ───────────────────────────────────────────────────
    @app.callback(
        Output('spectrum-stats-bar', 'children'),
        Input('store-spectra',            'data'),
        Input('store-selected-scan-idx',  'data'),
        Input('store-matched-ions',       'data'),
    )
    def update_spectrum_stats(spectra_data, scan_idx, ions_data):
        if not spectra_data:
            return []

        import numpy as np
        scan_idx = scan_idx or 0
        s_dict   = spectra_data[scan_idx] if scan_idx < len(spectra_data) else {}

        mz_arr  = np.asarray(s_dict.get('mz',        []), dtype=np.float64)
        int_arr = np.asarray(s_dict.get('intensity', []), dtype=np.float64)

        n_peaks   = mz_arr.size
        base_peak = float(mz_arr[np.argmax(int_arr)]) if n_peaks > 0 else 0.0
        dyn_range = (int_arr.max() / max(1.0, np.percentile(int_arr, 1))
                     if n_peaks > 1 else 1.0)

        # Approximate SNR: median of top 10% / median of bottom 50%
        if n_peaks > 10:
            top10  = np.percentile(int_arr, 90)
            bot50  = np.median(int_arr)
            snr    = top10 / max(1.0, bot50)
        else:
            snr = float('nan')

        # Matched ions
        ions_data = ions_data or []
        n_matched = sum(1 for d in ions_data if d.get('matched'))
        ppms      = [d.get('mass_error_ppm', 0.0) for d in ions_data if d.get('matched')]
        median_ppm = float(np.median(ppms)) if ppms else float('nan')

        ms_level  = int(s_dict.get('ms_level', 2))
        prec_mz   = float(s_dict.get('precursor_mz',    0.0))
        prec_z    = int(s_dict.get('precursor_charge', 0))
        rt        = float(s_dict.get('retention_time',  0.0))

        import math
        def _fmt(v, fmt='.1f'):
            return '—' if (isinstance(v, float) and math.isnan(v)) else f'{v:{fmt}}'

        ppm_color = '#2e7d32' if ppms and abs(median_ppm) <= 5 else '#e53935'
        dyn_exp   = int(np.log10(max(dyn_range, 1)))

        pills = [
            _metric_pill('RT',      f'{rt:.2f} min'),
            _metric_pill('Level',   f'MS{ms_level}'),
            _metric_pill('Peaks',   str(n_peaks)),
            _metric_pill('BPC m/z', f'{base_peak:.3f}'),
            _metric_pill('Dyn. range', f'10^{dyn_exp}'),
        ]
        if not math.isnan(snr):
            pills.append(_metric_pill('SNR', _fmt(snr, '.0f')))
        if prec_mz > 0:
            pills.append(_metric_pill('Prec. m/z', f'{prec_mz:.4f}'))
        if prec_z > 0:
            pills.append(_metric_pill('Prec. z', str(prec_z)))
        if n_matched:
            pills.append(_metric_pill('Matched ions', str(n_matched)))
        if ppms:
            pills.append(_metric_pill('Median Δm (ppm)', _fmt(median_ppm, '+.2f'), ppm_color))

        return html.Div(pills,
                        style={'display': 'flex', 'flexWrap': 'wrap',
                               'alignItems': 'center',
                               'background': '#f8f9fa',
                               'borderRadius': '6px',
                               'padding': '6px 12px',
                               'border': '1px solid #dadce0',
                               'gap': '0'})
