"""
Export modal and download callbacks.
"""
import json
import csv
import io
from dash import Input, Output, State, dcc, no_update


def register_callbacks(app):

    # Toggle modal open/close
    @app.callback(
        Output('export-modal', 'is_open'),
        Input('export-btn',      'n_clicks'),
        Input('export-close-btn','n_clicks'),
        State('export-modal',    'is_open'),
        prevent_initial_call=True,
    )
    def toggle_modal(open_click, close_click, is_open):
        return not is_open

    # Download
    @app.callback(
        Output('download-data', 'data'),
        Input('export-confirm-btn', 'n_clicks'),
        State('export-format',       'value'),
        State('store-search-results','data'),
        State('store-matched-ions',  'data'),
        State('store-features',      'data'),
        State('store-spectra',       'data'),
        State('store-selected-scan-idx', 'data'),
        prevent_initial_call=True,
    )
    def do_export(n_clicks, fmt, results_data, ions_data,
                  feats_data, spectra_data, scan_idx):
        if not fmt:
            return no_update

        if fmt == 'csv_results':
            if not results_data:
                return no_update
            buf = io.StringIO()
            fields = ['rank', 'protein_name', 'sequence', 'start_pos', 'end_pos',
                      'score', 'matched_ions', 'total_ions', 'sequence_coverage',
                      'theoretical_mass', 'observed_mass', 'mass_error_da',
                      'mass_error_ppm', 'modifications']
            writer = csv.DictWriter(buf, fieldnames=fields)
            writer.writeheader()
            for rank, r in enumerate(results_data, 1):
                pf   = r['proteoform']
                mods = '; '.join(f"{m['name']}@{m['position']}"
                                 for m in pf.get('modifications', []))
                writer.writerow({
                    'rank': rank,
                    'protein_name': pf.get('protein_name', ''),
                    'sequence': pf.get('sequence', ''),
                    'start_pos': pf.get('start_pos', ''),
                    'end_pos': pf.get('end_pos', ''),
                    'score': pf.get('score', ''),
                    'matched_ions': pf.get('matched_ions', ''),
                    'total_ions': pf.get('total_ions', ''),
                    'sequence_coverage': r.get('sequence_coverage', ''),
                    'theoretical_mass': pf.get('theoretical_mass', ''),
                    'observed_mass': pf.get('observed_mass', ''),
                    'mass_error_da': pf.get('mass_error_da', ''),
                    'mass_error_ppm': pf.get('mass_error_ppm', ''),
                    'modifications': mods,
                })
            return dcc.send_string(buf.getvalue(), 'proform_results.csv')

        elif fmt == 'csv_ions':
            if not ions_data:
                return no_update
            buf = io.StringIO()
            fields = ['ion_type', 'position', 'charge', 'th_mz', 'obs_mz',
                      'mass_error_da', 'mass_error_ppm', 'matched', 'sequence']
            writer = csv.DictWriter(buf, fieldnames=fields)
            writer.writeheader()
            for ion in ions_data:
                writer.writerow({
                    'ion_type': ion['ion_type'],
                    'position': ion['position'],
                    'charge':   ion['charge'],
                    'th_mz':    ion['mz'],
                    'obs_mz':   ion.get('observed_mz', ''),
                    'mass_error_da':  ion.get('mass_error_da', ''),
                    'mass_error_ppm': ion.get('mass_error_ppm', ''),
                    'matched':  ion.get('matched', False),
                    'sequence': ion.get('sequence', ''),
                })
            return dcc.send_string(buf.getvalue(), 'proform_ions.csv')

        elif fmt == 'csv_features':
            if not feats_data:
                return no_update
            buf = io.StringIO()
            if feats_data:
                fields = list(feats_data[0].keys())
                writer = csv.DictWriter(buf, fieldnames=fields)
                writer.writeheader()
                writer.writerows(feats_data)
            return dcc.send_string(buf.getvalue(), 'proform_features.csv')

        elif fmt == 'json_session':
            session = {
                'search_results': results_data or [],
                'matched_ions':   ions_data or [],
                'features':       feats_data or [],
            }
            return dcc.send_string(json.dumps(session, indent=2), 'proform_session.json')

        return no_update
