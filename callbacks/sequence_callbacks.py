"""
Sequence tab callback — updates the proteoform coverage viewer.
"""
from dash import Input, Output, State, html

from src.data.models import Proteoform, FragmentIon
from src.viz.sequence_plots import create_sequence_plot
from src.analysis.mass_utils import calc_sequence_mass
from src.data.amino_acids import AA_MASSES


def register_callbacks(app):

    @app.callback(
        Output('sequence-graph', 'figure'),
        Output('coverage-display', 'children'),
        Output('sequence-mass-header', 'children'),
        Input('store-selected-result', 'data'),
        Input('store-protein', 'data'),
        Input('show-cleavage', 'value'),
        State('store-matched-ions', 'data'),
    )
    def update_sequence(selected_result, prot_data, show_cleavage, ions_data):
        ions = ([FragmentIon.from_dict(d) for d in ions_data]
                if ions_data else [])

        # Use selected proteoform if available
        if selected_result:
            pf = Proteoform.from_dict(selected_result)
        elif prot_data and prot_data.get('sequence'):
            seq = ''.join(c for c in prot_data['sequence'] if c in AA_MASSES)
            pf  = Proteoform(
                sequence=seq,
                protein_name=prot_data.get('name', 'Protein'),
            )
        else:
            pf = Proteoform(sequence='', protein_name='')

        fig = create_sequence_plot(pf, ions, show_cleavage=bool(show_cleavage))

        # Coverage text — use vectorized coverage_count_map (numpy) instead of
        # coverage_map (Python loops), avoiding a redundant O(N×M) pass.
        if pf.sequence:
            from src.analysis.peak_matching import coverage_count_map
            cov_cnt_cb = coverage_count_map(pf.sequence, ions)
            n_cov = sum(1 for nn, nc in cov_cnt_cb.values() if nn + nc > 0)
            pct   = n_cov / len(pf.sequence) * 100
            cov_text = f"Sequence coverage: {n_cov}/{len(pf.sequence)} aa ({pct:.1f}%)"
        else:
            cov_text = ''

        # Sequence-view mass header
        if pf.sequence and (pf.theoretical_mass or pf.observed_mass):
            delta_da  = pf.observed_mass - pf.theoretical_mass if pf.theoretical_mass else 0.0
            name_part = pf.protein_name or 'Proteoform'
            mass_header = html.Span([
                html.Strong(f"{name_part}"),
                html.Span("  |  ", style={'color': '#aaaaaa'}),
                html.Span(f"Theoretical mass : {pf.theoretical_mass:.2f} Da",
                          style={'marginRight': '16px'}),
                html.Span(f"Observed mass : {pf.observed_mass:.2f} Da",
                          style={'marginRight': '16px'}),
                html.Span(
                    f"\u0394 Mass (Da) : {delta_da:+.2f}",
                    style={'color': '#EF5350' if abs(pf.mass_error_ppm) > 5 else '#2e7d32'},
                ),
            ])
        else:
            mass_header = ''

        return fig, cov_text, mass_header
