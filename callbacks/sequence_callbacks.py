"""
Sequence tab callback — updates the proteoform coverage viewer.
"""
from dash import Input, Output, State

from src.data.models import Proteoform, FragmentIon
from src.viz.sequence_plots import create_sequence_plot
from src.analysis.mass_utils import calc_sequence_mass
from src.data.amino_acids import AA_MASSES


def register_callbacks(app):

    @app.callback(
        Output('sequence-graph', 'figure'),
        Output('coverage-display', 'children'),
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

        # Coverage text
        from src.analysis.peak_matching import coverage_map
        if pf.sequence:
            cov = coverage_map(pf.sequence, ions)
            n_cov = sum(1 for v in cov.values() if v)
            pct   = n_cov / len(pf.sequence) * 100
            cov_text = f"Sequence coverage: {n_cov}/{len(pf.sequence)} aa ({pct:.1f}%)"
        else:
            cov_text = ''

        return fig, cov_text
