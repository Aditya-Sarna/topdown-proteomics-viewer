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
        Input('seq-residues-per-row', 'value'),
        Input('seq-focus-mass-btn', 'n_clicks'),
        State('store-matched-ions', 'data'),
        State('seq-focus-mass', 'value'),
    )
    def update_sequence(selected_result, prot_data, show_cleavage, residues_per_row,
                        _focus_clicks, ions_data, focus_mass):
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

        rpr = int(residues_per_row) if residues_per_row else 20
        fig = create_sequence_plot(pf, ions, show_cleavage=bool(show_cleavage),
                                   residues_per_row=rpr)

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

        # Focus mass: highlight the matching subsequence on the grid
        if focus_mass and pf.sequence:
            try:
                fm = float(focus_mass)
                if fm > 0:
                    from src.analysis.mass_utils import calc_sequence_mass
                    best_pos, best_err = None, 1e9
                    seq = pf.sequence
                    for i in range(len(seq)):
                        for j in range(i + 1, len(seq) + 1):
                            sub_mass = calc_sequence_mass(seq[i:j])
                            err = abs(sub_mass - fm) / fm * 1e6
                            if err < best_err:
                                best_err, best_pos = err, (i, j)
                    if best_pos is not None and best_err < 50:  # 50 ppm threshold
                        i0, i1 = best_pos
                        # Draw a highlight rectangle around each residue in the match
                        for k in range(i0, i1):
                            col = k % rpr
                            row = k // rpr
                            fig.add_shape(
                                type='rect',
                                x0=col - 0.45, x1=col + 0.45,
                                y0=-row - 0.45, y1=-row + 0.45,
                                line=dict(color='#e65100', width=3),
                                fillcolor='rgba(230, 81, 0, 0.12)',
                            )
                        # Label above the first residue of the match
                        first_col = i0 % rpr
                        first_row = i0 // rpr
                        fig.add_annotation(
                            x=first_col,
                            y=-first_row + 0.7,
                            text=f"<b>{fm:.1f} Da</b> ({seq[i0:i1][:12]}{'…' if i1 - i0 > 12 else ''}) Δ{best_err:.1f} ppm",
                            showarrow=True, arrowhead=2, arrowsize=0.8,
                            arrowcolor='#e65100', ax=0, ay=-25,
                            font=dict(color='#e65100', size=10, family='Arial'),
                            xanchor='center',
                            bgcolor='rgba(255,255,255,0.85)',
                            bordercolor='#e65100', borderwidth=1, borderpad=3,
                        )
            except (ValueError, TypeError):
                pass

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
