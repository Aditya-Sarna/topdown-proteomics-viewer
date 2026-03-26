"""
Mirror plot callback.
"""
from dash import Input, Output, State, no_update

from src.data.models import Spectrum, FragmentIon
from src.viz.mirror_plots import create_mirror_plot


def register_callbacks(app):

    @app.callback(
        Output('mirror-graph', 'figure'),
        Input('store-selected-scan-idx', 'data'),
        Input('store-matched-ions', 'data'),
        Input('mirror-mz-min', 'value'),
        Input('mirror-mz-max', 'value'),
        Input('mirror-annotate-btn', 'n_clicks'),
        State('store-spectra', 'data'),
        State('mirror-user-mass', 'value'),
        State('mirror-user-charge-max', 'value'),
    )
    def update_mirror(idx, ions_data, mz_min, mz_max, _annotate_clicks,
                      spectra_data, user_mass, user_charge_max):
        if not spectra_data:
            from src.viz.mirror_plots import go, DARK_BG, PLOT_BG
            fig = go.Figure()
            fig.update_layout(
                template='plotly_white',
                title='Load a spectrum first',
                paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
                font=dict(color='#555555'),
            )
            return fig

        idx      = idx or 0
        spectrum = Spectrum.from_dict(spectra_data[idx])
        ions     = ([FragmentIon.from_dict(d) for d in ions_data]
                    if ions_data else None)

        mz_range = (
            float(mz_min) if mz_min else 200,
            float(mz_max) if mz_max else 2000,
        )

        # Build user-defined neutral mass annotation (charge series)
        user_mass_annotation = None
        if user_mass and float(user_mass) > 0:
            from src.analysis.mass_utils import mass_to_mz, possible_charges
            mass_val   = float(user_mass)
            max_z      = int(user_charge_max or 15)
            lo, hi     = mz_range
            charges    = possible_charges(mass_val, min_mz=lo, max_mz=hi)
            charges    = [z for z in charges if z <= max_z]
            user_mass_annotation = {
                'mass': mass_val,
                'charge_mzs': [(z, mass_to_mz(mass_val, z)) for z in charges],
            }

        return create_mirror_plot(spectrum, ions, mz_range=mz_range,
                                   user_mass_annotation=user_mass_annotation)

    @app.callback(
        Output('mirror-user-mass-display', 'children'),
        Input('mirror-annotate-btn', 'n_clicks'),
        State('mirror-user-mass', 'value'),
        State('mirror-user-charge-max', 'value'),
        State('mirror-mz-min', 'value'),
        State('mirror-mz-max', 'value'),
        prevent_initial_call=True,
    )
    def update_mirror_mass_display(_n, user_mass, max_z, mz_min, mz_max):
        if not user_mass:
            return ''
        from src.analysis.mass_utils import mass_to_mz, possible_charges
        mass_val = float(user_mass)
        lo  = float(mz_min or 200)
        hi  = float(mz_max or 2000)
        charges = possible_charges(mass_val, min_mz=lo, max_mz=hi)
        charges = [z for z in charges if z <= int(max_z or 15)]
        if not charges:
            return f'Mass {mass_val:.2f} Da — no charge states in [{lo:.0f}, {hi:.0f}] m/z'
        parts = [f'z{z}: {mass_to_mz(mass_val, z):.3f}' for z in sorted(charges)]
        return (f'Mass {mass_val:.4f} Da — {len(charges)} charge state(s) in range: '
                + '  |  '.join(parts[:6]) + ('  …' if len(parts) > 6 else ''))
