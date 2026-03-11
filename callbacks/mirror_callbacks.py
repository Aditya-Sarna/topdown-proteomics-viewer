"""
Mirror plot callback.
"""
from dash import Input, Output, State

from src.data.models import Spectrum, FragmentIon
from src.viz.mirror_plots import create_mirror_plot


def register_callbacks(app):

    @app.callback(
        Output('mirror-graph', 'figure'),
        Input('store-selected-scan-idx', 'data'),
        Input('store-matched-ions', 'data'),
        Input('mirror-mz-min', 'value'),
        Input('mirror-mz-max', 'value'),
        State('store-spectra', 'data'),
    )
    def update_mirror(idx, ions_data, mz_min, mz_max, spectra_data):
        if not spectra_data:
            from src.viz.mirror_plots import go, DARK_BG, PLOT_BG
            fig = go.Figure()
            fig.update_layout(
                template='plotly_dark',
                title='Load a spectrum first',
                paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
                font=dict(color='#e0e0e0'),
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
        return create_mirror_plot(spectrum, ions, mz_range=mz_range)
