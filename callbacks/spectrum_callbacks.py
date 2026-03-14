"""
Spectrum tab callback — updates the spectrum graph when the scan or
matched-ions store changes.
"""
from dash import Input, Output, State, no_update

from src.data.models import Spectrum, FragmentIon
from src.viz.spectrum_plots import create_spectrum_plot, create_deconvolved_spectrum_plot


def register_callbacks(app):

    @app.callback(
        Output('spectrum-graph', 'figure'),
        Input('store-selected-scan-idx', 'data'),
        Input('store-matched-ions', 'data'),
        State('store-spectra', 'data'),
    )
    def update_spectrum(idx, ions_data, spectra_data):
        if not spectra_data:
            from src.viz.spectrum_plots import go, DARK_BG, PLOT_BG
            fig = go.Figure()
            fig.update_layout(
                template='plotly_dark',
                title='Load a spectrum or click "Load Demo" to begin',
                paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
                font=dict(color='#e0e0e0'),
            )
            return fig

        idx = idx or 0
        spectrum = Spectrum.from_dict(spectra_data[idx])
        ions = ([FragmentIon.from_dict(d) for d in ions_data]
                if ions_data else None)
        return create_spectrum_plot(spectrum, ions)

    @app.callback(
        Output('deconv-spectrum-graph', 'figure'),
        Input('store-selected-scan-idx', 'data'),
        State('store-spectra', 'data'),
    )
    def update_deconv_spectrum(idx, spectra_data):
        if not spectra_data:
            import plotly.graph_objects as go
            return go.Figure().update_layout(
                title='Load a spectrum to see deconvolved view',
                template='plotly_white', paper_bgcolor='#ffffff',
                font=dict(color='#aaaaaa'), height=360,
            )
        idx      = idx or 0
        spectrum = Spectrum.from_dict(spectra_data[idx])
        return create_deconvolved_spectrum_plot(spectrum)
