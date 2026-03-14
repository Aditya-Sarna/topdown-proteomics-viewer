"""
Heatmap callbacks — raw MS heatmap and deconvolved precursor-mass heatmap.
"""
from dash import Input, Output

from src.data.models import Spectrum
from src.viz.heatmap_plots import create_raw_heatmap, create_deconvolved_heatmap


def register_callbacks(app):

    @app.callback(
        Output('raw-heatmap-graph',    'figure'),
        Output('deconv-heatmap-graph', 'figure'),
        Input('store-spectra',      'data'),
        Input('heatmap-ms-level',   'value'),
    )
    def update_heatmaps(spectra_data, ms_level):
        spectra = [Spectrum.from_dict(d) for d in (spectra_data or [])]
        raw    = create_raw_heatmap(spectra, ms_level=int(ms_level or 2))
        deconv = create_deconvolved_heatmap(spectra)
        return raw, deconv
