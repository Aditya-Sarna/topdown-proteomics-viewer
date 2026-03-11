"""
Feature map and intensity trace callbacks.
"""
from dash import Input, Output, State

from src.data.models import Feature, Proteoform
from src.viz.feature_plots import create_feature_map, create_intensity_trace


def register_callbacks(app):

    # Feature map (main 2-D scatter)
    @app.callback(
        Output('feature-map-graph', 'figure'),
        Input('store-features', 'data'),
        Input('feature-color-by', 'value'),
        Input('feature-filter', 'value'),
        Input('th-mass-input', 'value'),
        Input('feature-map-graph', 'clickData'),
        State('store-search-results', 'data'),
    )
    def update_feature_map(feats_data, color_by, fid_filter, th_mass, click_data, results_data):
        if not feats_data:
            return create_feature_map([])

        features = [Feature.from_dict(d) for d in feats_data]
        selected_id = None

        if click_data:
            # Extract feature_id from customdata[0]
            try:
                selected_id = click_data['points'][0]['customdata'][0]
            except (KeyError, IndexError, TypeError):
                pass

        if fid_filter:
            selected_id = fid_filter.strip()

        # Build Proteoform overlays from search results
        proteoforms = []
        if th_mass:
            from src.data.models import Proteoform
            pf = Proteoform(sequence='', protein_name='Target',
                            theoretical_mass=float(th_mass))
            proteoforms = [pf]
        elif results_data:
            from src.data.models import SearchResult
            for r in results_data[:3]:
                sr = SearchResult.from_dict(r)
                proteoforms.append(sr.proteoform)

        return create_feature_map(features, proteoforms or None, selected_id, color_by)

    # Intensity trace (elution profile)
    @app.callback(
        Output('intensity-trace-graph', 'figure'),
        Input('feature-map-graph', 'clickData'),
        Input('feature-filter', 'value'),
        Input('th-mass-input', 'value'),
        Input('th-charge-input', 'value'),
        State('store-features', 'data'),
    )
    def update_trace(click_data, fid_filter, th_mass, th_charge, feats_data):
        if not feats_data:
            return create_intensity_trace([])

        features = [Feature.from_dict(d) for d in feats_data]
        selected_id = None

        if click_data:
            try:
                selected_id = click_data['points'][0]['customdata'][0]
            except (KeyError, IndexError, TypeError):
                pass
        if fid_filter:
            selected_id = fid_filter.strip()

        return create_intensity_trace(
            features,
            selected_id=selected_id,
            theoretical_mass=float(th_mass) if th_mass else 0.0,
            theoretical_charge=int(th_charge) if th_charge else 0,
        )
