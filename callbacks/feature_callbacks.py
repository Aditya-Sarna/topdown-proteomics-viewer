"""
Feature map and intensity trace callbacks.
"""
from dash import Input, Output, State

from src.data.models import Feature, Proteoform
from src.viz.feature_plots import create_feature_map, create_intensity_trace, create_feature_3d_plot, create_xic_plot


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

    # 3-D feature scatter
    @app.callback(
        Output('feature-3d-graph', 'figure'),
        Input('store-features', 'data'),
    )
    def update_feature_3d(feats_data):
        features = [Feature.from_dict(d) for d in (feats_data or [])]
        return create_feature_3d_plot(features)

    # XIC — Extracted Ion Chromatogram from real MS1 spectra
    @app.callback(
        Output('xic-graph', 'figure'),
        Output('xic-label', 'children'),
        Input('feature-map-graph', 'clickData'),
        Input('feature-filter', 'value'),
        State('store-features', 'data'),
        State('store-spectra',  'data'),
    )
    def update_xic(click_data, fid_filter, feats_data, spectra_data):
        selected_id = None
        if click_data:
            try:
                selected_id = click_data['points'][0]['customdata'][0]
            except (KeyError, IndexError, TypeError):
                pass
        if fid_filter:
            selected_id = fid_filter.strip()

        feature = None
        if selected_id and feats_data:
            for d in feats_data:
                if d.get('feature_id') == selected_id:
                    feature = Feature.from_dict(d)
                    break

        # If nothing selected, default to the highest-intensity feature
        if feature is None and feats_data:
            feature = Feature.from_dict(
                max(feats_data, key=lambda d: d.get('intensity', 0))
            )

        spectra_list = spectra_data or []
        fig = create_xic_plot(spectra_list, feature)

        label = (
            f'Extracted Ion Chromatogram — {feature.feature_id} '
            f'(m/z {feature.mz_apex:.3f}, z={feature.charge})'
            if feature else 'Extracted Ion Chromatogram (XIC)'
        )
        return fig, label
