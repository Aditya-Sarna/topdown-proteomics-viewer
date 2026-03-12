from dash import Output, Input
from callbacks.data_callbacks     import register_callbacks as _data
from callbacks.spectrum_callbacks import register_callbacks as _spectrum
from callbacks.sequence_callbacks import register_callbacks as _sequence
from callbacks.mirror_callbacks   import register_callbacks as _mirror
from callbacks.feature_callbacks  import register_callbacks as _feature
from callbacks.search_callbacks   import register_callbacks as _search
from callbacks.export_callbacks   import register_callbacks as _export


def register_all_callbacks(app):
    app.clientside_callback(
        """
        function(n) {
            if (n && n >= 1) {
                var el = document.getElementById('splash-screen');
                if (el) { el.classList.add('splash-hidden'); }
            }
            return '';
        }
        """,
        Output('splash-screen', 'className'),
        Input('splash-interval', 'n_intervals'),
        prevent_initial_call=True,
    )

    _data(app)
    _spectrum(app)
    _sequence(app)
    _mirror(app)
    _feature(app)
    _search(app)
    _export(app)
