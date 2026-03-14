from callbacks.data_callbacks     import register_callbacks as _data
from callbacks.spectrum_callbacks import register_callbacks as _spectrum
from callbacks.sequence_callbacks import register_callbacks as _sequence
from callbacks.mirror_callbacks   import register_callbacks as _mirror
from callbacks.feature_callbacks  import register_callbacks as _feature
from callbacks.search_callbacks   import register_callbacks as _search
from callbacks.export_callbacks   import register_callbacks as _export
from callbacks.heatmap_callbacks  import register_callbacks as _heatmap


def register_all_callbacks(app):
    _data(app)
    _spectrum(app)
    _sequence(app)
    _mirror(app)
    _feature(app)
    _search(app)
    _export(app)
    _heatmap(app)
