"""
ProForm Viewer — Top-Down Proteomics Analysis Tool
Run with:  python app.py
Then open: http://127.0.0.1:8050
"""
import dash
import dash_bootstrap_components as dbc

from layouts.main_layout import create_layout
from callbacks import register_all_callbacks

app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.FLATLY,
        'https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap',
    ],
    suppress_callback_exceptions=True,
    title='Top-Down Proteomics Viewer',
    update_title=None,
)

app.layout = create_layout()
register_all_callbacks(app)

server = app.server

if __name__ == '__main__':
    app.run(debug=True, host='127.0.0.1', port=8050)
