"""
Main Dash layout.
Structured as: NavBar | Sidebar (controls) | Main content area (5 tabs).
"""
import os
from dash import dcc, html, dash_table
import dash_bootstrap_components as dbc
from src.data.amino_acids import PTM_DATABASE

_DEMO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo_data')
_DEMO_FILES = sorted(f for f in os.listdir(_DEMO_DIR) if not f.startswith('.'))

# ──────────────────────────────────────────────────────────────────────────────
# Colour constants
# ──────────────────────────────────────────────────────────────────────────────
SIDEBAR_BG = '#f8f8f8'
CARD_BG    = '#ffffff'
ACCENT     = '#1a73e8'
TEXT_MUTED = '#555555'


# ──────────────────────────────────────────────────────────────────────────────
# Reusable components
# ──────────────────────────────────────────────────────────────────────────────

def _section(title, *children):
    return html.Div([
        html.H6(title, className='text-uppercase fw-bold mb-2',
                style={'color': TEXT_MUTED, 'fontSize': '0.7rem', 'letterSpacing': '0.08em'}),
        *children,
    ], className='mb-3')


def _label(text):
    return html.Small(text, className='text-muted d-block mb-1')


def _card(*children, **kwargs):
    return dbc.Card(dbc.CardBody(list(children), className='p-2'),
                    className='mb-2',
                    style={'background': CARD_BG, 'border': '1px solid #dddddd', **kwargs})


# ──────────────────────────────────────────────────────────────────────────────
# Navbar
# ──────────────────────────────────────────────────────────────────────────────

def _navbar():
    return dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col(html.Span([
                    html.Span("Top-Down Proteomics Viewer",
                              className='fw-bold',
                              style={'fontSize': '1.05rem', 'color': '#111111',
                                     'letterSpacing': '0.02em'}),
                    dbc.Badge("v1.0", color='primary', className='ms-2 align-middle'),
                ]), width='auto'),
                dbc.Col(html.Small("Single-spectrum proteoform analysis",
                                   style={'color': TEXT_MUTED}), width='auto'),
            ], align='center', className='g-2 flex-grow-1'),
            dbc.Row([
                dbc.Col(html.Div([
                    dcc.Dropdown(
                        id='demo-file-select',
                        options=[{'label': f, 'value': f} for f in _DEMO_FILES],
                        placeholder='Load Datasets…',
                        clearable=False,
                        style={'minWidth': '220px', 'fontSize': '0.8rem', 'display': 'inline-block'},
                    ),
                ], style={'display': 'inline-block', 'verticalAlign': 'middle'}), width='auto'),
                dbc.Col(dcc.Loading(
                    id='load-demo-loading',
                    type='circle',
                    color='#4CAF50',
                    children=dbc.Button('Load', id='load-demo-btn', color='success',
                                       size='sm', outline=True, className='ms-1 me-2'),
                ), width='auto'),
                dbc.Col(dbc.Button("Export", id='export-btn', color='secondary',
                                   size='sm', outline=True), width='auto'),
            ], align='center', className='g-1'),
        ], fluid=True),
        color='white', dark=False,
        style={'borderBottom': '1px solid #dddddd', 'boxShadow': '0 1px 3px rgba(0,0,0,0.08)'},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Sidebar
# ──────────────────────────────────────────────────────────────────────────────

def _sidebar():
    ptm_options = [{'label': name, 'value': name}
                   for name in sorted(PTM_DATABASE.keys())]

    return html.Div([

        # ── Data Input ─────────────────────────────────────────────────────
        _section("Data Input",
            _label("Upload Spectrum (.mzML / .pcml / peak CSV)"),
            dcc.Upload(
                id='upload-spectrum',
                children=html.Div(['Drag & drop or ', html.A('browse', style={'color': ACCENT}),
                                   html.Span(' (.mzML, .pcml, .csv)', style={'color': TEXT_MUTED, 'fontSize': '0.70rem'})]),
                style={
                    'borderWidth': '1px', 'borderStyle': 'dashed',
                    'borderRadius': '4px', 'borderColor': '#cccccc',
                    'textAlign': 'center', 'padding': '8px 4px',
                    'color': TEXT_MUTED, 'fontSize': '0.78rem', 'cursor': 'pointer',
                    'background': '#fafafa',
                },
                multiple=False,
            ),
            dcc.Loading(
                id='upload-spectrum-loading',
                type='circle',
                color='#2196F3',
                children=html.Div(id='upload-spectrum-status',
                                  className='mt-1',
                                  style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
            ),

            html.Div(className='mt-2'),
            _label("Upload Feature Table (.csv)"),
            dcc.Upload(
                id='upload-features',
                children=html.Div(['Drag & drop or ', html.A('browse', style={'color': ACCENT})]),
                style={
                    'borderWidth': '1px', 'borderStyle': 'dashed',
                    'borderRadius': '4px', 'borderColor': '#cccccc',
                    'textAlign': 'center', 'padding': '8px 4px',
                    'color': TEXT_MUTED, 'fontSize': '0.78rem', 'cursor': 'pointer',
                    'background': '#fafafa',
                },
                multiple=False,
            ),
            dcc.Loading(
                id='upload-features-loading',
                type='circle',
                color='#2196F3',
                children=html.Div(id='upload-features-status',
                                  className='mt-1',
                                  style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
            ),
        ),

        # ── Scan Selector ──────────────────────────────────────────────────
        _section("Scan",
            _label("Select Scan"),
            dcc.Dropdown(
                id='scan-selector',
                placeholder='Load a spectrum file first…',
                style={'fontSize': '0.8rem'},
                className='dark-dropdown',
            ),
            html.Div(id='scan-info',
                     className='mt-1',
                     style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
        ),

        # ── Protein ────────────────────────────────────────────────────────
        _section("Protein",
            _label("Protein Name"),
            dbc.Input(id='protein-name', placeholder='e.g. Ubiquitin',
                      size='sm', className='mb-1',
                      style={'background': '#ffffff', 'color': '#111111',
                             'border': '1px solid #cccccc'}),
            _label("Sequence (single-letter)"),
            dbc.Textarea(
                id='protein-sequence',
                placeholder='MQIFVKTLTGK…',
                rows=4,
                style={'background': '#ffffff', 'color': '#111111',
                       'fontSize': '0.78rem', 'fontFamily': 'Courier New',
                       'border': '1px solid #cccccc', 'resize': 'vertical'},
            ),
            html.Div(id='protein-mass-display',
                     className='mt-1',
                     style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
        ),

        # ── Search Parameters ──────────────────────────────────────────────
        _section("Search Parameters",
            dbc.Row([
                dbc.Col([
                    _label("Tolerance (ppm)"),
                    dbc.Input(id='tolerance-ppm', type='number', value=10, min=1, max=100,
                              size='sm',
                              style={'background': '#ffffff', 'color': '#111111',
                                     'border': '1px solid #cccccc'}),
                ], width=6),
                dbc.Col([
                    _label("Max z"),
                    dbc.Input(id='max-charge', type='number', value=4, min=1, max=20,
                              size='sm',
                              style={'background': '#ffffff', 'color': '#111111',
                                     'border': '1px solid #cccccc'}),
                ], width=6),
            ], className='g-1 mb-1'),

            _label("Ion Types"),
            dbc.Checklist(
                id='ion-types',
                options=[{'label': t.upper(), 'value': t}
                         for t in ['b', 'y', 'c', 'z', 'a']],
                value=['b', 'y', 'c', 'z'],
                inline=True,
                className='small',
                inputStyle={'marginRight': '3px'},
            ),
            html.Div(className='mt-1'),
            dbc.Row([
                dbc.Col([
                    _label("Truncations"),
                    dbc.Switch(id='search-truncations', value=True, label=''),
                ], width=6),
                dbc.Col([
                    _label("Modifications"),
                    dbc.Switch(id='search-mods', value=True, label=''),
                ], width=6),
            ], className='g-1 mb-1'),

            _label("Variable PTMs"),
            dcc.Dropdown(
                id='variable-mods',
                options=ptm_options,
                value=['Phosphorylation', 'Acetylation', 'Oxidation'],
                multi=True,
                placeholder='Select PTMs…',
                style={'fontSize': '0.75rem'},
                className='dark-dropdown',
            ),
        ),

        # ── Manual Mass Input ──────────────────────────────────────────────
        _section("Manual Mass Input",
            _label("Observed Neutral Mass (Da)"),
            dbc.InputGroup([
                dbc.Input(id='manual-mass', type='number', step=0.001,
                          placeholder='e.g. 8564.83',
                          size='sm',
                          style={'background': '#ffffff', 'color': '#111111',
                                 'border': '1px solid #cccccc'}),
                dbc.Button("Δ", id='calc-mass-diff-btn', color='secondary',
                           size='sm', outline=True),
            ], size='sm'),
            html.Div(id='mass-diff-display',
                     className='mt-1',
                     style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
        ),

        # ── Run Button ─────────────────────────────────────────────────────
        dcc.Loading(
            id='search-loading',
            type='circle',
            color='#2196F3',
            children=[
                dbc.Button('Run Search', id='run-search-btn',
                           color='primary', className='w-100 mb-2 fw-bold'),
                html.Div(id='search-status',
                         style={'fontSize': '0.78rem', 'color': TEXT_MUTED,
                                'textAlign': 'center'}),
            ],
        ),

        # ── Results Summary ─────────────────────────────────────────────────
        html.Hr(style={'borderColor': '#dddddd'}),
        _section("Top Hit Summary",
            html.Div(id='top-hit-summary',
                     style={'fontSize': '0.78rem', 'color': '#e0e0e0',
                            'lineHeight': '1.6'}),
        ),
        _section("Modification Suggestions",
            html.Div(id='mod-suggestions',
                     style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),
        ),

    ], style={
        'width': '270px', 'minWidth': '270px',
        'background': SIDEBAR_BG,
        'padding': '12px 12px',
        'overflowY': 'auto', 'height': 'calc(100vh - 56px)',
        'borderRight': '1px solid #dddddd',
        'flexShrink': '0',
    })


# ──────────────────────────────────────────────────────────────────────────────
# Tab contents
# ──────────────────────────────────────────────────────────────────────────────

def _tab_spectrum():
    return dbc.Tab(label='Spectrum', tab_id='tab-spectrum',
                   children=dbc.Card(dbc.CardBody([
                       dcc.Graph(id='spectrum-graph', config={'displayModeBar': True},
                                 style={'height': '520px'}),
                       html.Hr(style={'borderColor': '#dddddd'}),
                       html.Div(id='peak-click-info',
                                style={'fontSize': '0.78rem', 'color': TEXT_MUTED}),
                   ], className='p-2'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_sequence():
    return dbc.Tab(label='Sequence', tab_id='tab-sequence',
                   children=dbc.Card(dbc.CardBody([
                       dbc.Row([
                           dbc.Col(dbc.Switch(id='show-cleavage', value=True,
                                              label='Show cleavage ticks'), width='auto'),
                           dbc.Col(html.Div(id='coverage-display',
                                            style={'fontSize': '0.82rem', 'color': '#333333'}),
                                   width='auto'),
                       ], className='mb-2 align-items-center'),
                       dcc.Graph(id='sequence-graph', config={'displayModeBar': False},
                                 style={'minHeight': '300px'}),
                   ], className='p-2'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_mirror():
    return dbc.Tab(label='Mirror Plot', tab_id='tab-mirror',
                   children=dbc.Card(dbc.CardBody([
                       dbc.Row([
                           dbc.Col([
                               _label("m/z range — Min"),
                               dbc.Input(id='mirror-mz-min', type='number', value=200,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#111111',
                                                'border': '1px solid #cccccc'}),
                           ], width=2),
                           dbc.Col([
                               _label("m/z range — Max"),
                               dbc.Input(id='mirror-mz-max', type='number', value=2000,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#111111',
                                                'border': '1px solid #cccccc'}),
                           ], width=2),
                       ], className='mb-2'),
                       dcc.Graph(id='mirror-graph', config={'displayModeBar': True},
                                 style={'height': '560px'}),
                   ], className='p-2'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_features():
    return dbc.Tab(label='Feature Map', tab_id='tab-features',
                   children=dbc.Card(dbc.CardBody([
                       dbc.Row([
                           dbc.Col([
                               _label("Color by"),
                               dbc.RadioItems(
                                   id='feature-color-by',
                                   options=[{'label': 'Charge', 'value': 'charge'},
                                            {'label': 'Proteoform', 'value': 'proteoform'}],
                                   value='charge',
                                   inline=True, className='small',
                               ),
                           ], width='auto'),
                           dbc.Col([
                               _label("Feature ID filter"),
                               dbc.Input(id='feature-filter', placeholder='e.g. F001',
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#111111',
                                                'border': '1px solid #cccccc'}),
                           ], width=3),
                       ], className='mb-2 align-items-end'),
                       dcc.Graph(id='feature-map-graph', config={'displayModeBar': True},
                                 style={'height': '480px'}),
                       html.Hr(style={'borderColor': '#dddddd'}),
                       dbc.Row([
                           dbc.Col([
                               _label("Theoretical mass for overlay (Da)"),
                               dbc.Input(id='th-mass-input', type='number', step=0.001,
                                         placeholder='e.g. 8564.83',
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#111111',
                                                'border': '1px solid #cccccc'}),
                           ], width=3),
                           dbc.Col([
                               _label("Charge for trace"),
                               dbc.Input(id='th-charge-input', type='number', value=10, min=1,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#111111',
                                                'border': '1px solid #cccccc'}),
                           ], width=2),
                       ], className='mb-2'),
                       dcc.Graph(id='intensity-trace-graph', config={'displayModeBar': False},
                                 style={'height': '300px'}),
                   ], className='p-2'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_search():
    columns = [
        {'name': 'Rank',        'id': 'rank'},
        {'name': 'Sequence',    'id': 'sequence'},
        {'name': 'Range',       'id': 'range'},
        {'name': 'Score',       'id': 'score'},
        {'name': 'E-value',     'id': 'e_value'},
        {'name': 'q-value',     'id': 'q_value'},
        {'name': 'Matched',     'id': 'matched'},
        {'name': 'Coverage %',  'id': 'coverage'},
        {'name': 'Th. Mass',    'id': 'th_mass'},
        {'name': 'Obs. Mass',   'id': 'obs_mass'},
        {'name': 'Δmass (ppm)', 'id': 'ppm_error'},
        {'name': 'PTMs',        'id': 'ptms'},
    ]
    return dbc.Tab(label='Search Results', tab_id='tab-search',
                   children=dbc.Card(dbc.CardBody([
                       html.Div(id='search-result-summary',
                                className='mb-2',
                                style={'fontSize': '0.8rem', 'color': TEXT_MUTED}),
                       dash_table.DataTable(
                           id='results-table',
                           columns=columns,
                           data=[],
                           row_selectable='single',
                           selected_rows=[],
                           page_size=15,
                           style_table={'overflowX': 'auto'},
                           style_cell={
                               'backgroundColor': '#ffffff',
                               'color': '#111111',
                               'border': '1px solid #dddddd',
                               'fontSize': '0.78rem',
                               'padding': '5px 8px',
                               'whiteSpace': 'normal',
                               'maxWidth': '220px',
                               'overflow': 'hidden',
                               'textOverflow': 'ellipsis',
                           },
                           style_header={
                               'backgroundColor': '#f0f0f0',
                               'color': '#111111',
                               'fontWeight': 'bold',
                               'border': '1px solid #dddddd',
                               'fontSize': '0.75rem',
                           },
                           style_data_conditional=[
                               {'if': {'row_index': 0},
                                'backgroundColor': '#e8f0fe',
                                'color': '#1a73e8'},
                               {'if': {'state': 'selected'},
                                'backgroundColor': '#d2e3fc',
                                'border': '1px solid #1a73e8'},
                               # Highlight q-value column: green < 0.01, yellow < 0.05
                               {'if': {'filter_query': '{q_value} < 0.01', 'column_id': 'q_value'},
                                'color': '#2e7d32', 'fontWeight': 'bold'},
                               {'if': {'filter_query': '{q_value} >= 0.01 && {q_value} < 0.05', 'column_id': 'q_value'},
                                'color': '#f57f17'},
                           ],
                           tooltip_data=[],
                           tooltip_duration=None,
                           sort_action='native',
                           filter_action='native',
                       ),
                       html.Hr(style={'borderColor': '#dddddd'}),
                       html.H6("Selected Proteoform — Ion Table",
                               className='text-muted small mb-1'),
                       dash_table.DataTable(
                           id='ion-table',
                           columns=[
                               {'name': 'Ion',       'id': 'ion'},
                               {'name': 'Th. m/z',   'id': 'th_mz'},
                               {'name': 'Obs. m/z',  'id': 'obs_mz'},
                               {'name': 'Δ (ppm)',   'id': 'ppm'},
                               {'name': 'Charge',    'id': 'charge'},
                               {'name': 'Matched',   'id': 'matched'},
                               {'name': 'Sequence',  'id': 'seq'},
                           ],
                           data=[],
                           page_size=12,
                           style_table={'overflowX': 'auto'},
                           style_cell={
                               'backgroundColor': '#ffffff',
                               'color': '#111111',
                               'border': '1px solid #dddddd',
                               'fontSize': '0.75rem',
                               'padding': '4px 8px',
                           },
                           style_header={
                               'backgroundColor': '#f0f0f0',
                               'color': '#111111',
                               'fontWeight': 'bold',
                               'border': '1px solid #dddddd',
                               'fontSize': '0.73rem',
                           },
                           style_data_conditional=[
                               {'if': {'filter_query': '{matched} eq "✓"'},
                                'color': '#2e7d32'},
                           ],
                           sort_action='native',
                           filter_action='native',
                       ),
                   ], className='p-2'), style={'background': CARD_BG, 'border': 'none'}))


# ──────────────────────────────────────────────────────────────────────────────
# Export Modal
# ──────────────────────────────────────────────────────────────────────────────

def _export_modal():
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("Export Results")),
        dbc.ModalBody([
            _label("Export Format"),
            dbc.RadioItems(
                id='export-format',
                options=[
                    {'label': 'CSV — Search Results',   'value': 'csv_results'},
                    {'label': 'CSV — Ion Matches',       'value': 'csv_ions'},
                    {'label': 'CSV — Feature Table',     'value': 'csv_features'},
                    {'label': 'JSON — Full Session',     'value': 'json_session'},
                ],
                value='csv_results',
                className='small',
            ),
        ]),
        dbc.ModalFooter([
            dbc.Button("Download", id='export-confirm-btn', color='primary', size='sm'),
            dbc.Button("Close",    id='export-close-btn',  color='secondary',
                       size='sm', className='ms-2', outline=True),
        ]),
    ], id='export-modal', is_open=False)


# ──────────────────────────────────────────────────────────────────────────────
# Stores (client-side session state)
# ──────────────────────────────────────────────────────────────────────────────

def _stores():
    return html.Div([
        dcc.Store(id='store-spectra',            storage_type='memory', data=[]),
        dcc.Store(id='store-selected-scan-idx',  storage_type='memory', data=0),
        dcc.Store(id='store-features',           storage_type='memory', data=[]),
        dcc.Store(id='store-protein',            storage_type='memory', data={}),
        dcc.Store(id='store-search-results',     storage_type='memory', data=[]),
        dcc.Store(id='store-matched-ions',       storage_type='memory', data=[]),
        dcc.Store(id='store-selected-result',    storage_type='memory', data={}),
        dcc.Download(id='download-data'),
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Root layout
# ──────────────────────────────────────────────────────────────────────────────

def create_layout():
    tabs = dbc.Tabs([
        _tab_spectrum(),
        _tab_sequence(),
        _tab_mirror(),
        _tab_features(),
        _tab_search(),
    ], id='main-tabs', active_tab='tab-spectrum',
    className='mb-0',
    style={'background': '#ffffff'})

    main_content = html.Div([
        tabs,
    ], style={
        'flex': '1',
        'padding': '8px 12px',
        'overflowY': 'auto',
        'height': 'calc(100vh - 56px)',
    })

    body = html.Div([
        _sidebar(),
        main_content,
    ], style={
        'display': 'flex',
        'flexDirection': 'row',
        'height': 'calc(100vh - 56px)',
    })

    return html.Div([
        _stores(),
        _navbar(),
        body,
        _export_modal(),
        # Status toast
        dbc.Toast(
            id='status-toast',
            header="ProForm Viewer",
            is_open=False,
            dismissable=True,
            duration=4000,
            style={'position': 'fixed', 'bottom': '20px', 'right': '20px',
                   'background': CARD_BG, 'color': 'white', 'zIndex': 9999},
        ),
    ], style={'background': '#ffffff', 'minHeight': '100vh'})
