"""
Main Dash layout.
Structured as: NavBar | Sidebar (controls) | Main content area (5 tabs).
"""
import os
from dash import dcc, html, dash_table
import dash_bootstrap_components as dbc
from src.data.amino_acids import PTM_DATABASE

_DEMO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'demo_data')
# Clean protein name → filename mapping for the navbar dropdown
_DEMO_OPTIONS = [
    {'label': 'Hemoglobin Beta',   'value': 'hemoglobin_beta.mzML'},
    {'label': 'Hemoglobin Alpha',  'value': 'hemoglobin_alpha.mzML'},
    {'label': 'Insulin B Chain',   'value': 'insulin_b_chain.mzML'},
    {'label': 'Serum Albumin N49', 'value': 'serum_albumin_n49.mzML'},
    {'label': 'Ubiquitin',         'value': 'ubiquitin.mzML'},
    {'label': 'Cytochrome C',      'value': 'cytochrome_c.mzML'},
    {'label': 'Thioredoxin',       'value': 'thioredoxin.mzML'},
]

# ──────────────────────────────────────────────────────────────────────────────
# Colour constants
# ──────────────────────────────────────────────────────────────────────────────
SIDEBAR_BG = '#f8f9fa'
CARD_BG    = '#ffffff'
ACCENT     = '#1a73e8'
TEXT_MUTED = '#5f6368'
BORDER     = '#dadce0'

# Graph config shared across all graphs
_GRAPH_CONFIG = {'displayModeBar': True, 'displaylogo': False,
                 'modeBarButtonsToRemove': ['toImage'],
                 'toImageButtonOptions': {'format': 'svg', 'scale': 2}}


# ──────────────────────────────────────────────────────────────────────────────
# Reusable components
# ──────────────────────────────────────────────────────────────────────────────

def _section(title, *children):
    return html.Div([
        html.Div(title, className='sidebar-section-label mb-2'),
        *children,
    ], className='mb-3')


def _label(text):
    return html.Small(text, className='text-muted d-block mb-1',
                      style={'fontSize': '0.75rem'})


def _card(*children, **kwargs):
    return dbc.Card(dbc.CardBody(list(children), className='p-2'),
                    className='mb-2',
                    style={'background': CARD_BG, 'border': f'1px solid {BORDER}',
                           'borderRadius': '6px', **kwargs})


def _panel_header(text, subtitle=None):
    """Section header for tab panels."""
    return html.Div([
        html.H6(text, className='fw-semibold mb-0',
                style={'fontSize': '0.9rem', 'color': '#202124'}),
        html.Small(subtitle, className='text-muted', style={'fontSize': '0.75rem'})
        if subtitle else html.Span(),
    ], className='mb-2 mt-1')


# ──────────────────────────────────────────────────────────────────────────────
# Navbar
# ──────────────────────────────────────────────────────────────────────────────

def _navbar():
    return dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col(html.Div([
                    html.Span("ProForm Viewer",
                              className='fw-bold',
                              style={'fontSize': '1.1rem', 'color': '#202124',
                                     'letterSpacing': '-0.01em'}),
                    dbc.Badge("v1.0", color='primary', pill=True,
                              className='ms-2 align-middle',
                              style={'fontSize': '0.65rem', 'verticalAlign': 'middle'}),
                ]), width='auto'),
                dbc.Col(
                    html.Small("Top-Down Proteomics · Single-spectrum · One-protein Analysis",
                               style={'color': TEXT_MUTED, 'fontSize': '0.75rem',
                                      'fontWeight': '400'}),
                    width='auto'),
            ], align='center', className='g-2 flex-grow-1'),
            dbc.Row([
                dbc.Col(
                    dcc.Dropdown(
                        id='demo-file-select',
                        options=_DEMO_OPTIONS,
                        placeholder='Select demo dataset…',
                        clearable=False,
                        style={'minWidth': '235px', 'fontSize': '0.8rem'},
                    ),
                    width='auto'),
                dbc.Col(dcc.Loading(
                    id='load-demo-loading',
                    type='dot',
                    color='#4CAF50',
                    children=dbc.Button('Load', id='load-demo-btn', color='success',
                                        size='sm', outline=True,
                                        className='ms-1 me-1',
                                        style={'fontWeight': '600', 'paddingLeft': '14px',
                                               'paddingRight': '14px'}),
                ), width='auto'),
                dbc.Col(dbc.Button("⬇ Export", id='export-btn', color='secondary',
                                   size='sm', outline=True,
                                   style={'fontWeight': '500'}), width='auto'),
            ], align='center', className='g-1'),
        ], fluid=True),
        color='white', dark=False,
        style={'borderBottom': f'1px solid {BORDER}',
               'boxShadow': '0 1px 4px rgba(60,64,67,0.10)',
               'height': '58px'},
    )


# ──────────────────────────────────────────────────────────────────────────────
# Sidebar
# ──────────────────────────────────────────────────────────────────────────────

def _sidebar():
    ptm_options = [{'label': name, 'value': name}
                   for name in sorted(PTM_DATABASE.keys())]

    return html.Div([

        # ── Data Input (upload disabled — use navbar dropdown) ─────────────
        # Hidden placeholder components required by callbacks
        html.Div(style={'display': 'none'}, children=[
            dcc.Upload(id='upload-spectrum', children=[], multiple=False),
            html.Div(id='upload-spectrum-status'),
            html.Div(id='upload-spectrum-loading'),
            dcc.Upload(id='upload-features', children=[], multiple=False),
            html.Div(id='upload-features-status'),
            html.Div(id='upload-features-loading'),
        ]),

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
            _label("Search Mode"),
            dbc.RadioItems(
                id='search-mode',
                options=[
                    {'label': 'Targeted (single protein)', 'value': 'targeted'},
                    {'label': 'Database (FASTA)',           'value': 'database'},
                ],
                value='targeted',
                inline=False,
                className='small mb-2',
                inputStyle={'marginRight': '4px'},
            ),
            # Targeted controls (default visible)
            html.Div(id='targeted-protein-controls', children=[
                _label("Protein Name"),
                dbc.Input(id='protein-name', placeholder='e.g. Ubiquitin',
                          size='sm', className='mb-1', debounce=True,
                          style={'background': '#ffffff', 'color': '#111111',
                                 'border': '1px solid #cccccc'}),
                _label("Sequence (single-letter)"),
                dbc.Textarea(
                    id='protein-sequence',
                    placeholder='MQIFVKTLTGK…',
                    rows=4,
                    debounce=True,
                    style={'background': '#ffffff', 'color': '#111111',
                           'fontSize': '0.78rem', 'fontFamily': 'Courier New',
                           'border': '1px solid #cccccc', 'resize': 'vertical'},
                ),
                html.Div(id='protein-mass-display',
                         className='mt-1',
                         style={'fontSize': '0.72rem', 'color': TEXT_MUTED}),

                # ── N-terminal Modification quick-picker ──────────────────
                html.Hr(style={'margin': '8px 0', 'borderColor': '#e0e0e0'}),
                _label("N-terminal Modification"),
                html.Div([
                    # "No Mod" pill
                    dbc.Button('No Mod', id='nterm-mod-none',
                               size='sm', outline=True, color='secondary',
                               className='me-1 mb-1',
                               style={'fontSize': '0.7rem'}),
                    # Common mods
                    dbc.Button('Acetylation  +42.01', id='nterm-mod-ac',
                               size='sm', color='danger',
                               className='me-1 mb-1',
                               style={'fontSize': '0.7rem', 'background': '#e53935',
                                      'borderColor': '#e53935'}),
                    dbc.Button('Formylation  +27.99', id='nterm-mod-fo',
                               size='sm', color='warning',
                               className='me-1 mb-1',
                               style={'fontSize': '0.7rem', 'background': '#fb8c00',
                                      'borderColor': '#fb8c00'}),
                    dbc.Button('Trimethylation  +42.05', id='nterm-mod-tri',
                               size='sm', color='primary',
                               className='me-1 mb-1',
                               style={'fontSize': '0.7rem', 'background': '#8e24aa',
                                      'borderColor': '#8e24aa'}),
                    dbc.Button('Palmitate  +238.23', id='nterm-mod-palm',
                               size='sm', color='success',
                               className='me-1 mb-1',
                               style={'fontSize': '0.7rem', 'background': '#00897b',
                                      'borderColor': '#00897b'}),
                ], className='mb-1'),
                # Custom mass input
                dbc.InputGroup([
                    dbc.InputGroupText('Custom', style={'fontSize': '0.7rem'}),
                    dbc.Input(id='nterm-mod-custom-mass', type='number', step=0.0001,
                              placeholder='mass shift (Da)',
                              size='sm',
                              style={'background': '#fff', 'color': '#111',
                                     'border': '1px solid #ccc'}),
                    dbc.Button('Apply', id='nterm-mod-custom-btn',
                               size='sm', color='secondary', outline=True),
                ], size='sm', className='mb-1'),
                html.Div(id='nterm-mod-display',
                         style={'fontSize': '0.72rem', 'color': '#1a73e8',
                                'fontWeight': '600'}),
            ]),
            # Database / FASTA controls (hidden — disabled for now)
            html.Div(id='database-protein-controls', style={'display': 'none'}, children=[
                dbc.Button(
                    'Load Demo DB (Ubiquitin, HBB, BSA, HSA, Insulin)',
                    id='load-demo-db-btn',
                    color='success', outline=True, size='sm',
                    className='w-100 mb-2',
                    style={'fontSize': '0.72rem'},
                ),
                # FASTA upload disabled
                html.Div(style={'display': 'none'}, children=[
                    dcc.Upload(id='upload-fasta', children=[], multiple=False),
                    html.Div(id='upload-fasta-status'),
                    html.Div(id='upload-fasta-loading'),
                    html.Div(id='fasta-protein-count'),
                ]),
            ]),
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

            dbc.Row([
                dbc.Col([
                    _label("Deconvolute (OpenMS)"),
                    dbc.Switch(id='deconvolute-spectrum', value=False, label=''),
                ], width=12),
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
        'width': '285px', 'minWidth': '285px',
        'background': SIDEBAR_BG,
        'padding': '14px 14px',
        'overflowY': 'auto', 'height': 'calc(100vh - 58px)',
        'borderRight': f'1px solid {BORDER}',
        'flexShrink': '0',
    })


# ──────────────────────────────────────────────────────────────────────────────
# Tab contents
# ──────────────────────────────────────────────────────────────────────────────

def _tab_spectrum():
    return dbc.Tab(label='Spectrum', tab_id='tab-spectrum',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('MS Spectrum Viewer',
                                     'Interactive annotated spectrum with matched fragment ions'),
                       # ── Spectrum metrics pill bar ──────────────────────────────────
                       html.Div(id='spectrum-stats-bar', className='mb-2'),
                       dcc.Loading(
                           id='loading-spectrum-graph',
                           type='dot',
                           color='#1a73e8',
                           children=dcc.Graph(id='spectrum-graph',
                                             config=_GRAPH_CONFIG,
                                             style={'height': '620px'}),
                       ),
                       html.Div(id='peak-click-info',
                                className='mt-2 p-2',
                                style={'fontSize': '0.78rem', 'color': TEXT_MUTED,
                                       'background': '#f8f9fa', 'borderRadius': '4px',
                                       'minHeight': '28px'}),
                       html.Hr(style={'borderColor': BORDER, 'margin': '12px 0'}),
                       _panel_header('Charge-Deconvolved Spectrum',
                                     'Monoisotopic mass domain after isotope deconvolution'),
                       dcc.Loading(
                           id='loading-deconv-spectrum-graph',
                           type='dot',
                           color='#1a73e8',
                           children=dcc.Graph(id='deconv-spectrum-graph',
                                             config=_GRAPH_CONFIG,
                                             style={'height': '420px'}),
                       ),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_sequence():
    return dbc.Tab(label='Sequence', tab_id='tab-sequence',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('Proteoform Sequence Coverage',
                                     'ProSight Lite–style residue-level fragment ion coverage map'),
                       # Proteoform mass summary header
                       html.Div(id='sequence-mass-header',
                                className='mb-3 p-2',
                                style={'background': '#e8f0fe',
                                       'borderRadius': '6px',
                                       'fontSize': '0.83rem',
                                       'color': '#1557b0',
                                       'borderLeft': '4px solid #1a73e8',
                                       'minHeight': '30px'}),
                       dbc.Row([
                           dbc.Col([
                               dbc.Switch(id='show-cleavage', value=True,
                                          label='Cleavage ticks'),
                           ], width='auto'),
                           dbc.Col(html.Div(id='coverage-display',
                                            style={'fontSize': '0.83rem', 'color': '#202124',
                                                   'fontWeight': '500'}),
                                   width='auto'),
                           dbc.Col([
                               _label('Residues / row'),
                               dcc.Dropdown(
                                   id='seq-residues-per-row',
                                   options=[{'label': str(n), 'value': n}
                                            for n in [10, 15, 20, 25, 30]],
                                   value=20,
                                   clearable=False,
                                   style={'width': '80px', 'fontSize': '0.78rem'},
                               ),
                           ], width='auto'),
                           dbc.Col([
                               _label('Focus mass (Da)'),
                               dbc.InputGroup([
                                   dbc.Input(id='seq-focus-mass', type='number', step=0.001,
                                             placeholder='e.g. 8564.8',
                                             size='sm',
                                             style={'width': '130px', 'background': '#fff',
                                                    'color': '#202124',
                                                    'border': f'1px solid {BORDER}'}),
                                   dbc.Button('▶', id='seq-focus-mass-btn',
                                              size='sm', color='primary', outline=True),
                               ], size='sm'),
                           ], width='auto'),
                       ], className='mb-3 align-items-end'),
                       dcc.Loading(
                           id='loading-sequence-graph',
                           type='dot',
                           color='#1a73e8',
                           children=dcc.Graph(id='sequence-graph',
                                 config={
                                     'displayModeBar': True,
                                     'displaylogo': False,
                                     'modeBarButtonsToRemove': [
                                         'zoom2d', 'pan2d', 'select2d', 'lasso2d',
                                         'zoomIn2d', 'zoomOut2d', 'autoScale2d',
                                         'resetScale2d', 'toggleSpikelines',
                                         'hoverClosestCartesian', 'hoverCompareCartesian',
                                     ],
                                     'toImageButtonOptions': {
                                         'format': 'svg', 'scale': 2,
                                         'filename': 'sequence_coverage',
                                     },
                                 },
                                 style={'minHeight': '340px'}),
                       ),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_mirror():
    return dbc.Tab(label='Mirror Plot', tab_id='tab-mirror',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('Mirror Plot — Theoretical vs Experimental',
                                     'Compare predicted fragment ions against the measured spectrum'),
                       dbc.Row([
                           dbc.Col([
                               _label('m/z range — Min'),
                               dbc.Input(id='mirror-mz-min', type='number', value=200,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=2),
                           dbc.Col([
                               _label('m/z range — Max'),
                               dbc.Input(id='mirror-mz-max', type='number', value=2000,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=2),
                           dbc.Col([
                               _label('User Neutral Mass (Da)'),
                               dbc.InputGroup([
                                   dbc.Input(id='mirror-user-mass', type='number', step=0.001,
                                             placeholder='e.g. 8564.83',
                                             size='sm',
                                             style={'background': '#fff', 'color': '#202124',
                                                    'border': f'1px solid {BORDER}'}),
                                   dbc.Button('▶', id='mirror-annotate-btn',
                                              size='sm', color='primary', outline=True,
                                              title='Annotate charge series'),
                               ], size='sm'),
                           ], width=3),
                           dbc.Col([
                               _label('Max z'),
                               dbc.Input(id='mirror-user-charge-max', type='number', value=15,
                                         min=1, max=30, size='sm',
                                         style={'background': '#fff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=1),
                           dbc.Col([
                               _label(' '),
                               html.Div(id='mirror-user-mass-display',
                                        style={'fontSize': '0.75rem', 'color': TEXT_MUTED,
                                               'paddingTop': '4px'}),
                           ], width=4),
                       ], className='mb-3 align-items-end g-2'),
                       dcc.Loading(
                           id='loading-mirror-graph',
                           type='dot', color='#1a73e8',
                           children=dcc.Graph(id='mirror-graph', config=_GRAPH_CONFIG,
                                             style={'height': '700px'}),
                       ),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_features():
    return dbc.Tab(label='Feature Map', tab_id='tab-features',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('Feature Map — m/z × Retention Time',
                                     'Quantitative feature visualization with theoretical overlay'),
                       dbc.Row([
                           dbc.Col([
                               _label("Color by"),
                               dbc.RadioItems(
                                   id='feature-color-by',
                                   options=[{'label': 'Charge state', 'value': 'charge'},
                                            {'label': 'Proteoform ID', 'value': 'proteoform'}],
                                   value='charge',
                                   inline=True, className='small',
                               ),
                           ], width='auto'),
                           dbc.Col([
                               _label("Filter by Feature ID"),
                               dbc.Input(id='feature-filter', placeholder='e.g. F001',
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=3),
                           dbc.Col([
                               _label("Theoretical mass (Da)"),
                               dbc.Input(id='th-mass-input', type='number', step=0.001,
                                         placeholder='e.g. 8564.83',
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=3),
                           dbc.Col([
                               _label("Charge for trace"),
                               dbc.Input(id='th-charge-input', type='number', value=10, min=1,
                                         size='sm',
                                         style={'background': '#ffffff', 'color': '#202124',
                                                'border': f'1px solid {BORDER}'}),
                           ], width=2),
                       ], className='mb-3 align-items-end g-2'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='feature-map-graph', config=_GRAPH_CONFIG,
                                     style={'height': '560px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Extracted Ion Chromatogram (XIC)',
                                     'XIC extracted from MS1 scans for selected feature'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='xic-graph', config=_GRAPH_CONFIG,
                                     style={'height': '300px'})
                       ),
                       html.Div(id='xic-label',
                                className='mb-1 mt-2',
                                style={'fontSize': '0.8rem', 'color': TEXT_MUTED,
                                       'fontWeight': '500'}),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Elution Profile \u2014 Intensity Trace',
                                     'Gaussian-approximated elution trace for selected feature with theoretical overlay'),
                       dcc.Graph(id='intensity-trace-graph', config={'displayModeBar': False},
                                 style={'height': '320px'}),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Theoretical vs Observed — Isotope Envelope & Elution',
                                     'Isotope pattern mirror (top) and normalised elution profile overlay (bottom)'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='comparison-panel-graph', config=_GRAPH_CONFIG,
                                     style={'height': '600px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Mass Accuracy Map',
                                     'PPM error distribution of observed feature masses vs theoretical'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='mass-accuracy-graph', config=_GRAPH_CONFIG,
                                     style={'height': '320px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Feature 3D Plot — RT × Charge × Intensity',
                                     '3D scatter for multi-dimensional quantitative overview'),
                       dcc.Graph(id='feature-3d-graph', config=_GRAPH_CONFIG,
                                 style={'height': '580px'}),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_heatmap():
    return dbc.Tab(label='MS Heatmap', tab_id='tab-heatmap',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('MS Chromatographic Heatmaps',
                                     'RT × m/z (raw) and RT × neutral mass (deconvolved) intensity maps'),
                       dbc.Row([
                           dbc.Col([
                               _label('MS Level filter'),
                               dbc.RadioItems(
                                   id='heatmap-ms-level',
                                   options=[
                                       {'label': 'All', 'value': 0},
                                       {'label': 'MS1', 'value': 1},
                                       {'label': 'MS2', 'value': 2},
                                   ],
                                   value=2, inline=True, className='small',
                               ),
                           ], width='auto'),
                       ], className='mb-3'),
                       html.H6('Raw MS Heatmap (m/z × RT)',
                               className='text-muted fw-semibold small mb-1'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='raw-heatmap-graph', config=_GRAPH_CONFIG,
                                     style={'height': '520px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       html.H6('Deconvolved Precursor-Mass Heatmap (Neutral Mass × RT)',
                               className='text-muted fw-semibold small mt-1 mb-1'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='deconv-heatmap-graph', config=_GRAPH_CONFIG,
                                     style={'height': '520px'})
                       ),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


def _tab_search():
    columns = [
        {'name': 'Rank',        'id': 'rank'},
        {'name': 'Protein',     'id': 'protein'},
        {'name': 'Sequence',    'id': 'sequence'},
        {'name': 'Range',       'id': 'range'},
        {'name': '#AA',         'id': 'matched_aa'},
        {'name': 'Coverage %',  'id': 'coverage'},
        {'name': '# Tags',      'id': 'n_tags'},
        {'name': 'Score',       'id': 'score'},
        {'name': 'E-value',     'id': 'e_value'},
        {'name': 'q-value',     'id': 'q_value'},
        {'name': 'Matched',     'id': 'matched'},
        {'name': 'Th. Mass',    'id': 'th_mass'},
        {'name': 'Obs. Mass',   'id': 'obs_mass'},
        {'name': 'Δmass (ppm)', 'id': 'ppm_error'},
        {'name': 'PTMs',        'id': 'ptms'},
    ]
    return dbc.Tab(label='Search Results', tab_id='tab-search',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('Proteoform Search Results',
                                     'Ranked candidates with mass accuracy, coverage, and PTM annotations'),
                       html.Div(id='search-result-summary',
                                className='mb-2',
                                style={'fontSize': '0.82rem', 'color': TEXT_MUTED}),
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
                               'color': '#202124',
                               'border': f'1px solid {BORDER}',
                               'fontSize': '0.78rem',
                               'padding': '6px 10px',
                               'whiteSpace': 'normal',
                               'maxWidth': '220px',
                               'overflow': 'hidden',
                               'textOverflow': 'ellipsis',
                           },
                           style_header={
                               'backgroundColor': '#f8f9fa',
                               'color': '#202124',
                               'fontWeight': '600',
                               'border': f'1px solid {BORDER}',
                               'fontSize': '0.75rem',
                               'borderBottom': f'2px solid {BORDER}',
                           },
                           style_data_conditional=[
                               {'if': {'row_index': 0},
                                'backgroundColor': '#e8f0fe',
                                'color': '#1a73e8'},
                               {'if': {'state': 'selected'},
                                'backgroundColor': '#d2e3fc',
                                'border': '1px solid #1a73e8'},
                               {'if': {'filter_query': '{q_value} < 0.01', 'column_id': 'q_value'},
                                'color': '#2e7d32', 'fontWeight': 'bold'},
                               {'if': {'filter_query': '{q_value} >= 0.01 && {q_value} < 0.05',
                                       'column_id': 'q_value'},
                                'color': '#f57f17'},
                           ],
                           tooltip_data=[],
                           tooltip_duration=None,
                           sort_action='native',
                           filter_action='native',
                       ),
                       html.Div(
                           id='result-row-detail',
                           style={
                               'display': 'none',
                               'backgroundColor': '#e8f0fe',
                               'borderLeft': '3px solid #1a73e8',
                               'borderBottom': f'1px solid #b0c4de',
                               'padding': '8px 16px',
                               'marginBottom': '8px',
                               'fontSize': '0.78rem',
                               'borderRadius': '0 0 6px 0',
                           },
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       # ── FLASHTnT-style Truncation & Modification Ladder ────────
                       _panel_header('FLASHTnT — Truncation & Modification Ladder',
                                     'One-spectrum, one-protein: all candidate truncations and modifications ranked by tag score'),
                       html.Small(
                           'Each bar spans the matched sequence range (start → end positions). '
                           'Green = N-terminal truncation, Red = C-terminal truncation, '
                           'Purple = internal fragment, Gold star = PTM site. '
                           'Select a row above to highlight the corresponding candidate.',
                           className='text-muted d-block mb-2',
                           style={'fontSize': '0.75rem', 'lineHeight': '1.5'},
                       ),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='truncation-ladder-graph', config=_GRAPH_CONFIG,
                                     style={'height': '420px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       html.Div(id='fragment-stats-header', className='mb-1'),
                       dash_table.DataTable(
                           id='ion-table',
                           columns=[
                               {'name': 'Name',           'id': 'name'},
                               {'name': 'Ion type',       'id': 'ion_type'},
                               {'name': 'Ion number',     'id': 'ion_num'},
                               {'name': 'Th. mass (Da)',  'id': 'th_mass'},
                               {'name': 'Obs. mass (Da)', 'id': 'obs_mass'},
                               {'name': 'Δ mass (Da)',    'id': 'da_err'},
                               {'name': 'Δ mass (ppm)',   'id': 'ppm_err'},
                           ],
                           data=[],
                           page_size=12,
                           style_table={'overflowX': 'auto'},
                           style_cell={
                               'backgroundColor': '#ffffff',
                               'color': '#202124',
                               'border': f'1px solid {BORDER}',
                               'fontSize': '0.75rem',
                               'padding': '5px 10px',
                           },
                           style_header={
                               'backgroundColor': '#f8f9fa',
                               'color': '#202124',
                               'fontWeight': '600',
                               'border': f'1px solid {BORDER}',
                               'fontSize': '0.73rem',
                           },
                           style_data_conditional=[
                               {'if': {'filter_query': '{matched} eq "1"'},
                                'color': '#2e7d32'},
                           ],
                           sort_action='native',
                           filter_action='native',
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Score Distribution',
                                     'Histogram of all candidate scores'),
                       dcc.Graph(id='score-dist-graph', config={'displayModeBar': False},
                                 style={'height': '260px'}),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Target-Decoy Score Distribution & FDR Curve',
                                     'Score histogram (target vs decoy) + cumulative q-value; '
                                     'shaded region = identifications at 1 % FDR'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='fdr-curve-graph', config=_GRAPH_CONFIG,
                                     style={'height': '320px'})
                       ),
                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),
                       _panel_header('Internal Fragment Map',
                                     '2D dot map of matched internal fragment ions (start × end position)'),
                       dcc.Graph(id='internal-frag-graph', config=_GRAPH_CONFIG,
                                 style={'height': '500px'}),
                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


# ──────────────────────────────────────────────────────────────────────────────
# Diagnostics Tab
# ──────────────────────────────────────────────────────────────────────────────

def _tab_diagnostics():
    return dbc.Tab(label='Diagnostics', tab_id='tab-diagnostics',
                   children=dbc.Card(dbc.CardBody([
                       _panel_header('Diagnostic Dashboard',
                                     'Spectrum QC, TIC browser, precursor envelope, '
                                     'ion breakdown, sequence tag map'),

                       # ── TIC ───────────────────────────────────────────────
                       _panel_header('Total Ion Chromatogram',
                                     'Click any point to select that scan across all tabs'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='tic-graph', config=_GRAPH_CONFIG,
                                     style={'height': '240px'})
                       ),

                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),

                       # ── QC + precursor envelope (side-by-side) ────────────
                       dbc.Row([
                           dbc.Col([
                               _panel_header('Spectral QC Panel',
                                             'Dynamic range · matched-ion ppm errors · '
                                             'Δm/z spacing · intensity waterfall'),
                               dcc.Loading(type='dot', color='#1a73e8', children=
                                   dcc.Graph(id='diag-qc-graph', config=_GRAPH_CONFIG,
                                             style={'height': '500px'})
                               ),
                           ], md=8),
                           dbc.Col([
                               _panel_header('Precursor Isotope Envelope',
                                             'Zoom on precursor cluster with theoretical annotations'),
                               dcc.Loading(type='dot', color='#1a73e8', children=
                                   dcc.Graph(id='diag-precursor-envelope',
                                             config=_GRAPH_CONFIG,
                                             style={'height': '270px'})
                               ),
                               html.Hr(style={'borderColor': BORDER, 'margin': '12px 0'}),
                               _panel_header('Fragment Ion Type Breakdown',
                                             'Matched-ion count by type and charge state'),
                               dcc.Loading(type='dot', color='#1a73e8', children=
                                   dcc.Graph(id='diag-ion-breakdown',
                                             config={'displayModeBar': False},
                                             style={'height': '210px'})
                               ),
                           ], md=4),
                       ], className='mb-2 g-3'),

                       html.Hr(style={'borderColor': BORDER, 'margin': '14px 0'}),

                       # ── Sequence tag map ──────────────────────────────────
                       _panel_header('Sequence Tag Map',
                                     'Consecutive matched-ion runs (b/c lane & y/z lane) '
                                     '— longer unbroken bars = higher quality identifications'),
                       dcc.Loading(type='dot', color='#1a73e8', children=
                           dcc.Graph(id='diag-seq-tag-map', config=_GRAPH_CONFIG,
                                     style={'height': '220px'})
                       ),

                   ], className='p-3'), style={'background': CARD_BG, 'border': 'none'}))


# ──────────────────────────────────────────────────────────────────────────────
# Export Modal
# ──────────────────────────────────────────────────────────────────────────────

def _export_modal():
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("Export Results"),
                        style={'borderBottom': f'1px solid {BORDER}'}),
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
                inputStyle={'marginRight': '6px'},
            ),
        ], style={'padding': '20px'}),
        dbc.ModalFooter([
            dbc.Button("Download", id='export-confirm-btn', color='primary', size='sm',
                       style={'fontWeight': '600', 'paddingLeft': '18px', 'paddingRight': '18px'}),
            dbc.Button("Close",    id='export-close-btn',  color='secondary',
                       size='sm', className='ms-2', outline=True),
        ]),
    ], id='export-modal', is_open=False, size='md')


# ──────────────────────────────────────────────────────────────────────────────
# Stores (client-side session state)
# ──────────────────────────────────────────────────────────────────────────────

def _stores():
    return html.Div([
        dcc.Store(id='store-spectra',            storage_type='session', data=[]),
        dcc.Store(id='store-selected-scan-idx',  storage_type='session', data=0),
        dcc.Store(id='store-features',           storage_type='session', data=[]),
        dcc.Store(id='store-protein',            storage_type='session', data={}),
        dcc.Store(id='store-fasta-proteins',     storage_type='session', data=[]),
        dcc.Store(id='store-search-results',     storage_type='session', data=[]),
        dcc.Store(id='store-matched-ions',       storage_type='memory',  data=[]),
        dcc.Store(id='store-selected-result',    storage_type='session', data={}),
        dcc.Store(id='store-nterm-mod',          storage_type='session',
                  data={'name': 'No Modification', 'mass_shift': 0.0}),
        dcc.Download(id='download-data'),
        # Hidden labels needed by callbacks
        html.Div(id='comparison-panel-label', style={'display': 'none'}),
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Root layout
# ──────────────────────────────────────────────────────────────────────────────

def create_layout():
    tabs = dbc.Tabs([
        _tab_spectrum(),
        _tab_sequence(),
        _tab_heatmap(),
        _tab_mirror(),
        _tab_features(),
        _tab_search(),
        _tab_diagnostics(),
    ], id='main-tabs', active_tab='tab-spectrum',
    className='mb-0',
    style={'background': '#ffffff'})

    main_content = html.Div([
        tabs,
    ], style={
        'flex': '1',
        'padding': '10px 14px',
        'overflowY': 'auto',
        'height': 'calc(100vh - 58px)',
        'background': '#f4f5f7',
    })

    body = html.Div([
        _sidebar(),
        main_content,
    ], style={
        'display': 'flex',
        'flexDirection': 'row',
        'height': 'calc(100vh - 58px)',
    })

    return html.Div([
        _stores(),
        _navbar(),
        body,
        _export_modal(),
        dbc.Toast(
            id='status-toast',
            header='ProForm Viewer',
            is_open=False,
            dismissable=True,
            duration=4500,
            style={'position': 'fixed', 'bottom': '24px', 'right': '24px',
                   'background': '#ffffff', 'zIndex': 9999,
                   'boxShadow': '0 4px 12px rgba(60,64,67,0.20)',
                   'borderRadius': '8px', 'minWidth': '280px'},
        ),
    ], style={'background': '#f4f5f7', 'minHeight': '100vh'})
