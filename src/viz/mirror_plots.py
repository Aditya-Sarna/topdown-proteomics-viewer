"""
Mirror plot: theoretical spectrum (top, pointing up) vs
            experimental spectrum (bottom, pointing down).
Matching peaks are coloured by ion type; non-matching experimental peaks are grey.
"""
from typing import List, Optional
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..data.models import Spectrum, FragmentIon

ION_COLORS = {
    'b': '#2196F3', 'y': '#F44336',
    'c': '#4CAF50', 'z': '#FF9800', 'a': '#9C27B0',
}
UNMATCHED_COLOR = '#aaaaaa'
DARK_BG = '#ffffff'
PLOT_BG = '#ffffff'


def _sticks(mz_list, height=100.0):
    xs, ys = [], []
    for m in mz_list:
        xs += [m, m, None]
        ys += [0, height, None]
    return xs, ys


def create_mirror_plot(spectrum: Spectrum,
                        matched_ions: Optional[List[FragmentIon]] = None,
                        mz_range: Optional[tuple] = None) -> go.Figure:
    """
    Mirror plot with shared x-axis.
    Top panel  → theoretical fragment ions (stick, coloured by type).
    Bottom panel → experimental spectrum (inverted sticks, coloured where matched).
    """
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.04,
        subplot_titles=('<b>Theoretical</b>', '<b>Experimental</b>'),
        row_heights=[0.5, 0.5],
    )

    obs_mz  = spectrum.mz_array
    obs_int = spectrum.normalize()

    if mz_range is None:
        if len(obs_mz) > 0:
            mz_range = (max(50, obs_mz.min() - 50), obs_mz.max() + 50)
        else:
            mz_range = (200, 2000)

    # ---- Theoretical (top) --------------------------------
    if matched_ions:
        by_type: dict = {}
        for ion in matched_ions:
            # Always draw theoretical peaks for all calculated ions (matched or not)
            by_type.setdefault(ion.ion_type, []).append((ion.mz, ion.matched))

        shown_legend = set()
        for itype, peaks in by_type.items():
            matched_mz = [p[0] for p in peaks if p[1] and mz_range[0] <= p[0] <= mz_range[1]]
            unmatched_mz = [p[0] for p in peaks if not p[1] and mz_range[0] <= p[0] <= mz_range[1]]
            color = ION_COLORS.get(itype, '#666666')

            if matched_mz:
                xs, ys = _sticks(matched_mz, 100.0)
                fig.add_trace(go.Scatter(
                    x=xs, y=ys, mode='lines',
                    line=dict(color=color, width=2),
                    name=f"{itype.upper()}-ions",
                    legendgroup=itype,
                    showlegend=(itype not in shown_legend),
                ), row=1, col=1)
                shown_legend.add(itype)

            if unmatched_mz:
                xs, ys = _sticks(unmatched_mz, 100.0)
                fig.add_trace(go.Scatter(
                    x=xs, y=ys, mode='lines',
                    line=dict(color=color, width=1, dash='dot'),
                    opacity=0.3,
                    legendgroup=itype,
                    showlegend=False,
                ), row=1, col=1)

            if matched_mz:
                candidate_labels = [
                    (ion.mz, ion.label())
                    for ion in matched_ions
                    if ion.matched and ion.ion_type == itype
                    and mz_range[0] <= ion.mz <= mz_range[1]
                ]
                span = mz_range[1] - mz_range[0] if mz_range else 1800.0
                min_gap = span * 0.018
                placed: list = []
                filtered_lx, filtered_ly, filtered_lt = [], [], []
                for lm, ll in sorted(candidate_labels, key=lambda x: x[0]):
                    if any(abs(lm - pm) < min_gap for pm in placed):
                        continue
                    placed.append(lm)
                    filtered_lx.append(lm)
                    filtered_ly.append(106)
                    filtered_lt.append(ll)
                if filtered_lx:
                    fig.add_trace(go.Scatter(
                        x=filtered_lx, y=filtered_ly,
                        mode='text',
                        text=filtered_lt,
                        textfont=dict(size=8, color=color),
                        showlegend=False, legendgroup=itype,
                        hoverinfo='skip',
                    ), row=1, col=1)
    else:
        fig.add_annotation(
            text="Run a search to see matched theoretical ions",
            xref='paper', yref='paper', x=0.5, y=0.75,
            showarrow=False, font=dict(size=13, color='#666666'),
        )

    # ---- Experimental (bottom) ----------------------------
    if len(obs_mz) > 0:
        mask = (obs_mz >= mz_range[0]) & (obs_mz <= mz_range[1])
        obs_mz_f, obs_int_f = obs_mz[mask], obs_int[mask]

        # Colour each experimental peak if it matches a theoretical one
        matched_set = {}
        if matched_ions:
            for ion in matched_ions:
                if ion.matched:
                    matched_set[round(ion.observed_mz, 2)] = ion.ion_type

        # Group into: matched (per type) and unmatched
        by_type_exp: dict = {}
        unk_exp = []
        for m, i in zip(obs_mz_f, obs_int_f):
            key = round(m, 2)
            nearest_key = min(matched_set.keys(), key=lambda k: abs(k - key), default=None)
            if nearest_key and abs(nearest_key - key) < 0.05:
                itype = matched_set[nearest_key]
                by_type_exp.setdefault(itype, {'mz': [], 'int': []})
                by_type_exp[itype]['mz'].append(m)
                by_type_exp[itype]['int'].append(i)
            else:
                unk_exp.append((m, i))

        # Unmatched experimental
        if unk_exp:
            ux = [p[0] for p in unk_exp]
            ui = [-p[1] for p in unk_exp]
            xs, ys = _sticks([], 0)
            xs2, ys2 = [], []
            for m, i in zip(ux, ui):
                xs2 += [m, m, None]
                ys2 += [0, i, None]
            fig.add_trace(go.Scatter(
                x=xs2, y=ys2, mode='lines',
                line=dict(color=UNMATCHED_COLOR, width=1.2),
                name='Experimental', showlegend=True,
            ), row=2, col=1)

        # Matched experimental (coloured)
        for itype, data in by_type_exp.items():
            color = ION_COLORS.get(itype, UNMATCHED_COLOR)
            xs2, ys2 = [], []
            for m, i in zip(data['mz'], data['int']):
                xs2 += [m, m, None]
                ys2 += [0, -i, None]
            fig.add_trace(go.Scatter(
                x=xs2, y=ys2, mode='lines',
                line=dict(color=color, width=2),
                legendgroup=itype, showlegend=False,
                name=f"{itype.upper()}-matched",
            ), row=2, col=1)

    # Zero divider line
    fig.add_shape(type='line', x0=mz_range[0], x1=mz_range[1],
                  y0=0, y1=0, yref='y2',
                  line=dict(color='rgba(0,0,0,0.15)', width=1))

    # Axes
    fig.update_yaxes(title_text='Rel. Intensity (%)', row=1, col=1,
                     range=[-5, 115],
                     gridcolor='rgba(0,0,0,0.08)')
    fig.update_yaxes(title_text='Rel. Intensity (%)',
                     tickvals=[-100, -75, -50, -25, 0],
                     ticktext=['100', '75', '50', '25', '0'],
                     row=2, col=1,
                     gridcolor='rgba(0,0,0,0.08)')
    fig.update_xaxes(title_text='m/z', row=2, col=1,
                     gridcolor='rgba(0,0,0,0.08)')

    fig.update_layout(
        template='plotly_white',
        title='<b>Mirror Plot</b> — Theoretical (top) vs Experimental (bottom)',
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        height=580,
        hovermode='x unified',
        legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0,
                    font=dict(size=11)),
        margin=dict(l=60, r=20, t=80, b=50),
    )
    return fig
