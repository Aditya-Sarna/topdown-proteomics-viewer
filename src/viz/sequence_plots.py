"""
ProSight Lite–inspired sequence coverage viewer.

Each amino acid is rendered as a coloured square on a grid (20 residues/row).
Colours encode fragment-ion coverage type:
  Blue   → N-terminal (b/c/a)
  Red    → C-terminal (y/z)
  Purple → Both
  Dark   → Uncovered
PTMs are shown as gold asterisks above the residue.
Cleavage positions are shown as tick marks between residues.
"""
from typing import List, Optional
import numpy as np
import plotly.graph_objects as go

from ..data.models import Proteoform, FragmentIon, Modification
from ..analysis.peak_matching import coverage_map

RESIDUES_PER_ROW = 20
DARK_BG  = '#ffffff'
PLOT_BG  = '#ffffff'

_COLORS = {
    'n_term':    '#1976D2',  # blue
    'c_term':    '#D32F2F',  # red
    'both':      '#7B1FA2',  # purple
    'uncovered': '#e8e8e8',  # light gray
    'ptm':       '#FFD54F',  # gold
}

_LEGEND = [
    ('N-terminal (b/c/a)', _COLORS['n_term']),
    ('C-terminal (y/z)',   _COLORS['c_term']),
    ('Both ends',          _COLORS['both']),
    ('Uncovered',          _COLORS['uncovered']),
]


def create_sequence_plot(proteoform: Proteoform,
                          matched_ions: Optional[List[FragmentIon]] = None,
                          show_cleavage: bool = True) -> go.Figure:
    """
    Draw the sequence coverage grid plot.
    """
    fig = go.Figure()
    seq = proteoform.sequence

    if not seq:
        fig.update_layout(
            template='plotly_white',
            title='Enter protein sequence and run search to see coverage',
            paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
            font=dict(color='#111111'),
        )
        return fig

    n_aa     = len(seq)
    n_rows   = (n_aa + RESIDUES_PER_ROW - 1) // RESIDUES_PER_ROW
    cov      = coverage_map(seq, matched_ions or [])

    # Modification lookup: 0-based position -> list of mod names
    mod_lookup: dict = {}
    for m in proteoform.modifications:
        p = m.position - 1
        mod_lookup.setdefault(p, []).append(m.name)

    # Build grid coordinates
    xs, ys, colors, labels, hovers = [], [], [], [], []
    for i, aa in enumerate(seq):
        col = i % RESIDUES_PER_ROW
        row = i // RESIDUES_PER_ROW
        xs.append(col)
        ys.append(-row)
        labels.append(aa)

        covered = cov.get(i, [])
        has_n = any(t in ('b', 'c', 'a') for t in covered)
        has_c = any(t in ('y', 'z') for t in covered)

        if has_n and has_c:
            colors.append(_COLORS['both'])
        elif has_n:
            colors.append(_COLORS['n_term'])
        elif has_c:
            colors.append(_COLORS['c_term'])
        else:
            colors.append(_COLORS['uncovered'])

        mnames = mod_lookup.get(i, [])
        hov = f"<b>{aa}{i+1}</b>"
        if mnames:
            hov += f"<br>PTMs: {', '.join(mnames)}"
        hov += f"<br>Ions: {', '.join(covered) or 'none'}"
        hovers.append(hov)

    # Residue squares
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        mode='markers+text',
        marker=dict(symbol='square', size=26, color=colors,
                    line=dict(color='rgba(0,0,0,0.15)', width=1)),
        text=labels,
        textfont=dict(color='#111111', size=11, family='Courier New'),
        textposition='middle center',
        hovertext=hovers, hoverinfo='text',
        showlegend=False,
    ))

    # PTM asterisk markers
    ptm_xs, ptm_ys, ptm_hov = [], [], []
    for idx, mnames in mod_lookup.items():
        if 0 <= idx < n_aa:
            ptm_xs.append(idx % RESIDUES_PER_ROW)
            ptm_ys.append(-(idx // RESIDUES_PER_ROW) + 0.65)
            ptm_hov.append('; '.join(mnames))
    if ptm_xs:
        fig.add_trace(go.Scatter(
            x=ptm_xs, y=ptm_ys,
            mode='markers+text',
            marker=dict(size=0, opacity=0),
            text=['*'] * len(ptm_xs),
            textfont=dict(color=_COLORS['ptm'], size=16, family='Arial Black'),
            textposition='middle center',
            hovertext=ptm_hov, hoverinfo='text',
            showlegend=False, name='PTM',
        ))

    # Cleavage site ticks between residues
    if show_cleavage and matched_ions:
        cleavage_b, cleavage_y = set(), set()
        for ion in (matched_ions or []):
            if not ion.matched:
                continue
            if ion.ion_type in ('b', 'c', 'a') and ion.charge == 1:
                cleavage_b.add(ion.position)     # cleavage after position ion.position
            elif ion.ion_type in ('y', 'z') and ion.charge == 1:
                cleavage_y.add(n_aa - ion.position)  # 0-based left boundary

        for pos in cleavage_b:
            if 0 < pos < n_aa:
                col = (pos - 1) % RESIDUES_PER_ROW
                row = (pos - 1) // RESIDUES_PER_ROW
                fig.add_shape(type='line',
                    x0=col + 0.48, x1=col + 0.48,
                    y0=-row - 0.42, y1=-row + 0.42,
                    line=dict(color=_COLORS['n_term'], width=2))

        for pos in cleavage_y:
            if 0 < pos < n_aa:
                col = (pos - 1) % RESIDUES_PER_ROW
                row = (pos - 1) // RESIDUES_PER_ROW
                fig.add_shape(type='line',
                    x0=col + 0.48, x1=col + 0.48,
                    y0=-row - 0.42, y1=-row + 0.42,
                    line=dict(color=_COLORS['c_term'], width=2, dash='dot'))

    # Row position labels (left margin)
    for row_i in range(n_rows):
        fig.add_annotation(
            x=-1.2, y=-row_i,
            text=str(row_i * RESIDUES_PER_ROW + 1),
            showarrow=False,
            font=dict(color='#888888', size=9),
            xanchor='right',
        )

    # Legend traces
    for lbl, col in _LEGEND:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(symbol='square', size=12, color=col),
            name=lbl, showlegend=True,
        ))

    # Build title
    n_covered = sum(1 for v in coverage_map(seq, matched_ions or []).values() if v)
    cov_pct   = n_covered / n_aa * 100 if n_aa else 0
    n_mods    = len(proteoform.modifications)
    title = (f"<b>{proteoform.protein_name or 'Protein'}</b>  "
             f"{proteoform.start_pos}–{proteoform.end_pos}  |  "
             f"{n_aa} aa  |  Coverage: {cov_pct:.1f}%  |  "
             f"PTMs: {n_mods}")
    if proteoform.mass_error_ppm != 0:
        title += f"  |  Δmass: {proteoform.mass_error_ppm:.1f} ppm"

    height = max(300, n_rows * 52 + 120)

    fig.update_layout(
        template='plotly_white',
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-2, RESIDUES_PER_ROW + 0.5]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-(n_rows - 0.5), 1.2]),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        height=height,
        legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0,
                    font=dict(size=11)),
        margin=dict(l=50, r=20, t=80, b=20),
    )
    return fig
