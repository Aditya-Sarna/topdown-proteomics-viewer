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
from ..analysis.peak_matching import coverage_map, coverage_count_map

RESIDUES_PER_ROW = 20
DARK_BG  = '#ffffff'
PLOT_BG  = '#ffffff'

# Gradient colour palettes: index = coverage depth [0=none, 1=low, 2=mid, 3=high]
_N_FILL   = ['#F0F0F0', '#E3F2FD', '#90CAF9', '#42A5F5']   # N-terminal (b/c/a)
_N_BORD   = ['#cccccc', '#90CAF9', '#1976D2', '#1565C0']
_C_FILL   = ['#F0F0F0', '#FCE4EC', '#EF9A9A', '#EF5350']   # C-terminal (y/z)
_C_BORD   = ['#cccccc', '#EF9A9A', '#E53935', '#C62828']
_B_FILL   = ['#F0F0F0', '#F3E5F5', '#CE93D8', '#AB47BC']   # Both
_B_BORD   = ['#cccccc', '#CE93D8', '#8E24AA', '#7B1FA2']
_PTM_GOLD = '#F9A825'

# Legend uses the "mid" (index-2) shades
_LEGEND = [
    ('N-terminal (b/c/a)', _N_FILL[2]),
    ('C-terminal (y/z)',   _C_FILL[2]),
    ('Both ends',          _B_FILL[2]),
    ('Uncovered',          _N_FILL[0]),
]


def _depth_idx(n: int) -> int:
    if n == 0: return 0
    if n == 1: return 1
    if n <= 3: return 2
    return 3


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
    cov_cnt  = coverage_count_map(seq, matched_ions or [])

    # Modification lookup: 0-based position -> list of (name, mass_shift)
    mod_lookup: dict = {}
    for m in proteoform.modifications:
        p = m.position - 1
        mod_lookup.setdefault(p, []).append((m.name, m.mass_shift))

    # Build grid coordinates
    xs, ys, colors, border_colors, labels, hovers = [], [], [], [], [], []
    for i, aa in enumerate(seq):
        col = i % RESIDUES_PER_ROW
        row = i // RESIDUES_PER_ROW
        xs.append(col)
        ys.append(-row)
        labels.append(aa)

        nn, nc = cov_cnt.get(i, (0, 0))
        if nn > 0 and nc > 0:
            idx = _depth_idx(nn + nc)
            colors.append(_B_FILL[idx])
            border_colors.append(_B_BORD[idx])
        elif nn > 0:
            idx = _depth_idx(nn)
            colors.append(_N_FILL[idx])
            border_colors.append(_N_BORD[idx])
        elif nc > 0:
            idx = _depth_idx(nc)
            colors.append(_C_FILL[idx])
            border_colors.append(_C_BORD[idx])
        else:
            colors.append(_N_FILL[0])
            border_colors.append('#cccccc')

        # Gold border on PTM-modified cells
        if i in mod_lookup:
            border_colors[-1] = _PTM_GOLD

        mods_list = mod_lookup.get(i, [])
        hov = f"<b>{aa}{i+1}</b>"
        if mods_list:
            hov += '<br>PTMs: ' + ', '.join(f"{n} ({s:+.3f})" for n, s in mods_list)
        hov += f"<br>Ions: {', '.join(cov.get(i, [])) or 'none'}"
        hov += f"<br>Depth: N={nn} C={nc}"
        hovers.append(hov)

    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        mode='markers+text',
        marker=dict(symbol='square', size=26, color=colors,
                    line=dict(color=border_colors, width=1.5)),
        text=labels,
        textfont=dict(color='#111111', size=11, family='Courier New'),
        textposition='middle center',
        hovertext=hovers, hoverinfo='text',
        showlegend=False,
    ))

    # PTM mass-shift labels (replace asterisk with delta-mass text)
    ptm_xs, ptm_ys, ptm_texts, ptm_hov = [], [], [], []
    for idx, mods_list in mod_lookup.items():
        if 0 <= idx < n_aa:
            ptm_xs.append(idx % RESIDUES_PER_ROW)
            ptm_ys.append(-(idx // RESIDUES_PER_ROW) + 0.68)
            total_shift = sum(ms for _, ms in mods_list)
            # Compact format: +114.9 or +79.97 depending on magnitude
            label = f'{total_shift:+.1f}' if abs(total_shift) >= 10 else f'{total_shift:+.2f}'
            ptm_texts.append(label)
            ptm_hov.append(' | '.join(f"{nm} ({ms:+.4f})" for nm, ms in mods_list))
    if ptm_xs:
        fig.add_trace(go.Scatter(
            x=ptm_xs, y=ptm_ys,
            mode='markers+text',
            marker=dict(size=0, opacity=0),
            text=ptm_texts,
            textfont=dict(color=_PTM_GOLD, size=7, family='Arial'),
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
                    line=dict(color=_N_BORD[2], width=2))

        for pos in cleavage_y:
            if 0 < pos < n_aa:
                col = (pos - 1) % RESIDUES_PER_ROW
                row = (pos - 1) // RESIDUES_PER_ROW
                fig.add_shape(type='line',
                    x0=col + 0.48, x1=col + 0.48,
                    y0=-row - 0.42, y1=-row + 0.42,
                    line=dict(color=_C_BORD[2], width=2, dash='dot'))

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
    title_line1 = (f"<b>{proteoform.protein_name or 'Protein'}</b>  "
                   f"{proteoform.start_pos}–{proteoform.end_pos}  |  {n_aa} aa")
    title_line2 = f"Coverage: {cov_pct:.1f}%  |  PTMs: {n_mods}"
    if proteoform.mass_error_ppm != 0:
        title_line2 += f"  |  Δmass: {proteoform.mass_error_ppm:.1f} ppm"
    title = title_line1 + f"<br><sup>{title_line2}</sup>"

    height = max(300, n_rows * 52 + 140)

    fig.update_layout(
        template='plotly_white',
        title=dict(text=title, font=dict(size=12), pad=dict(b=6)),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-2, RESIDUES_PER_ROW + 0.5]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-(n_rows - 0.5), 1.2]),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        height=height,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
        margin=dict(l=50, r=20, t=100, b=20),
    )
    return fig


def create_internal_fragment_map(
        sequence: str,
        internal_ions,
) -> go.Figure:
    """
    2-D dot map of internal fragment ions.
    x = N-terminal cleavage start position (1-based)
    y = C-terminal end position (1-based)
    Green dot = matched, grey dot = unmatched.
    """
    fig = go.Figure()
    if not sequence or not internal_ions:
        fig.update_layout(
            title='Run search and select a result to view internal fragment map',
            template='plotly_white',
            paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
            font=dict(color='#aaaaaa'),
            height=420,
        )
        return fig

    xs_m, ys_m, hov_m = [], [], []  # matched
    xs_u, ys_u, hov_u = [], [], []  # unmatched

    for (i, j, mz_th, matched, obs_mz, ppm) in internal_ions:
        seq_frag = sequence[i:j]
        label    = f'{seq_frag} [{i+1}–{j}]<br>Th: {mz_th:.4f}<br>'
        if matched:
            label += f'Obs: {obs_mz:.4f}<br>Δ: {ppm:.1f} ppm'
            xs_m.append(i + 1)
            ys_m.append(j)
            hov_m.append(label)
        else:
            label += 'Not matched'
            xs_u.append(i + 1)
            ys_u.append(j)
            hov_u.append(label)

    if xs_u:
        fig.add_trace(go.Scatter(
            x=xs_u, y=ys_u, mode='markers',
            marker=dict(size=6, color='#cccccc', opacity=0.5),
            hovertext=hov_u, hoverinfo='text',
            name='Unmatched',
        ))
    if xs_m:
        fig.add_trace(go.Scatter(
            x=xs_m, y=ys_m, mode='markers',
            marker=dict(size=9, color='#4CAF50',
                        line=dict(color='#2e7d32', width=1)),
            hovertext=hov_m, hoverinfo='text',
            name='Matched',
        ))

    n = len(sequence)
    n_matched = len(xs_m)
    n_total   = len(xs_m) + len(xs_u)
    fig.update_layout(
        template='plotly_white',
        title=dict(
            text=(f'Internal Fragment Map — {n_matched}/{n_total} matched '
                  f'({100*n_matched/n_total:.0f}%)' if n_total else 'Internal Fragment Map'),
            font=dict(size=12),
        ),
        xaxis=dict(title='N-terminal start position', range=[0, n + 1],
                   showgrid=True, gridcolor='rgba(0,0,0,0.06)', automargin=True),
        yaxis=dict(title='C-terminal end position',   range=[0, n + 1],
                   showgrid=True, gridcolor='rgba(0,0,0,0.06)', automargin=True),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
        margin=dict(l=65, r=20, t=75, b=55),
        height=420,
    )
    return fig
