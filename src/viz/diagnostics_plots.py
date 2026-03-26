"""
Diagnostic / QC visualizations for ProForm Viewer.

  create_tic_plot          — Total Ion Chromatogram with current scan highlighted
  create_spectrum_qc       — 2×2 QC panel (dynamic range, ppm errors, Δm/z spacing, rank)
  create_ion_breakdown     — horizontal bar chart of matched-ion counts by type
  create_fdr_curve         — target-decoy score distribution + cumulative FDR (q-value) curve
  create_precursor_envelope — zoomed isotope cluster around precursor m/z
  create_sequence_tag_stats — consecutive matched-ion tag coverage visualisation
"""
from typing import List, Optional

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

_BG  = '#ffffff'
_PBG = '#ffffff'
_PROTON = 1.007276


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Total Ion Chromatogram
# ─────────────────────────────────────────────────────────────────────────────

def create_tic_plot(spectra: list,
                    selected_rt: Optional[float] = None,
                    selected_idx: Optional[int] = None) -> go.Figure:
    """
    Total Ion Chromatogram: sum of all peak intensities per scan vs RT.
    The currently selected scan is highlighted with a red dashed vline.
    Clicking on the TIC fires a `clickData` event used to update the scan selection.
    """
    fig = go.Figure()

    if not spectra:
        fig.update_layout(
            title='Load a spectrum file to see the Total Ion Chromatogram',
            template='plotly_white', height=230,
            paper_bgcolor=_BG, font=dict(color='#aaaaaa'),
            margin=dict(l=60, r=20, t=48, b=42),
        )
        return fig

    rts, tics, ms_levels = [], [], []
    for s in spectra:
        rt  = float(s.get('retention_time', 0.0))
        raw = s.get('intensity', [])
        tic = float(np.sum(raw)) if raw else 0.0
        rts.append(rt)
        tics.append(tic)
        ms_levels.append(int(s.get('ms_level', 2)))

    rts       = np.array(rts)
    tics      = np.array(tics)
    ms_levels = np.array(ms_levels)
    order     = np.argsort(rts)
    rts       = rts[order]
    tics      = tics[order]
    ms_levels = ms_levels[order]

    # MS1 base-peak chromatogram (blue)
    ms1_mask = ms_levels == 1
    if ms1_mask.any():
        fig.add_trace(go.Scatter(
            x=rts[ms1_mask], y=tics[ms1_mask],
            mode='lines',
            line=dict(color='#1a73e8', width=1.6),
            fill='tozeroy', fillcolor='rgba(26,115,232,0.08)',
            name='MS1',
            hovertemplate='RT: %{x:.2f} min<br>TIC: %{y:.2e}<extra>MS1</extra>',
        ))

    # MS2 events (orange ticks at bottom)
    ms2_mask = ms_levels == 2
    if ms2_mask.any():
        max_tic = tics.max() if tics.max() > 0 else 1.0
        tick_h  = max_tic * 0.06
        tick_xs, tick_ys = [], []
        for rt in rts[ms2_mask]:
            tick_xs += [rt, rt, None]
            tick_ys += [0, tick_h, None]
        fig.add_trace(go.Scatter(
            x=tick_xs, y=tick_ys,
            mode='lines',
            line=dict(color='rgba(230,81,0,0.50)', width=1.0),
            name='MS2',
            hoverinfo='skip',
        ))

    if selected_rt is not None:
        fig.add_vline(
            x=selected_rt,
            line_dash='dash', line_color='#E53935', line_width=1.5,
            annotation_text=f' RT {selected_rt:.2f}',
            annotation_font_color='#E53935',
            annotation_font_size=10,
        )

    fig.update_layout(
        template='plotly_white',
        title=dict(text='<b>Total Ion Chromatogram (TIC)</b> — click to select scan',
                   font=dict(size=12)),
        xaxis=dict(title='Retention Time (min)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.07)', automargin=True),
        yaxis=dict(title='Σ Intensity', showgrid=True,
                   gridcolor='rgba(0,0,0,0.07)', automargin=True),
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111'),
        height=240,
        margin=dict(l=65, r=25, t=55, b=45),
        legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
        hovermode='x unified',
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Spectrum QC — 2 × 2 panel
# ─────────────────────────────────────────────────────────────────────────────

def create_spectrum_qc(spectrum_dict: dict,
                       matched_ions_dicts: Optional[list] = None) -> go.Figure:
    """
    2×2 spectral quality-control figure:
      (1,1) Peak intensity distribution  (log₁₀ histogram → dynamic range)
      (1,2) Matched-ion mass error histogram (ppm) — spot systematic shifts
      (2,1) Adjacent Δm/z histogram — reveals dominant charge-state spacing
      (2,2) Intensity rank / waterfall — log-log, shows signal-to-noise profile
    """
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            'Peak Intensity Distribution',
            'Matched-Ion Mass Error (ppm)',
            'Adjacent Δm/z Histogram',
            'Intensity Rank (waterfall)',
        ],
        vertical_spacing=0.22,
        horizontal_spacing=0.13,
    )

    mz_arr  = np.asarray(spectrum_dict.get('mz',        []), dtype=np.float64)
    int_arr = np.asarray(spectrum_dict.get('intensity', []), dtype=np.float64)

    # Panel 1 — intensity log-histogram
    if int_arr.size > 0 and int_arr.max() > 0:
        log_int = np.log10(np.clip(int_arr, 1, None))
        fig.add_trace(go.Histogram(
            x=log_int, nbinsx=30,
            marker_color='#1a73e8', opacity=0.80,
            name='Peaks',
            hovertemplate='log₁₀(I): %{x:.1f}<br>Count: %{y}<extra></extra>',
        ), row=1, col=1)

    # Panel 2 — ppm error histogram for matched ions
    matched_ions_dicts = matched_ions_dicts or []
    ppms = [d.get('mass_error_ppm', 0.0)
            for d in matched_ions_dicts if d.get('matched')]
    if ppms:
        fig.add_trace(go.Histogram(
            x=ppms, nbinsx=min(30, max(5, len(ppms) // 2)),
            marker_color='#43A047', opacity=0.80,
            name='Matched',
            hovertemplate='ppm: %{x:.2f}<br>Count: %{y}<extra></extra>',
        ), row=1, col=2)
        fig.add_vline(x=0,  line_dash='dash', line_color='rgba(0,0,0,0.35)', row=1, col=2)
        median_ppm = float(np.median(ppms))
        rmse_ppm   = float(np.sqrt(np.mean(np.array(ppms) ** 2)))
        fig.add_annotation(
            xref='x2 domain', yref='y2 domain', x=0.97, y=0.95,
            text=f'Median: {median_ppm:+.2f} ppm<br>RMSE: {rmse_ppm:.2f} ppm',
            showarrow=False, font=dict(size=9, color='#333333'),
            xanchor='right', bgcolor='rgba(255,255,255,0.85)',
            bordercolor='rgba(0,0,0,0.12)', borderwidth=1, borderpad=4,
            align='left',
        )
    else:
        fig.add_annotation(
            xref='x2 domain', yref='y2 domain', x=0.5, y=0.5,
            text='Run a search to see mass errors',
            showarrow=False, font=dict(size=10, color='#aaaaaa'),
            xanchor='center',
        )

    # Panel 3 — adjacent Δm/z histogram (reveals isotope / charge spacing)
    if mz_arr.size > 2:
        diffs = np.diff(np.sort(mz_arr))
        diffs = diffs[(diffs > 0) & (diffs < 2.5)]
        if diffs.size > 0:
            fig.add_trace(go.Histogram(
                x=diffs, nbinsx=60,
                marker_color='#FB8C00', opacity=0.80,
                name='Δm/z',
                hovertemplate='Δm/z: %{x:.3f}<br>Count: %{y}<extra></extra>',
            ), row=2, col=1)
            # Annotate expected 1/z lines
            for z in [1, 2, 3, 4, 5]:
                exp = 1.003355 / z
                fig.add_vline(x=exp, line_dash='dot',
                              line_color=f'rgba(180,0,0,{0.5 - z*0.05:.2f})',
                              annotation_text=f'z={z}',
                              annotation_font_size=8,
                              annotation_font_color='#c0392b',
                              row=2, col=1)

    # Panel 4 — intensity rank (waterfall / Pareto)
    if int_arr.size > 0:
        sorted_int = np.sort(int_arr)[::-1]
        ranks      = np.arange(1, len(sorted_int) + 1)
        fig.add_trace(go.Scatter(
            x=ranks, y=sorted_int,
            mode='lines',
            line=dict(color='#8E24AA', width=1.5),
            name='Rank',
            hovertemplate='Rank: %{x}<br>I: %{y:.2e}<extra></extra>',
        ), row=2, col=2)

    fig.update_layout(
        template='plotly_white',
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111', size=10),
        height=500,
        margin=dict(l=55, r=20, t=80, b=50),
        showlegend=False,
    )
    fig.update_xaxes(title_text='log₁₀(Intensity)',    row=1, col=1)
    fig.update_xaxes(title_text='Mass Error (ppm)',     row=1, col=2)
    fig.update_xaxes(title_text='Δm/z to next peak',   row=2, col=1)
    fig.update_xaxes(title_text='Rank',                 row=2, col=2)
    fig.update_yaxes(title_text='Peak count',           row=1, col=1)
    fig.update_yaxes(title_text='Ion count',            row=1, col=2)
    fig.update_yaxes(title_text='Peak count',           row=2, col=1)
    fig.update_yaxes(title_text='Intensity',            type='log', row=2, col=2)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Matched-Ion Type Breakdown
# ─────────────────────────────────────────────────────────────────────────────

def create_ion_breakdown(matched_ions_dicts: list) -> go.Figure:
    """
    Horizontal bar chart: count of matched fragment ions split by type and charge.
    Shows per-type subtotals and per-charge breakdown as stacked bars.
    """
    fig = go.Figure()

    palette = {
        'b': '#1a73e8', 'y': '#E53935', 'c': '#43A047',
        'z': '#FB8C00', 'a': '#8E24AA', 'x': '#00ACC1', '?': '#9E9E9E',
    }

    if not matched_ions_dicts:
        fig.update_layout(
            title='No matched ions — run a search first',
            template='plotly_white', height=220,
            paper_bgcolor=_BG, font=dict(color='#aaaaaa'),
            margin=dict(l=75, r=30, t=48, b=35),
        )
        return fig

    from collections import defaultdict
    counts     = defaultdict(int)
    z_breakdown = defaultdict(lambda: defaultdict(int))

    for d in matched_ions_dicts:
        if not d.get('matched'):
            continue
        t = d.get('ion_type', '?')
        z = int(d.get('charge', 1))
        counts[t] += 1
        z_breakdown[t][z] += 1

    ion_order = [t for t in ['b', 'y', 'c', 'z', 'a', 'x']
                 if t in counts] + [t for t in counts if t not in 'byczax']

    if not ion_order:
        fig.update_layout(
            title='No matched ions found',
            template='plotly_white', height=220,
            paper_bgcolor=_BG, font=dict(color='#aaaaaa'),
            margin=dict(l=75, r=30, t=48, b=35),
        )
        return fig

    all_charges = sorted({z for t in ion_order for z in z_breakdown[t]})

    for z in all_charges:
        vals = [z_breakdown[t].get(z, 0) for t in ion_order]
        fig.add_trace(go.Bar(
            y=ion_order,
            x=vals,
            orientation='h',
            name=f'z={z}',
            marker_color=px.colors.qualitative.Plotly[z % len(px.colors.qualitative.Plotly)],
            opacity=0.85,
            hovertemplate='%{y}-ions z=%{name}: %{x}<extra></extra>',
        ))

    # Total labels at end of each bar
    totals = [counts[t] for t in ion_order]
    fig.add_trace(go.Bar(
        y=ion_order, x=[0] * len(ion_order),
        orientation='h',
        text=[str(v) for v in totals],
        textposition='outside',
        marker_color='rgba(0,0,0,0)',
        showlegend=False,
        hoverinfo='skip',
    ))

    fig.update_layout(
        template='plotly_white',
        title=dict(text='<b>Matched Fragment Ion Breakdown</b> by Type & Charge',
                   font=dict(size=12)),
        xaxis=dict(title='Ion count', automargin=True,
                   showgrid=True, gridcolor='rgba(0,0,0,0.07)'),
        yaxis=dict(title='Ion type', autorange='reversed'),
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111'),
        barmode='stack',
        height=max(220, len(ion_order) * 42 + 100),
        margin=dict(l=55, r=55, t=55, b=45),
        legend=dict(title=dict(text='Charge', font=dict(size=10)),
                    font=dict(size=10), orientation='v',
                    bgcolor='rgba(255,255,255,0.7)'),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Target-Decoy FDR Curve
# ─────────────────────────────────────────────────────────────────────────────

def create_fdr_curve(results_data: list) -> go.Figure:
    """
    Two-panel FDR figure:
      Left  — Target vs Decoy score distribution (overlapping histograms)
      Right — Cumulative q-value curve: rank vs q-value with 1 % / 5 % FDR guides
    """
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            'Score Distribution — Target vs Decoy',
            'Cumulative FDR (q-value) Curve',
        ],
        horizontal_spacing=0.14,
    )

    if not results_data:
        fig.update_layout(
            title='Run a search to see the FDR curve',
            template='plotly_white', paper_bgcolor=_BG,
            font=dict(color='#aaaaaa'), height=300,
            margin=dict(l=55, r=20, t=55, b=45),
        )
        return fig

    targets = [r['proteoform']['score'] for r in results_data
               if '[DECOY]' not in r['proteoform'].get('protein_name', '')]
    decoys  = [r['proteoform']['score'] for r in results_data
               if '[DECOY]' in r['proteoform'].get('protein_name', '')]

    if targets:
        fig.add_trace(go.Histogram(
            x=targets, nbinsx=max(10, len(targets) // 3),
            name='Target',
            marker_color='rgba(26,115,232,0.70)',
            opacity=0.80,
            hovertemplate='Score: %{x:.2f}<br>Count: %{y}<extra>Target</extra>',
        ), row=1, col=1)

    if decoys:
        fig.add_trace(go.Histogram(
            x=decoys, nbinsx=max(6, len(decoys) // 2),
            name='Decoy',
            marker_color='rgba(229,57,53,0.65)',
            opacity=0.80,
            hovertemplate='Score: %{x:.2f}<br>Count: %{y}<extra>Decoy</extra>',
        ), row=1, col=1)

    # FDR curve — rank vs q-value
    pairs = [(r['proteoform']['score'], r.get('q_value', 1.0))
             for r in results_data]
    pairs.sort(key=lambda x: -x[0])
    if pairs:
        sc_vals, q_vals = zip(*pairs)
        cum_targets = list(range(1, len(sc_vals) + 1))
        fig.add_trace(go.Scatter(
            x=cum_targets,
            y=list(q_vals),
            mode='lines',
            line=dict(color='#1a73e8', width=2.2),
            name='q-value',
            hovertemplate='Rank: %{x}<br>q: %{y:.4f}<extra></extra>',
        ), row=1, col=2)

        # 1 % and 5 % FDR reference lines
        fig.add_hline(y=0.01, line_dash='dash', line_color='#43A047',
                      annotation_text='q = 0.01 (1 % FDR)',
                      annotation_font_size=9,
                      annotation_font_color='#2e7d32',
                      row=1, col=2)
        fig.add_hline(y=0.05, line_dash='dot', line_color='#FB8C00',
                      annotation_text='q = 0.05 (5 % FDR)',
                      annotation_font_size=9,
                      annotation_font_color='#e65100',
                      row=1, col=2)

        # Shade passing at 1 % FDR
        n_pass_1pct = sum(1 for q in q_vals if q <= 0.01)
        if n_pass_1pct:
            fig.add_vrect(
                x0=1, x1=n_pass_1pct,
                fillcolor='rgba(67,160,71,0.08)',
                line_width=0,
                annotation_text=f'{n_pass_1pct} @ 1 %',
                annotation_font_size=9,
                annotation_font_color='#2e7d32',
                row=1, col=2,
            )

    fig.update_layout(
        template='plotly_white',
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111', size=11),
        barmode='overlay',
        height=320,
        margin=dict(l=55, r=20, t=70, b=50),
        legend=dict(orientation='h', yanchor='bottom', y=1.08, x=0,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
    )
    fig.update_xaxes(title_text='Score',  row=1, col=1)
    fig.update_yaxes(title_text='Count',  row=1, col=1)
    fig.update_xaxes(title_text='Rank',   row=1, col=2)
    fig.update_yaxes(title_text='q-value (FDR)', row=1, col=2)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Precursor Isotope Envelope Viewer
# ─────────────────────────────────────────────────────────────────────────────

def create_precursor_envelope(spectrum_dict: dict,
                               precursor_mz: float = 0.0,
                               precursor_z:  int   = 0,
                               window_da:    float = 8.0) -> go.Figure:
    """
    Zoom-in on the precursor isotope cluster.
    Observed sticks are blue; theoretical isotope positions (if z known) are
    annotated with red dotted vlines labelled M, M+1, …
    """
    fig = go.Figure()
    empty = dict(template='plotly_white', height=270, paper_bgcolor=_BG,
                 font=dict(color='#aaaaaa'),
                 margin=dict(l=60, r=20, t=50, b=42))

    mz_arr  = np.asarray(spectrum_dict.get('mz',        []), dtype=np.float64)
    int_arr = np.asarray(spectrum_dict.get('intensity', []), dtype=np.float64)

    if mz_arr.size == 0:
        fig.update_layout(title='No spectrum loaded', **empty)
        return fig

    # Set window
    if precursor_mz > 0:
        lo = precursor_mz - window_da / 2
        hi = precursor_mz + window_da / 2
    else:
        apex = mz_arr[np.argmax(int_arr)]
        lo, hi = apex - 4.0, apex + 4.0

    mask   = (mz_arr >= lo) & (mz_arr <= hi)
    zm_mz  = mz_arr[mask]
    zm_int = int_arr[mask]

    if zm_mz.size == 0:
        apex = mz_arr[np.argmax(int_arr)]
        mask   = (mz_arr >= apex - 4) & (mz_arr <= apex + 4)
        zm_mz  = mz_arr[mask]
        zm_int = int_arr[mask]

    norm = zm_int / zm_int.max() * 100.0 if zm_int.size > 0 and zm_int.max() > 0 else zm_int

    # Stick traces
    xs, ys = [], []
    for m, h in zip(zm_mz, norm):
        xs += [float(m), float(m), None]
        ys += [0.0, float(h), None]

    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='lines',
        line=dict(color='#1a73e8', width=1.8),
        name='Observed',
        hoverinfo='skip',
    ))
    fig.add_trace(go.Scatter(
        x=zm_mz, y=norm, mode='markers',
        marker=dict(color='#1a73e8', size=5, symbol='circle'),
        showlegend=False,
        customdata=zm_int,
        hovertemplate='m/z: %{x:.4f}<br>Rel.Int: %{y:.1f}%<br>Abs: %{customdata:.2e}<extra></extra>',
    ))

    # Theoretical isotope lines (if precursor_z known)
    if precursor_mz > 0 and precursor_z > 0:
        spacing = 1.003355 / precursor_z
        for i in range(-1, 8):
            th = precursor_mz + i * spacing
            if lo <= th <= hi:
                label = f'M{i:+d}' if i != 0 else 'M₀'
                fig.add_vline(
                    x=th, line_dash='dot',
                    line_color='rgba(229,57,53,0.45)',
                    annotation_text=label,
                    annotation_font_size=9,
                    annotation_font_color='#c62828',
                    annotation_position='top right',
                )

    # Precursor marker
    if precursor_mz > 0 and lo <= precursor_mz <= hi:
        fig.add_vline(
            x=precursor_mz, line_dash='dash',
            line_color='rgba(0,0,0,0.30)',
            annotation_text=f'prec {precursor_mz:.3f}',
            annotation_font_size=9,
        )

    z_info = f'  z = {precursor_z}' if precursor_z > 0 else ''
    title_text = (
        f'<b>Precursor Isotope Envelope</b>'
        + (f' — m/z {precursor_mz:.4f}{z_info}' if precursor_mz > 0 else '')
    )
    fig.update_layout(
        template='plotly_white',
        title=dict(text=title_text, font=dict(size=12)),
        xaxis=dict(title='m/z', showgrid=True, gridcolor='rgba(0,0,0,0.08)',
                   automargin=True),
        yaxis=dict(title='Rel. Intensity (%)', range=[-3, 112],
                   showgrid=True, gridcolor='rgba(0,0,0,0.08)'),
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111'),
        height=270,
        margin=dict(l=60, r=20, t=60, b=45),
        showlegend=False,
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Sequence Tag Statistics (consecutive matched-ion blocks)
# ─────────────────────────────────────────────────────────────────────────────

def create_sequence_tag_stats(sequence: str,
                               matched_ions_dicts: list) -> go.Figure:
    """
    Show matched-ion runs (consecutive positions) as horizontal coloured bars
    on top of the sequence.

    A 'tag' is a maximal run of consecutive matched b/c/y/z ions.
    Colours encode ion type (blue=b/c, red=y/z).
    The chart shows the sequence along x and ion positions along y,
    making sequence coverage gaps immediately visible.
    """
    fig = go.Figure()
    empty = dict(template='plotly_white', height=220, paper_bgcolor=_BG,
                 font=dict(color='#aaaaaa'),
                 margin=dict(l=60, r=20, t=50, b=42))

    if not sequence or not matched_ions_dicts:
        fig.update_layout(title='Run a search to see sequence tag statistics', **empty)
        return fig

    n = len(sequence)
    matched_b  = set()
    matched_y  = set()

    for d in matched_ions_dicts:
        if not d.get('matched'):
            continue
        pos  = int(d.get('position', 0))
        itype = d.get('ion_type', '')
        if itype in ('b', 'c', 'a') and 0 < pos <= n:
            matched_b.add(pos)
        elif itype in ('y', 'z') and 0 < pos <= n:
            matched_y.add(n - pos)   # convert to N→C position

    def _find_runs(positions, n_total):
        runs = []
        if not positions:
            return runs
        sorted_p = sorted(positions)
        start = sorted_p[0]
        end   = sorted_p[0]
        for p in sorted_p[1:]:
            if p == end + 1:
                end = p
            else:
                runs.append((start, end))
                start = end = p
        runs.append((start, end))
        return runs

    b_runs = _find_runs(matched_b, n)
    y_runs = _find_runs(matched_y, n)

    LANE_B = 1.3
    LANE_Y = 0.7

    # b/c runs (blue, top lane)
    for start, end in b_runs:
        length = end - start + 1
        fig.add_trace(go.Bar(
            y=[LANE_B], x=[length], base=[start - 1],
            orientation='h',
            marker_color='rgba(26,115,232,0.75)',
            marker_line=dict(color='rgba(26,115,232,0.9)', width=1),
            name='b/c tag' if start == b_runs[0][0] else '',
            showlegend=(start == b_runs[0][0]),
            legendgroup='b',
            hovertemplate=f'b/c tag: pos {start}–{end} ({length} aa)<extra></extra>',
            width=0.4,
        ))

    # y/z runs (red, bottom lane)
    for start, end in y_runs:
        length = end - start + 1
        fig.add_trace(go.Bar(
            y=[LANE_Y], x=[length], base=[start - 1],
            orientation='h',
            marker_color='rgba(229,57,53,0.70)',
            marker_line=dict(color='rgba(229,57,53,0.9)', width=1),
            name='y/z tag' if start == y_runs[0][0] else '',
            showlegend=(start == y_runs[0][0]),
            legendgroup='y',
            hovertemplate=f'y/z tag: pos {start}–{end} ({length} aa)<extra></extra>',
            width=0.4,
        ))

    # Sequence x-axis ruler
    fig.add_trace(go.Scatter(
        x=list(range(1, n + 1)), y=[1.0] * n,
        mode='markers',
        marker=dict(size=4, color='rgba(0,0,0,0.12)', symbol='line-ns-open'),
        showlegend=False,
        hoverinfo='skip',
    ))

    n_bc = len(matched_b)
    n_yz = len(matched_y)
    longest_b = max((e - s + 1 for s, e in b_runs), default=0)
    longest_y = max((e - s + 1 for s, e in y_runs), default=0)

    fig.update_layout(
        template='plotly_white',
        title=dict(
            text=(f'<b>Sequence Tag Map</b>  —  '
                  f'b/c: {n_bc} pos, longest tag {longest_b} aa  |  '
                  f'y/z: {n_yz} pos, longest tag {longest_y} aa'),
            font=dict(size=11),
        ),
        xaxis=dict(title='Sequence position (N→C)', range=[0, n + 1],
                   showgrid=True, gridcolor='rgba(0,0,0,0.06)', automargin=True),
        yaxis=dict(
            tickvals=[LANE_Y, LANE_B],
            ticktext=['y/z ions', 'b/c ions'],
            range=[0.2, 2.0],
            showgrid=False,
            automargin=True,
        ),
        barmode='overlay',
        paper_bgcolor=_BG, plot_bgcolor=_PBG,
        font=dict(color='#111111'),
        height=220,
        margin=dict(l=75, r=20, t=70, b=48),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
    )
    return fig
