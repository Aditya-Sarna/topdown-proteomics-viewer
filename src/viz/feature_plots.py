"""
Feature map visualizations:
  1. 2-D scatter: m/z vs retention time (dot area ∝ intensity)
  2. Intensity trace: Gaussian-approximated elution profile for selected feature
  3. Comparison panel: theoretical vs observed isotope envelope + elution profile
"""
from typing import List, Optional
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.interpolate import interp1d

from ..data.models import Feature, Proteoform
from ..analysis.mass_utils import mass_to_mz, possible_charges, generate_isotope_envelope
from ..data.amino_acids import PROTON

DARK_BG = '#ffffff'
PLOT_BG = '#ffffff'


def create_feature_map(features: List[Feature],
                        proteoforms: Optional[List[Proteoform]] = None,
                        selected_id: Optional[str] = None,
                        color_by: str = 'charge') -> go.Figure:
    """
    2-D feature map: m/z vs retention time.
    Dot size  → log10(intensity)
    Dot color → charge state OR proteoform identity
    """
    fig = go.Figure()

    if not features:
        fig.update_layout(
            template='plotly_white',
            title='No features loaded — upload a feature table or click Load Demo',
            xaxis_title='m/z', yaxis_title='Retention Time (min)',
            paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
            font=dict(color='#111111'),
        )
        return fig

    rows = []
    for f in features:
        rows.append({
            'feature_id': f.feature_id,
            'mz':   f.mz_apex,
            'rt':   f.rt_apex,
            'intensity': f.intensity,
            'charge': str(f.charge),
            'mass':   f.monoisotopic_mass,
            'proteoform_id': f.proteoform_id or 'Unknown',
            'seq_short': (f.sequence[:15] + '…') if len(f.sequence) > 15 else f.sequence,
        })
    df = pd.DataFrame(rows)

    max_int = df['intensity'].max()
    df['size'] = np.log10(df['intensity'].clip(1) + 1) / np.log10(max_int + 2) * 22 + 6
    df['alpha'] = df['feature_id'].apply(
        lambda fid: 1.0 if selected_id is None or fid == selected_id else 0.35)

    palette = px.colors.qualitative.Plotly
    group_col = 'charge' if color_by == 'charge' else 'proteoform_id'
    groups = sorted(df[group_col].unique())

    for gi, grp in enumerate(groups):
        sub = df[df[group_col] == grp]
        col = palette[gi % len(palette)]
        fig.add_trace(go.Scatter(
            x=sub['mz'], y=sub['rt'],
            mode='markers',
            marker=dict(
                size=sub['size'],
                color=col,
                opacity=sub['alpha'].tolist(),
                line=dict(color='rgba(0,0,0,0.3)', width=0.8),
            ),
            name=f'z={grp}' if color_by == 'charge' else str(grp),
            customdata=sub[['feature_id', 'mass', 'intensity', 'charge']].values,
            hovertemplate=(
                '<b>%{customdata[0]}</b><br>'
                'm/z: %{x:.4f}<br>'
                'RT: %{y:.2f} min<br>'
                'Mass: %{customdata[1]:.3f} Da<br>'
                'z: %{customdata[3]}<br>'
                'Intensity: %{customdata[2]:.2e}<extra></extra>'
            ),
        ))

    # Overlay theoretical m/z lines for input proteoforms
    if proteoforms:
        for pf in proteoforms:
            if pf.theoretical_mass <= 0:
                continue
            for z in possible_charges(pf.theoretical_mass,
                                       min_mz=df['mz'].min() - 20,
                                       max_mz=df['mz'].max() + 20):
                th_mz = mass_to_mz(pf.theoretical_mass, z)
                fig.add_vline(
                    x=th_mz, line_dash='dash',
                    line_color='rgba(255,165,0,0.5)',
                    annotation_text=f"{pf.protein_name or 'Target'} z={z}",
                    annotation_font_size=9,
                    annotation_font_color='#E65100',
                )

    # Highlight selected feature with ring
    if selected_id:
        sel = df[df['feature_id'] == selected_id]
        if not sel.empty:
            fig.add_trace(go.Scatter(
                x=sel['mz'], y=sel['rt'],
                mode='markers',
                marker=dict(symbol='circle-open', size=sel['size'] + 8,
                            color='white', line=dict(width=2, color='white')),
                showlegend=False, hoverinfo='skip',
                name='selected',
            ))

    fig.update_layout(
        template='plotly_white',
        title=dict(text='<b>Feature Map</b> — m/z vs Retention Time',
                   font=dict(size=12)),
        xaxis=dict(title='m/z', showgrid=True, gridcolor='rgba(0,0,0,0.08)', automargin=True),
        yaxis=dict(title='Retention Time (min)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.08)', automargin=True),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        hovermode='closest',
        legend=dict(
            title=dict(text='Charge' if color_by == 'charge' else 'Proteoform',
                       font=dict(size=10)),
            font=dict(size=10),
            orientation='h', yanchor='bottom', y=1.02, x=0,
            bgcolor='rgba(255,255,255,0.7)',
        ),
        height=560,
        margin=dict(l=65, r=20, t=90, b=55),
    )
    return fig


def create_intensity_trace(features: List[Feature],
                             selected_id: Optional[str] = None,
                             theoretical_mass: float = 0.0,
                             theoretical_charge: int = 0) -> go.Figure:
    """
    Simulated Gaussian elution trace for selected / all features.
    Overlays theoretical expected m/z (as annotation).
    """
    fig = go.Figure()

    show_feats = [f for f in features
                  if selected_id is None or f.feature_id == selected_id]

    if not show_feats:
        fig.update_layout(template='plotly_white', title='No feature selected',
                          paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
                          font=dict(color='#111111'))
        return fig

    palette = px.colors.qualitative.Plotly
    for gi, feat in enumerate(show_feats[:10]):
        rt_pts = np.linspace(feat.rt_start, feat.rt_end, 80)
        sigma  = (feat.rt_end - feat.rt_start) / 4.0
        sigma  = max(sigma, 0.05)
        trace  = feat.intensity * np.exp(-0.5 * ((rt_pts - feat.rt_apex) / sigma) ** 2)
        col    = palette[gi % len(palette)]

        fig.add_trace(go.Scatter(
            x=rt_pts, y=trace,
            mode='lines',
            name=f"{feat.feature_id} (z={feat.charge}, m/z={feat.mz_apex:.3f})",
            line=dict(color=col, width=2),
            fill='tozeroy', fillcolor=col.replace(')', ',0.10)').replace('rgb', 'rgba'),
        ))
        fig.add_vline(x=feat.rt_apex, line_dash='dot', line_color=col,
                      annotation_text=f"apex {feat.rt_apex:.2f}", annotation_font_size=9)

    if theoretical_mass > 0 and theoretical_charge > 0:
        th_mz = mass_to_mz(theoretical_mass, theoretical_charge)
        fig.add_annotation(
            xref='paper', yref='paper', x=0.98, y=0.98,
            text=f"Theoretical: {th_mz:.4f} m/z (z={theoretical_charge})",
            showarrow=False, font=dict(color='#333333', size=11),
            xanchor='right', bgcolor='rgba(255,255,255,0.85)'
        )

    fig.update_layout(
        template='plotly_white',
        title=dict(text='<b>Elution Trace</b> — Intensity vs RT', font=dict(size=12)),
        xaxis=dict(title='Retention Time (min)', automargin=True),
        yaxis=dict(title='Intensity', automargin=True),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        height=320,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                    font=dict(size=9), bgcolor='rgba(255,255,255,0.7)'),
        margin=dict(l=65, r=20, t=80, b=45),
    )
    return fig


def create_feature_3d_plot(features: List[Feature]) -> go.Figure:
    """
    3-D scatter: Retention Time (x) × Charge (y) × log₁₀(Intensity) (z).
    Each proteoform identity gets its own colour.
    """
    fig = go.Figure()
    if not features:
        fig.update_layout(
            title='No features loaded',
            paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
            font=dict(color='#111111'),
            height=520,
        )
        return fig

    palette = px.colors.qualitative.Plotly
    prot_ids = sorted(set(f.proteoform_id or 'Unknown' for f in features))

    for gi, prot in enumerate(prot_ids):
        sub = [f for f in features if (f.proteoform_id or 'Unknown') == prot]
        col = palette[gi % len(palette)]
        fig.add_trace(go.Scatter3d(
            x=[f.rt_apex        for f in sub],
            y=[f.charge         for f in sub],
            z=[np.log10(max(f.intensity, 1)) for f in sub],
            mode='markers',
            marker=dict(size=5, color=col, opacity=0.85,
                        line=dict(color='rgba(0,0,0,0.2)', width=0.5)),
            name=str(prot),
            customdata=[[f.feature_id, f.monoisotopic_mass, f.intensity] for f in sub],
            hovertemplate=(
                '<b>%{customdata[0]}</b><br>'
                'RT: %{x:.2f} min<br>'
                'Charge: %{y}<br>'
                'log₁₀(I): %{z:.2f}<extra></extra>'
            ),
        ))

    fig.update_layout(
        title='Feature 3D Plot — RT × Charge × Intensity',
        scene=dict(
            xaxis=dict(title='RT (min)', backgroundcolor='#f8f8f8',
                       gridcolor='#cccccc', showbackground=True),
            yaxis=dict(title='Charge', backgroundcolor='#f8f8f8',
                       gridcolor='#cccccc', showbackground=True),
            zaxis=dict(title='log₁₀(I)', backgroundcolor='#f8f8f8',
                       gridcolor='#cccccc', showbackground=True),
        ),
        paper_bgcolor=DARK_BG,
        font=dict(color='#111111'),
        height=580,
        legend=dict(font=dict(size=10)),
        margin=dict(l=0, r=0, t=50, b=0),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for the comparison panel
# ─────────────────────────────────────────────────────────────────────────────

def _make_sticks(mzs, heights):
    """(x, y) arrays for vertical stick traces, None-separated."""
    xs, ys = [], []
    for m, h in zip(mzs, heights):
        xs += [float(m), float(m), None]
        ys += [0.0, float(h), None]
    return xs, ys


def _fwhm_from_trace(rt: np.ndarray, intensity: np.ndarray) -> float:
    """Estimate FWHM of an elution peak via half-max crossing."""
    if len(intensity) < 3 or intensity.max() == 0:
        return 0.0
    half_max = intensity.max() / 2.0
    above = rt[intensity >= half_max]
    return float(above[-1] - above[0]) if len(above) >= 2 else 0.0


def create_comparison_panel(
    feature: Optional[Feature],
    spectra: list,
    theoretical_mass: float = 0.0,
    ppm_tolerance: float = 10.0,
) -> go.Figure:
    """
    Two-row comparison panel for a selected feature.

    Row 1 — Isotope Envelope (m/z domain)
        Blue sticks (↑): theoretical averagine isotope pattern
        Red sticks  (↓): observed peaks from the nearest MS1 scan, mirrored

    Row 2 — Elution Profile (RT domain)
        Blue dotted line: theoretical symmetric Gaussian (from feature RT bounds)
        Orange solid line: observed XIC extracted from MS1 spectra, normalised 0–100 %

    Annotated metrics:
        Row 1 — mass accuracy (ppm), isotopes matched, avg isotope Δm/z
        Row 2 — Pearson R², FWHM observed vs theoretical, apex intensity
    """
    _empty = dict(
        template='plotly_white',
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#aaaaaa'),
        height=600,
        margin=dict(l=55, r=20, t=55, b=40),
    )

    if feature is None:
        fig = go.Figure()
        fig.update_layout(
            title='Select a feature on the Feature Map to see the comparison panel',
            **_empty,
        )
        return fig

    mass   = theoretical_mass if theoretical_mass > 0 else feature.monoisotopic_mass
    charge = feature.charge if feature.charge > 0 else 1

    if mass <= 0:
        fig = go.Figure()
        fig.update_layout(
            title='No theoretical mass — enter a value in the "Theoretical mass" field above',
            **_empty,
        )
        return fig

    # ── Theoretical isotope envelope ─────────────────────────────────────────
    n_iso = min(10, max(4, int(mass / 3000) + 4))
    th_mzs, th_ints = generate_isotope_envelope(mass, charge, n_isotopes=n_iso)

    # ── Nearest MS1 scan at feature RT apex ──────────────────────────────────
    ms1 = [sd for sd in spectra
           if sd.get('ms_level', 2) == 1 and len(sd.get('mz', [])) > 0]

    obs_mzs_env = np.array([], dtype=np.float64)
    obs_ints_env = np.array([], dtype=np.float64)
    scan_rt_used: Optional[float] = None

    if ms1:
        scan_rts = np.array([sd.get('retention_time', 0.0) for sd in ms1])
        ni = int(np.argmin(np.abs(scan_rts - feature.rt_apex)))
        scan_rt_used = float(scan_rts[ni])
        s = ms1[ni]
        raw_mz  = np.asarray(s.get('mz',        []), dtype=np.float64)
        raw_int = np.asarray(s.get('intensity',  []), dtype=np.float64)
        # envelope window: monoisotopic to last isotope, with half-spacing margin
        half_sp  = 0.5 * (1.003355 / charge)
        env_lo   = th_mzs[0]  - half_sp
        env_hi   = th_mzs[-1] + half_sp
        mask     = (raw_mz >= env_lo) & (raw_mz <= env_hi)
        obs_mzs_env  = raw_mz[mask]
        obs_ints_env = raw_int[mask]

    # Normalise observed to 0–100
    obs_norm_env = np.zeros_like(obs_ints_env)
    if obs_ints_env.size > 0 and obs_ints_env.max() > 0:
        obs_norm_env = obs_ints_env / obs_ints_env.max() * 100.0

    # ── Isotope matching metrics ──────────────────────────────────────────────
    matched_count = 0
    ppm_errs: list = []
    mono_obs_mz: Optional[float] = None

    for i, th_mz in enumerate(th_mzs):
        if obs_mzs_env.size == 0:
            break
        tol_da = th_mz * ppm_tolerance / 1e6
        dists   = np.abs(obs_mzs_env - th_mz)
        best    = int(np.argmin(dists))
        if dists[best] <= tol_da:
            matched_count += 1
            ppm_errs.append((obs_mzs_env[best] - th_mz) / th_mz * 1e6)
            if i == 0:
                mono_obs_mz = float(obs_mzs_env[best])

    avg_iso_ppm = float(np.mean(np.abs(ppm_errs))) if ppm_errs else float('nan')

    # Mass accuracy from observed monoisotopic peak (or from stored feature mass)
    if mono_obs_mz is not None:
        obs_mass = mono_obs_mz * charge - charge * PROTON
        mass_acc_ppm = (obs_mass - mass) / mass * 1e6
    elif feature.monoisotopic_mass > 0 and theoretical_mass > 0:
        mass_acc_ppm = (feature.monoisotopic_mass - theoretical_mass) / theoretical_mass * 1e6
    else:
        mass_acc_ppm = float('nan')

    # ── XIC extraction ────────────────────────────────────────────────────────
    if feature.mz_start > 0 and feature.mz_end > feature.mz_start:
        xic_mz_lo, xic_mz_hi = feature.mz_start, feature.mz_end
    else:
        half_sp    = 0.5 * (1.003355 / charge)
        xic_mz_lo  = th_mzs[0]  - half_sp
        xic_mz_hi  = th_mzs[-1] + half_sp

    xic_rts_l, xic_vals_l = [], []
    for sd in spectra:
        if sd.get('ms_level', 2) != 1:
            continue
        rt   = sd.get('retention_time', 0.0)
        mzs  = np.asarray(sd.get('mz',        []), dtype=np.float64)
        ints = np.asarray(sd.get('intensity',  []), dtype=np.float64)
        if not mzs.size:
            continue
        mask = (mzs >= xic_mz_lo) & (mzs <= xic_mz_hi)
        xic_rts_l.append(rt)
        xic_vals_l.append(float(ints[mask].sum()))

    xic_rts  = np.array(xic_rts_l)
    xic_vals = np.array(xic_vals_l)
    if xic_rts.size:
        order    = np.argsort(xic_rts)
        xic_rts  = xic_rts[order]
        xic_vals = xic_vals[order]

    # ── Theoretical Gaussian ──────────────────────────────────────────────────
    rt_start = feature.rt_start if feature.rt_start > 0 else feature.rt_apex - 1.0
    rt_end   = feature.rt_end   if feature.rt_end   > 0 else feature.rt_apex + 1.0
    sigma    = max((rt_end - rt_start) / 4.0, 0.05)
    th_fwhm  = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma

    th_rts   = np.linspace(max(rt_start - 0.5, 0.0), rt_end + 0.5, 300)
    th_trace = np.exp(-0.5 * ((th_rts - feature.rt_apex) / sigma) ** 2) * 100.0

    # Normalise observed XIC to 0–100 over the feature window
    xic_norm = np.zeros_like(xic_vals)
    if xic_vals.size:
        win_mask = (xic_rts >= rt_start - 0.2) & (xic_rts <= rt_end + 0.2)
        win_max  = xic_vals[win_mask].max() if win_mask.any() else xic_vals.max()
        if win_max > 0:
            xic_norm = np.clip(xic_vals / win_max * 100.0, 0, 110)

    # ── Elution shape metrics (R², FWHM) ─────────────────────────────────────
    r2: Optional[float] = None
    obs_fwhm: Optional[float] = None

    if xic_rts.size >= 4:
        cmp_lo = max(rt_start - 1.0, float(xic_rts[0]))
        cmp_hi = min(rt_end   + 1.0, float(xic_rts[-1]))
        th_cmp = np.linspace(cmp_lo, cmp_hi, 100)
        th_ref = np.exp(-0.5 * ((th_cmp - feature.rt_apex) / sigma) ** 2)

        f_obs  = interp1d(xic_rts, xic_norm, bounds_error=False, fill_value=0.0)
        obs_rs = f_obs(th_cmp)
        if obs_rs.max() > 0:
            obs_rs = obs_rs / obs_rs.max()

        corr = np.corrcoef(th_ref, obs_rs)[0, 1]
        r2   = float(corr ** 2) if np.isfinite(corr) else None

        win_mask_fwhm = (xic_rts >= rt_start - 0.2) & (xic_rts <= rt_end + 0.2)
        if win_mask_fwhm.sum() >= 3:
            obs_fwhm = _fwhm_from_trace(xic_rts[win_mask_fwhm], xic_vals[win_mask_fwhm])

    # ── Build figure ──────────────────────────────────────────────────────────
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=False,
        subplot_titles=[
            f'Isotope Envelope  —  theoretical m/z {mass_to_mz(mass, charge):.4f}  z={charge}',
            'Elution Profile  —  Normalised Intensity vs. Retention Time',
        ],
        row_heights=[0.42, 0.58],
        vertical_spacing=0.14,
    )

    # ── Row 1: isotope mirror ─────────────────────────────────────────────────
    th_sx, th_sy = _make_sticks(th_mzs, th_ints)
    fig.add_trace(go.Scatter(
        x=th_sx, y=th_sy, mode='lines',
        line=dict(color='#1a73e8', width=2.5),
        name='Theoretical', legendgroup='iso',
        hoverinfo='skip',
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=th_mzs, y=th_ints, mode='markers',
        marker=dict(color='#1a73e8', size=7, symbol='circle'),
        showlegend=False, legendgroup='iso',
        hovertemplate='m/z: %{x:.4f}<br>Rel. Int: %{y:.1f} %<extra>Theoretical</extra>',
    ), row=1, col=1)

    if obs_mzs_env.size > 0:
        ob_sx, ob_sy = _make_sticks(obs_mzs_env, -obs_norm_env)
        rt_label = f'Observed (RT {scan_rt_used:.2f} min)' if scan_rt_used else 'Observed'
        fig.add_trace(go.Scatter(
            x=ob_sx, y=ob_sy, mode='lines',
            line=dict(color='#d32f2f', width=2.0),
            name=rt_label, legendgroup='iso',
            hoverinfo='skip',
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=obs_mzs_env, y=-obs_norm_env, mode='markers',
            marker=dict(color='#d32f2f', size=5, symbol='circle'),
            showlegend=False, legendgroup='iso',
            customdata=obs_norm_env,
            hovertemplate='m/z: %{x:.4f}<br>Rel. Int: %{customdata:.1f} %<extra>Observed</extra>',
        ), row=1, col=1)

    # Zero-line separator
    mz_x0 = float(th_mzs[0]) - 0.3
    mz_x1 = float(th_mzs[-1]) + 0.3
    fig.add_trace(go.Scatter(
        x=[mz_x0, mz_x1], y=[0.0, 0.0], mode='lines',
        line=dict(color='rgba(0,0,0,0.25)', width=1),
        showlegend=False, hoverinfo='skip',
    ), row=1, col=1)

    # ── Row 2: elution profile ────────────────────────────────────────────────
    fig.add_trace(go.Scatter(
        x=th_rts, y=th_trace, mode='lines',
        line=dict(color='#1a73e8', width=2.5, dash='dot'),
        name='Theoretical (Gaussian)', legendgroup='elut',
        fill='tozeroy', fillcolor='rgba(26,115,232,0.08)',
        hovertemplate='RT: %{x:.3f} min<br>Norm. Int: %{y:.1f} %<extra>Theoretical</extra>',
    ), row=2, col=1)

    if xic_rts.size > 0:
        fig.add_trace(go.Scatter(
            x=xic_rts, y=xic_norm, mode='lines',
            line=dict(color='#e65100', width=2.0),
            name='Observed (XIC)', legendgroup='elut',
            fill='tozeroy', fillcolor='rgba(230,81,0,0.08)',
            hovertemplate='RT: %{x:.3f} min<br>Norm. Int: %{y:.1f} %<extra>Observed XIC</extra>',
        ), row=2, col=1)

    # RT apex vertical line
    fig.add_trace(go.Scatter(
        x=[feature.rt_apex, feature.rt_apex], y=[0, 105],
        mode='lines',
        line=dict(color='rgba(0,0,0,0.30)', width=1, dash='dash'),
        showlegend=False, hoverinfo='skip',
    ), row=2, col=1)

    # ── Annotations: isotope metrics (row 1) ─────────────────────────────────
    iso_lines = []
    if np.isfinite(mass_acc_ppm):
        col_m = '#1b8f1b' if abs(mass_acc_ppm) <= 5 else '#d32f2f'
        iso_lines.append(
            f"<span style='color:{col_m}'><b>Mass accuracy:</b> {mass_acc_ppm:+.2f} ppm</span>"
        )
    iso_lines.append(f"<b>Isotopes matched:</b> {matched_count} / {n_iso}")
    if ppm_errs:
        iso_lines.append(f"<b>Avg Δm/z:</b> {avg_iso_ppm:.2f} ppm")

    fig.add_annotation(
        xref='x domain', yref='y domain',
        x=0.98, y=0.97,
        text='<br>'.join(iso_lines),
        showarrow=False,
        font=dict(size=10, color='#333333'),
        xanchor='right', yanchor='top',
        bgcolor='rgba(255,255,255,0.88)',
        bordercolor='rgba(0,0,0,0.15)',
        borderwidth=1, borderpad=5,
        align='left',
        row=1, col=1,
    )

    # ── Annotations: elution metrics (row 2) ─────────────────────────────────
    elu_lines = []
    if r2 is not None:
        r2_col = '#1b8f1b' if r2 >= 0.90 else ('#e65100' if r2 >= 0.70 else '#d32f2f')
        elu_lines.append(
            f"<span style='color:{r2_col}'><b>Shape R²:</b> {r2:.3f}</span>"
        )
    elu_lines.append(f"<b>FWHM (theoretical):</b> {th_fwhm:.3f} min")
    if obs_fwhm is not None and obs_fwhm > 0:
        elu_lines.append(f"<b>FWHM (observed):</b>  {obs_fwhm:.3f} min")
    if feature.intensity > 0:
        elu_lines.append(f"<b>Apex intensity:</b>   {feature.intensity:.3e}")

    fig.add_annotation(
        xref='x2 domain', yref='y2 domain',
        x=0.98, y=0.97,
        text='<br>'.join(elu_lines),
        showarrow=False,
        font=dict(size=10, color='#333333'),
        xanchor='right', yanchor='top',
        bgcolor='rgba(255,255,255,0.88)',
        bordercolor='rgba(0,0,0,0.15)',
        borderwidth=1, borderpad=5,
        align='left',
    )

    # ── Layout ────────────────────────────────────────────────────────────────
    fig.update_layout(
        template='plotly_white',
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111', size=11),
        height=600,
        margin=dict(l=65, r=20, t=90, b=50),
        hovermode='closest',
        legend=dict(
            orientation='h', yanchor='bottom', y=1.04, x=0,
            font=dict(size=9), bgcolor='rgba(255,255,255,0.7)',
        ),
        title=dict(
            text=(
                f'<b>Theoretical vs. Observed</b> — {feature.feature_id}<br>'
                f'<sup>Theoretical mass: {mass:.4f} Da  |  z = {charge}  |  '
                f'Observed (feature): {feature.monoisotopic_mass:.4f} Da</sup>'
            ),
            font=dict(size=12),
        ),
    )
    fig.update_yaxes(title_text='Rel. Intensity (%)', tickformat='.0f', row=1, col=1)
    fig.update_xaxes(title_text='m/z', row=1, col=1)
    fig.update_yaxes(title_text='Norm. Intensity (%)', range=[-5, 118], row=2, col=1)
    fig.update_xaxes(title_text='Retention Time (min)', row=2, col=1)

    return fig


def create_xic_plot(spectra: list, feature: Optional[Feature] = None,
                    mz_tolerance_ppm: float = 10.0) -> go.Figure:
    """
    Extracted Ion Chromatogram (XIC) computed from real MS1 spectra.

    For every MS1 scan in *spectra*, sums the intensity of all peaks that
    fall within the feature's m/z window (mz_start..mz_end, or
    mz_apex ± mz_tolerance_ppm if the window is unknown).

    Parameters
    ----------
    spectra : list of Spectrum.to_dict() dicts (from store-spectra)
    feature : selected Feature whose m/z window defines the extraction
    mz_tolerance_ppm : fallback ppm window when mz_start/mz_end are 0
    """
    fig = go.Figure()
    empty_layout = dict(
        template='plotly_white',
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#aaaaaa'),
        height=280,
        margin=dict(l=55, r=20, t=55, b=40),
    )

    if not spectra:
        fig.update_layout(title='Load a spectrum file to see XIC', **empty_layout)
        return fig

    if feature is None:
        fig.update_layout(title='Click a feature on the map to see its XIC', **empty_layout)
        return fig

    # Determine m/z extraction window
    mz_lo = feature.mz_start
    mz_hi = feature.mz_end
    if mz_lo <= 0 or mz_hi <= 0 or mz_hi <= mz_lo:
        tol = feature.mz_apex * mz_tolerance_ppm / 1e6
        mz_lo = feature.mz_apex - tol
        mz_hi = feature.mz_apex + tol

    # Extract XIC across MS1 scans
    rts, xic, all_rts, all_xic = [], [], [], []
    for sd in spectra:
        if sd.get('ms_level', 2) != 1:
            continue
        rt  = sd.get('retention_time', 0.0)
        mzs = np.asarray(sd.get('mz', []), dtype=np.float64)
        ints = np.asarray(sd.get('intensity', []), dtype=np.float64)
        if len(mzs) == 0:
            continue
        mask = (mzs >= mz_lo) & (mzs <= mz_hi)
        xic_val = float(ints[mask].sum())
        all_rts.append(rt)
        all_xic.append(xic_val)

    if not all_rts:
        fig.update_layout(title='No MS1 scans in loaded data — XIC unavailable',
                          **empty_layout)
        return fig

    all_rts  = np.array(all_rts)
    all_xic  = np.array(all_xic)
    order    = np.argsort(all_rts)
    all_rts  = all_rts[order]
    all_xic  = all_xic[order]

    # Background TIC-like trace (faint)
    fig.add_trace(go.Scatter(
        x=all_rts, y=all_xic,
        mode='lines',
        line=dict(color='rgba(0,0,0,0.12)', width=1),
        fill='tozeroy',
        fillcolor='rgba(0,0,0,0.04)',
        name='XIC (full RT)',
        showlegend=True,
        hovertemplate='RT: %{x:.3f} min<br>Intensity: %{y:.3e}<extra></extra>',
    ))

    # Highlight the feature's RT window
    rt_lo = feature.rt_start if feature.rt_start > 0 else (feature.rt_apex - 0.5)
    rt_hi = feature.rt_end   if feature.rt_end   > 0 else (feature.rt_apex + 0.5)
    in_window = (all_rts >= rt_lo) & (all_rts <= rt_hi)
    if in_window.any():
        fig.add_trace(go.Scatter(
            x=all_rts[in_window], y=all_xic[in_window],
            mode='lines',
            line=dict(color='#1a73e8', width=2.5),
            fill='tozeroy',
            fillcolor='rgba(26,115,232,0.18)',
            name=f"{feature.feature_id} window",
            hovertemplate='RT: %{x:.3f} min<br>Intensity: %{y:.3e}<extra></extra>',
        ))

    # Apex vertical line
    fig.add_vline(
        x=feature.rt_apex,
        line_dash='dash', line_color='#1a73e8', line_width=1.5,
        annotation_text=f"apex {feature.rt_apex:.2f} min",
        annotation_font_size=9, annotation_font_color='#1a73e8',
    )

    # Feature RT window shaded region
    fig.add_vrect(
        x0=rt_lo, x1=rt_hi,
        fillcolor='rgba(26,115,232,0.06)',
        layer='below', line_width=0,
    )

    # Peak intensity annotation
    if in_window.any():
        peak_rt  = float(all_rts[in_window][np.argmax(all_xic[in_window])])
        peak_int = float(all_xic[in_window].max())
        fig.add_annotation(
            x=peak_rt, y=peak_int,
            text=f"{peak_int:.2e}",
            showarrow=True, arrowhead=2, arrowsize=0.8,
            arrowcolor='#1a73e8', font=dict(size=9, color='#1a73e8'),
            ay=-28, ax=0,
        )

    fig.update_layout(
        template='plotly_white',
        title=dict(
            text=(f'<b>XIC</b> — {feature.feature_id}<br>'
                  f'<sup>m/z {mz_lo:.3f}–{mz_hi:.3f}  |  z={feature.charge}</sup>'),
            font=dict(size=12),
        ),
        xaxis=dict(title='Retention Time (min)', automargin=True),
        yaxis=dict(title='Summed Intensity', automargin=True),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        height=300,
        legend=dict(orientation='h', yanchor='bottom', y=1.02,
                    font=dict(size=10), bgcolor='rgba(255,255,255,0.7)'),
        margin=dict(l=65, r=20, t=85, b=45),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Mass Accuracy Plot — ppm error per feature vs m/z (or index)
# ─────────────────────────────────────────────────────────────────────────────

def create_mass_accuracy_plot(features: List[Feature],
                               theoretical_mass: float = 0.0,
                               ppm_tolerance: float = 10.0) -> go.Figure:
    """
    Scatter of mass accuracy (ppm) for loaded features.

    x-axis  → m/z apex of each feature
    y-axis  → mass accuracy in ppm (observed_mass − theoretical) / theoretical × 10⁶
              If theoretical_mass is not provided, uses the first feature's mass
              as the reference.
    Colour  → charge state
    Size    → log₁₀(intensity)
    Reference lines at ±ppm_tolerance (dashed grey) and ±5 ppm (solid blue).
    """
    fig = go.Figure()
    _empty = dict(
        template='plotly_white',
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#aaaaaa'),
        height=320,
        margin=dict(l=60, r=20, t=55, b=45),
    )
    if not features:
        fig.update_layout(title='No features loaded — mass accuracy unavailable', **_empty)
        return fig

    # Determine reference mass
    ref_mass = theoretical_mass
    if ref_mass <= 0:
        candidates = [f.monoisotopic_mass for f in features if f.monoisotopic_mass > 0]
        ref_mass = float(np.median(candidates)) if candidates else 0.0

    if ref_mass <= 0:
        fig.update_layout(title='Enter a theoretical mass to compute mass accuracy', **_empty)
        return fig

    palette = px.colors.qualitative.Plotly
    charges = sorted(set(f.charge for f in features))
    charge_color = {z: palette[i % len(palette)] for i, z in enumerate(charges)}

    for z in charges:
        sub = [f for f in features if f.charge == z and f.monoisotopic_mass > 0]
        if not sub:
            continue
        xs   = [f.mz_apex for f in sub]
        ppms = [(f.monoisotopic_mass - ref_mass) / ref_mass * 1e6 for f in sub]
        sizes = [max(5, np.log10(max(f.intensity, 1)) * 3) for f in sub]
        fids  = [f.feature_id for f in sub]
        masses = [f.monoisotopic_mass for f in sub]

        fig.add_trace(go.Scatter(
            x=xs, y=ppms,
            mode='markers',
            marker=dict(
                size=sizes,
                color=charge_color[z],
                opacity=0.80,
                line=dict(color='rgba(0,0,0,0.2)', width=0.7),
            ),
            name=f'z={z}',
            customdata=list(zip(fids, masses, ppms)),
            hovertemplate=(
                '<b>%{customdata[0]}</b><br>'
                'm/z: %{x:.4f}<br>'
                'Mass: %{customdata[1]:.4f} Da<br>'
                'Error: %{customdata[2]:+.2f} ppm<extra></extra>'
            ),
        ))

    # Reference lines
    for ppm_val, style, colour, label in [
        (0,             'solid', 'rgba(0,0,0,0.35)',  '0 ppm'),
        (5,             'dash',  'rgba(26,115,232,0.5)', '+5 ppm'),
        (-5,            'dash',  'rgba(26,115,232,0.5)', '-5 ppm'),
        (ppm_tolerance, 'dot',   'rgba(180,0,0,0.4)',   f'+{ppm_tolerance:.0f} ppm'),
        (-ppm_tolerance,'dot',   'rgba(180,0,0,0.4)',   f'−{ppm_tolerance:.0f} ppm'),
    ]:
        fig.add_hline(y=ppm_val, line_dash=style, line_color=colour,
                      annotation_text=label,
                      annotation_font_size=9,
                      annotation_font_color=colour,
                      annotation_position='right')

    fig.update_layout(
        template='plotly_white',
        title=dict(text=f'<b>Mass Accuracy</b> — ref {ref_mass:.4f} Da',
                   font=dict(size=12)),
        xaxis=dict(title='m/z', showgrid=True, gridcolor='rgba(0,0,0,0.07)',
                   automargin=True),
        yaxis=dict(title='Mass Error (ppm)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.07)', zeroline=True,
                   zerolinewidth=1.5, zerolinecolor='rgba(0,0,0,0.3)',
                   automargin=True),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        hovermode='closest',
        legend=dict(
            title=dict(text='Charge', font=dict(size=10)),
            font=dict(size=10),
            orientation='h', yanchor='bottom', y=1.02, x=0,
            bgcolor='rgba(255,255,255,0.7)',
        ),
        height=320,
        margin=dict(l=65, r=80, t=70, b=50),
    )
    return fig
