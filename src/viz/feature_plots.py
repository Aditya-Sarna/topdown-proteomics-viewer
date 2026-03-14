"""
Feature map visualizations:
  1. 2-D scatter: m/z vs retention time (dot area ∝ intensity)
  2. Intensity trace: Gaussian-approximated elution profile for selected feature
"""
from typing import List, Optional
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from ..data.models import Feature, Proteoform
from ..analysis.mass_utils import mass_to_mz, possible_charges

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
        height=520,
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
        height=520,
        legend=dict(font=dict(size=10)),
        margin=dict(l=0, r=0, t=50, b=0),
    )
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
