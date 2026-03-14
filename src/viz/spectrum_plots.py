"""
Spectrum viewer plot — efficient vertical-line rendering with per-ion-type colouring.
"""
from typing import List, Optional
import numpy as np
import plotly.graph_objects as go

from ..data.models import Spectrum, FragmentIon

ION_COLORS = {
    'b': '#2196F3',
    'y': '#F44336',
    'c': '#4CAF50',
    'z': '#FF9800',
    'a': '#9C27B0',
}
UNMATCHED_COLOR = '#aaaaaa'
DARK_BG = '#ffffff'
PLOT_BG = '#ffffff'


def _lines(mz_arr, int_arr):
    """Convert arrays to None-separated x/y for vertical line rendering."""
    xs, ys = [], []
    for m, i in zip(mz_arr, int_arr):
        xs += [m, m, None]
        ys += [0, i, None]
    return xs, ys


def create_spectrum_plot(spectrum: Spectrum,
                          matched_ions: Optional[List[FragmentIon]] = None,
                          normalize: bool = True,
                          annotate_top_n: int = 30) -> go.Figure:
    """
    Interactive spectrum viewer.
    - Unmatched peaks: grey sticks (one trace).
    - Matched peaks: one coloured stick trace per ion type.
    - Top-N matched peaks receive text labels.
    """
    fig = go.Figure()

    mz  = spectrum.mz_array
    raw = spectrum.intensity_array.copy()

    if len(mz) == 0:
        fig.update_layout(template='plotly_white',
                          title='No spectrum loaded',
                          xaxis_title='m/z', yaxis_title='Rel. Intensity (%)',
                          paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG)
        return fig

    intensity = raw / raw.max() * 100.0 if normalize and raw.max() > 0 else raw

    # Build lookup: observed_mz (rounded) -> FragmentIon
    match_lookup: dict = {}
    if matched_ions:
        for ion in matched_ions:
            if ion.matched:
                key = round(ion.observed_mz, 2)
                if key not in match_lookup or ion.charge < match_lookup[key].charge:
                    match_lookup[key] = ion

    # Separate peaks by ion type
    by_type: dict = {t: {'mz': [], 'int': [], 'labels': [], 'hov': []} for t in ION_COLORS}
    unk = {'mz': [], 'int': [], 'hov': []}

    for m, i in zip(mz, intensity):
        key = round(m, 2)
        if key in match_lookup:
            ion  = match_lookup[key]
            itype = ion.ion_type
            lbl  = ion.label()
            hov  = (f"<b>{lbl}</b><br>m/z: {m:.4f}<br>"
                    f"Intensity: {i:.1f}%<br>"
                    f"Δ: {ion.mass_error_ppm:.1f} ppm")
            by_type[itype]['mz'].append(m)
            by_type[itype]['int'].append(i)
            by_type[itype]['labels'].append(lbl)
            by_type[itype]['hov'].append(hov)
        else:
            unk['mz'].append(m)
            unk['int'].append(i)
            unk['hov'].append(f"m/z: {m:.4f}<br>Intensity: {i:.1f}%")

    # --- Unmatched peaks ---
    if unk['mz']:
        xs, ys = _lines(unk['mz'], unk['int'])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode='lines',
            line=dict(color=UNMATCHED_COLOR, width=1.2),
            name='Unmatched', showlegend=True,
            hoverinfo='skip',
        ))
        fig.add_trace(go.Scatter(
            x=unk['mz'], y=unk['int'], mode='markers',
            marker=dict(size=4, color=UNMATCHED_COLOR, opacity=0.0),
            hovertext=unk['hov'], hoverinfo='text',
            showlegend=False, name='Unmatched peaks',
        ))

    # --- Matched peaks per ion type ---
    for itype, data in by_type.items():
        if not data['mz']:
            continue
        color = ION_COLORS[itype]
        xs, ys = _lines(data['mz'], data['int'])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode='lines',
            line=dict(color=color, width=2),
            name=f"{itype.upper()}-ions", legendgroup=itype, showlegend=True,
            hoverinfo='skip',
        ))
        fig.add_trace(go.Scatter(
            x=data['mz'], y=data['int'], mode='markers',
            marker=dict(size=5, color=color, opacity=0.0),
            hovertext=data['hov'], hoverinfo='text',
            showlegend=False, legendgroup=itype,
        ))

    if matched_ions:
        all_matched = [(ion.observed_mz, intensity[np.argmin(np.abs(mz - ion.observed_mz))], ion.label())
                       for ion in matched_ions if ion.matched]
        all_matched.sort(key=lambda x: -x[1])
        mz_span = (mz.max() - mz.min()) if len(mz) > 1 else 1000.0
        min_gap = mz_span * 0.012
        placed_mz = []
        annotations = []
        for obs_m, obs_i, lbl in all_matched[:annotate_top_n]:
            if any(abs(obs_m - pm) < min_gap for pm in placed_mz):
                continue
            placed_mz.append(obs_m)
            itype = lbl[0] if lbl[0] in ION_COLORS else 'b'
            annotations.append(dict(
                x=obs_m, y=obs_i + 2,
                text=f"<b>{lbl}</b>",
                showarrow=False, yshift=4,
                font=dict(size=9, color=ION_COLORS.get(itype, '#444444')),
                xanchor='center',
            ))
        fig.update_layout(annotations=annotations)

    prec_info = (f"Precursor: {spectrum.precursor_mz:.3f} m/z  "
                 f"z={spectrum.precursor_charge}  "
                 f"Mass={spectrum.precursor_mass:.2f} Da") if spectrum.precursor_mz else ""
    n_matched_total = sum(1 for ion in (matched_ions or []) if ion.matched)
    title = (f"<b>{spectrum.scan_id}</b>  |  RT {spectrum.retention_time:.2f} min  "
             f"|  {prec_info}  |  Matched: {n_matched_total}")

    fig.update_layout(
        template='plotly_white',
        title=dict(text=title, font=dict(size=12)),
        xaxis=dict(title='m/z', showgrid=True, gridcolor='rgba(0,0,0,0.08)'),
        yaxis_color='#111111', xaxis_color='#111111',
        yaxis=dict(title='Rel. Intensity (%)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.08)', range=[-3, 115]),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        hovermode='closest',
        legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0,
                    font=dict(size=11)),
        margin=dict(l=50, r=20, t=80, b=50),
    )
    return fig


def create_deconvolved_spectrum_plot(spectrum: Spectrum) -> go.Figure:
    """
    Plot the charge-deconvolved spectrum: intensity vs neutral monoisotopic mass.
    Runs pyopenms Deisotoper when available; falls back to the original spectrum.
    """
    from ..analysis.deconvolution import deconvolute_spectrum

    fig = go.Figure()
    dr   = deconvolute_spectrum(spectrum)
    spec = dr.spectrum
    masses = spec.mz_array
    ints   = spec.intensity_array

    if len(masses) == 0:
        fig.update_layout(title='No deconvolved peaks', template='plotly_white',
                          paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG)
        return fig

    norm = ints / ints.max() * 100.0 if ints.max() > 0 else ints
    xs, ys = _lines(masses, norm)
    note   = 'OpenMS' if dr.used_openms else 'raw (pyopenms unavailable)'

    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='lines',
        line=dict(color='#4CAF50', width=1.2),
        name=f'Deconvolved ({note})',
        hoverinfo='skip',
    ))

    # Annotate top 15 peaks by intensity
    top_idx = np.argsort(norm)[-15:]
    for i in top_idx:
        fig.add_annotation(
            x=masses[i], y=norm[i],
            text=f'{masses[i]:.1f}',
            showarrow=False, yshift=8,
            font=dict(size=8, color='#388E3C'),
        )

    fig.update_layout(
        template='plotly_white',
        title=dict(
            text=(f'Deconvolved Spectrum — {dr.n_original_peaks}→'
                  f'{dr.n_deconvoluted_peaks} peaks ({note})'),
            font=dict(size=12),
        ),
        xaxis=dict(title='Monoisotopic Mass (Da, z=1 equiv)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.08)'),
        yaxis=dict(title='Rel. Intensity (%)', showgrid=True,
                   gridcolor='rgba(0,0,0,0.08)', range=[-3, 115]),
        paper_bgcolor=DARK_BG, plot_bgcolor=PLOT_BG,
        font=dict(color='#111111'),
        showlegend=False,
        margin=dict(l=50, r=20, t=60, b=50),
        height=360,
    )
    return fig
