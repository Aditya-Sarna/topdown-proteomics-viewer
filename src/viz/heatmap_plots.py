"""
2-D heatmaps for raw and deconvolved MS data:
  - Raw heatmap:        m/z (y) × RT (x), log-intensity as colour
  - Deconvolved heatmap: neutral precursor mass (y) × RT (x)
"""
from typing import List
import numpy as np
import plotly.graph_objects as go

from ..data.models import Spectrum

DARK_BG = '#ffffff'
PLOT_BG = '#ffffff'
_LAYOUT_BASE = dict(
    paper_bgcolor=DARK_BG,
    plot_bgcolor=PLOT_BG,
    font=dict(color='#111111'),
    template='plotly_white',
    margin=dict(t=50, b=50, l=70, r=20),
)


def create_raw_heatmap(
        spectra: List[Spectrum],
        mz_min: float = 200.0,
        mz_max: float = 2000.0,
        n_mz_bins: int = 400,
        ms_level: int = 2,
) -> go.Figure:
    """
    Raw MS heatmap: RT (x) × m/z (y), log10(max intensity per cell) as colour.
    ms_level: 0=all, 1=MS1 only, 2=MS2 only
    """
    fig = go.Figure()
    if not spectra:
        fig.update_layout(title='Load spectra to view heatmap', height=450, **_LAYOUT_BASE)
        return fig

    pool = [s for s in spectra if ms_level == 0 or s.ms_level == ms_level]
    if not pool:
        pool = spectra

    rts = sorted(set(round(s.retention_time, 4) for s in pool))
    if len(rts) == 1:
        rts = [rts[0] - 0.001, rts[0], rts[0] + 0.001]
    rt_to_idx = {rt: i for i, rt in enumerate(rts)}
    n_rt = len(rts)

    mz_step = (mz_max - mz_min) / n_mz_bins
    matrix  = np.zeros((n_mz_bins, n_rt), dtype=np.float64)

    for spec in pool:
        rkey = round(spec.retention_time, 4)
        ri   = rt_to_idx.get(rkey, 0)
        mask = (spec.mz_array >= mz_min) & (spec.mz_array <= mz_max)
        for mz, inty in zip(spec.mz_array[mask], spec.intensity_array[mask]):
            mi = int((mz - mz_min) / mz_step)
            mi = min(mi, n_mz_bins - 1)
            if matrix[mi, ri] < inty:
                matrix[mi, ri] = inty

    log_mat    = np.log10(matrix + 1.0)
    mz_edges   = np.linspace(mz_min, mz_max, n_mz_bins + 1)
    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2.0

    fig.add_trace(go.Heatmap(
        z=log_mat,
        x=rts,
        y=mz_centers,
        colorscale='Viridis',
        colorbar=dict(title='log₁₀(I+1)', len=0.75, tickfont=dict(size=9)),
        hovertemplate='RT: %{x:.2f} min<br>m/z: %{y:.1f}<br>log₁₀(I): %{z:.2f}<extra></extra>',
    ))
    lvl = {0: 'All', 1: 'MS1', 2: 'MS2'}.get(ms_level, 'MS2')
    fig.update_layout(
        title=f'Raw MS Heatmap ({lvl}) — {len(pool)} scans',
        xaxis_title='Retention Time (min)',
        yaxis_title='m/z',
        height=450,
        **_LAYOUT_BASE,
    )
    return fig


def create_deconvolved_heatmap(spectra: List[Spectrum]) -> go.Figure:
    """
    Deconvolved heatmap: RT (x) × neutral monoisotopic mass (y).
    Uses precursor_mass from MS2 scans; estimates from mz×z for any scan that
    has precursor info but no explicit neutral mass.
    """
    PROTON = 1.007276
    fig = go.Figure()
    if not spectra:
        fig.update_layout(title='No spectra loaded', height=450, **_LAYOUT_BASE)
        return fig

    rts, masses, ints = [], [], []
    for spec in spectra:
        mass = spec.precursor_mass
        if mass <= 0 and spec.precursor_mz > 0 and spec.precursor_charge > 0:
            mass = spec.precursor_mz * spec.precursor_charge - spec.precursor_charge * PROTON
        if mass <= 0:
            continue
        peak_int = float(spec.intensity_array.max()) if len(spec.intensity_array) > 0 else 1.0
        rts.append(spec.retention_time)
        masses.append(mass)
        ints.append(peak_int)

    if not rts:
        fig.update_layout(
            title='No precursor mass data — load an mzML file with MS2 precursor information',
            height=450, **_LAYOUT_BASE,
        )
        return fig

    rts_arr    = np.array(rts)
    masses_arr = np.array(masses)
    ints_arr   = np.array(ints)

    mass_min  = max(100.0, masses_arr.min() * 0.95)
    mass_max  = masses_arr.max() * 1.05
    n_mass    = 300
    mass_step = (mass_max - mass_min) / n_mass

    rt_min, rt_max = rts_arr.min(), rts_arr.max()
    if rt_max == rt_min:
        rt_max = rt_min + 1.0
    n_rt    = min(200, max(10, len(rts)))
    rt_step = (rt_max - rt_min) / n_rt

    matrix = np.zeros((n_mass, n_rt), dtype=np.float64)
    for rt, mass, inty in zip(rts_arr, masses_arr, ints_arr):
        ri = int((rt   - rt_min)   / rt_step);   ri = min(ri, n_rt   - 1)
        mi = int((mass - mass_min) / mass_step);  mi = min(mi, n_mass - 1)
        if mi >= 0 and matrix[mi, ri] < inty:
            matrix[mi, ri] = inty

    log_mat      = np.log10(matrix + 1.0)
    mass_edges   = np.linspace(mass_min, mass_max, n_mass + 1)
    mass_centers = (mass_edges[:-1] + mass_edges[1:]) / 2.0
    rt_centers   = rt_min + (np.arange(n_rt) + 0.5) * rt_step

    fig.add_trace(go.Heatmap(
        z=log_mat,
        x=rt_centers,
        y=mass_centers,
        colorscale='Plasma',
        colorbar=dict(title='log₁₀(I+1)', len=0.75, tickfont=dict(size=9)),
        hovertemplate='RT: %{x:.2f} min<br>Mass: %{y:.1f} Da<br>log₁₀(I): %{z:.2f}<extra></extra>',
    ))
    fig.update_layout(
        title=f'Deconvolved MS2 Heatmap — {len(rts)} precursor masses',
        xaxis_title='Retention Time (min)',
        yaxis_title='Neutral Monoisotopic Mass (Da)',
        height=450,
        **_LAYOUT_BASE,
    )
    return fig
