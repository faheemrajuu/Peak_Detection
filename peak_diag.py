# peak_diag.py
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from scipy.signal import find_peaks, peak_widths

# ---- FONT STACKS (bold) ----
BOLD_FONT = "Arial Black, Arial Bold, Helvetica Bold, Arial, sans-serif"
REG_FONT  = "Arial, sans-serif"

# ==========================================================
# CONFIG
# ==========================================================
@dataclass
class PeakDiagConfig:
    # constants
    wavelength: float = 0.7863  # Å
    q_split: float = 1.0        # SAXS/WAXS boundary

    # behavior
    save_output: bool = True
    use_manual_peaks: bool = True

    # peak settings
    saxs_prominence_frac: float = 6e-3
    waxs_prominence_frac: float = 0.02

    # manual peaks
    manual_q_saxs: Tuple[float, ...] = (0.041, 0.047, 0.072)
    manual_q_waxs: Tuple[float, ...] = ()

    # styling / export
    width: int = 900
    height: int = 600
    export_width: int = 1200
    export_height: int = 800
    export_scale: int = 2

# ==========================================================
# IO
# ==========================================================
def load_dat(path: Path, skiprows: int = 2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Expects whitespace-delimited columns:
      col0=q, col1=I, col2=dI
    """
    df = pd.read_csv(path, sep=r"\s+", skiprows=skiprows)
    q = df.iloc[:, 0].to_numpy(dtype=float)
    I = df.iloc[:, 1].to_numpy(dtype=float)
    dI = df.iloc[:, 2].to_numpy(dtype=float)
    return q, I, dI

# ==========================================================
# PEAK DIAGNOSTICS
# ==========================================================
def peak_diagnostics(d: np.ndarray, y: np.ndarray, prominence_frac: float, wavelength: float) -> List[Dict[str, Any]]:
    if len(y) == 0:
        return []

    prom_abs = prominence_frac * np.nanmax(y)
    peaks, props = find_peaks(y, prominence=prom_abs)

    if len(peaks) == 0:
        return []

    widths = peak_widths(y, peaks, rel_height=0.5)
    dx = np.mean(np.diff(d))
    fwhm_d = widths[0] * dx

    results: List[Dict[str, Any]] = []
    for i, p in enumerate(peaks):
        d_peak = d[p]
        q_peak = 2 * np.pi / d_peak
        arg = q_peak * wavelength / (4 * np.pi)

        # numerical safety
        if not np.isfinite(q_peak):
            continue
        if arg <= 0 or arg > 1:
            continue

        theta = np.arcsin(arg)
        two_t = np.rad2deg(2 * theta)

        results.append(
            dict(
                d_peak=float(d_peak),
                q_peak=float(q_peak),
                **{"2theta_peak": float(two_t)},
                I_peak=float(y[p]),
                prominence=float(props["prominences"][i]),
                FWHM_d=float(fwhm_d[i]),
            )
        )
    return results


def peak_diagnostics_q(q: np.ndarray, y: np.ndarray, prominence_frac: float, wavelength: float) -> List[Dict[str, Any]]:
    if len(q) == 0:
        return []
    d = 2 * np.pi / q
    idx = np.argsort(d)  # ascending d
    return peak_diagnostics(d[idx], y[idx], prominence_frac, wavelength=wavelength)


def add_manual_peaks(q_manual: Tuple[float, ...], q: np.ndarray, y: np.ndarray) -> List[Dict[str, Any]]:
    peaks: List[Dict[str, Any]] = []
    if not q_manual:
        return peaks

    for qm in q_manual:
        idx = int(np.argmin(np.abs(q - qm)))
        d_peak = 2 * np.pi / q[idx]
        peaks.append(
            dict(
                d_peak=float(d_peak),
                q_peak=float(q[idx]),
                **{"2theta_peak": np.nan},
                I_peak=float(y[idx]),
                prominence=np.nan,
                FWHM_d=np.nan,
                source="manual",
            )
        )
    return peaks

# ==========================================================
# PLOTTING
# ==========================================================

def style_common(fig: go.Figure, width: int, height: int) -> None:
    fig.update_layout(
        # Global font (affects many things)
        font=dict(family=BOLD_FONT, size=28, color="black"),
        xaxis_title=dict(
            text=r"$\mathbf{q\ (Å^{-1})}$",
            font=dict(size=38, family=BOLD_FONT, color="black"),
        ),
        yaxis_title=dict(
            text="Intensity (a.u.)",
            font=dict(size=38, family=BOLD_FONT, color="black"),
        ),
        
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        template="plotly_white",
        width=width,
        height=height,
        margin=dict(l=80, r=40, t=40, b=80),
        showlegend=True,
        legend=dict(
            x=0.6, y=0.95,
            bgcolor="rgba(0,0,0,0)",
            font=dict(size=28, family=BOLD_FONT, color="black"),
        ),
    )

def make_saxs_figure(
    q_saxs: np.ndarray,
    I_saxs: np.ndarray,
    dI_saxs: np.ndarray,
    saxs_peaks: List[Dict[str, Any]],
    sample_name: str,
    cfg: PeakDiagConfig,
) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=q_saxs,
            y=I_saxs,
            mode="lines",
            name=f"{sample_name}_SAXS",
            line=dict(width=4, color="#0072B2"),
            error_y=dict(type="data", array=dI_saxs, visible=True, thickness=2, width=0),
        )
    )

    for p in saxs_peaks:
        fig.add_trace(
            go.Scatter(
                x=[p["q_peak"]],
                y=[p["I_peak"]],
                mode="markers",
                marker=dict(symbol="triangle-down", size=14, color="black"),
                showlegend=False,
            )
        )

    style_common(fig, width=cfg.width, height=cfg.height)

    # log axes with safe ranges
    xmin = np.min(q_saxs[q_saxs > 0]) if np.any(q_saxs > 0) else np.min(q_saxs)
    xmax = np.max(q_saxs)

    ymin = np.min(I_saxs[I_saxs > 0]) if np.any(I_saxs > 0) else np.min(I_saxs)
    ymax = np.max(I_saxs)

    fig.update_xaxes(
        type="log",
        showline=True,
        mirror=True,
        linewidth=3,
        linecolor="black",
        tickfont=dict(size=38, family=BOLD_FONT, color="black"),
        range=[np.log10(xmin), np.log10(xmax)],
    )
    fig.update_yaxes(
        type="log",
        showline=True,
        mirror=True,
        linewidth=3,
        linecolor="black",
        tickfont=dict(size=38, family=BOLD_FONT, color="black"),
        showticklabels=False,  # (you had this)
        range=[np.log10(ymin), np.log10(ymax * 1.05)],
    )
    return fig

def make_waxs_figure(
    q_waxs: np.ndarray,
    I_waxs: np.ndarray,
    waxs_peaks: List[Dict[str, Any]],
    sample_name: str,
    cfg: PeakDiagConfig,
) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=q_waxs,
            y=I_waxs,
            mode="lines",
            name=f"{sample_name}_WAXS",
            line=dict(width=4, color="#D55E00"),
        )
    )

    for p in waxs_peaks:
        fig.add_trace(
            go.Scatter(
                x=[p["q_peak"]],
                y=[p["I_peak"]],
                mode="markers",
                marker=dict(symbol="triangle-down", size=14, color="black"),
                showlegend=False,
            )
        )

    style_common(fig, width=cfg.width, height=cfg.height)

    fig.update_xaxes(
        type="linear",
        showline=True,
        mirror=True,
        linewidth=3,
        linecolor="black",
        tickfont=dict(size=38, family=BOLD_FONT, color="black"),
    )

    fig.update_yaxes(
        type="linear",
        showline=True,
        mirror=True,
        linewidth=3,
        linecolor="black",
        tickfont=dict(size=38, family=BOLD_FONT, color="black"),
    )
    
    return fig

# ==========================================================
# MAIN PIPELINE
# ==========================================================
def run_peak_diag(
    data_dir: Path,
    filename: str,
    cfg: PeakDiagConfig = PeakDiagConfig(),
    out_subdir: str = "peak_diag_output",
    skiprows: int = 2,
    show: bool = True,
) -> dict[str, Any]:
    data_path = data_dir / filename
    sample_name = Path(filename).stem

    out_dir = data_dir / out_subdir
    out_dir.mkdir(exist_ok=True)

    # load
    q, I, dI = load_dat(data_path, skiprows=skiprows)

    # split
    mask_saxs = q < cfg.q_split
    mask_waxs = q >= cfg.q_split

    q_saxs, I_saxs, dI_saxs = q[mask_saxs], I[mask_saxs], dI[mask_saxs]
    q_waxs, I_waxs = q[mask_waxs], I[mask_waxs]

    # sort by q descending (as you had)
    idx_saxs = np.argsort(q_saxs)[::-1]
    idx_waxs = np.argsort(q_waxs)[::-1]

    q_saxs, I_saxs, dI_saxs = q_saxs[idx_saxs], I_saxs[idx_saxs], dI_saxs[idx_saxs]
    q_waxs, I_waxs = q_waxs[idx_waxs], I_waxs[idx_waxs]

    # peaks: auto
    saxs_auto = peak_diagnostics_q(q_saxs, I_saxs, cfg.saxs_prominence_frac, wavelength=cfg.wavelength)
    waxs_auto = peak_diagnostics_q(q_waxs, I_waxs, cfg.waxs_prominence_frac, wavelength=cfg.wavelength)
    for p in saxs_auto:
        p["source"] = "auto"
    for p in waxs_auto:
        p["source"] = "auto"

    # peaks: manual
    saxs_manual = add_manual_peaks(cfg.manual_q_saxs, q_saxs, I_saxs) if cfg.use_manual_peaks else []
    waxs_manual = add_manual_peaks(cfg.manual_q_waxs, q_waxs, I_waxs) if cfg.use_manual_peaks else []

    saxs_peaks = saxs_auto + saxs_manual
    waxs_peaks = waxs_auto + waxs_manual

    # save tables
    if cfg.save_output:
        if saxs_peaks:
            pd.DataFrame(saxs_peaks).sort_values("q_peak").to_csv(
                out_dir / f"{sample_name}_SAXS_peaks.csv", index=False, float_format="%.3f"
            )
        if waxs_peaks:
            pd.DataFrame(waxs_peaks).sort_values("q_peak").to_csv(
                out_dir / f"{sample_name}_WAXS_peaks.csv", index=False, float_format="%.3f"
            )

    # figures
    fig_saxs = make_saxs_figure(q_saxs, I_saxs, dI_saxs, saxs_peaks, sample_name, cfg)
    fig_waxs = make_waxs_figure(q_waxs, I_waxs, waxs_peaks, sample_name, cfg)

    # save figures (needs kaleido installed)
    if cfg.save_output:
        fig_saxs.write_image(out_dir / f"{sample_name}_SAXS.png", width=cfg.export_width, height=cfg.export_height, scale=cfg.export_scale)
        fig_waxs.write_image(out_dir / f"{sample_name}_WAXS.png", width=cfg.export_width, height=cfg.export_height, scale=cfg.export_scale)
        print(f"✅ Saved tables + figures in: {out_dir}")

    # print summary
    print(f"SAXS peaks found: {len(saxs_peaks)}")
    print(f"WAXS peaks found: {len(waxs_peaks)}")

    if show:
        fig_saxs.show()
        fig_waxs.show()

    return dict(
        sample_name=sample_name,
        out_dir=out_dir,
        saxs_peaks=saxs_peaks,
        waxs_peaks=waxs_peaks,
        fig_saxs=fig_saxs,
        fig_waxs=fig_waxs,
    )

# ==========================================================
# CLI
# ==========================================================
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Peak diagnostics for SAXS/WAXS 1D data (.dat).")
    p.add_argument("--data_dir", type=str, required=True, help="Folder containing the .dat file")
    p.add_argument("--file", type=str, required=True, help="Filename, e.g. F33_T35.0C.dat")
    p.add_argument("--no_save", action="store_true", help="Disable saving CSV/PNG outputs")
    p.add_argument("--no_show", action="store_true", help="Do not display interactive figures")
    p.add_argument("--skiprows", type=int, default=2, help="Number of header rows to skip")
    p.add_argument("--q_split", type=float, default=1.0, help="SAXS/WAXS boundary in q (1/Å)")
    p.add_argument("--lambda_A", type=float, default=0.7863, help="Wavelength in Å")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    cfg = PeakDiagConfig(
        wavelength=args.lambda_A,
        q_split=args.q_split,
        save_output=(not args.no_save),
    )

    run_peak_diag(
        data_dir=Path(args.data_dir),
        filename=args.file,
        cfg=cfg,
        skiprows=args.skiprows,
        show=(not args.no_show),
    )


if __name__ == "__main__":
    main()
