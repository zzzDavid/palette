#!/usr/bin/env python3
"""Extract a publication-friendly color palette from a painting.

Cluster image pixels in CIELab, lightly retouch each centroid so it sits in
a print-safe band (L*, C* clipped, min pairwise ΔE enforced), then name each
color via nearest-neighbor in Lab against the XKCD color survey.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from colorspacious import cspace_convert
from sklearn.cluster import KMeans

# Publication-safe ranges (CIELab / CIELCh). See README for sources.
L_MIN, L_MAX = 25.0, 80.0       # avoid press-black and paper-white
C_MAX = 65.0                    # avoid garish, out-of-CMYK chroma
C_NEUTRAL = 6.0                 # below this we treat the color as a neutral
MIN_PAIR_DELTAE = 15.0          # min Lab Euclidean distance between any two palette colors
CVD_WARN_DELTAE = 12.0          # warn if any pair drops below this under deuteranopia


# ---------- image / clustering ----------

def load_image(path: str | Path, max_pixels: int = 400_000) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    w, h = img.size
    if w * h > max_pixels:
        scale = (max_pixels / (w * h)) ** 0.5
        img = img.resize((max(1, int(w * scale)), max(1, int(h * scale))), Image.LANCZOS)
    return np.asarray(img, dtype=np.float32) / 255.0


def extract_centers_lab(img_rgb: np.ndarray, k: int, seed: int = 0) -> np.ndarray:
    pixels_rgb = img_rgb.reshape(-1, 3)
    pixels_lab = cspace_convert(pixels_rgb, "sRGB1", "CIELab")
    km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(pixels_lab)
    return km.cluster_centers_


# ---------- Lab <-> LCh helpers ----------

def lab_to_lch(lab: np.ndarray) -> np.ndarray:
    L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]
    C = np.hypot(a, b)
    h = np.degrees(np.arctan2(b, a)) % 360.0
    return np.stack([L, C, h], axis=-1)


def lch_to_lab(lch: np.ndarray) -> np.ndarray:
    L, C, h = lch[..., 0], lch[..., 1], lch[..., 2]
    a = C * np.cos(np.radians(h))
    b = C * np.sin(np.radians(h))
    return np.stack([L, a, b], axis=-1)


def lab_to_srgb(lab: np.ndarray) -> np.ndarray:
    return np.clip(cspace_convert(lab, "CIELab", "sRGB1"), 0.0, 1.0)


def lab_to_hex(lab: np.ndarray) -> str:
    rgb255 = np.round(lab_to_srgb(np.asarray(lab)) * 255).astype(int)
    return "#{:02X}{:02X}{:02X}".format(*rgb255)


# ---------- treatment ----------

def treat_color(lab: np.ndarray) -> tuple[np.ndarray, bool]:
    """Pull a single Lab color into the publication-safe band.

    Lightness clipped to [L_MIN, L_MAX]; chroma capped at C_MAX (but very low
    chroma is preserved so intentional neutrals stay neutral).
    """
    L, C, h = lab_to_lch(lab)
    L_new = float(np.clip(L, L_MIN, L_MAX))
    if C <= C_NEUTRAL:
        C_new = float(C)
    else:
        C_new = float(min(C, C_MAX))
    changed = not (np.isclose(L_new, L) and np.isclose(C_new, C))
    return lch_to_lab(np.array([L_new, C_new, h])), changed


def enforce_min_distance(labs: np.ndarray,
                         min_d: float = MIN_PAIR_DELTAE,
                         max_iter: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Rotate hue of the offending member of any too-close pair until ΔE ≥ min_d."""
    labs = labs.astype(np.float64).copy()
    changed = np.zeros(len(labs), dtype=bool)
    for _ in range(max_iter):
        violated = False
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                if np.linalg.norm(labs[i] - labs[j]) < min_d:
                    violated = True
                    # rotate the dimmer (lower L*) one — keeps perceptual order stable
                    t = i if labs[i, 0] < labs[j, 0] else j
                    lch = lab_to_lch(labs[t])
                    lch[2] = (lch[2] + 20.0) % 360.0
                    labs[t] = lch_to_lab(lch)
                    changed[t] = True
        if not violated:
            break
    return labs, changed


# ---------- naming ----------

def build_name_table(use_xkcd: bool = True) -> list[tuple[str, np.ndarray]]:
    src = mcolors.XKCD_COLORS if use_xkcd else mcolors.CSS4_COLORS
    table = []
    for raw_name, hex_str in src.items():
        rgb = np.array(mcolors.to_rgb(hex_str), dtype=np.float32)
        lab = cspace_convert(rgb, "sRGB1", "CIELab")
        clean = raw_name.split(":", 1)[-1]
        table.append((clean, lab))
    return table


def closest_name(lab: np.ndarray, table: list[tuple[str, np.ndarray]]) -> str:
    arr = np.stack([row[1] for row in table])
    dists = np.linalg.norm(arr - lab, axis=1)
    return table[int(np.argmin(dists))][0]


# ---------- CVD check ----------

def cvd_min_pair_deltaE(labs: np.ndarray, kind: str = "deuteranomaly", severity: int = 100) -> float:
    rgb = lab_to_srgb(labs)
    cvd_space = {"name": "sRGB1+CVD", "cvd_type": kind, "severity": severity}
    rgb_cvd = cspace_convert(rgb, cvd_space, "sRGB1")
    rgb_cvd = np.clip(rgb_cvd, 0.0, 1.0)
    lab_cvd = cspace_convert(rgb_cvd, "sRGB1", "CIELab")
    n = len(lab_cvd)
    best = np.inf
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(lab_cvd[i] - lab_cvd[j]))
            if d < best:
                best = d
    return best


# ---------- plotting ----------

def _swatch(ax, hex_code: str) -> None:
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=hex_code, linewidth=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def make_figure(img_rgb: np.ndarray, palette: list[dict], out_path: Path) -> None:
    h, w = img_rgb.shape[:2]
    aspect = w / h
    n = len(palette)

    if aspect >= 1.35:
        # wide painting → stack: painting on top, palette row underneath
        fig_w = 11.0
        img_h = fig_w / aspect
        swatch_h = 1.6
        fig = plt.figure(figsize=(fig_w, img_h + swatch_h + 0.4),
                         constrained_layout=True)
        gs = fig.add_gridspec(2, n, height_ratios=[img_h, swatch_h])
        ax_img = fig.add_subplot(gs[0, :])
        ax_img.imshow(img_rgb)
        ax_img.set_axis_off()
        for i, c in enumerate(palette):
            ax = fig.add_subplot(gs[1, i])
            _swatch(ax, c["hex"])
            label = f"{c['name']}\n{c['hex']}"
            if c["treated"]:
                label += "  *"
            ax.set_xlabel(label, fontsize=10, family="DejaVu Sans Mono", labelpad=8)
    else:
        # square / portrait → side by side: painting left, palette column right
        fig_h = 7.5
        img_w = fig_h * aspect
        text_w = 3.2
        fig = plt.figure(figsize=(img_w + text_w + 0.4, fig_h),
                         constrained_layout=True)
        gs = fig.add_gridspec(n, 2, width_ratios=[img_w, text_w])
        ax_img = fig.add_subplot(gs[:, 0])
        ax_img.imshow(img_rgb)
        ax_img.set_axis_off()
        for i, c in enumerate(palette):
            ax = fig.add_subplot(gs[i, 1])
            _swatch(ax, c["hex"])
            label = f"  {c['name']}   {c['hex']}"
            if c["treated"]:
                label += "  *"
            ax.text(1.08, 0.5, label, transform=ax.transAxes,
                    va="center", ha="left",
                    fontsize=11, family="DejaVu Sans Mono")

    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------- pipeline ----------

def pick_palette(img_path: str | Path,
                 k: int = 5,
                 seed: int = 0,
                 out_path: str | Path | None = None,
                 use_xkcd: bool = True,
                 verbose: bool = True) -> list[dict]:
    img_path = Path(img_path)
    img_rgb = load_image(img_path)

    centers_lab = extract_centers_lab(img_rgb, k, seed=seed)
    centers_lab = centers_lab[np.argsort(centers_lab[:, 0])]   # darkest → lightest

    treated_labs = []
    treated_flags = []
    for lab in centers_lab:
        new_lab, changed = treat_color(lab)
        treated_labs.append(new_lab)
        treated_flags.append(changed)
    treated_labs = np.array(treated_labs)

    before = treated_labs.copy()
    treated_labs, hue_changed = enforce_min_distance(treated_labs)
    treated_flags = [bool(treated_flags[i] or hue_changed[i] or not np.allclose(treated_labs[i], before[i]))
                     for i in range(len(treated_labs))]

    name_table = build_name_table(use_xkcd=use_xkcd)
    palette = []
    for lab, flag in zip(treated_labs, treated_flags):
        palette.append({
            "hex": lab_to_hex(lab),
            "name": closest_name(lab, name_table),
            "lab": lab,
            "treated": flag,
        })

    if verbose:
        print(f"\nExtracted {k} colors from: {img_path}")
        print(f"{'#':>2}  {'Name':<24}  {'Hex':<8}  Treated?")
        print("-" * 52)
        for i, c in enumerate(palette, 1):
            print(f"{i:>2}  {c['name']:<24}  {c['hex']:<8}  {'yes' if c['treated'] else 'no'}")

        for kind in ("deuteranomaly", "protanomaly", "tritanomaly"):
            d = cvd_min_pair_deltaE(treated_labs, kind=kind)
            tag = "ok " if d >= CVD_WARN_DELTAE else "low"
            print(f"  min ΔE under {kind:<14}: {d:5.1f}  [{tag}]")

    if out_path is None:
        out_path = Path("outputs") / (img_path.stem + "_palette.png")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    make_figure(img_rgb, palette, out_path)
    if verbose:
        print(f"\nFigure saved to: {out_path}")
    return palette


def main():
    ap = argparse.ArgumentParser(description="Extract a publication-friendly color palette from a painting.")
    ap.add_argument("image", help="path to input painting image")
    ap.add_argument("-k", "--num-colors", type=int, default=5,
                    help="number of colors to extract (default: 5)")
    ap.add_argument("-o", "--output", default=None,
                    help="output figure path (default: outputs/<image>_palette.png)")
    ap.add_argument("--seed", type=int, default=0, help="K-means random seed")
    ap.add_argument("--css", action="store_true",
                    help="use the 148 CSS4 color names instead of the ~950 XKCD survey names")
    args = ap.parse_args()
    pick_palette(args.image, k=args.num_colors, seed=args.seed,
                 out_path=args.output, use_xkcd=not args.css)


if __name__ == "__main__":
    main()
