#!/usr/bin/env python3
"""Animated versions of the dual-unitarity circle (thesis fig:circle) — three candidates.

Data: mu0_ladders.json, dumped from data/local/extended_rungs.jld2 — the SAME
mu_0(T) rungs the thesis Eq.(3) fits and the res_circle figure use, at
p = 0, 0.1, 0.3, 0.5 (windows T<=24, 17, 10, 6).

Candidates (all share one slow clock in T, so the p-dependence is visible):
  A  circle_A_panels.gif     2x2 panels, one coupling each, REAL data points
                             appearing at their measured T with a fading trail.
  B  circle_B_overlay.gif    one panel, all four trajectories overlaid.
  C  circle_C_synthetic.gif  2x2 panels, smooth synthetic trajectory per
                             coupling (radius = mean |mu0| of the data; phase
                             from a fit of theta(T) = a0*T + b + C/T to the
                             unwrapped data phases) — continuous winding
                             instead of discrete rungs.

What they demonstrate: the winding rate per unit T is a0 + O(1/T^2) — it slows
towards a constant as T grows, and the constant a0 (= -f v) grows with p, so
the p = 0.5 panel visibly races ahead of p = 0 on the same clock.

Usage:  python3 generate_circle_gifs.py [A|B|C|all]     (default: all)
Requires matplotlib + pillow.
"""

import json
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

HERE = pathlib.Path(__file__).parent
DATA = json.load(open(HERE / "mu0_ladders.json"))

PS = ["0.0", "0.1", "0.3", "0.5"]
COL = {"0.0": "#1F77B4", "0.1": "#C13639", "0.3": "#2E7D46", "0.5": "#E8862A"}
MRK = {"0.0": "o", "0.1": "s", "0.3": "D", "0.5": "^"}
INK = "#333A44"
REF = "#9AA3AD"

T_RATE = 0.55   # clock speed in T-units per second (small = slow, as requested)
FPS = 12
T_MIN, T_MAX = 2.0, 24.0
LIM = 2.15

series = {}
for p in PS:
    d = DATA[p]
    T = np.array(d["T"]); z = np.array(d["re"]) + 1j * np.array(d["im"])
    phi = np.array(d["phi"])          # unwrapped phase from the thesis analysis
    order = np.argsort(T)
    series[p] = (T[order], z[order], phi[order])


def path_upto(p, Tnow, n=400):
    """Trail along the circle: interpolate the unwrapped phase and the modulus
    in T between the measured rungs, so the drawn path follows the actual
    winding instead of cutting chords across the circle."""
    T, z, phi = series[p]
    Tend = min(Tnow, T[-1])
    if Tend <= T[0]:
        return np.array([]), np.array([])
    Td = np.linspace(T[0], Tend, n)
    ph = np.interp(Td, T, phi)
    rr = np.interp(Td, T, np.abs(z))
    # phi is the unwrapped phase of -mu (the analysis works with log(-mu)),
    # so the trajectory of mu itself is at phi + pi: without this shift the
    # drawn path is the data rotated by pi and misses every marker.
    zz = rr * np.exp(1j * (ph + np.pi))
    return zz, np.array([Td, ph])

n_frames = int((T_MAX - T_MIN) / T_RATE * FPS) + FPS  # + 1 s hold at the end
clock = np.minimum(T_MIN + np.arange(n_frames) / FPS * T_RATE, T_MAX)


def style(ax):
    ax.set_aspect("equal")
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM)
    ax.set_xticks([-2, -1, 0, 1, 2]); ax.set_yticks([-2, -1, 0, 1, 2])
    ax.grid(True, color=REF, alpha=0.22, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelsize=7)


def ref_circle(ax, r):
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(r * np.cos(th), r * np.sin(th), color=REF, ls="--", lw=1.0, zorder=1)


def trail_artists(ax, p):
    line, = ax.plot([], [], color=COL[p], lw=1.2, alpha=0.6, zorder=2)
    dots = ax.scatter([], [], s=14, marker=MRK[p], color=COL[p], zorder=3)
    head = ax.scatter([], [], s=55, marker=MRK[p], color=COL[p],
                      edgecolors="white", linewidths=0.8, zorder=4)
    return line, dots, head


def update_series(artists, p, Tnow):
    T, z, _phi = series[p]
    m = T <= Tnow + 1e-9
    line, dots, head = artists
    zz, _ = path_upto(p, Tnow)
    if zz.size:
        line.set_data(zz.real, zz.imag)
        head.set_offsets([[zz.real[-1], zz.imag[-1]]])
    dots.set_offsets(np.c_[z.real[m], z.imag[m]] if m.any() else np.empty((0, 2)))
    return artists


# --------------------------------------------------------------- candidate A
def make_A(out):
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 8.0), dpi=90)
    fig.subplots_adjust(hspace=0.28, wspace=0.24, top=0.90, bottom=0.06)
    clock_txt = fig.suptitle("", fontsize=13, color=INK)
    art = {}
    for ax, p in zip(axes.flat, PS):
        style(ax)
        T, z, _ = series[p]
        r = np.mean(np.abs(z))
        ref_circle(ax, r)
        ax.set_title(rf"p = {p}   $|\mu_0|$ = {r:.3f}   (T $\leq$ {T.max():.0f})",
                     fontsize=9, color=INK)
        art[p] = trail_artists(ax, p)

    def update(i):
        Tnow = clock[i]
        clock_txt.set_text(rf"$T = {Tnow:.1f}$")
        for p in PS:
            update_series(art[p], p, Tnow)
        return []

    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / FPS)
    anim.save(out, writer=PillowWriter(fps=FPS))
    plt.close(fig)


# --------------------------------------------------------------- candidate B
def build_B():
    """Figure + artists for candidate B, styled like the thesis res_circle:
    Re/Im axis labels and the thesis legend format."""
    fig, ax = plt.subplots(figsize=(7.6, 6.4), dpi=95)
    fig.subplots_adjust(top=0.90, right=0.76, left=0.11, bottom=0.10)
    style(ax)
    ax.set_xlabel(r"Re $\mu_0$", fontsize=10, color=INK)
    ax.set_ylabel(r"Im $\mu_0$", fontsize=10, color=INK)
    clock_txt = fig.suptitle("", fontsize=13, color=INK)
    art = {}
    for p in PS:
        T, z, _ = series[p]
        r = np.mean(np.abs(z))
        ref_circle(ax, r)
        art[p] = trail_artists(ax, p)
        art[p][1].set_label(rf"p={p},  $|\mu_0|$={r:.3f}")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8,
              frameon=False)

    def update(i):
        Tnow = clock[i]
        clock_txt.set_text(rf"$T = {Tnow:.1f}$")
        for p in PS:
            update_series(art[p], p, Tnow)
        return []

    return fig, update


def make_B(out):
    fig, update = build_B()
    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / FPS)
    anim.save(out, writer=PillowWriter(fps=FPS))
    plt.close(fig)


def export_frames_B(outdir, step=2, dpi=72):
    """PNG frame stack for embedding in the Beamer deck via \animategraphics.
    Every `step`-th frame; play at FPS/step to keep the real-time pacing.
    Frames are palette-quantized (the plots are thin lines on white), which
    roughly halves the stack size at no visible cost."""
    from PIL import Image
    outdir = pathlib.Path(outdir); outdir.mkdir(exist_ok=True)
    fig, update = build_B()
    k = 0
    for i in range(0, n_frames, step):
        update(i)
        p = outdir / f"frame_{k:03d}.png"
        fig.savefig(p, dpi=dpi)
        im = Image.open(p).convert("RGB").quantize(colors=128, method=Image.MEDIANCUT)
        im.save(p, optimize=True)
        k += 1
    plt.close(fig)
    print(f"wrote {k} frames to {outdir} (play at {FPS/step:.0f} fps)")


# --------------------------------------------------------------- candidate C
def synth(p):
    """Smooth trajectory: radius = data mean; phase from theta = a0*T + b + C/T fit."""
    T, z, _ = series[p]
    phi = series[p][2]
    A = np.c_[T, np.ones_like(T), 1.0 / T]
    coef, *_ = np.linalg.lstsq(A, phi, rcond=None)
    r = np.mean(np.abs(z))
    Td = np.linspace(T.min(), T.max(), 600)
    ph = coef[0] * Td + coef[1] + coef[2] / Td
    # same +pi as path_upto: phi is the phase of -mu
    return Td, r * np.exp(1j * (ph + np.pi)), r


def make_C(out):
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 8.0), dpi=90)
    fig.subplots_adjust(hspace=0.28, wspace=0.24, top=0.90, bottom=0.06)
    clock_txt = fig.suptitle("", fontsize=13, color=INK)
    art, dense = {}, {}
    for ax, p in zip(axes.flat, PS):
        style(ax)
        Td, zd, r = synth(p)
        dense[p] = (Td, zd)
        ref_circle(ax, r)
        ax.set_title(rf"p = {p}   $|\mu_0|$ = {r:.3f}   (T $\leq$ {Td.max():.0f})",
                     fontsize=9, color=INK)
        line, = ax.plot([], [], color=COL[p], lw=1.4, alpha=0.75, zorder=2)
        head = ax.scatter([], [], s=55, marker=MRK[p], color=COL[p],
                          edgecolors="white", linewidths=0.8, zorder=4)
        art[p] = (line, head)

    def update(i):
        Tnow = clock[i]
        clock_txt.set_text(rf"$T = {Tnow:.1f}$")
        for p in PS:
            Td, zd = dense[p]
            m = Td <= Tnow + 1e-9
            line, head = art[p]
            line.set_data(zd.real[m], zd.imag[m])
            if m.any():
                head.set_offsets([[zd.real[m][-1], zd.imag[m][-1]]])
        return []

    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / FPS)
    anim.save(out, writer=PillowWriter(fps=FPS))
    plt.close(fig)


if __name__ == "__main__":
    which = (sys.argv[1] if len(sys.argv) > 1 else "all").upper()
    if which in ("A", "ALL"):
        make_A(HERE / "circle_A_panels.gif"); print("wrote circle_A_panels.gif")
    if which in ("B", "ALL"):
        make_B(HERE / "circle_B_overlay.gif"); print("wrote circle_B_overlay.gif")
    if which in ("FRAMES",):
        export_frames_B(HERE / "imgs" / "circle_frames")
    if which in ("C", "ALL"):
        make_C(HERE / "circle_C_synthetic.gif"); print("wrote circle_C_synthetic.gif")
