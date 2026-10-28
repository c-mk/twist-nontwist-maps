# program description:

'''
    repeated iteration of the standard map (sometimes called the Chirikov-Taylor map)
    to reproduce figures 

    a map satisfies the twist condition (is a twist map) if points higher up 
    along a vertical line make larger jumps in the horzontal direction. 
    Thus, we expect the rotation number (when it exist) to 
    increase monotically in the vertical direction
    '''
"""
phase_space_explorer.py
========================
One script, two maps: Chirikov standard (twist) map and the standard
nontwist map (SNM). Static phase portraits with labeled structures, or
a parameter-scan animation that SLOWS DOWN through the interesting
mixed regular/chaotic regime and speeds through the boring parts.

MAPS
----
Twist (standard) map:
    y_{n+1} = y_n - (K/2pi) sin(2 pi x_n)
    x_{n+1} = x_n + y_{n+1}                    (mod 1)

Nontwist map (SNM):
    y_{n+1} = y_n + b sin(2 pi x_n)
    x_{n+1} = x_n + a (1 - y_{n+1}^2)           (mod 1)

WHAT "X" AND "O" POINTS ARE  (this is what your PI means)
-----------------------------------------------------------
Both maps are symplectic (area-preserving, det[Jacobian]=1 everywhere
-- true regardless of K, a, b; you can check this analytically). A
period-q point p satisfies T^q(p) = p. Linearize the q-fold composed
map at p to get a 2x2 matrix M = DT^q. Because det(M)=1, its type is
set entirely by the trace, via Greene's residue:

    R = (2 - Tr(M)) / 4

    0 < R < 1   -> ELLIPTIC  -> "O-point"  (center of an island;
                                 nearby orbits circle it, like the O)
    R < 0 or R>1-> HYPERBOLIC-> "X-point"  (saddle; nearby orbits are
                                 stretched/squeezed along crossing
                                 directions, like the X)

O-points anchor island chains (regions of regular, bounded motion).
X-points sit on separatrices; the stable/unstable manifolds emanating
from an X-point intersect transversally (homoclinic tangle), and THAT
is the local seed of chaos -- as K (or b) grows, neighboring X-point
tangles widen and merge (resonance overlap), turning isolated chaotic
threads into one connected chaotic sea.

Both maps have UNIVERSAL, parameter-independent period-1 points, so
this script can always find and label them automatically:
    twist:    O at (0,0),  X at (0.5,0)      -- for any K
    nontwist: O at (0,1),  X at (0.5,1)      -- for any a,b
              (the other root branch is at y=-1; also always present)

PHASE-SPACE STRUCTURES YOU WILL SEE (labeled on the static plots)
-------------------------------------------------------------------
- Rotational invariant curve (RIC): a smooth curve an orbit traces out
  forever without closing -- the survivor of a KAM torus with
  irrational rotation number. Acts as a hard transport barrier: orbits
  above it can never cross to below it.
- Island (resonance chain): closed loops surrounding an O-point,
  formed where the local rotation number is rational (Poincare-
  Birkhoff theorem). Motion inside is regular and bounded.
- Chaotic sea: scattered, area-filling points with no visible curve
  structure -- orbits here wander without bound, sensitive to initial
  conditions (positive Lyapunov exponent), seeded by X-point tangles.
- (Nontwist map only) Shearless curve: y=0, where the twist condition
  omega'(y)=0 is violated (a single simple extremum of omega(y) =
  a(1-y^2)). This curve is often the LAST one to break as b increases,
  making it usually the most robust transport barrier in the system.


USAGE
-----
Static portrait, twist map:
    python3 phase_space_explorer.py --mode static --map twist --k 0.9

Static portrait, nontwist map:
    python3 phase_space_explorer.py --mode static --map nontwist --a 0.615 --b 0.30

Slow-through-the-interesting-part animation (twist map, sweeps K):
    python3 phase_space_explorer.py --mode animate --map twist

Slow-through-the-interesting-part animation (nontwist map, fixed a, sweeps b):
    python3 phase_space_explorer.py --mode animate --map nontwist --a 0.615

Custom slow-motion window and speeds:
    python3 phase_space_explorer.py --mode animate --map twist \\
        --slow-range 0.85 1.05 --slow-ms 300 --fast-ms 60

Output path (png for static, gif for animate) via --output.
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from PIL import Image
import io


# ----------------------------------------------------------------------
# 1. THE MAPS AND THEIR (UNWRAPPED, for Newton) STEP FUNCTIONS
# ----------------------------------------------------------------------
def twist_map(x, y, k):
    y_new = y - (k / (2.0 * np.pi)) * np.sin(2.0 * np.pi * x)
    x_new = (x + y_new) % 1.0
    return x_new, y_new


def twist_step_unwrapped(x, y, params):
    (k,) = params
    y_new = y - (k / (2.0 * np.pi)) * np.sin(2.0 * np.pi * x)
    x_new = x + y_new
    return x_new, y_new


def twist_jac(x, y, k):
    c = np.cos(2.0 * np.pi * x)
    dyn_dx = -k * c
    dyn_dy = 1.0
    dxn_dx = 1.0 + dyn_dx
    dxn_dy = 1.0
    return np.array([[dxn_dx, dxn_dy], [dyn_dx, dyn_dy]])


def nontwist_map(x, y, a, b):
    y_new = y + b * np.sin(2.0 * np.pi * x)
    x_new = (x + a * (1.0 - y_new**2)) % 1.0
    return x_new, y_new


def nontwist_step_unwrapped(x, y, params):
    a, b = params
    y_new = y + b * np.sin(2.0 * np.pi * x)
    x_new = x + a * (1.0 - y_new**2)
    return x_new, y_new


def nontwist_jac(x, y, a, b):
    c = np.cos(2.0 * np.pi * x)
    y_new = y + b * np.sin(2.0 * np.pi * x)
    dyn_dx = 2.0 * np.pi * b * c
    dyn_dy = 1.0
    dxn_dx = 1.0 - 2.0 * a * y_new * dyn_dx
    dxn_dy = -2.0 * a * y_new * dyn_dy
    return np.array([[dxn_dx, dxn_dy], [dyn_dx, dyn_dy]])


MAPS = {
    'twist':    dict(step=twist_map,               unwrapped=twist_step_unwrapped,
                      jac=twist_jac,
                      period1_guesses=[(0.0, 0.0), (0.5, 0.0)],
                      default_params=dict(k=0.9), ylim=(0, 1),
                      param_name='k'),
    'nontwist': dict(step=nontwist_map,             unwrapped=nontwist_step_unwrapped,
                      jac=nontwist_jac,
                      period1_guesses=[(0.0, 1.0), (0.5, 1.0), (0.0, -1.0), (0.5, -1.0)],
                      default_params=dict(a=0.615, b=0.30), ylim=(-2, 2),
                      param_name='b'),
}


# ----------------------------------------------------------------------
# 2. PERIODIC POINT SEARCH + CLASSIFICATION (Greene's residue)
# ----------------------------------------------------------------------
def iterate_unwrapped(step_func, jac_func, x, y, q, params):
    X, Y = x, y
    M = np.eye(2)
    for _ in range(q):
        J = jac_func(X, Y, *params)
        M = J @ M
        X, Y = step_func(X, Y, params)
    return X, Y, M


def periodic_residual(guess, step_func, jac_func, q, params, p):
    x, y = guess
    Xf, Yf, _ = iterate_unwrapped(step_func, jac_func, x, y, q, params)
    return [Xf - x - p, Yf - y]


def find_periodic_point(step_func, jac_func, q, params, guess, p=0):
    sol = fsolve(periodic_residual, guess,
                  args=(step_func, jac_func, q, params, p), xtol=1e-13)
    x0, y0 = sol[0] % 1.0, sol[1]
    _, _, M = iterate_unwrapped(step_func, jac_func, x0, y0, q, params)
    return x0, y0, M


def classify(M):
    tr = np.trace(M)
    R = (2.0 - tr) / 4.0
    kind = 'O' if 0.0 < R < 1.0 else 'X'
    return kind, R, tr


# ----------------------------------------------------------------------
# 3. STATIC PHASE PORTRAIT WITH X/O + STRUCTURE LABELS
# ----------------------------------------------------------------------
def draw_portrait(ax, map_name, params, orbit_y0s, x0=0.5, iterations=1500,
                   annotate_xo=True, annotate_structures=True):
    spec = MAPS[map_name]
    step, jac = spec['step'], spec['jac']

    cmap = plt.get_cmap('viridis')
    for i, y0 in enumerate(orbit_y0s):
        x, y = x0, y0
        xs, ys = [x], [y]
        for _ in range(iterations):
            x, y = step(x, y, *params)
            xs.append(x); ys.append(y)
        ax.plot(xs, ys, '.', ms=0.4, color=cmap(i / max(len(orbit_y0s)-1, 1)), alpha=0.6)

    if annotate_xo:
        uf, jf = spec['unwrapped'], spec['jac']
        for guess in spec['period1_guesses']:
            x0f, y0f, M = find_periodic_point(uf, jf, 1, params, guess, p=0)
            kind, R, tr = classify(M)
            marker = 'o' if kind == 'O' else 'x'
            color = 'lime' if kind == 'O' else 'red'
            ax.plot(x0f, y0f, marker, ms=11, mew=2.2, color=color, zorder=5)
            ax.annotate(f'{kind} (R={R:.2f})', (x0f, y0f),
                        textcoords='offset points', xytext=(8, 8),
                        fontsize=8, color=color, weight='bold')

    if annotate_structures:
        _add_structure_labels(ax, map_name, params)

    ax.set_xlim(0, 1)
    ax.set_ylim(*spec['ylim'])
    ax.set_xlabel('x'); ax.set_ylabel('y')
    return ax


def _add_structure_labels(ax, map_name, params):
    """Hand-placed labels calibrated to this script's DEFAULT parameters
    (twist: k=0.9 ; nontwist: a=0.615, b=0.30). If you pass very
    different parameters the islands/chaos move, so treat these as a
    starting template -- reposition by eye for your own parameter
    choice."""
    if map_name == 'twist':
        ax.annotate('island\n(regular)', xy=(0.45, 0.55), xytext=(0.15, 0.85),
                    fontsize=9, color='navy',
                    arrowprops=dict(arrowstyle='->', color='navy'))
        ax.annotate('chaotic sea\n(near X-point)', xy=(0.5, 0.05), xytext=(0.65, 0.15),
                    fontsize=9, color='darkred',
                    arrowprops=dict(arrowstyle='->', color='darkred'))
        ax.annotate('RIC\n(survives)', xy=(0.5, 0.72), xytext=(0.75, 0.9),
                    fontsize=9, color='darkgreen',
                    arrowprops=dict(arrowstyle='->', color='darkgreen'))
    elif map_name == 'nontwist':
        ax.axhline(0, color='crimson', ls='--', lw=0.8, alpha=0.6)
        ax.annotate('shearless curve (y=0)\nlast barrier to break',
                    xy=(0.9, 0.0), xytext=(0.55, -0.9),
                    fontsize=9, color='crimson',
                    arrowprops=dict(arrowstyle='->', color='crimson'))
        ax.annotate('island', xy=(0.5, 0.45), xytext=(0.15, 1.2),
                    fontsize=9, color='navy',
                    arrowprops=dict(arrowstyle='->', color='navy'))
        ax.annotate('chaotic sea', xy=(0.5, -1.3), xytext=(0.7, -1.7),
                    fontsize=9, color='darkred',
                    arrowprops=dict(arrowstyle='->', color='darkred'))


# ----------------------------------------------------------------------
# 4. SLOW-THROUGH-THE-INTERESTING-PART ANIMATION
#    Built by rendering frames to PIL images with PER-FRAME duration:
#    slow (long duration) inside the mixed regular/chaotic window,
#    fast (short duration) outside it.
# ----------------------------------------------------------------------
def build_schedule(p_min, p_max, slow_range, n_slow, n_fast_each, slow_ms, fast_ms):
    lo, hi = slow_range
    before = np.linspace(p_min, lo, n_fast_each, endpoint=False)
    middle = np.linspace(lo, hi, n_slow)
    after = np.linspace(hi, p_max, n_fast_each)[1:]
    values = np.concatenate([before, middle, after])
    durations = ([fast_ms] * len(before) + [slow_ms] * len(middle) +
                 [fast_ms] * len(after))
    return values, durations


def render_frame(map_name, params, orbit_y0s, x0, iterations, ylim, label):
    fig, ax = plt.subplots(figsize=(6, 6), dpi=110)
    step = MAPS[map_name]['step']
    for y0 in orbit_y0s:
        x, y = x0, y0
        xs, ys = [x], [y]
        for _ in range(iterations):
            x, y = step(x, y, *params)
            xs.append(x); ys.append(y)
        ax.plot(xs, ys, '.', ms=0.3, color='black', alpha=0.6)
    ax.set_xlim(0, 1); ax.set_ylim(*ylim)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(label)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert('RGB')


def animate(map_name, fixed_params, sweep_range, slow_range, output,
            num_orbits=40, iterations=350, n_slow=50, n_fast_each=8,
            slow_ms=300, fast_ms=70):
    spec = MAPS[map_name]
    ylim = spec['ylim']
    y0s = np.linspace(ylim[0]*0.95, ylim[1]*0.95, num_orbits)
    x0 = 0.5

    values, durations = build_schedule(sweep_range[0], sweep_range[1], slow_range,
                                        n_slow, n_fast_each, slow_ms, fast_ms)

    frames = []
    for v in values:
        if map_name == 'twist':
            params = (v,)
            label = f'Twist (standard) map, K={v:.4f}'
        else:
            a = fixed_params.get('a', 0.615)
            params = (a, v)
            label = f'Nontwist map, a={a}, b={v:.4f}'
        frames.append(render_frame(map_name, params, y0s, x0, iterations, ylim, label))

    frames[0].save(output, save_all=True, append_images=frames[1:],
                    duration=durations, loop=0)
    print(f"saved {output}  ({len(frames)} frames, "
          f"slow window [{slow_range[0]},{slow_range[1]}] at {slow_ms}ms/frame, "
          f"elsewhere {fast_ms}ms/frame)")


# ----------------------------------------------------------------------
# 5. CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=['static', 'animate'], default='static')
    ap.add_argument('--map', choices=['twist', 'nontwist'], default='twist')
    ap.add_argument('--k', type=float, default=None, help='twist map parameter')
    ap.add_argument('--a', type=float, default=None, help='nontwist map parameter a')
    ap.add_argument('--b', type=float, default=None, help='nontwist map parameter b')
    ap.add_argument('--num-orbits', type=int, default=40)
    ap.add_argument('--iterations', type=int, default=1500)
    ap.add_argument('--output', type=str, default=None)
    # animation-only options
    ap.add_argument('--slow-range', type=float, nargs=2, default=None,
                     help='parameter window to slow down through, e.g. --slow-range 0.85 1.05')
    ap.add_argument('--sweep-range', type=float, nargs=2, default=None,
                     help='full parameter sweep range')
    ap.add_argument('--slow-ms', type=int, default=300, help='ms per frame inside slow window')
    ap.add_argument('--fast-ms', type=int, default=70, help='ms per frame outside slow window')
    ap.add_argument('--n-slow', type=int, default=50, help='number of frames inside slow window')
    args = ap.parse_args()

    spec = MAPS[args.map]

    if args.mode == 'static':
        if args.map == 'twist':
            k = args.k if args.k is not None else spec['default_params']['k']
            params = (k,)
            title = f'Twist (standard) map, K={k}'
        else:
            a = args.a if args.a is not None else spec['default_params']['a']
            b = args.b if args.b is not None else spec['default_params']['b']
            params = (a, b)
            title = f'Nontwist map, a={a}, b={b}'

        orbit_y0s = (np.linspace(0, 1, args.num_orbits) if args.map == 'twist'
                     else np.concatenate([np.linspace(-1.8, -0.05, args.num_orbits//2),
                                           np.linspace(0.05, 1.8, args.num_orbits//2)]))

        fig, ax = plt.subplots(figsize=(8, 8))
        draw_portrait(ax, args.map, params, orbit_y0s, x0=0.5,
                      iterations=args.iterations)
        ax.set_title(title)
        plt.tight_layout()
        out = args.output or f'{args.map}_map_portrait.png'
        plt.savefig(out, dpi=200)
        plt.close(fig)
        print(f"saved {out}")

    else:  # animate
        if args.map == 'twist':
            sweep = args.sweep_range or (0.0, 5.0)
            slow = args.slow_range or (0.85, 1.05)   # brackets K_c = 0.9716...
            fixed = {}
            out = args.output or 'twist_map_scan.gif'
        else:
            a = args.a if args.a is not None else spec['default_params']['a']
            sweep = args.sweep_range or (0.0, 0.65)
            slow = args.slow_range or (0.15, 0.45)   # typical reconnection/mixed window
            fixed = {'a': a}
            out = args.output or 'nontwist_map_scan.gif'

        animate(args.map, fixed, sweep, slow, out,
                num_orbits=args.num_orbits, iterations=350,
                n_slow=args.n_slow, slow_ms=args.slow_ms, fast_ms=args.fast_ms)


if __name__ == '__main__':
    main()