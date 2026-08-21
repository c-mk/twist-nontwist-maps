"""
Standard Nontwist Map (SNM)
============================
    y_{n+1} = y_n + b*sin(2*pi*x_n)
    x_{n+1} = x_n + a*(1 - y_{n+1}**2)      (mod 1)

Reference: del-Castillo-Negrete, Greene & Morrison,
"Area preserving nontwist maps: periodic orbits and
transition to chaos", Physica D 91 (1996) 1-23.

This script:
  1. Implements the map with the CORRECT sign/order (y updates first
     from the OLD x, then x updates using the NEW y).
  2. Verifies the map is symplectic (det[Jacobian] = 1) analytically,
     as a built-in regression test -- this is the check that would
     have caught the sign bug in the earlier version.
  3. Plots the twist (rotation) function omega(y) = a*(1 - y**2) and
     marks the shearless curve y = 0, i.e. the single simple extremum
     that violates the twist condition domega/dy != 0.
  4. Generates and plots the phase-space portrait (Figs. 3-4 style).
"""

import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# 1. THE MAP
# ----------------------------------------------------------------------
def nontwist_map(x, y, a, b):
    """One iterate of the standard nontwist map (SNM).

    Order matters: y_{n+1} is computed from x_n FIRST, then x_{n+1}
    is computed from the NEW y_{n+1}. This ordering is what makes the
    map symplectic (see jacobian() / check_symplectic() below).
    """
    y_new = y + b * np.sin(2.0 * np.pi * x)
    x_new = (x + a * (1.0 - y_new**2)) % 1.0
    return x_new, y_new


def jacobian(x, y, a, b):
    """Analytic Jacobian matrix of the map at (x, y)."""
    c = np.cos(2.0 * np.pi * x)
    y_new = y + b * np.sin(2.0 * np.pi * x)

    dyn_dy = 1.0
    dyn_dx = 2.0 * np.pi * b * c

    dxn_dx = 1.0 - 2.0 * a * y_new * dyn_dx
    dxn_dy = -2.0 * a * y_new * dyn_dy

    return np.array([[dxn_dx, dxn_dy],
                      [dyn_dx, dyn_dy]])


def check_symplectic(a, b, n_tests=8):
    """Sanity check: det(Jacobian) should equal 1 everywhere. This is
    the regression test that catches sign/ordering bugs immediately."""
    rng = np.random.default_rng(0)
    for _ in range(n_tests):
        x0, y0 = rng.random(), rng.uniform(-1.5, 1.5)
        detJ = np.linalg.det(jacobian(x0, y0, a, b))
        assert np.isclose(detJ, 1.0), f"Map is NOT symplectic: det J = {detJ}"
    print(f"[OK] map is symplectic for a={a}, b={b}  (det J = 1)")


# ----------------------------------------------------------------------
# 2. TWIST FUNCTION  omega(y) = a (1 - y^2)
# ----------------------------------------------------------------------
def plot_twist_function(a, ax=None):
    y = np.linspace(-1.5, 1.5, 400)
    omega = a * (1.0 - y**2)
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    ax.plot(y, omega, color='navy')
    ax.axvline(0, color='gray', ls='--', lw=1, label='shearless curve, y=0')
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel('y')
    ax.set_ylabel(r'$\omega(y)=a(1-y^2)$')
    ax.set_title('Twist function: single extremum at y=0')
    ax.legend(fontsize=8)
    return ax


# ----------------------------------------------------------------------
# 3. ORBITS / PHASE PORTRAIT
# ----------------------------------------------------------------------
def generate_orbits(a, b, y0_values, x0=0.5, iterations=2000):
    """Iterate the map for a family of initial conditions.

    y0_values : 1D array of initial y's, scanned across the
                shearless region so both twist-violating families
                (y>0 and y<0 branches with the same omega) show up.
    x0        : common starting x (0.5 is unremarkable / generic).
    """
    orbits = []
    for y0 in y0_values:
        x, y = x0, y0
        xs = np.empty(iterations + 1)
        ys = np.empty(iterations + 1)
        xs[0], ys[0] = x, y
        for n in range(iterations):
            x, y = nontwist_map(x, y, a, b)
            xs[n + 1], ys[n + 1] = x, y
        orbits.append((xs, ys))
    return orbits


def plot_phase_space(a, b, y0_values=None, x0=0.5, iterations=2000,
                      ylim=(-2.0, 2.0), ax=None):
    if y0_values is None:
        # denser sampling near the shearless curve, where the
        # interesting reconnection/meander structure lives
        y0_values = np.concatenate([
            np.linspace(-1.8, -0.2, 25),
            np.linspace(-0.2, 0.2, 15),
            np.linspace(0.2, 1.8, 25),
        ])

    orbits = generate_orbits(a, b, y0_values, x0=x0, iterations=iterations)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 7))

    cmap = plt.get_cmap('viridis')
    for i, (xs, ys) in enumerate(orbits):
        color = cmap(i / max(len(orbits) - 1, 1))
        ax.plot(xs, ys, '.', ms=0.4, color=color, alpha=0.7)

    ax.set_xlim(0, 1)
    ax.set_ylim(*ylim)
    ax.set_xlabel('x (mod 1)')
    ax.set_ylabel('y')
    ax.set_title(f'Standard Nontwist Map: a={a}, b={b}')
    return ax


# ----------------------------------------------------------------------
# 4. RUN  -- change a, b here to explore the map
# ----------------------------------------------------------------------
if __name__ == '__main__':
    a, b = 0.615, 0.4

    check_symplectic(a, b)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6),
                              gridspec_kw={'width_ratios': [1, 2]})
    plot_twist_function(a, ax=axes[0])
    plot_phase_space(a, b, ax=axes[1])
    plt.tight_layout()
    plt.savefig('nontwist_map.png', dpi=200)
    plt.show()