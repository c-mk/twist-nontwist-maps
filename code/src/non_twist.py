"""
Standard Twist and Nontwist Maps
================================
Generates and saves phase space portraits for Hamiltonian map systems.
Plots are strictly single-panel phase portraits in grayscale/black.
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. MAP DEFINITIONS
# ----------------------------------------------------------------------
def nontwist_map(x, y, a, b):
    """Standard Nontwist Map (SNM)"""
    y_new = y + b * np.sin(2.0 * np.pi * x)
    x_new = (x + a * (1.0 - y_new**2)) % 1.0
    return x_new, y_new

def twist_map(x, y, K):
    """Standard Twist Map (Chirikov-Taylor)"""
    y_new = (y + (K / (2.0 * np.pi)) * np.sin(2.0 * np.pi * x)) % 1.0
    x_new = (x + y_new) % 1.0
    return x_new, y_new

# ----------------------------------------------------------------------
# 2. ORBIT GENERATION
# ----------------------------------------------------------------------
def generate_orbits(map_func, map_args, y0_values, x0_values=None, iterations=2000):
    """Iterate a generic map for a family of initial conditions."""
    orbits = []
    
    # If no specific x0 array is provided, start all at x=0.5
    if x0_values is None:
        x0_values = [0.5] * len(y0_values)
        
    for x0, y0 in zip(x0_values, y0_values):
        x, y = x0, y0
        xs = np.empty(iterations + 1)
        ys = np.empty(iterations + 1)
        xs[0], ys[0] = x, y
        for n in range(iterations):
            x, y = map_func(x, y, *map_args)
            xs[n + 1], ys[n + 1] = x, y
        orbits.append((xs, ys))
    return orbits

# ----------------------------------------------------------------------
# 3. PLOTTING PIPELINE
# ----------------------------------------------------------------------
def plot_phase_portrait(map_func, map_args, y0_values, filename, title, 
                        param_text, ylim, x0_values=None):
    """Generates a single-panel phase portrait and saves it."""
    fig, ax = plt.subplots(figsize=(8, 7))
    
    orbits = generate_orbits(map_func, map_args, y0_values, x0_values, iterations=2000)
    
    # Rendering in black with alpha for depth
    for xs, ys in orbits:
        ax.plot(xs, ys, '.', ms=0.8, color='black', alpha=0.5)

    ax.set_xlim(0, 1)
    ax.set_ylim(*ylim)
    ax.set_xlabel('x (mod 1)')
    ax.set_ylabel('y')
    ax.set_title(title)
    
    # Add parameter text box on the chart
    bbox_props = dict(boxstyle="round,pad=0.4", fc="white", ec="black", lw=1.5)
    ax.text(0.05, 0.95, param_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=bbox_props, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    print(f"Saved: {filename}")
    plt.close()

# ----------------------------------------------------------------------
# 4. EXECUTION
# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # --- Generate Nontwist Map ---
    a_nt, b_nt = 0.3, 0.715
    y0_nt = np.concatenate([
        np.linspace(-1.8, -0.2, 25),
        np.linspace(-0.2, 0.2, 15),
        np.linspace(0.2, 1.8, 25),
    ])
    ''
    plot_phase_portrait(
        map_func=nontwist_map,
        map_args=(a_nt, b_nt),
        y0_values=y0_nt,
        filename='nontwist_map.png',
        title='Standard Nontwist Map',
        param_text=f"a = {a_nt}\nb = {b_nt}",
        ylim=(-2.0, 2.0)
    )
    
    # --- Generate Standard Twist Map ---
    K_t = 0.97
    
    # Using a 2D grid of initial conditions ensures the twist map 
    # phase space is fully populated and doesn't look blank.
    grid_size = 15
    y0_grid, x0_grid = np.meshgrid(np.linspace(0, 1, grid_size), 
                                   np.linspace(0, 1, grid_size))
    
    plot_phase_portrait(
        map_func=twist_map,
        map_args=(K_t,),
        y0_values=y0_grid.flatten(),
        x0_values=x0_grid.flatten(),
        filename='standard_twist_map.png',
        title='Standard Twist Map',
        param_text=f"K = {K_t}", # The Standard Twist Map only uses 'K'
        ylim=(0.0, 1.0)
    )