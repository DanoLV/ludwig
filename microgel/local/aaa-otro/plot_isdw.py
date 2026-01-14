#!/usr/bin/env python3
"""
Plot the d_isdw (Inverse Squared Distance Weighting) interpolation weights
This script visualizes how the ISDW function distributes weights among neighboring nodes
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches

def d_isdw(r0, node_i, node_j, node_k, i_min, i_max, j_min, j_max, k_min, k_max):
    """
    Calculate the ISDW weight for a specific node.

    Parameters:
    -----------
    r0 : array-like [x, y, z]
        Particle position in lattice coordinates
    node_i, node_j, node_k : int
        Coordinates of the node for which we calculate the weight
    i_min, i_max, j_min, j_max, k_min, k_max : int
        Range of neighboring nodes to consider

    Returns:
    --------
    weight : float
        Normalized weight for the given node
    """
    epsilon = 1e-15  # Small value to avoid division by zero

    # Position of current node relative to particle
    r_current = np.array([r0[0] - node_i, r0[1] - node_j, r0[2] - node_k])

    # Distance from particle to current node
    dist_current = np.sqrt(np.sum(r_current**2))

    # Compute sum of inverse squared distances to all neighboring nodes
    inv_dist_sum = 0.0
    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            for k in range(k_min, k_max + 1):
                r_node = np.array([r0[0] - i, r0[1] - j, r0[2] - k])
                dist_node = np.sqrt(np.sum(r_node**2))
                inv_dist_sum += 1.0 / ((dist_node + epsilon)**2)

    # Calculate normalized weight for current node
    inv_dist_current = 1.0 / ((dist_current + epsilon)**2)

    if inv_dist_sum > epsilon:
        weight = inv_dist_current / inv_dist_sum
    else:
        weight = 0.0

    return weight


def plot_isdw_1d(drange=2.0):
    """Plot ISDW weights along a 1D line"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Test different particle positions
    positions = [0.0, 0.25, 0.5, 0.75]

    for idx, x0 in enumerate(positions):
        ax = axes[idx // 2, idx % 2]

        # Define range of nodes
        i_min = int(np.floor(x0 - drange))
        i_max = int(np.ceil(x0 + drange))

        # Calculate weights for each node
        nodes = np.arange(i_min, i_max + 1)
        weights = []

        for node in nodes:
            weight = d_isdw([x0, 0.5, 0.5], node, 0, 0,
                          i_min, i_max, 0, 0, 0, 0)
            weights.append(weight)

        # Plot
        ax.bar(nodes, weights, width=0.8, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax.axvline(x0, color='red', linestyle='--', linewidth=2, label=f'Particle at x={x0}')
        ax.set_xlabel('Node index i', fontsize=11)
        ax.set_ylabel('Weight', fontsize=11)
        ax.set_title(f'ISDW weights (particle at x={x0})', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()

        # Add weight values on bars
        for i, (node, weight) in enumerate(zip(nodes, weights)):
            if weight > 0.01:  # Only show significant weights
                ax.text(node, weight, f'{weight:.3f}',
                       ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    return fig


def plot_isdw_2d(x0=0.3, y0=0.7, drange=2.0):
    """Plot ISDW weights on a 2D grid"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Define range of nodes
    i_min = int(np.floor(x0 - drange))
    i_max = int(np.ceil(x0 + drange))
    j_min = int(np.floor(y0 - drange))
    j_max = int(np.ceil(y0 + drange))

    # Calculate weights for each node
    nodes_i = np.arange(i_min, i_max + 1)
    nodes_j = np.arange(j_min, j_max + 1)

    weights = np.zeros((len(nodes_j), len(nodes_i)))

    for ii, i in enumerate(nodes_i):
        for jj, j in enumerate(nodes_j):
            weight = d_isdw([x0, y0, 0.5], i, j, 0,
                          i_min, i_max, j_min, j_max, 0, 0)
            weights[jj, ii] = weight

    # Plot 1: Heatmap
    im = ax1.imshow(weights, origin='lower', cmap='YlOrRd', aspect='auto',
                    extent=[i_min-0.5, i_max+0.5, j_min-0.5, j_max+0.5])

    # Add grid lines
    for i in range(i_min, i_max + 2):
        ax1.axvline(i - 0.5, color='gray', linewidth=0.5, alpha=0.5)
    for j in range(j_min, j_max + 2):
        ax1.axhline(j - 0.5, color='gray', linewidth=0.5, alpha=0.5)

    # Add weight values
    for ii, i in enumerate(nodes_i):
        for jj, j in enumerate(nodes_j):
            if weights[jj, ii] > 0.01:
                ax1.text(i, j, f'{weights[jj, ii]:.3f}',
                        ha='center', va='center', fontsize=9,
                        color='black' if weights[jj, ii] < 0.5 else 'white',
                        fontweight='bold')

    # Mark particle position
    ax1.plot(x0, y0, 'b*', markersize=20, markeredgecolor='white', markeredgewidth=2,
            label=f'Particle at ({x0:.2f}, {y0:.2f})')

    ax1.set_xlabel('Node index i', fontsize=12)
    ax1.set_ylabel('Node index j', fontsize=12)
    ax1.set_title('ISDW weights (2D heatmap)', fontsize=13, fontweight='bold')
    ax1.legend(loc='upper right')

    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('Weight', fontsize=11)

    # Plot 2: 3D bar chart
    ax2 = fig.add_subplot(122, projection='3d')

    xpos, ypos = np.meshgrid(nodes_i, nodes_j)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    dx = dy = 0.8
    dz = weights.flatten()

    colors = plt.cm.YlOrRd(dz / (dz.max() if dz.max() > 0 else 1))

    ax2.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, shade=True, alpha=0.8)
    ax2.scatter([x0], [y0], [0], color='blue', s=200, marker='*',
               edgecolors='white', linewidths=2, depthshade=False)

    ax2.set_xlabel('Node index i', fontsize=11)
    ax2.set_ylabel('Node index j', fontsize=11)
    ax2.set_zlabel('Weight', fontsize=11)
    ax2.set_title('ISDW weights (3D view)', fontsize=13, fontweight='bold')

    plt.tight_layout()
    return fig


def plot_isdw_3d(x0=0.3, y0=0.5, z0=0.7, drange=2.0):
    """Plot ISDW weights for 3D stencil"""
    fig = plt.figure(figsize=(16, 6))

    # Define range of nodes
    i_min = int(np.floor(x0 - drange))
    i_max = int(np.ceil(x0 + drange))
    j_min = int(np.floor(y0 - drange))
    j_max = int(np.ceil(y0 + drange))
    k_min = int(np.floor(z0 - drange))
    k_max = int(np.ceil(z0 + drange))

    # Calculate weights for each node
    nodes_data = []

    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            for k in range(k_min, k_max + 1):
                weight = d_isdw([x0, y0, z0], i, j, k,
                              i_min, i_max, j_min, j_max, k_min, k_max)
                if weight > 0.001:  # Only store significant weights
                    nodes_data.append((i, j, k, weight))

    # Sort by weight
    nodes_data.sort(key=lambda x: x[3], reverse=True)

    # Plot 1: 3D scatter plot
    ax1 = fig.add_subplot(131, projection='3d')

    for i, j, k, weight in nodes_data:
        size = weight * 1000  # Scale for visibility
        color = plt.cm.YlOrRd(weight / nodes_data[0][3])  # Normalize to max weight
        ax1.scatter([i], [j], [k], s=size, c=[color], alpha=0.7, edgecolors='black')

    # Plot particle position
    ax1.scatter([x0], [y0], [z0], s=300, c='blue', marker='*',
               edgecolors='white', linewidths=2)

    ax1.set_xlabel('i', fontsize=11)
    ax1.set_ylabel('j', fontsize=11)
    ax1.set_zlabel('k', fontsize=11)
    ax1.set_title(f'ISDW 3D weights\nParticle at ({x0:.2f}, {y0:.2f}, {z0:.2f})',
                 fontsize=12, fontweight='bold')

    # Plot 2: XY projection (sum over z)
    ax2 = fig.add_subplot(132)

    xy_weights = {}
    for i, j, k, weight in nodes_data:
        key = (i, j)
        xy_weights[key] = xy_weights.get(key, 0) + weight

    for (i, j), weight in xy_weights.items():
        size = weight * 500
        color = plt.cm.YlOrRd(weight / max(xy_weights.values()))
        ax2.scatter([i], [j], s=size, c=[color], alpha=0.7, edgecolors='black')
        ax2.text(i, j, f'{weight:.2f}', ha='center', va='center', fontsize=8)

    ax2.plot(x0, y0, 'b*', markersize=15, markeredgecolor='white', markeredgewidth=1.5)
    ax2.set_xlabel('i', fontsize=11)
    ax2.set_ylabel('j', fontsize=11)
    ax2.set_title('XY projection (summed over k)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')

    # Plot 3: Weight distribution table
    ax3 = fig.add_subplot(133)
    ax3.axis('off')

    # Create table with top weights
    table_data = [['Node (i,j,k)', 'Weight', '% of total']]
    total_weight = sum([w for _, _, _, w in nodes_data])

    for idx, (i, j, k, weight) in enumerate(nodes_data[:10]):  # Top 10
        percentage = (weight / total_weight) * 100
        table_data.append([f'({i},{j},{k})', f'{weight:.4f}', f'{percentage:.1f}%'])

    table = ax3.table(cellText=table_data, cellLoc='left', loc='center',
                     colWidths=[0.4, 0.3, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header row
    for i in range(3):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')

    # Alternate row colors
    for i in range(1, len(table_data)):
        for j in range(3):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f0f0f0')

    ax3.set_title('Top contributing nodes', fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()
    return fig


def plot_weight_vs_distance(drange=2.0):
    """Plot how weight varies with distance from particle"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Test different ranges
    distances = np.linspace(0.01, drange, 100)

    for x0 in [0.0, 0.25, 0.5]:
        weights = []

        for d in distances:
            # Particle at x0, node at x0 + d
            node_pos = int(np.round(x0 + d))
            i_min = int(np.floor(x0 - drange))
            i_max = int(np.ceil(x0 + drange))

            # For this plot, we calculate weight for a node at distance d
            weight = d_isdw([x0, 0.5, 0.5], node_pos, 0, 0,
                          i_min, i_max, 0, 0, 0, 0)
            weights.append(weight)

        ax1.plot(distances, weights, linewidth=2, label=f'Particle at x={x0}', marker='o',
                markersize=3, markevery=10)

    ax1.set_xlabel('Distance from particle', fontsize=12)
    ax1.set_ylabel('Weight', fontsize=12)
    ax1.set_title('ISDW: Weight vs Distance', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Log-log plot
    ax2.loglog(distances, 1/distances**2, 'k--', linewidth=2, alpha=0.5,
              label='$1/r^2$ (unnormalized)')

    for x0 in [0.0, 0.25, 0.5]:
        weights = []

        for d in distances:
            node_pos = int(np.round(x0 + d))
            i_min = int(np.floor(x0 - drange))
            i_max = int(np.ceil(x0 + drange))

            weight = d_isdw([x0, 0.5, 0.5], node_pos, 0, 0,
                          i_min, i_max, 0, 0, 0, 0)
            weights.append(weight)

        ax2.loglog(distances, weights, linewidth=2, label=f'ISDW at x={x0}', marker='o',
                  markersize=3, markevery=10)

    ax2.set_xlabel('Distance from particle', fontsize=12)
    ax2.set_ylabel('Weight (log scale)', fontsize=12)
    ax2.set_title('ISDW: Log-log plot', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend()

    plt.tight_layout()
    return fig


def main():
    print("=" * 70)
    print("ISDW (Inverse Squared Distance Weighting) Interpolation Plotter")
    print("=" * 70)
    print()

    # 1D plots
    print("Generating 1D plots...")
    fig1 = plot_isdw_1d(drange=2.0)
    fig1.savefig('isdw_1d.png', dpi=150, bbox_inches='tight')
    print("  Saved: isdw_1d.png")

    # 2D plots
    print("Generating 2D plots...")
    fig2 = plot_isdw_2d(x0=0.3, y0=0.7, drange=2.0)
    fig2.savefig('isdw_2d.png', dpi=150, bbox_inches='tight')
    print("  Saved: isdw_2d.png")

    # 3D plots
    print("Generating 3D plots...")
    fig3 = plot_isdw_3d(x0=0.3, y0=0.5, z0=0.7, drange=2.0)
    fig3.savefig('isdw_3d.png', dpi=150, bbox_inches='tight')
    print("  Saved: isdw_3d.png")

    # Weight vs distance
    print("Generating weight vs distance plots...")
    fig4 = plot_weight_vs_distance(drange=2.0)
    fig4.savefig('isdw_weight_vs_distance.png', dpi=150, bbox_inches='tight')
    print("  Saved: isdw_weight_vs_distance.png")

    print()
    print("=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  - isdw_1d.png: 1D weights for different particle positions")
    print("  - isdw_2d.png: 2D heatmap and 3D bar chart")
    print("  - isdw_3d.png: Full 3D stencil with projections and table")
    print("  - isdw_weight_vs_distance.png: Weight decay with distance")

    plt.show()


if __name__ == '__main__':
    main()
