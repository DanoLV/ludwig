#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def peskin_delta(r):
    """
    Peskin delta function implementation
    """
    rmod = np.abs(r)
    delta = np.zeros_like(rmod)

    # Case 1: rmod <= 1.0
    mask1 = rmod <= 1.0
    delta[mask1] = 0.125 * (3.0 - 2.0 * rmod[mask1] +
                            np.sqrt(1.0 + 4.0 * rmod[mask1] - 4.0 * rmod[mask1]**2))

    # Case 2: 1.0 < rmod <= 2.0
    mask2 = (rmod > 1.0) & (rmod <= 2.0)
    delta[mask2] = 0.125 * (5.0 - 2.0 * rmod[mask2] -
                            np.sqrt(-7.0 + 12.0 * rmod[mask2] - 4.0 * rmod[mask2]**2))

    # For rmod > 2.0, delta remains 0

    return delta

# Generate r values from -5 to 5
r = np.linspace(-5, 5, 1000)

# Calculate delta values
delta = peskin_delta(r)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(r, delta, 'b-', linewidth=2, label='Peskin delta function')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=-2, color='r', linestyle=':', alpha=0.5, label='Support boundaries (±2)')
plt.axvline(x=2, color='r', linestyle=':', alpha=0.5)
plt.grid(True, alpha=0.3)
plt.xlabel('r', fontsize=12)
plt.ylabel('δ(r)', fontsize=12)
plt.title('Peskin Delta Function', fontsize=14, fontweight='bold')
plt.legend(fontsize=10)
plt.xlim(-2.1, 2.1)
plt.ylim(-0.05, 0.55)

# Add annotations for key points
plt.annotate('δ(0) ≈ 0.5', xy=(0, peskin_delta(np.array([0]))[0]),
             xytext=(0.5, 0.52), fontsize=10,
             arrowprops=dict(arrowstyle='->', color='black', alpha=0.7))

plt.tight_layout()
plt.savefig('peskin_delta_function.png', dpi=150, bbox_inches='tight')
print("Plot saved as 'peskin_delta_function.png'")
plt.show()

# Print some key values
print("\nKey values:")
print(f"δ(0) = {peskin_delta(np.array([0.0]))[0]:.6f}")
print(f"δ(1) = {peskin_delta(np.array([1.0]))[0]:.6f}")
print(f"δ(2) = {peskin_delta(np.array([2.0]))[0]:.6f}")
print(f"δ(2.5) = {peskin_delta(np.array([2.5]))[0]:.6f}")
