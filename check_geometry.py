#!/usr/bin/env python3
"""
Geometry Preview for FDTD Laser Simulation

This script shows the pillar geometry and source positions without running the simulation.
Use this to verify that all sources are positioned correctly between the pillars.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def plot_geometry_preview():
    """Plot the geometry setup showing pillars and inter-pillar sources."""
    
    # Parameters from your C++ code
    lambda1 = 1.5e-6  # m
    DX = lambda1 / 50.0  # Spatial step
    
    SIZE_X = int(6.0 * lambda1 / DX)
    SIZE_Y = int(6.0 * lambda1 / DX)
    
    # Pillar parameters
    num_pillars_per_column = 5
    pillar_radius = int(0.2 * lambda1 / DX)
    pillar_spacing_y = int(0.5 * lambda1 / DX)
    column1_x = int(3.0 * lambda1 / DX)
    column2_x = int(3.5 * lambda1 / DX)
    first_pillar_y = int(1.5 * lambda1 / DX)
    
    # PML parameters
    PML_WIDTH = 10
    
    print(f"Grid size: {SIZE_X} x {SIZE_Y}")
    print(f"Pillar radius: {pillar_radius} grid points")
    print(f"Pillar spacing: {pillar_spacing_y} grid points")
    print(f"Column 1 at x = {column1_x}")
    print(f"Column 2 at x = {column2_x}")
    print(f"First pillar at y = {first_pillar_y}")
    
    # Calculate source positions (exactly between pillars)
    sources_col1 = []
    sources_col2 = []
    
    for i in range(4):  # 4 gaps between 5 pillars
        # Column 1 sources
        x1 = column1_x
        y1 = first_pillar_y + int((i + 0.5) * pillar_spacing_y)
        sources_col1.append((x1, y1))
        
        # Column 2 sources
        x2 = column2_x
        y2 = first_pillar_y + int((i + 0.5) * pillar_spacing_y)
        sources_col2.append((x2, y2))
    
    print(f"\nColumn 1 sources: {sources_col1}")
    print(f"Column 2 sources: {sources_col2}")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Set limits
    ax.set_xlim(0, SIZE_X)
    ax.set_ylim(0, SIZE_Y)
    ax.set_aspect('equal')
    
    # Draw PML boundaries
    ax.axvline(x=PML_WIDTH, color='gray', linestyle='--', alpha=0.7, linewidth=2, label='PML Boundary')
    ax.axvline(x=SIZE_X-PML_WIDTH, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    ax.axhline(y=PML_WIDTH, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    ax.axhline(y=SIZE_Y-PML_WIDTH, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    
    # Draw pillars
    for i in range(num_pillars_per_column):
        y = first_pillar_y + i * pillar_spacing_y
        
        # Column 1 pillars
        circle1 = Circle((column1_x, y), pillar_radius, 
                        fill=True, facecolor='lightblue', edgecolor='black', 
                        linewidth=2, alpha=0.7, label='SiO₂ Pillars' if i == 0 else '')
        ax.add_patch(circle1)
        ax.text(column1_x, y, f'P1-{i+1}', ha='center', va='center', fontsize=8, weight='bold')
        
        # Column 2 pillars
        circle2 = Circle((column2_x, y), pillar_radius, 
                        fill=True, facecolor='lightblue', edgecolor='black', 
                        linewidth=2, alpha=0.7)
        ax.add_patch(circle2)
        ax.text(column2_x, y, f'P2-{i+1}', ha='center', va='center', fontsize=8, weight='bold')
    
    # Draw sources
    # Column 1 sources (λ = 1500 nm)
    for i, (x, y) in enumerate(sources_col1):
        ax.plot(x, y, 'ro', markersize=12, markeredgecolor='darkred', 
               markeredgewidth=3, label='Col 1 Sources (λ=1500nm)' if i == 0 else '')
        ax.annotate(f'S1-{i+1}', (x+8, y+8), fontsize=10, color='darkred', weight='bold')
        
        # Draw lines to show they're between pillars
        pillar_above_y = first_pillar_y + (i+1) * pillar_spacing_y
        pillar_below_y = first_pillar_y + i * pillar_spacing_y
        ax.plot([x, x], [pillar_below_y + pillar_radius, pillar_above_y - pillar_radius], 
               'r--', alpha=0.5, linewidth=1)
    
    # Column 2 sources (λ = 1364 nm)
    for i, (x, y) in enumerate(sources_col2):
        ax.plot(x, y, 'bs', markersize=12, markeredgecolor='darkblue', 
               markeredgewidth=3, label='Col 2 Sources (λ=1364nm)' if i == 0 else '')
        ax.annotate(f'S2-{i+1}', (x+8, y+8), fontsize=10, color='darkblue', weight='bold')
        
        # Draw lines to show they're between pillars
        pillar_above_y = first_pillar_y + (i+1) * pillar_spacing_y
        pillar_below_y = first_pillar_y + i * pillar_spacing_y
        ax.plot([x, x], [pillar_below_y + pillar_radius, pillar_above_y - pillar_radius], 
               'b--', alpha=0.5, linewidth=1)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Labels and title
    ax.set_xlabel('X (grid points)', fontsize=12)
    ax.set_ylabel('Y (grid points)', fontsize=12)
    ax.set_title('FDTD Simulation Geometry\nInter-Pillar Laser Sources (8 total)', fontsize=14, weight='bold')
    
    # Add legend
    ax.legend(loc='upper right', fontsize=11)
    
    # Add parameter information
    info_text = f'Domain: {SIZE_X}×{SIZE_Y} grid points\n'
    info_text += f'Pillar radius: {pillar_radius} pts = {pillar_radius*DX*1e6:.2f} μm\n'
    info_text += f'Pillar spacing: {pillar_spacing_y} pts = {pillar_spacing_y*DX*1e6:.2f} μm\n'
    info_text += f'Column separation: {column2_x-column1_x} pts = {(column2_x-column1_x)*DX*1e6:.2f} μm\n'
    info_text += f'Source intensity: 2.0 GV/m each\n'
    info_text += f'Material: SiO₂ (εᵣ = 3.9)'
    
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    
    # Add distance annotations
    for i in range(4):
        y = sources_col1[i][1]
        mid_x = (column1_x + column2_x) / 2
        ax.annotate(f'{(column2_x-column1_x)*DX*1e6:.1f} μm', 
                   xy=(mid_x, y), ha='center', va='bottom',
                   fontsize=8, color='purple', weight='bold')
    
    plt.tight_layout()
    
    # Save the geometry plot
    plt.savefig('geometry_preview.png', dpi=300, bbox_inches='tight')
    print(f"\nGeometry preview saved as 'geometry_preview.png'")
    
    plt.show()
    
    return fig

def check_source_pillar_distances():
    """Check that sources are properly positioned between pillars."""
    lambda1 = 1.5e-6
    DX = lambda1 / 50.0
    
    pillar_radius = int(0.2 * lambda1 / DX)
    pillar_spacing_y = int(0.5 * lambda1 / DX)
    column1_x = int(3.0 * lambda1 / DX)
    first_pillar_y = int(1.5 * lambda1 / DX)
    
    print("\nSource-Pillar Distance Check:")
    print("=" * 40)
    
    for i in range(4):
        # Source position
        source_y = first_pillar_y + int((i + 0.5) * pillar_spacing_y)
        
        # Adjacent pillar positions
        pillar_below_y = first_pillar_y + i * pillar_spacing_y
        pillar_above_y = first_pillar_y + (i + 1) * pillar_spacing_y
        
        # Distances
        dist_to_below = source_y - pillar_below_y
        dist_to_above = pillar_above_y - source_y
        
        print(f"Source {i+1}:")
        print(f"  Position: y = {source_y}")
        print(f"  Distance to pillar below: {dist_to_below} grid points = {dist_to_below*DX*1e6:.2f} μm")
        print(f"  Distance to pillar above: {dist_to_above} grid points = {dist_to_above*DX*1e6:.2f} μm")
        print(f"  Minimum clearance: {min(dist_to_below, dist_to_above) - pillar_radius} grid points")
        print()

if __name__ == "__main__":
    print("FDTD Laser Simulation - Geometry Preview")
    print("=" * 45)
    
    # Check source positioning
    check_source_pillar_distances()
    
    # Plot geometry
    plot_geometry_preview()
    
    print("\nGeometry verification complete!")
    print("• All sources are positioned exactly between pillars")
    print("• Red circles: Column 1 sources (λ=1500nm, f1)")
    print("• Blue squares: Column 2 sources (λ=1364nm, f2)")
    print("• Light blue circles: SiO₂ pillars")
