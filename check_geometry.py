#!/usr/bin/env python3
"""
Geometry Check for Bar and Teeth Structure

This script helps visualize the new dielectric bar and teeth geometry
to verify the structure is correct before running the simulation.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches

def visualize_bar_geometry():
    """Create a visualization of the bar and teeth geometry."""
    
    # Parameters from C++ code (approximate values)
    lambda1 = 2e-6  # wavelength
    DX = lambda1 / 50.0  # grid spacing
    
    # Convert to grid units
    SIZE_X = int(6.0 * lambda1 / DX)
    SIZE_Y = int(6.0 * lambda1 / DX)
    PML_WIDTH = 10
    
    # Bar parameters
    bar_width = int(0.15 * lambda1 / DX)
    bar1_x = int(2.8 * lambda1 / DX)
    bar2_x = int(3.7 * lambda1 / DX)
    
    # Teeth parameters
    num_teeth = 6
    tooth_width = int(0.3 * lambda1 / DX)
    tooth_height = int(0.4 * lambda1 / DX)
    tooth_spacing = int(0.5 * lambda1 / DX)
    first_tooth_y = int(1.5 * lambda1 / DX)
    
    # Bar length
    bar_start_y = PML_WIDTH + 10
    bar_end_y = SIZE_Y - PML_WIDTH - 10
    
    # Source positions
    gap_center_x = (bar1_x + bar2_x) // 2
    
    print("FDTD Laser Simulation - Bar and Teeth Geometry Check")
    print("=" * 55)
    print(f"Grid size: {SIZE_X} x {SIZE_Y}")
    print(f"Bar width: {bar_width} grid points")
    print(f"Bar 1 center at x = {bar1_x}")
    print(f"Bar 2 center at x = {bar2_x}")
    print(f"Gap width: {bar2_x - bar1_x} grid points")
    print(f"Number of teeth per bar: {num_teeth}")
    print(f"Tooth dimensions: {tooth_height} x {tooth_width} grid points")
    print(f"Tooth spacing: {tooth_spacing} grid points")
    
    # Calculate source positions (between teeth in gap)
    sources = []
    print(f"\nInter-teeth Source Positions:")
    print("=" * 35)
    
    for i in range(5):  # 5 sources between 6 teeth
        x = gap_center_x
        y = first_tooth_y + int((i + 0.5) * tooth_spacing)
        sources.append((x, y))
        print(f"Source {i+1}: x={x}, y={y} (gap center, between teeth {i+1} and {i+2})")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Draw PML boundaries
    pml_color = 'lightgray'
    ax.axvline(x=PML_WIDTH, color=pml_color, linestyle='--', alpha=0.7, label='PML Boundary')
    ax.axvline(x=SIZE_X-PML_WIDTH, color=pml_color, linestyle='--', alpha=0.7)
    ax.axhline(y=PML_WIDTH, color=pml_color, linestyle='--', alpha=0.7)
    ax.axhline(y=SIZE_Y-PML_WIDTH, color=pml_color, linestyle='--', alpha=0.7)
    
    # Draw Bar 1 (left bar)
    bar1_rect = Rectangle((bar1_x - bar_width//2, bar_start_y), 
                         bar_width, bar_end_y - bar_start_y,
                         facecolor='lightblue', edgecolor='blue', 
                         alpha=0.7, label='Dielectric Bars')
    ax.add_patch(bar1_rect)
    
    # Draw Bar 2 (right bar)
    bar2_rect = Rectangle((bar2_x - bar_width//2, bar_start_y), 
                         bar_width, bar_end_y - bar_start_y,
                         facecolor='lightblue', edgecolor='blue', alpha=0.7)
    ax.add_patch(bar2_rect)
    
    # Draw teeth extending from both bars
    for t in range(num_teeth):
        tooth_center_y = first_tooth_y + t * tooth_spacing
        
        # Tooth from Bar 1 (extending right)
        tooth1_rect = Rectangle((bar1_x + bar_width//2, tooth_center_y - tooth_width//2),
                               tooth_height, tooth_width,
                               facecolor='lightblue', edgecolor='blue', alpha=0.7)
        ax.add_patch(tooth1_rect)
        
        # Tooth from Bar 2 (extending left)
        tooth2_rect = Rectangle((bar2_x - bar_width//2 - tooth_height, tooth_center_y - tooth_width//2),
                               tooth_height, tooth_width,
                               facecolor='lightblue', edgecolor='blue', alpha=0.7)
        ax.add_patch(tooth2_rect)
        
        # Label teeth
        ax.text(bar1_x + bar_width//2 + tooth_height//2, tooth_center_y, f'T{t+1}', 
               ha='center', va='center', fontsize=8, color='darkblue')
    
    # Draw laser sources
    for i, (x, y) in enumerate(sources):
        ax.plot(x, y, 'ro', markersize=10, 
               markeredgecolor='darkred', markeredgewidth=2,
               label='Laser Sources' if i == 0 else '')
        ax.text(x + 5, y, f'S{i+1}', 
               fontsize=8, color='darkred')
    
    # Add dimensions and annotations
    ax.annotate('', xy=(bar1_x - bar_width//2, SIZE_Y - 20), 
               xytext=(bar1_x + bar_width//2, SIZE_Y - 20),
               arrowprops=dict(arrowstyle='<->', color='black'))
    ax.text(bar1_x, SIZE_Y - 15, f'Bar Width\n{bar_width} pts', 
           ha='center', va='bottom', fontsize=8)
    
    ax.annotate('', xy=(bar1_x + bar_width//2, first_tooth_y), 
               xytext=(bar1_x + bar_width//2 + tooth_height, first_tooth_y),
               arrowprops=dict(arrowstyle='<->', color='green'))
    ax.text(bar1_x + bar_width//2 + tooth_height//2, first_tooth_y - 10, 
           f'Tooth Height\n{tooth_height} pts', ha='center', va='top', fontsize=8)
    
    ax.annotate('', xy=(bar1_x, bar_end_y + 5), 
               xytext=(bar2_x, bar_end_y + 5),
               arrowprops=dict(arrowstyle='<->', color='purple'))
    ax.text((bar1_x + bar2_x)//2, bar_end_y + 15, 
           f'Gap: {bar2_x - bar1_x} pts', ha='center', va='bottom', fontsize=8)
    
    # Labels and formatting
    ax.set_xlim(0, SIZE_X)
    ax.set_ylim(0, SIZE_Y)
    ax.set_xlabel('X (grid points)', fontsize=12)
    ax.set_ylabel('Y (grid points)', fontsize=12)
    ax.set_title('Dielectric Laser Accelerator - Bar and Teeth Geometry', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Add parameter text
    param_text = f'Parameters:\n'
    param_text += f'• {num_teeth} teeth per bar\n'
    param_text += f'• Tooth spacing: {tooth_spacing} pts\n'
    param_text += f'• Bar width: {bar_width} pts\n'
    param_text += f'• Tooth height: {tooth_height} pts\n'
    param_text += f'• 5 laser sources in gap'
    
    ax.text(0.02, 0.98, param_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('geometry_check.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ Geometry visualization saved as 'geometry_check.png'")
    print(f"📐 Structure: {num_teeth} teeth per bar, {bar2_x - bar1_x} point gap")
    print(f"🔍 Grid size: {SIZE_X} x {SIZE_Y} points")
    print(f"📡 Sources: 5 laser sources in the gap between teeth")

if __name__ == "__main__":
    print("🔧 Checking Bar and Teeth Geometry...")
    visualize_bar_geometry()
    
    # Calculate source positions (between teeth in gap)
    sources = []
    
    for i in range(5):  # 5 sources between 6 teeth
        x = gap_center_x
        y = first_tooth_y + int((i + 0.5) * tooth_spacing)
        sources.append((x, y))
        
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
