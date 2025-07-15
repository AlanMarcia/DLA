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

def load_geometry_params():
    """Load geometry parameters from the C++ simulation output file."""
    params = {}
    try:
        with open('geometry_params.txt', 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    key = parts[0]
                    try:
                        value = float(parts[1])
                        params[key] = value
                    except ValueError:
                        params[key] = parts[1]
        print(f"✅ Loaded {len(params)} geometry parameters from file")
        return params
    except FileNotFoundError:
        print("⚠️ Warning: geometry_params.txt not found. Using default values.")
        return None

def visualize_bar_geometry():
    """Create a visualization of the bar and teeth geometry."""
    
    # Load parameters from file first
    loaded_params = load_geometry_params()
    
    if loaded_params:
        # Use parameters from the C++ simulation
        SIZE_X = int(loaded_params.get('SIZE_X', 600))
        SIZE_Y = int(loaded_params.get('SIZE_Y', 600))
        lambda1 = loaded_params.get('lambda1', 2e-6)
        DX = loaded_params.get('DX', 20e-9)
        PML_WIDTH = int(loaded_params.get('PML_WIDTH', 10))
        
        # Bar parameters (from SI units in the file)
        bar_width = int(loaded_params.get('bar_width', 15))
        bar1_x = int(loaded_params.get('bar1_x', 200))
        bar2_x = int(loaded_params.get('bar2_x', 400))
        bar_start_y = int(loaded_params.get('bar_start_y', 50))
        bar_end_y = int(loaded_params.get('bar_end_y', 550))
        
        # Teeth parameters (wavelength-based)
        num_teeth = int(loaded_params.get('num_teeth', 10))
        tooth_width = int(loaded_params.get('tooth_width', 30))
        tooth_height = int(loaded_params.get('tooth_height', 40))
        tooth_spacing = int(loaded_params.get('tooth_spacing', 50))
        first_tooth_y = int(loaded_params.get('first_tooth_y', 100))
        gap_center_x = int(loaded_params.get('gap_center_x', 300))
        
        print("📁 Using parameters from geometry_params.txt")
        print(f"Grid spacing: {DX*1e9:.1f} nm")
        print(f"Bar width: {loaded_params.get('bar_width_nm', 'N/A')} nm")
        print(f"Gap width: {loaded_params.get('gap_between_bars_um', 'N/A')} μm")
    else:
        # Fallback to hardcoded parameters from C++ code (approximate values)
        lambda1 = 2e-6  # wavelength
        DX = 20e-9  # grid spacing = 20 nm
        
        # Convert to grid units
        SIZE_X = 600
        SIZE_Y = 600
        PML_WIDTH = 10
        
        # Bar parameters (SI units converted to grid points)
        bar_width_m = 300e-9  # 300 nm
        bar1_x_m = 4e-6  # 4 μm
        bar2_x_m = 8e-6  # 8 μm
        bar_start_y_m = 1e-6  # 1 μm
        bar_end_y_m = 11e-6  # 11 μm
        
        bar_width = int(bar_width_m / DX)
        bar1_x = int(bar1_x_m / DX)
        bar2_x = int(bar2_x_m / DX)
        bar_start_y = int(bar_start_y_m / DX)
        bar_end_y = int(bar_end_y_m / DX)
        
        # Teeth parameters (wavelength-based)
        num_teeth = 10
        tooth_width = int(0.3 * lambda1 / DX)
        tooth_height = int(0.4 * lambda1 / DX)
        tooth_spacing = int(0.5 * lambda1 / DX)
        first_tooth_y = int(2e-6 / DX)  # 2 μm in SI units
        gap_center_x = (bar1_x + bar2_x) // 2
        
        print("🔧 Using fallback hardcoded parameters")
    
    # Bar length (same for both cases)
    # Note: bar_start_y and bar_end_y are already set above in both branches
    
    print("FDTD Laser Simulation - Bar and Teeth Geometry Check")
    print("=" * 55)
    print(f"Grid size: {SIZE_X} x {SIZE_Y} points")
    print(f"Domain size: {SIZE_X * DX * 1e6:.1f} x {SIZE_Y * DX * 1e6:.1f} μm")
    print(f"Grid spacing: {DX * 1e9:.1f} nm")
    print(f"Bar width: {bar_width} grid points ({bar_width * DX * 1e9:.0f} nm)")
    print(f"Bar 1 center at x = {bar1_x} pts ({bar1_x * DX * 1e6:.1f} μm)")
    print(f"Bar 2 center at x = {bar2_x} pts ({bar2_x * DX * 1e6:.1f} μm)")
    print(f"Gap width: {bar2_x - bar1_x} grid points ({(bar2_x - bar1_x) * DX * 1e6:.1f} μm)")
    print(f"Number of teeth per bar: {num_teeth}")
    print(f"Tooth dimensions: {tooth_height} x {tooth_width} pts ({tooth_height * DX * 1e9:.0f} x {tooth_width * DX * 1e9:.0f} nm)")
    print(f"Tooth spacing: {tooth_spacing} grid points ({tooth_spacing * DX * 1e9:.0f} nm)")
    
    # Calculate external source positions (left and right sides, between teeth)
    left_sources = []
    right_sources = []
    
    # Load external source positions from parameters or use defaults
    if loaded_params:
        left_source_x = int(loaded_params.get('left_source_x', 10))
        right_source_x = int(loaded_params.get('right_source_x', SIZE_X - 10))
    else:
        left_source_x = int(200e-9 / DX)  # 200 nm from left edge
        right_source_x = int(11.8e-6 / DX)  # 11.8 μm from left edge
    
    print(f"\nExternal Source Positions:")
    print("=" * 35)
    print("Left Side Sources:")
    
    for i in range(9):  # 5 sources on left side, between teeth
        y = first_tooth_y + int((i + 0.5) * tooth_spacing)
        left_sources.append((left_source_x, y))
        print(f"Left S{i+1}: x={left_source_x} pts ({left_source_x * DX * 1e6:.1f} μm), y={y} pts ({y * DX * 1e6:.1f} μm)")
    
    print("Right Side Sources:")
    for i in range(9):  # 5 sources on right side, between teeth
        y = first_tooth_y + int((i + 0.5) * tooth_spacing)
        right_sources.append((right_source_x, y))
        print(f"Right S{i+1}: x={right_source_x} pts ({right_source_x * DX * 1e6:.1f} μm), y={y} pts ({y * DX * 1e6:.1f} μm)")
    
    # Create figure with physical units
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Convert all coordinates to micrometers for plotting
    def to_um(grid_points):
        return grid_points * DX * 1e6
    
    # Draw PML boundaries
    pml_color = 'lightgray'
    ax.axvline(x=to_um(PML_WIDTH), color=pml_color, linestyle='--', alpha=0.7, label='PML Boundary')
    ax.axvline(x=to_um(SIZE_X-PML_WIDTH), color=pml_color, linestyle='--', alpha=0.7)
    ax.axhline(y=to_um(PML_WIDTH), color=pml_color, linestyle='--', alpha=0.7)
    ax.axhline(y=to_um(SIZE_Y-PML_WIDTH), color=pml_color, linestyle='--', alpha=0.7)
    
    # Draw Bar 1 (left bar)
    bar1_rect = Rectangle((to_um(bar1_x - bar_width//2), to_um(bar_start_y)), 
                         to_um(bar_width), to_um(bar_end_y - bar_start_y),
                         facecolor='lightblue', edgecolor='blue', 
                         alpha=0.7, label='Dielectric Bars')
    ax.add_patch(bar1_rect)
    
    # Draw Bar 2 (right bar)
    bar2_rect = Rectangle((to_um(bar2_x - bar_width//2), to_um(bar_start_y)), 
                         to_um(bar_width), to_um(bar_end_y - bar_start_y),
                         facecolor='lightblue', edgecolor='blue', alpha=0.7)
    ax.add_patch(bar2_rect)
    
    # Draw teeth extending from both bars
    for t in range(num_teeth):
        tooth_center_y = first_tooth_y + t * tooth_spacing
        
        # Tooth from Bar 1 (extending right)
        tooth1_rect = Rectangle((to_um(bar1_x + bar_width//2), to_um(tooth_center_y - tooth_width//2)),
                               to_um(tooth_height), to_um(tooth_width),
                               facecolor='lightblue', edgecolor='blue', alpha=0.7)
        ax.add_patch(tooth1_rect)
        
        # Tooth from Bar 2 (extending left)
        tooth2_rect = Rectangle((to_um(bar2_x - bar_width//2 - tooth_height), to_um(tooth_center_y - tooth_width//2)),
                               to_um(tooth_height), to_um(tooth_width),
                               facecolor='lightblue', edgecolor='blue', alpha=0.7)
        ax.add_patch(tooth2_rect)
        
        # Label teeth
        ax.text(to_um(bar1_x + bar_width//2 + tooth_height//2), to_um(tooth_center_y), f'T{t+1}', 
               ha='center', va='center', fontsize=8, color='darkblue')
    
    # Draw external laser sources
    # Left side sources
    for i, (x, y) in enumerate(left_sources):
        ax.plot(to_um(x), to_um(y), 'ro', markersize=10, 
               markeredgecolor='darkred', markeredgewidth=2,
               label='Left Sources' if i == 0 else '')
        ax.text(to_um(x) - 0.3, to_um(y), f'L{i+1}', 
               fontsize=8, color='darkred', ha='right')
    
    # Right side sources
    for i, (x, y) in enumerate(right_sources):
        ax.plot(to_um(x), to_um(y), 'go', markersize=10, 
               markeredgecolor='darkgreen', markeredgewidth=2,
               label='Right Sources' if i == 0 else '')
        ax.text(to_um(x) + 0.3, to_um(y), f'R{i+1}', 
               fontsize=8, color='darkgreen', ha='left')
    
    # Add dimensions and annotations
    ax.annotate('', xy=(to_um(bar1_x - bar_width//2), to_um(SIZE_Y - 20)), 
               xytext=(to_um(bar1_x + bar_width//2), to_um(SIZE_Y - 20)),
               arrowprops=dict(arrowstyle='<->', color='black'))
    ax.text(to_um(bar1_x), to_um(SIZE_Y - 15), f'Bar Width\n{bar_width * DX * 1e9:.0f} nm', 
           ha='center', va='bottom', fontsize=8)
    
    ax.annotate('', xy=(to_um(bar1_x + bar_width//2), to_um(first_tooth_y)), 
               xytext=(to_um(bar1_x + bar_width//2 + tooth_height), to_um(first_tooth_y)),
               arrowprops=dict(arrowstyle='<->', color='green'))
    ax.text(to_um(bar1_x + bar_width//2 + tooth_height//2), to_um(first_tooth_y) - 0.2, 
           f'Tooth Height\n{tooth_height * DX * 1e9:.0f} nm', ha='center', va='top', fontsize=8)
    
    ax.annotate('', xy=(to_um(bar1_x), to_um(bar_end_y) + 0.1), 
               xytext=(to_um(bar2_x), to_um(bar_end_y) + 0.1),
               arrowprops=dict(arrowstyle='<->', color='purple'))
    ax.text(to_um((bar1_x + bar2_x)//2), to_um(bar_end_y) + 0.3, 
           f'Gap: {(bar2_x - bar1_x) * DX * 1e6:.1f} μm', ha='center', va='bottom', fontsize=8)
    
    # Labels and formatting
    ax.set_xlim(0, to_um(SIZE_X))
    ax.set_ylim(0, to_um(SIZE_Y))
    ax.set_xlabel('X (μm)', fontsize=12)
    ax.set_ylabel('Y (μm)', fontsize=12)
    ax.set_title('Dielectric Laser Accelerator - Bar and Teeth Geometry', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Add parameter text with loaded values
    if loaded_params:
        param_text = f'Parameters (from file):\n'
        param_text += f'• Domain: {SIZE_X * DX * 1e6:.1f} × {SIZE_Y * DX * 1e6:.1f} μm\n'
        param_text += f'• λ: {loaded_params.get("lambda1_nm", "N/A"):.0f} nm\n'
        param_text += f'• Grid: {DX * 1e9:.1f} nm spacing\n'
        param_text += f'• {num_teeth} teeth per bar\n'
        param_text += f'• Gap: {(bar2_x - bar1_x) * DX * 1e6:.1f} μm\n'
        param_text += f'• 10 external sources (5L + 5R)'
    else:
        param_text = f'Parameters (hardcoded):\n'
        param_text += f'• Domain: {SIZE_X * DX * 1e6:.1f} × {SIZE_Y * DX * 1e6:.1f} μm\n'
        param_text += f'• {num_teeth} teeth per bar\n'
        param_text += f'• Tooth spacing: {tooth_spacing * DX * 1e9:.0f} nm\n'
        param_text += f'• Bar width: {bar_width * DX * 1e9:.0f} nm\n'
        param_text += f'• Tooth height: {tooth_height * DX * 1e9:.0f} nm\n'
        param_text += f'• 10 external sources (5L + 5R)'
    
    ax.text(0.02, 0.98, param_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('geometry_check.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✅ Geometry visualization saved as 'geometry_check.png'")
    print(f"📐 Structure: {num_teeth} teeth per bar, {(bar2_x - bar1_x) * DX * 1e6:.1f} μm gap")
    print(f"🔍 Domain size: {SIZE_X * DX * 1e6:.1f} x {SIZE_Y * DX * 1e6:.1f} μm")
    print(f"📡 Sources: 10 external laser sources (5 left + 5 right) positioned between teeth")

if __name__ == "__main__":
    print("🔧 Checking Bar and Teeth Geometry...")
    visualize_bar_geometry()
