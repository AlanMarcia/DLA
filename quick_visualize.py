#!/usr/bin/env python3
"""
Quick FDTD Results Viewer with Source Overlay

This script quickly displays FDTD simulation results with all laser sources marked.
Run this after your C++ simulation to see the results.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import glob
import os
import re

def load_geometry_params(filename="geometry_params.txt"):
    """Load geometry parameters from file."""
    params = {}
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        params[parts[0]] = float(parts[1])
                    except ValueError:
                        params[parts[0]] = parts[1]
    return params

def calculate_source_positions(params):
    """Calculate all 8 inter-pillar source positions."""
    # Get parameters
    column1_x = int(params.get('column1_x', 150))
    column2_x = int(params.get('column2_x', 175))
    first_pillar_y = int(params.get('first_pillar_y', 75))
    pillar_spacing_y = int(params.get('pillar_spacing_y', 50))
    
    # Sources are now centered on column x-positions and exactly between pillars
    # Column 1 sources (f1 = 1500 nm)
    sources_col1 = []
    for i in range(4):
        x = column1_x  # Centered on pillar column
        y = first_pillar_y + int((i + 0.5) * pillar_spacing_y)  # Exactly between pillars
        sources_col1.append((x, y))
    
    # Column 2 sources (f2 = 1364 nm)
    sources_col2 = []
    for i in range(4):
        x = column2_x  # Centered on pillar column
        y = first_pillar_y + int((i + 0.5) * pillar_spacing_y)  # Exactly between pillars
        sources_col2.append((x, y))
    
    return sources_col1, sources_col2

def plot_field_with_sources(field_data, timestep, params, sources_col1, sources_col2):
    """Plot field data with all sources and geometry."""
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Plot the field
    im = ax.imshow(field_data.T, extent=[0, field_data.shape[0], 0, field_data.shape[1]], 
                  origin='lower', cmap='RdBu_r', aspect='equal')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Electric Field Ey (V/m)', fontsize=12)
    
    # Draw pillars
    pillar_radius = params.get('pillar_radius', 20)
    column1_x = int(params.get('column1_x', 150))
    column2_x = int(params.get('column2_x', 175))
    first_pillar_y = int(params.get('first_pillar_y', 75))
    pillar_spacing_y = int(params.get('pillar_spacing_y', 50))
    num_pillars = int(params.get('num_pillars_per_column', 5))
    
    for i in range(num_pillars):
        y = first_pillar_y + i * pillar_spacing_y
        # Column 1 pillars
        circle1 = Circle((column1_x, y), pillar_radius, 
                        fill=False, edgecolor='black', linewidth=2, alpha=0.7)
        ax.add_patch(circle1)
        # Column 2 pillars
        circle2 = Circle((column2_x, y), pillar_radius, 
                        fill=False, edgecolor='black', linewidth=2, alpha=0.7)
        ax.add_patch(circle2)
    
    # Draw PML boundaries
    pml_width = params.get('PML_WIDTH', 10)
    size_x = params.get('SIZE_X', 300)
    size_y = params.get('SIZE_Y', 300)
    
    ax.axvline(x=pml_width, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=size_x-pml_width, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=pml_width, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=size_y-pml_width, color='gray', linestyle='--', alpha=0.5)
    
    # Draw ALL sources with clear labels
    # Column 1 sources (frequency f1 = 1500 nm) - Red circles
    for i, (x, y) in enumerate(sources_col1):
        ax.plot(x, y, 'ro', markersize=10, markeredgecolor='darkred', 
               markeredgewidth=2, label='Column 1 Sources (λ=1500nm)' if i == 0 else '')
        ax.annotate(f'S1-{i+1}', (x+3, y+3), fontsize=10, color='darkred', weight='bold')
    
    # Column 2 sources (frequency f2 = 1364 nm) - Blue squares
    for i, (x, y) in enumerate(sources_col2):
        ax.plot(x, y, 'bs', markersize=10, markeredgecolor='darkblue', 
               markeredgewidth=2, label='Column 2 Sources (λ=1364nm)' if i == 0 else '')
        ax.annotate(f'S2-{i+1}', (x+3, y+3), fontsize=10, color='darkblue', weight='bold')
    
    # Labels and title
    ax.set_xlabel('X (grid points)', fontsize=12)
    ax.set_ylabel('Y (grid points)', fontsize=12)
    ax.set_title(f'FDTD Simulation - Timestep {timestep}\n'
                f'8 Inter-pillar Laser Sources (2 GV/m each)', fontsize=14, weight='bold')
    
    # Add legend
    ax.legend(loc='upper right', fontsize=11)
    
    # Add parameter information
    lambda1_nm = params.get('lambda1_nm', 1500)
    lambda2_nm = params.get('lambda2_nm', 1364)
    intensity_gv = params.get('inter_pillar_intensity_GVm', 2.0)
    
    info_text = f'Wavelengths: λ₁={lambda1_nm:.0f}nm, λ₂={lambda2_nm:.0f}nm\n'
    info_text += f'Source Intensity: {intensity_gv:.1f} GV/m\n'
    info_text += f'Total Sources: 8 (4 per column)\n'
    info_text += f'Pillar Material: SiO₂ (εᵣ=3.9)'
    
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    plt.tight_layout()
    return fig

def main():
    """Main visualization function."""
    print("FDTD Laser Simulation Results Viewer")
    print("=" * 40)
    
    # Load parameters
    params = load_geometry_params()
    if not params:
        print("Warning: No geometry_params.txt found. Using default values.")
        params = {
            'SIZE_X': 300, 'SIZE_Y': 300, 'column1_x': 150, 'column2_x': 175,
            'first_pillar_y': 75, 'pillar_spacing_y': 50, 'pillar_radius': 20,
            'lambda1': 1.5e-6, 'DX': 3e-8, 'lambda1_nm': 1500, 'lambda2_nm': 1364,
            'inter_pillar_intensity_GVm': 2.0, 'PML_WIDTH': 10, 'num_pillars_per_column': 5
        }
    
    # Calculate source positions
    sources_col1, sources_col2 = calculate_source_positions(params)
    print(f"Column 1 sources (λ=1500nm): {sources_col1}")
    print(f"Column 2 sources (λ=1364nm): {sources_col2}")
    
    # Find field data files
    field_files = sorted(glob.glob("Ey_field_*.dat"), 
                        key=lambda x: int(re.search(r'(\d+)', x).group()))
    
    if not field_files:
        print("No field data files found. Please run the C++ simulation first.")
        print("Expected files: Ey_field_0.dat, Ey_field_50.dat, etc.")
        return
    
    print(f"Found {len(field_files)} field data files")
    
    # Show first few timesteps
    show_count = min(4, len(field_files))
    for i in range(show_count):
        file_path = field_files[i]
        timestep = int(re.search(r'(\d+)', file_path).group())
        
        try:
            # Load field data
            field_data = np.loadtxt(file_path)
            print(f"Plotting timestep {timestep}...")
            
            # Create plot
            fig = plot_field_with_sources(field_data, timestep, params, sources_col1, sources_col2)
            
            # Save the plot
            output_name = f"field_with_sources_t{timestep:04d}.png"
            plt.savefig(output_name, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_name}")
            
            # Show the plot
            plt.show()
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    # Ask if user wants to process all files
    if len(field_files) > show_count:
        response = input(f"\nProcess all {len(field_files)} timesteps? (y/n): ").strip().lower()
        if response == 'y':
            for file_path in field_files[show_count:]:
                timestep = int(re.search(r'(\d+)', file_path).group())
                try:
                    field_data = np.loadtxt(file_path)
                    fig = plot_field_with_sources(field_data, timestep, params, sources_col1, sources_col2)
                    
                    output_name = f"field_with_sources_t{timestep:04d}.png"
                    plt.savefig(output_name, dpi=300, bbox_inches='tight')
                    plt.close()  # Close to save memory
                    print(f"Processed timestep {timestep}")
                    
                except Exception as e:
                    print(f"Error processing timestep {timestep}: {e}")
    
    print("\nVisualization complete!")
    print("All plots show:")
    print("  • Red circles: Column 1 sources (λ=1500nm)")
    print("  • Blue squares: Column 2 sources (λ=1364nm)")
    print("  • Black circles: SiO₂ pillars")
    print("  • Gray dashed lines: PML boundaries")

if __name__ == "__main__":
    main()
