#!/usr/bin/env python3
"""
FDTD Laser Simulation Visualization

This script creates comprehensive visualizations of the FDTD simulation results,
showing the electric field evolution with all laser sources clearly marked.

Features:
- Shows all 8 inter-pillar sources in every plot
- Displays pillar geometry
- Animates field evolution
- Creates individual timestep plots
- Saves high-quality figures
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import glob
import os
import re
from pathlib import Path

class FDTDVisualizerWithSources:
    def __init__(self, data_directory="."):
        """Initialize the visualizer with simulation data directory."""
        self.data_dir = Path(data_directory)
        self.geometry_params = {}
        self.field_files = []
        self.field_data = {}
        
        # Load geometry parameters
        self.load_geometry_params()
        
        # Find and load field data files
        self.load_field_data()
        
        # Calculate source positions
        self.calculate_source_positions()
        
    def load_geometry_params(self):
        """Load geometry parameters from the parameter file."""
        param_file = self.data_dir / "geometry_params.txt"
        
        if param_file.exists():
            with open(param_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        key = parts[0]
                        try:
                            value = float(parts[1])
                            self.geometry_params[key] = value
                        except ValueError:
                            self.geometry_params[key] = parts[1]
            print(f"Loaded {len(self.geometry_params)} geometry parameters")
        else:
            print("Warning: geometry_params.txt not found. Using default values.")
            self._set_default_params()
    
    def _set_default_params(self):
        """Set default parameters if geometry file is not found."""
        # From your C++ code constants
        self.geometry_params = {
            'SIZE_X': 300,
            'SIZE_Y': 300,
            'DX': 3e-8,  # 30 nm
            'pillar_radius': 20,
            'pillar_spacing_y': 50,
            'column1_x': 150,
            'column2_x': 175,
            'first_pillar_y': 75,
            'num_pillars_per_column': 5,
            'PML_WIDTH': 10
        }
    
    def calculate_source_positions(self):
        """Calculate all source positions based on geometry parameters."""
        # Extract parameters
        column1_x = int(self.geometry_params.get('column1_x', 150))
        column2_x = int(self.geometry_params.get('column2_x', 175))
        first_pillar_y = int(self.geometry_params.get('first_pillar_y', 75))
        pillar_spacing_y = int(self.geometry_params.get('pillar_spacing_y', 50))
        
        # Inter-pillar sources in column 1 - exactly between pillars
        self.sources_column1 = []
        for i in range(4):  # 4 inter-pillar sources per column
            x = column1_x  # Centered on pillar column
            y = first_pillar_y + int((i + 0.5) * pillar_spacing_y)  # Exactly between pillars
            self.sources_column1.append((x, y))
        
        # Inter-pillar sources in column 2 - exactly between pillars
        self.sources_column2 = []
        for i in range(4):  # 4 inter-pillar sources per column
            x = column2_x  # Centered on pillar column
            y = first_pillar_y + int((i + 0.5) * pillar_spacing_y)  # Exactly between pillars
            self.sources_column2.append((x, y))
        
        # All sources combined
        self.all_sources = self.sources_column1 + self.sources_column2
        
        print(f"Column 1 sources: {self.sources_column1}")
        print(f"Column 2 sources: {self.sources_column2}")
        print(f"Total sources: {len(self.all_sources)}")

    def load_field_data(self):
        """Load all field data files."""
        # Find all Ey field files
        field_pattern = self.data_dir / "Ey_field_*.dat"
        self.field_files = sorted(glob.glob(str(field_pattern)), 
                                key=lambda x: int(re.search(r'(\d+)', x).group()))
        
        print(f"Found {len(self.field_files)} field data files")
        
        # Load each field file
        for file_path in self.field_files:
            timestep = int(re.search(r'(\d+)', Path(file_path).name).group())
            try:
                data = np.loadtxt(file_path)
                self.field_data[timestep] = data
                print(f"Loaded timestep {timestep}: {data.shape}")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
    
    def get_pillar_positions(self):
        """Calculate pillar center positions."""
        column1_x = int(self.geometry_params.get('column1_x', 150))
        column2_x = int(self.geometry_params.get('column2_x', 175))
        first_pillar_y = int(self.geometry_params.get('first_pillar_y', 75))
        pillar_spacing_y = int(self.geometry_params.get('pillar_spacing_y', 50))
        num_pillars = int(self.geometry_params.get('num_pillars_per_column', 5))
        
        pillars = []
        for i in range(num_pillars):
            y = first_pillar_y + i * pillar_spacing_y
            pillars.append((column1_x, y))  # Column 1
            pillars.append((column2_x, y))  # Column 2
        
        return pillars
    
    def create_base_plot(self, ax, show_geometry=True):
        """Create base plot with geometry and sources."""
        if show_geometry:
            # Draw pillars
            pillars = self.get_pillar_positions()
            pillar_radius = self.geometry_params.get('pillar_radius', 20)
            
            for x, y in pillars:
                circle = Circle((x, y), pillar_radius, 
                              fill=False, edgecolor='black', linewidth=2, alpha=0.7)
                ax.add_patch(circle)
            
            # Draw PML boundaries
            pml_width = self.geometry_params.get('PML_WIDTH', 10)
            size_x = self.geometry_params.get('SIZE_X', 300)
            size_y = self.geometry_params.get('SIZE_Y', 300)
            
            # PML boundary rectangles
            ax.axvline(x=pml_width, color='gray', linestyle='--', alpha=0.5, label='PML boundary')
            ax.axvline(x=size_x-pml_width, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=pml_width, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=size_y-pml_width, color='gray', linestyle='--', alpha=0.5)
        
        # Draw all sources with different colors and markers
        # Column 1 sources (frequency f1) - Red circles
        for i, (x, y) in enumerate(self.sources_column1):
            ax.plot(x, y, 'ro', markersize=8, markeredgecolor='darkred', 
                   markeredgewidth=2, label='Col1 Sources (f1)' if i == 0 else '')
            ax.annotate(f'S1-{i+1}', (x, y), xytext=(5, 5), 
                       textcoords='offset points', fontsize=8, color='darkred')
        
        # Column 2 sources (frequency f2) - Blue squares
        for i, (x, y) in enumerate(self.sources_column2):
            ax.plot(x, y, 'bs', markersize=8, markeredgecolor='darkblue', 
                   markeredgewidth=2, label='Col2 Sources (f2)' if i == 0 else '')
            ax.annotate(f'S2-{i+1}', (x, y), xytext=(5, 5), 
                       textcoords='offset points', fontsize=8, color='darkblue')
    
    def plot_single_timestep(self, timestep, save_figure=True, show_plot=True):
        """Plot a single timestep with all sources marked."""
        if timestep not in self.field_data:
            print(f"No data for timestep {timestep}")
            return None
        
        data = self.field_data[timestep]
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot field data
        im = ax.imshow(data.T, extent=[0, data.shape[0], 0, data.shape[1]], 
                      origin='lower', cmap='RdBu_r', aspect='equal')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources
        self.create_base_plot(ax)

        # Labels and title
        ax.set_xlabel('X (grid points)', fontsize=12)
        ax.set_ylabel('Y (grid points)', fontsize=12)
        ax.set_title(f'FDTD Simulation - Timestep {timestep}\n'
                    f'Inter-pillar Laser Sources (8 total)', fontsize=14)
        
        # Add legend
        ax.legend(loc='upper right', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        plt.tight_layout()
        
        if save_figure:
            output_file = self.data_dir / f"field_plot_t{timestep:04d}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return fig
    
    def add_parameter_text(self, ax):
        """Add key parameters as text on the plot."""
        # Get key parameters
        lambda1_nm = self.geometry_params.get('lambda1_nm', 1500)
        lambda2_nm = self.geometry_params.get('lambda2_nm', 1364)
        intensity_gv = self.geometry_params.get('inter_pillar_intensity_GVm', 2.0)
        pillar_radius_um = self.geometry_params.get('pillar_radius_um', 0.6)
        
        param_text = f'λ₁ = {lambda1_nm:.0f} nm\n'
        param_text += f'λ₂ = {lambda2_nm:.0f} nm\n'
        param_text += f'Intensity = {intensity_gv:.1f} GV/m\n'
        param_text += f'Pillar R = {pillar_radius_um:.2f} μm'
        
        ax.text(0.02, 0.98, param_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    def plot_all_timesteps(self):
        """Plot all available timesteps."""
        print(f"Plotting {len(self.field_data)} timesteps...")
        
        for timestep in sorted(self.field_data.keys()):
            self.plot_single_timestep(timestep, save_figure=True, show_plot=False)
        
        print("All plots saved!")
    
    def create_animation(self, interval=400, save_animation=True):
        """Create a slow animation of the field evolution and save as GIF."""
        if not self.field_data:
            print("No field data available for animation")
            return None
        
        timesteps = sorted(self.field_data.keys())
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Get field data range for consistent colormap
        all_data = np.concatenate([data.flatten() for data in self.field_data.values()])
        vmin, vmax = np.percentile(all_data, [5, 95])  # Use 5-95 percentile for better contrast
        
        # Initialize with first frame
        first_data = self.field_data[timesteps[0]]
        im = ax.imshow(first_data.T, extent=[0, first_data.shape[0], 0, first_data.shape[1]], 
                      origin='lower', cmap='RdBu_r', aspect='equal', vmin=vmin, vmax=vmax)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources (they don't change)
        self.create_base_plot(ax)
        
        # Labels
        ax.set_xlabel('X (grid points)', fontsize=12)
        ax.set_ylabel('Y (grid points)', fontsize=12)
        title = ax.set_title('FDTD Simulation Animation\nInter-pillar Laser Sources', fontsize=14)
        
        # Add legend
        ax.legend(loc='upper right', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        # Prepare animation elements
        def animate(frame):
            timestep = timesteps[frame]
            data = self.field_data[timestep]
            im.set_array(data.T)
            title.set_text(f'FDTD Simulation - Timestep {timestep}\n'
                          f'Inter-pillar Laser Sources (8 total)')
            
            return [im, title]
        
        # Create animation with a longer interval for a slower speed
        anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                                     interval=interval, blit=False, repeat=False)
        
        if save_animation:
            output_file = self.data_dir / "field_animation.gif"
            print(f"Saving animation to {output_file}...")
            # Use a low fps for a slower GIF
            anim.save(str(output_file), writer='pillow', fps=2, dpi=120)
            print(f"Animation saved: {output_file}")

        plt.tight_layout()
        plt.show()
        
        return anim
    
    def plot_source_overview(self):
        """Create an overview plot showing just the geometry and sources."""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Set up empty field for background
        size_x = int(self.geometry_params.get('SIZE_X', 300))
        size_y = int(self.geometry_params.get('SIZE_Y', 300))
        
        # Create empty background
        ax.set_xlim(0, size_x)
        ax.set_ylim(0, size_y)
        ax.set_aspect('equal')
        
        # Add geometry and sources
        self.create_base_plot(ax)
        
        # Labels and title
        ax.set_xlabel('X (grid points)', fontsize=12)
        ax.set_ylabel('Y (grid points)', fontsize=12)
        ax.set_title('FDTD Simulation Setup\nPillar Geometry and Laser Sources', fontsize=14)
        
        # Add legend
        ax.legend(loc='upper right', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save the overview
        output_file = self.data_dir / "simulation_overview.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Overview saved: {output_file}")
        
        plt.show()
        
        return fig

def main():
    """Main function to run the visualization."""
    print("FDTD Laser Simulation Visualizer")
    print("=" * 35)
    
    # Initialize visualizer
    visualizer = FDTDVisualizerWithSources(".")
    
    if not visualizer.field_data:
        print("No field data found. Please run the C++ simulation first.")
        return
    
    # Create overview plot
    print("\n1. Creating simulation overview...")
    visualizer.plot_source_overview()
    
    # Plot a few key timesteps
    print("\n2. Plotting key timesteps...")
    timesteps = sorted(visualizer.field_data.keys())
    key_timesteps = [timesteps[0], timesteps[len(timesteps)//3], 
                    timesteps[2*len(timesteps)//3], timesteps[-1]]
    
    for ts in key_timesteps:
        if ts in visualizer.field_data:
            visualizer.plot_single_timestep(ts, save_figure=True, show_plot=True)
    
    # Ask user about full processing
    response = input("\n3. Do you want to plot ALL timesteps? (y/n): ").strip().lower()
    if response == 'y':
        visualizer.plot_all_timesteps()
    
    # Ask about animation
    response = input("\n4. Do you want to create an animation? (y/n): ").strip().lower()
    if response == 'y':
        visualizer.create_animation(interval=100) # Speed up animation a bit
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    main()
