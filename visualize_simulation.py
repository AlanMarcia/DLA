#!/usr/bin/env python3
"""
FDTD Laser Simulation Visualization

This script creates comprehensive visualizations of the FDTD simulation results,
showing the electric field evolution with all laser sources clearly marked.

Features:
- Shows all 5 inter-teeth sources in every plot
- Displays bar and teeth geometry
- Animates field evolution
- Creates individual timestep plots
- Saves high-quality figures
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import glob
import os
import re
from pathlib import Path

# Optional GUI imports
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    print("Warning: tkinter not available. GUI folder selection disabled.")

class FDTDVisualizerWithSources:
    def __init__(self, data_directory=".", output_directory=None):
        """Initialize the visualizer with simulation data directory.
        
        Parameters:
        -----------
        data_directory : str
            Directory containing simulation data files
        output_directory : str, optional
            Directory where output files will be saved. If None, uses data_directory
        """
        self.data_dir = Path(data_directory)
        self.output_dir = Path(output_directory) if output_directory else self.data_dir
        self.geometry_params = {}
        self.field_files = []
        self.field_data = {}
        self.potential_files = []
        self.potential_data = {}
        
        # Load geometry parameters
        self.load_geometry_params()
        
        # Find and load field data files
        self.load_field_data()
        
        # Find and load potential data files
        self.load_potential_data()
        
        # Calculate source positions
        self.calculate_source_positions()
        
        # Create output directory if it doesn't exist
        self.ensure_output_directory()
        
    def ensure_output_directory(self):
        """Create output directory if it doesn't exist."""
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory: {self.output_dir}")
        
    def load_geometry_params(self):
        """Load geometry parameters from the parameter file."""
        param_file = self.data_dir / "geometry_params.txt"
        
        print(f"Looking for geometry parameters in: {param_file}")
        
        if param_file.exists():
            print(f"Found geometry_params.txt, loading parameters...")
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
            
            # Print key parameters for verification
            dx_value = self.geometry_params.get('DX', 'NOT_FOUND')
            print(f"DX from file: {dx_value}")
            grid_origin = self.geometry_params.get('grid_origin_x_m', 'NOT_FOUND')
            print(f"Grid origin from file: {grid_origin}")
        else:
            print("Warning: geometry_params.txt not found. Using default values.")
            print("Make sure to run the C++ simulation first to generate geometry_params.txt")
            self._set_default_params()
    
    def _set_default_params(self):
        """Set default parameters if geometry file is not found."""
        # Updated from actual C++ simulation with finer grid (10 nm)
        self.geometry_params = {
            'SIZE_X': 570,  # Updated to match actual simulation
            'SIZE_Y': 600,  # Updated to match actual simulation
            'lambda1': 2e-6,
            'DX': 10e-9,  # Fine grid (10 nm)
            'bar_width': 50,  # Updated to match actual simulation (50 grid points = 500 nm)
            'bar1_x': 185,  # Updated to match actual simulation
            'bar2_x': 384,  # Updated to match actual simulation
            'num_teeth': 33,  # Updated to match actual simulation (calculated from bar length)
            'tooth_width': 20,  # Updated for finer grid
            'tooth_height': 30,  # Updated for finer grid
            'tooth_spacing': 50,  # Updated for finer grid
            'first_tooth_y': 150,  # Updated for finer grid and larger domain
            'left_source_x': 45,   # Updated for new domain
            'right_source_x': 525, # Updated for new domain
            'PML_WIDTH': 10,
            'grid_origin_x_m': 3.15e-6  # Updated to match actual simulation
        }
    
    def calculate_source_positions(self):
        """Calculate all source positions based on geometry parameters."""
        # Get grid parameters for coordinate conversion
        self.DX = self.geometry_params.get('DX', 10e-9)  # Grid spacing in meters
        self.grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        
        # Get time parameters for timestep to time conversion
        self.DT = self.geometry_params.get('DT', 1e-15)  # Time step in seconds
        self.f1 = self.geometry_params.get('f1', 1.5e14)  # Frequency in Hz
        
        # Extract parameters for external sources (in grid points) - directly from file
        first_tooth_y = int(self.geometry_params.get('first_tooth_y', 150))
        tooth_spacing = int(self.geometry_params.get('tooth_spacing', 100))
        left_source_x = int(self.geometry_params.get('left_source_x', 134))
        right_source_x = int(self.geometry_params.get('right_source_x', 435))
        num_teeth = int(self.geometry_params.get('num_teeth', 10))  # Read dynamic number of teeth
        
        print(f"Using source parameters from file:")
        print(f"  first_tooth_y = {first_tooth_y} (grid points)")
        print(f"  tooth_spacing = {tooth_spacing} (grid points)")
        print(f"  left_source_x = {left_source_x} (grid points)")
        print(f"  right_source_x = {right_source_x} (grid points)")
        print(f"  num_teeth = {num_teeth}")
        
        # Left side sources - now including source before first tooth + num_teeth between teeth
        self.left_sources = []
        self.left_sources_physical = []  # Physical coordinates in micrometers
        
        # Source 0: before first tooth
        left_source0_y = int(self.geometry_params.get('left_source0_y', first_tooth_y - int(0.5 * tooth_spacing)))
        x_grid = left_source_x
        y_grid = left_source0_y
        
        # Convert to physical coordinates (micrometers)
        x_phys = (self.grid_origin_x + x_grid * self.DX) * 1e6  # μm
        y_phys = (y_grid * self.DX) * 1e6  # μm
        
        self.left_sources.append((x_grid, y_grid))  # Grid coordinates
        self.left_sources_physical.append((x_phys, y_phys))  # Physical coordinates
        
        # Sources 1-num_teeth: between teeth (dynamic based on num_teeth)
        for i in range(num_teeth):
            x_grid = left_source_x
            y_grid = first_tooth_y + int((i + 0.5) * tooth_spacing)
            
            # Convert to physical coordinates (micrometers)
            x_phys = (self.grid_origin_x + x_grid * self.DX) * 1e6  # μm
            y_phys = (y_grid * self.DX) * 1e6  # μm
            
            self.left_sources.append((x_grid, y_grid))  # Grid coordinates
            self.left_sources_physical.append((x_phys, y_phys))  # Physical coordinates
        
        # Right side sources - now including source before first tooth + num_teeth between teeth
        self.right_sources = []
        self.right_sources_physical = []  # Physical coordinates in micrometers
        
        # Source 0: before first tooth
        right_source0_y = int(self.geometry_params.get('right_source0_y', first_tooth_y - int(0.5 * tooth_spacing)))
        x_grid = right_source_x
        y_grid = right_source0_y
        
        # Convert to physical coordinates (micrometers)
        x_phys = (self.grid_origin_x + x_grid * self.DX) * 1e6  # μm
        y_phys = (y_grid * self.DX) * 1e6  # μm
        
        self.right_sources.append((x_grid, y_grid))  # Grid coordinates
        self.right_sources_physical.append((x_phys, y_phys))  # Physical coordinates
        
        # Sources 1-num_teeth: between teeth (dynamic based on num_teeth)
        for i in range(num_teeth):
            x_grid = right_source_x
            y_grid = first_tooth_y + int((i + 0.5) * tooth_spacing)
            
            # Convert to physical coordinates (micrometers)
            x_phys = (self.grid_origin_x + x_grid * self.DX) * 1e6  # μm
            y_phys = (y_grid * self.DX) * 1e6  # μm
            
            self.right_sources.append((x_grid, y_grid))  # Grid coordinates
            self.right_sources_physical.append((x_phys, y_phys))  # Physical coordinates
        
        # All sources combined
        self.all_sources = self.left_sources + self.right_sources
        self.all_sources_physical = self.left_sources_physical + self.right_sources_physical
        
        print(f"Left external sources (physical): {self.left_sources_physical}")
        print(f"Right external sources (physical): {self.right_sources_physical}")
        print(f"Total sources: {len(self.all_sources)}")

    def timestep_to_time_info(self, timestep):
        """Convert timestep to physical time information."""
        current_time_s = timestep * self.DT  # Physical time in seconds
        current_time_fs = current_time_s * 1e15  # Physical time in femtoseconds
        current_time_periods = current_time_s * self.f1  # Time in periods of wavelength 1
        
        return {
            'time_s': current_time_s,
            'time_fs': current_time_fs,
            'periods': current_time_periods
        }

    def format_time_string(self, timestep):
        """Format a nice time string for plot titles."""
        time_info = self.timestep_to_time_info(timestep)
        return f"t={timestep} ({time_info['time_fs']:.1f} fs, {time_info['periods']:.2f} periods)"

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
    
    def load_potential_data(self):
        """Load all potential data files."""
        # Find all potential files
        potential_pattern = self.data_dir / "potential_*.dat"
        self.potential_files = sorted(glob.glob(str(potential_pattern)), 
                                    key=lambda x: int(re.search(r'(\d+)', x).group()))
        
        print(f"Found {len(self.potential_files)} potential data files")
        
        # Load each potential file
        for file_path in self.potential_files:
            timestep = int(re.search(r'(\d+)', Path(file_path).name).group())
            try:
                data = np.loadtxt(file_path)
                self.potential_data[timestep] = data
                print(f"Loaded potential timestep {timestep}: {data.shape}")
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
    
    def get_bar_and_teeth_geometry(self):
        """Calculate bar and teeth geometry for visualization in physical units."""
        # Get grid parameters for coordinate conversion
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm grid spacing
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        
        # Extract parameters (grid points)
        bar1_x = int(self.geometry_params.get('bar1_x', 140))
        bar2_x = int(self.geometry_params.get('bar2_x', 185))
        bar_width = int(self.geometry_params.get('bar_width', 6))
        num_teeth = int(self.geometry_params.get('num_teeth', 6))
        tooth_width = int(self.geometry_params.get('tooth_width', 12))
        tooth_height = int(self.geometry_params.get('tooth_height', 20))
        tooth_spacing = int(self.geometry_params.get('tooth_spacing', 25))
        first_tooth_y = int(self.geometry_params.get('first_tooth_y', 75))
        PML_WIDTH = int(self.geometry_params.get('PML_WIDTH', 10))
        SIZE_Y = int(self.geometry_params.get('SIZE_Y', 300))
        
        # Bar dimensions (grid points)
        bar_start_y = PML_WIDTH + 10
        bar_end_y = SIZE_Y - PML_WIDTH - 10
        
        # Convert to physical coordinates (micrometers)
        geometry = {
            'bar1_x': (grid_origin_x + bar1_x * DX) * 1e6,  # μm
            'bar2_x': (grid_origin_x + bar2_x * DX) * 1e6,  # μm
            'bar_width': bar_width * DX * 1e6,  # μm
            'bar_start_y': bar_start_y * DX * 1e6,  # μm
            'bar_end_y': bar_end_y * DX * 1e6,  # μm
            'num_teeth': num_teeth,
            'tooth_width': tooth_width * DX * 1e6,  # μm
            'tooth_height': tooth_height * DX * 1e6,  # μm
            'tooth_spacing': tooth_spacing * DX * 1e6,  # μm
            'first_tooth_y': first_tooth_y * DX * 1e6  # μm
        }
        
        return geometry
    
    def create_base_plot(self, ax, show_geometry=True):
        """Create base plot with geometry and sources in physical units."""
        if show_geometry:
            # Get bar and teeth geometry (now in micrometers)
            geom = self.get_bar_and_teeth_geometry()
            
            # Draw Bar 1 (left bar)
            bar1_rect = Rectangle((geom['bar1_x'] - geom['bar_width']/2, geom['bar_start_y']), 
                                 geom['bar_width'], geom['bar_end_y'] - geom['bar_start_y'],
                                 facecolor='lightblue', edgecolor='blue', 
                                 alpha=0.7, linewidth=2, label='Dielectric Bars')
            ax.add_patch(bar1_rect)
            
            # Draw Bar 2 (right bar)
            bar2_rect = Rectangle((geom['bar2_x'] - geom['bar_width']/2, geom['bar_start_y']), 
                                 geom['bar_width'], geom['bar_end_y'] - geom['bar_start_y'],
                                 facecolor='lightblue', edgecolor='blue', alpha=0.7, linewidth=2)
            ax.add_patch(bar2_rect)
            
            # Draw teeth extending from both bars
            for t in range(geom['num_teeth']):
                tooth_center_y = geom['first_tooth_y'] + t * geom['tooth_spacing']
                
                # Tooth from Bar 1 (extending right)
                tooth1_rect = Rectangle((geom['bar1_x'] + geom['bar_width']/2, 
                                       tooth_center_y - geom['tooth_width']/2),
                                       geom['tooth_height'], geom['tooth_width'],
                                       facecolor='lightblue', edgecolor='blue', 
                                       alpha=0.7, linewidth=1)
                ax.add_patch(tooth1_rect)
                
                # Tooth from Bar 2 (extending left)
                tooth2_rect = Rectangle((geom['bar2_x'] - geom['bar_width']/2 - geom['tooth_height'], 
                                       tooth_center_y - geom['tooth_width']/2),
                                       geom['tooth_height'], geom['tooth_width'],
                                       facecolor='lightblue', edgecolor='blue', 
                                       alpha=0.7, linewidth=1)
                ax.add_patch(tooth2_rect)
                
                # Label teeth
                ax.text(geom['bar1_x'] + geom['bar_width']/2 + geom['tooth_height']/2, 
                       tooth_center_y, f'T{t+1}', 
                       ha='center', va='center', fontsize=7, color='darkblue')
            
            # Draw PML boundaries (convert to physical units)
            DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
            grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)
            pml_width = self.geometry_params.get('PML_WIDTH', 10)
            size_x = self.geometry_params.get('SIZE_X', 300)
            size_y = self.geometry_params.get('SIZE_Y', 300)
            
            pml_left = (grid_origin_x + pml_width * DX) * 1e6  # μm
            pml_right = (grid_origin_x + (size_x - pml_width) * DX) * 1e6  # μm
            pml_bottom = pml_width * DX * 1e6  # μm
            pml_top = (size_y - pml_width) * DX * 1e6  # μm
            
            ax.axvline(x=pml_left, color='gray', linestyle='--', alpha=0.5, label='PML boundary')
            ax.axvline(x=pml_right, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=pml_bottom, color='gray', linestyle='--', alpha=0.5)
            ax.axhline(y=pml_top, color='gray', linestyle='--', alpha=0.5)
        
        # Draw external sources (using physical coordinates)
        # Left side sources
        for i, (x, y) in enumerate(self.left_sources_physical):
            ax.plot(x, y, 'ro', markersize=10, markeredgecolor='darkred', 
                   markeredgewidth=2, label='Left Sources' if i == 0 else '')
            ax.annotate(f'L{i}', (x, y), xytext=(-10, 5), 
                       textcoords='offset points', fontsize=8, color='darkred')
        
        # Right side sources
        for i, (x, y) in enumerate(self.right_sources_physical):
            ax.plot(x, y, 'go', markersize=10, markeredgecolor='darkgreen', 
                   markeredgewidth=2, label='Right Sources' if i == 0 else '')
            ax.annotate(f'R{i}', (x, y), xytext=(10, 5), 
                       textcoords='offset points', fontsize=8, color='darkgreen')
    
    def plot_single_timestep(self, timestep, save_figure=True, show_plot=True):
        """Plot a single timestep with all sources marked in physical units."""
        if timestep not in self.field_data:
            print(f"No data for timestep {timestep}")
            return None
        
        data = self.field_data[timestep]
        
        # Calculate physical extent in micrometers
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + data.shape[0] * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = data.shape[1] * DX * 1e6  # μm
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot field data with physical extent
        im = ax.imshow(data.T, extent=[x_min, x_max, y_min, y_max], 
                      origin='lower', cmap='RdBu_r', aspect='equal')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources
        self.create_base_plot(ax)

        # Labels and title
        ax.set_xlabel('X (μm)', fontsize=12)
        ax.set_ylabel('Y (μm)', fontsize=12)
        time_str = self.format_time_string(timestep)
        ax.set_title(f'FDTD Simulation - {time_str}\n'
                    f'External Laser Sources (22 total: 11L + 11R)', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax, timestep)
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        
        if save_figure:
            output_file = self.output_dir / f"field_plot_t{timestep:04d}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return fig
    
    def add_parameter_text(self, ax, timestep=None):
        """Add key parameters as text on the plot."""
        # Get key parameters
        lambda1_nm = self.geometry_params.get('lambda1_nm', 2000)
        lambda2_nm = self.geometry_params.get('lambda2_nm', 2000)
        intensity_gv = self.geometry_params.get('external_source_intensity_GVm', 2.0)
        num_teeth = self.geometry_params.get('num_teeth', 6)
        gap_width_um = self.geometry_params.get('gap_between_bars_um', 2.0)
        domain_width_um = self.geometry_params.get('total_domain_x_um', 5.0)
        margin_um = self.geometry_params.get('margin_x_um', 1.0)
        
        param_text = f'λ₁ = {lambda1_nm:.0f} nm\n'
        param_text += f'λ₂ = {lambda2_nm:.0f} nm\n'
        param_text += f'Intensity = {intensity_gv:.1f} GV/m\n'
        param_text += f'Teeth per bar = {num_teeth}\n'
        param_text += f'Gap = {gap_width_um:.1f} μm\n'
        param_text += f'Domain = {domain_width_um:.1f} μm\n'
        param_text += f'Margin = {margin_um:.1f} μm'
        
        # Add time information if timestep is provided
        if timestep is not None:
            time_info = self.timestep_to_time_info(timestep)
            param_text += f'\nTime: {time_info["time_fs"]:.1f} fs'
            param_text += f'\nPeriods: {time_info["periods"]:.2f}'
        
        # Position text outside the plot area on the right bottom
        ax.text(1.07, 0.40, param_text, transform=ax.transAxes, 
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9))
    
    def plot_all_timesteps(self):
        """Plot all available timesteps."""
        print(f"Plotting {len(self.field_data)} timesteps...")
        
        for timestep in sorted(self.field_data.keys()):
            self.plot_single_timestep(timestep, save_figure=True, show_plot=False)
        
        print("All plots saved!")
    
    def create_animation(self, interval=1000, save_animation=True):
        """Create a slow animation of the field evolution and save as GIF."""
        if not self.field_data:
            print("No field data available for animation")
            return None
        
        timesteps = sorted(self.field_data.keys())
        
        # Calculate physical extent in micrometers
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        first_data = self.field_data[timesteps[0]]
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + first_data.shape[0] * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = first_data.shape[1] * DX * 1e6  # μm
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Get field data range for consistent colormap
        all_data = np.concatenate([data.flatten() for data in self.field_data.values()])
        vmin, vmax = np.percentile(all_data, [5, 95])  # Use 5-95 percentile for better contrast
        
        # Initialize with first frame
        im = ax.imshow(first_data.T, extent=[x_min, x_max, y_min, y_max], 
                      origin='lower', cmap='RdBu_r', aspect='equal', vmin=vmin, vmax=vmax)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources (they don't change)
        self.create_base_plot(ax)
        
        # Labels
        ax.set_xlabel('X (μm)', fontsize=12)
        ax.set_ylabel('Y (μm)', fontsize=12)
        title = ax.set_title('FDTD Simulation Animation\nExternal Laser Sources (22 total)', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text (no timestep for animation overview)
        self.add_parameter_text(ax)
        
        # Prepare animation elements
        def animate(frame):
            timestep = timesteps[frame]
            data = self.field_data[timestep]
            im.set_array(data.T)
            time_str = self.format_time_string(timestep)
            title.set_text(f'FDTD Simulation - {time_str}\n'
                          f'External Laser Sources (22 total: 11L + 11R)')
            
            return [im, title]
        
        # Create animation with a longer interval for a slower speed
        anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                                     interval=interval, blit=False, repeat=False)
        
        if save_animation:
            output_file = self.output_dir / "field_animation.gif"
            print(f"Saving animation to {output_file}...")
            # Use a low fps for a slower GIF
            anim.save(str(output_file), writer='pillow', fps=0.01, dpi=300)
            print(f"Animation saved: {output_file}")

        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        plt.show()
        
        return anim
    
    def create_potential_animation(self, interval=1000, save_animation=True):
        """Create animation of the potential evolution and save as GIF."""
        if not self.potential_data:
            print("No potential data available for animation")
            return None
        
        timesteps = sorted(self.potential_data.keys())
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Get potential data range for consistent colormap
        all_data = np.concatenate([data.flatten() for data in self.potential_data.values()])
        vmin, vmax = np.percentile(all_data, [5, 95])  # Use 5-95 percentile for better contrast
        
        # Initialize with first frame
        first_data = self.potential_data[timesteps[0]]
        im = ax.imshow(first_data.T, extent=[0, first_data.shape[0], 0, first_data.shape[1]], 
                      origin='lower', cmap='viridis', aspect='equal', vmin=vmin, vmax=vmax)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Potential (V)', fontsize=12)
        
        # Add geometry and sources (they don't change)
        self.create_base_plot(ax)
        
        # Labels
        ax.set_xlabel('X (grid points)', fontsize=12)
        ax.set_ylabel('Y (grid points)', fontsize=12)
        title = ax.set_title('Potential Animation\nExternal Laser Sources (10 total)', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        # Prepare animation elements
        def animate(frame):
            timestep = timesteps[frame]
            data = self.potential_data[timestep]
            im.set_array(data.T)
            title.set_text(f'Electric Potential - Timestep {timestep}\n'
                          f'External Laser Sources (10 total: 5L + 5R)')
            
            return [im, title]
        
        # Create animation
        anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                                     interval=interval, blit=False, repeat=False)
        
        if save_animation:
            output_file = self.output_dir / "potential_animation.gif"
            print(f"Saving potential animation to {output_file}...")
            anim.save(str(output_file), writer='pillow', fps=0.01, dpi=120)
            print(f"Potential animation saved: {output_file}")

        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        plt.show()
        
        return anim
    
    def plot_source_overview(self):
        """Create an overview plot showing just the geometry and sources in physical units."""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate physical domain extent
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        size_x = int(self.geometry_params.get('SIZE_X', 300))
        size_y = int(self.geometry_params.get('SIZE_Y', 300))
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + size_x * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = size_y * DX * 1e6  # μm
        
        # Set up domain
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect('equal')
        
        # Add geometry and sources
        self.create_base_plot(ax)
        
        # Labels and title
        ax.set_xlabel('X (μm)', fontsize=12)
        ax.set_ylabel('Y (μm)', fontsize=12)
        ax.set_title('FDTD Simulation Setup\nBar and Teeth Geometry with External Laser Sources', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        
        # Save the overview
        output_file = self.output_dir / "simulation_overview.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Overview saved: {output_file}")
        
        plt.show()
        
        return fig
    
    def plot_field_and_potential(self, timestep, save_figure=True, show_plot=True):
        """Plot electric field and potential side by side for a single timestep in physical units."""
        if timestep not in self.field_data:
            print(f"No field data for timestep {timestep}")
            return None
        
        field_data = self.field_data[timestep]
        has_potential = timestep in self.potential_data
        
        # Calculate physical extent in micrometers
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + field_data.shape[0] * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = field_data.shape[1] * DX * 1e6  # μm
        
        if has_potential:
            potential_data = self.potential_data[timestep]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        else:
            fig, ax1 = plt.subplots(figsize=(12, 10))
        
        # Plot electric field
        im1 = ax1.imshow(field_data.T, extent=[x_min, x_max, y_min, y_max], 
                        origin='lower', cmap='RdBu_r', aspect='equal')
        
        # Add colorbar for field
        cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
        cbar1.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources to field plot
        self.create_base_plot(ax1)
        
        # Labels for field plot
        ax1.set_xlabel('X (μm)', fontsize=12)
        ax1.set_ylabel('Y (μm)', fontsize=12)
        time_str = self.format_time_string(timestep)
        ax1.set_title(f'Electric Field - {time_str}', fontsize=14)
        ax1.legend(loc='upper right', fontsize=10)
        
        if has_potential:
            # Plot electric potential
            im2 = ax2.imshow(potential_data.T, extent=[x_min, x_max, y_min, y_max], 
                           origin='lower', cmap='viridis', aspect='equal')
            
            # Add colorbar for potential
            cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
            cbar2.set_label('Electric Potential (V)', fontsize=12)
            
            # Add geometry and sources to potential plot
            self.create_base_plot(ax2)
            
            # Labels for potential plot
            ax2.set_xlabel('X (μm)', fontsize=12)
            ax2.set_ylabel('Y (μm)', fontsize=12)
            ax2.set_title(f'Electric Potential - {time_str}', fontsize=14)
            ax2.legend(loc='upper right', fontsize=10)
            
            # Add parameter text to potential plot
            self.add_parameter_text(ax2)
        else:
            # Add parameter text to field plot if no potential
            self.add_parameter_text(ax1)
        
        plt.tight_layout()
        
        if save_figure:
            if has_potential:
                output_file = self.output_dir / f"field_potential_t{timestep:04d}.png"
            else:
                output_file = self.output_dir / f"field_only_t{timestep:04d}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return fig
    
    def plot_potential_only(self, timestep, save_figure=True, show_plot=True):
        """Plot only the electric potential for a single timestep in physical units."""
        if timestep not in self.potential_data:
            print(f"No potential data for timestep {timestep}")
            return None
        
        potential_data = self.potential_data[timestep]
        
        # Calculate physical extent in micrometers
        DX = self.geometry_params.get('DX', 10e-9)  # Updated default to 10 nm
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)  # Grid origin in meters
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + potential_data.shape[0] * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = potential_data.shape[1] * DX * 1e6  # μm
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot potential data
        im = ax.imshow(potential_data.T, extent=[x_min, x_max, y_min, y_max], 
                      origin='lower', cmap='viridis', aspect='equal')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Electric Potential (V)', fontsize=12)
        
        # Add geometry and sources
        self.create_base_plot(ax)

        # Labels and title
        ax.set_xlabel('X (μm)', fontsize=12)
        ax.set_ylabel('Y (μm)', fontsize=12)
        time_str = self.format_time_string(timestep)
        ax.set_title(f'Electric Potential - {time_str}\n'
                    f'External Laser Sources (22 total: 11L + 11R)', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        
        if save_figure:
            output_file = self.output_dir / f"potential_t{timestep:04d}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return fig

    def create_geometry_only_plot(self):
        """Create a plot showing only the geometry and sources for debugging."""
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Print debug information
        print("\n=== GEOMETRY DEBUG INFORMATION ===")
        print(f"DX from file: {self.geometry_params.get('DX', 'NOT_FOUND')}")
        print(f"Grid origin: {self.geometry_params.get('grid_origin_x_m', 'NOT_FOUND')}")
        print(f"SIZE_X: {self.geometry_params.get('SIZE_X', 'NOT_FOUND')}")
        print(f"SIZE_Y: {self.geometry_params.get('SIZE_Y', 'NOT_FOUND')}")
        print(f"Bar1_x (grid): {self.geometry_params.get('bar1_x', 'NOT_FOUND')}")
        print(f"Bar2_x (grid): {self.geometry_params.get('bar2_x', 'NOT_FOUND')}")
        print(f"Left source x (grid): {self.geometry_params.get('left_source_x', 'NOT_FOUND')}")
        print(f"Right source x (grid): {self.geometry_params.get('right_source_x', 'NOT_FOUND')}")
        
        # Calculate physical domain extent
        DX = self.geometry_params.get('DX', 10e-9)
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)
        size_x = int(self.geometry_params.get('SIZE_X', 300))
        size_y = int(self.geometry_params.get('SIZE_Y', 300))
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + size_x * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = size_y * DX * 1e6  # μm
        
        print(f"Physical domain: x=[{x_min:.2f}, {x_max:.2f}] μm, y=[{y_min:.2f}, {y_max:.2f}] μm")
        print("=====================================\n")
        
        # Set up domain
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect('equal')
        
        # Add geometry and sources
        self.create_base_plot(ax)
        
        # Labels and title
        ax.set_xlabel('X (μm)', fontsize=12)
        ax.set_ylabel('Y (μm)', fontsize=12)
        ax.set_title('FDTD Simulation Geometry Debug\nFiner Grid (10 nm) - Check Positions', fontsize=14)
        
        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add parameter text
        self.add_parameter_text(ax)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        
        # Save the debug plot
        output_file = self.output_dir / "geometry_debug.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Debug geometry plot saved: {output_file}")
        
        plt.show()
        
        return fig

    def plot_field_profile_y(self, timestep, save_figure=True, show_plot=True):
        """Plot electric field profile along y-direction at the center of the structure."""
        if timestep not in self.field_data:
            print(f"No field data for timestep {timestep}")
            return None
        
        field_data = self.field_data[timestep]
        
        # Get geometry parameters
        DX = self.geometry_params.get('DX', 10e-9)
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)
        bar1_x = int(self.geometry_params.get('bar1_x', 140))
        bar2_x = int(self.geometry_params.get('bar2_x', 185))
        
        # Calculate center position between the bars (in grid points)
        center_x_grid = (bar1_x + bar2_x) // 2
        
        # Extract field profile along y-direction at center
        if center_x_grid < field_data.shape[0]:
            field_profile = field_data[center_x_grid, :]
            
            # Create y-coordinate array in micrometers
            y_coords = np.arange(len(field_profile)) * DX * 1e6  # Convert to μm
            
            # Create figure
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot field profile
            ax.plot(y_coords, field_profile, 'b-', linewidth=2, label=f'Ey at center (x={center_x_grid} grid pts)')
            
            # Mark the tooth positions
            first_tooth_y = int(self.geometry_params.get('first_tooth_y', 150))
            tooth_spacing = int(self.geometry_params.get('tooth_spacing', 100))
            num_teeth = int(self.geometry_params.get('num_teeth', 10))
            
            for i in range(num_teeth):
                tooth_center_y = first_tooth_y + i * tooth_spacing
                tooth_y_phys = tooth_center_y * DX * 1e6  # Convert to μm
                ax.axvline(x=tooth_y_phys, color='red', linestyle='--', alpha=0.5, 
                          label='Teeth positions' if i == 0 else '')
            
            # Mark source positions
            for i, (x_phys, y_phys) in enumerate(self.left_sources_physical):
                ax.axvline(x=y_phys, color='green', linestyle=':', alpha=0.7,
                          label='Source positions' if i == 0 else '')
            
            # Labels and formatting
            ax.set_xlabel('Y Position (μm)', fontsize=12)
            ax.set_ylabel('Electric Field Ey (V/m)', fontsize=12)
            time_str = self.format_time_string(timestep)
            ax.set_title(f'Electric Field Profile at Center - {time_str}\n'
                        f'Cross-section at x = {(grid_origin_x + center_x_grid * DX) * 1e6:.2f} μm', fontsize=14)
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            
            # Add parameter text outside the plot area at bottom right
            center_x_phys = (grid_origin_x + center_x_grid * DX) * 1e6
            time_info = self.timestep_to_time_info(timestep)
            param_text = f'Profile Position: x = {center_x_phys:.2f} μm\n'
            param_text += f'Between bars at gap center\n'
            param_text += f'Time: {time_info["time_fs"]:.1f} fs ({time_info["periods"]:.2f} periods)'
            
            ax.text(1.07, 0.30, param_text, transform=ax.transAxes, 
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9))
            
            plt.tight_layout()
            plt.subplots_adjust(right=0.85)  # Make room for legend
            
            if save_figure:
                output_file = self.output_dir / f"field_profile_y_t{timestep:04d}.png"
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"Saved field profile: {output_file}")
            
            if show_plot:
                plt.show()
            else:
                plt.close()
            
            return fig
        else:
            print(f"Center position {center_x_grid} is outside field data bounds")
            return None

    def plot_field_and_profile(self, timestep, save_figure=True, show_plot=True):
        """Plot electric field 2D map and 1D profile side by side."""
        if timestep not in self.field_data:
            print(f"No field data for timestep {timestep}")
            return None
        
        field_data = self.field_data[timestep]
        
        # Calculate physical extent in micrometers
        DX = self.geometry_params.get('DX', 10e-9)
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)
        
        x_min = grid_origin_x * 1e6  # μm
        x_max = (grid_origin_x + field_data.shape[0] * DX) * 1e6  # μm
        y_min = 0  # μm
        y_max = field_data.shape[1] * DX * 1e6  # μm
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
        # Left plot: 2D field map
        im1 = ax1.imshow(field_data.T, extent=[x_min, x_max, y_min, y_max], 
                        origin='lower', cmap='RdBu_r', aspect='equal')
        
        # Add colorbar for field
        cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
        cbar1.set_label('Electric Field Ey (V/m)', fontsize=12)
        
        # Add geometry and sources to field plot
        self.create_base_plot(ax1)
        
        # Mark the center line
        bar1_x = int(self.geometry_params.get('bar1_x', 140))
        bar2_x = int(self.geometry_params.get('bar2_x', 185))
        center_x_grid = (bar1_x + bar2_x) // 2
        center_x_phys = (grid_origin_x + center_x_grid * DX) * 1e6
        ax1.axvline(x=center_x_phys, color='yellow', linestyle='-', linewidth=3, 
                   alpha=0.8, label='Profile line')
        
        # Labels for field plot
        ax1.set_xlabel('X (μm)', fontsize=12)
        ax1.set_ylabel('Y (μm)', fontsize=12)
        time_str = self.format_time_string(timestep)
        ax1.set_title(f'Electric Field 2D - {time_str}', fontsize=14)
        ax1.legend(loc='upper right', fontsize=10)
        
        # Right plot: 1D profile
        field_profile = field_data[center_x_grid, :]
        y_coords = np.arange(len(field_profile)) * DX * 1e6  # Convert to μm
        
        ax2.plot(y_coords, field_profile, 'b-', linewidth=2, label=f'Ey at x={center_x_phys:.2f} μm')
        
        # Mark the tooth positions
        first_tooth_y = int(self.geometry_params.get('first_tooth_y', 150))
        tooth_spacing = int(self.geometry_params.get('tooth_spacing', 100))
        num_teeth = int(self.geometry_params.get('num_teeth', 10))
        
        for i in range(num_teeth):
            tooth_center_y = first_tooth_y + i * tooth_spacing
            tooth_y_phys = tooth_center_y * DX * 1e6  # Convert to μm
            ax2.axvline(x=tooth_y_phys, color='red', linestyle='--', alpha=0.5, 
                       label='Teeth positions' if i == 0 else '')
        
        # Mark source positions
        for i, (x_phys, y_phys) in enumerate(self.left_sources_physical):
            ax2.axvline(x=y_phys, color='green', linestyle=':', alpha=0.7,
                       label='Source positions' if i == 0 else '')
        
        # Labels for profile plot
        ax2.set_xlabel('Y Position (μm)', fontsize=12)
        ax2.set_ylabel('Electric Field Ey (V/m)', fontsize=12)
        ax2.set_title(f'Field Profile at Gap Center', fontsize=14)
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        
        if save_figure:
            output_file = self.output_dir / f"field_2d_and_profile_t{timestep:04d}.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {output_file}")
        
        if show_plot:
            plt.show()
        else:
            plt.close()
        
        return fig

    def create_profile_animation(self, interval=1000, save_animation=True):
        """Create animation of the field profile evolution along y-direction."""
        if not self.field_data:
            print("No field data available for animation")
            return None
        
        timesteps = sorted(self.field_data.keys())
        
        # Get geometry parameters
        DX = self.geometry_params.get('DX', 10e-9)
        grid_origin_x = self.geometry_params.get('grid_origin_x_m', 0.0)
        bar1_x = int(self.geometry_params.get('bar1_x', 140))
        bar2_x = int(self.geometry_params.get('bar2_x', 185))
        
        # Calculate center position between the bars (in grid points)
        center_x_grid = (bar1_x + bar2_x) // 2
        center_x_phys = (grid_origin_x + center_x_grid * DX) * 1e6
        
        # Get y-coordinate array
        first_data = self.field_data[timesteps[0]]
        y_coords = np.arange(first_data.shape[1]) * DX * 1e6  # Convert to μm
        
        # Get field data range for consistent y-axis
        all_profiles = []
        for timestep in timesteps:
            if center_x_grid < self.field_data[timestep].shape[0]:
                all_profiles.append(self.field_data[timestep][center_x_grid, :])
        
        if not all_profiles:
            print("No valid profile data found")
            return None
        
        all_data = np.concatenate(all_profiles)
        y_min, y_max = np.percentile(all_data, [1, 99])  # Use 1-99 percentile for range
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Initialize with first frame
        first_profile = all_profiles[0]
        line, = ax.plot(y_coords, first_profile, 'b-', linewidth=2, 
                       label=f'Ey at x={center_x_phys:.2f} μm')
        
        # Mark the tooth positions (static)
        first_tooth_y = int(self.geometry_params.get('first_tooth_y', 150))
        tooth_spacing = int(self.geometry_params.get('tooth_spacing', 100))
        num_teeth = int(self.geometry_params.get('num_teeth', 10))
        
        for i in range(num_teeth):
            tooth_center_y = first_tooth_y + i * tooth_spacing
            tooth_y_phys = tooth_center_y * DX * 1e6  # Convert to μm
            ax.axvline(x=tooth_y_phys, color='red', linestyle='--', alpha=0.5, 
                      label='Teeth positions' if i == 0 else '')
        
        # Mark source positions (static)
        for i, (x_phys, y_phys) in enumerate(self.left_sources_physical):
            ax.axvline(x=y_phys, color='green', linestyle=':', alpha=0.7,
                      label='Source positions' if i == 0 else '')
        
        # Set up plot formatting
        ax.set_xlim(y_coords[0], y_coords[-1])
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel('Y Position (μm)', fontsize=12)
        ax.set_ylabel('Electric Field Ey (V/m)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Title that will be updated
        title = ax.set_title(f'Electric Field Profile Animation\nCross-section at x = {center_x_phys:.2f} μm', fontsize=14)
        
        # Animation function
        def animate(frame):
            timestep = timesteps[frame]
            if center_x_grid < self.field_data[timestep].shape[0]:
                profile = self.field_data[timestep][center_x_grid, :]
                line.set_ydata(profile)
                time_str = self.format_time_string(timestep)
                title.set_text(f'Electric Field Profile - {time_str}\n'
                              f'Cross-section at x = {center_x_phys:.2f} μm')
            return [line, title]
        
        # Create animation
        anim = animation.FuncAnimation(fig, animate, frames=len(timesteps), 
                                     interval=interval, blit=False, repeat=False)
        
        if save_animation:
            output_file = self.output_dir / "field_profile_animation.gif"
            print(f"Saving profile animation to {output_file}...")
            anim.save(str(output_file), writer='pillow', fps=1, dpi=150)
            print(f"Profile animation saved: {output_file}")

        plt.tight_layout()
        plt.subplots_adjust(right=0.85)  # Make room for legend
        plt.show()
        
        return anim

def select_directories_gui():
    """GUI function to select data and output directories."""
    if not GUI_AVAILABLE:
        print("GUI not available. Use command line interface.")
        return None, None
    
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    
    # Select data directory
    messagebox.showinfo("Select Directories", "First, select the directory containing simulation data files")
    data_dir = filedialog.askdirectory(title="Select Data Directory")
    
    if not data_dir:
        return None, None
    
    # Ask if user wants custom output directory
    use_custom = messagebox.askyesno("Output Directory", 
                                   "Do you want to specify a custom output directory?\n\n"
                                   "Choose 'No' to save results in the same directory as data files.")
    
    output_dir = None
    if use_custom:
        output_dir = filedialog.askdirectory(title="Select Output Directory")
        if not output_dir:
            output_dir = data_dir  # Fall back to data directory
    
    root.destroy()
    return data_dir, output_dir

def main_gui():
    """Main function with GUI directory selection."""
    print("FDTD Laser Simulation Visualizer (GUI Mode)")
    print("=" * 45)
    
    if not GUI_AVAILABLE:
        print("GUI not available. Falling back to command line mode.")
        main()
        return
    
    # Select directories using GUI
    data_dir, output_dir = select_directories_gui()
    
    if not data_dir:
        print("No data directory selected. Exiting.")
        return
    
    # Initialize visualizer
    visualizer = FDTDVisualizerWithSources(data_dir, output_dir)
    
    print(f"Data directory: {visualizer.data_dir}")
    print(f"Output directory: {visualizer.output_dir}")
    
    # Continue with normal visualization workflow
    run_visualization_workflow(visualizer)

def run_visualization_workflow(visualizer):
    """Run the main visualization workflow."""
    # First, create a geometry debug plot to check if everything is positioned correctly
    print("\n0. Creating geometry debug plot...")
    visualizer.create_geometry_only_plot()
    
    if not visualizer.field_data:
        print("No field data found. Please run the C++ simulation first.")
        print("But you can see the geometry setup above.")
        return
    
    # Create overview plot
    print("\n1. Creating simulation overview...")
    visualizer.plot_source_overview()
    
    # Check if potential data is available
    has_potential = bool(visualizer.potential_data)
    if has_potential:
        print(f"Found potential data for {len(visualizer.potential_data)} timesteps")
    else:
        print("No potential data found - field plots only")
    
    # Plot a few key timesteps
    print("\n2. Plotting key timesteps...")
    timesteps = sorted(visualizer.field_data.keys())
    # first 5 timesteps or all if less than 5
    key_timesteps = timesteps[:5] if len(timesteps) >= 5 else timesteps
    
    for ts in key_timesteps:
        if ts in visualizer.field_data:
            # Plot 2D field and 1D profile side by side
            visualizer.plot_field_and_profile(ts, save_figure=True, show_plot=True)
    
    # Ask user about full processing
    response = input("\n3. Do you want to plot ALL timesteps? (y/n): ").strip().lower()
    if response == 'y':
        plot_type = input("Plot (f)ield 2D only, (p)rofile 1D only, or (b)oth side-by-side? (f/p/b): ").strip().lower()
        if plot_type == 'p':
            print("Plotting field profiles for all timesteps...")
            for timestep in sorted(visualizer.field_data.keys()):
                visualizer.plot_field_profile_y(timestep, save_figure=True, show_plot=False)
        elif plot_type == 'b':
            print("Plotting 2D field and 1D profile for all timesteps...")
            for timestep in sorted(visualizer.field_data.keys()):
                visualizer.plot_field_and_profile(timestep, save_figure=True, show_plot=False)
        else:
            print("Plotting 2D field for all timesteps...")
            visualizer.plot_all_timesteps()
    
    # Ask about animation
    response = input("\n4. Do you want to create animations? (y/n): ").strip().lower()
    if response == 'y':
        anim_type = input("Create (f)ield 2D animation or (p)rofile 1D animation? (f/p): ").strip().lower()
        if anim_type == 'p':
            visualizer.create_profile_animation(interval=1000)
        else:
            visualizer.create_animation(interval=1000)
    
    print(f"\nVisualization complete! Results saved in: {visualizer.output_dir}")

def main():
    """Main function to run the visualization."""
    print("FDTD Laser Simulation Visualizer")
    print("=" * 35)
    
    # Ask user for output directory
    output_choice = input("\nChoose output directory:\n"
                         "1. Same as data directory (default)\n"
                         "2. Specify custom directory\n"
                         "3. Use GUI to select directories\n"
                         "Enter choice (1/2/3): ").strip()
    
    if output_choice == "3" and GUI_AVAILABLE:
        main_gui()
        return
    elif output_choice == "3":
        print("GUI not available. Using command line interface.")
    
    output_dir = None
    if output_choice == "2":
        output_dir = input("Enter output directory path: ").strip()
        if output_dir and not Path(output_dir).exists():
            create_dir = input(f"Directory '{output_dir}' doesn't exist. Create it? (y/n): ").strip().lower()
            if create_dir != 'y':
                print("Using default data directory for output.")
                output_dir = None
    
    # Initialize visualizer
    visualizer = FDTDVisualizerWithSources(".", output_dir)
    
    print(f"Data directory: {visualizer.data_dir}")
    print(f"Output directory: {visualizer.output_dir}")
    
    # Run the visualization workflow
    run_visualization_workflow(visualizer)

if __name__ == "__main__":
    main()
