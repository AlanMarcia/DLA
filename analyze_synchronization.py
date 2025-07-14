#!/usr/bin/env python3
"""
Field Analysis Tool for FDTD Laser Simulation

This script helps analyze the electromagnetic field evolution in the dielectric laser accelerator structure.
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from pathlib import Path

def load_field_data(filename):
    """Load electric field data from file"""
    try:
        data = np.loadtxt(filename)
        return data
    except:
        print(f"Could not load {filename}")
        return None

def load_electron_trajectory():
    """Load electron trajectory data - currently not used"""
    print("Electron trajectory analysis not available in this version")
    return None

def analyze_field_evolution():
    """Analyze how the field evolves over time at key locations"""
    field_files = sorted(glob.glob('Ey_field_*.dat'))
    
    if not field_files:
        print("No field files found")
        return
    
    times = []
    field_at_monitor = []
    
    # Parameters (should match C++ code)
    DX = 1e-9  # 1 nm grid spacing
    monitor_x_grid = 162  # Grid point to monitor
    monitor_y_grid = 200  # Grid point to monitor
    
    for filename in field_files:
        # Extract time from filename
        match = re.search(r'Ey_field_(\d+)\.dat', filename)
        if match:
            time_step = int(match.group(1))
            times.append(time_step)
            
            # Load field data
            field_data = load_field_data(filename)
            if field_data is not None:
                # Extract field at monitor point
                field_value = field_data[monitor_y_grid, monitor_x_grid]
                field_at_monitor.append(field_value)
            else:
                field_at_monitor.append(0)
    
    return np.array(times), np.array(field_at_monitor)

def find_field_peaks(times, field_values, threshold_fraction=0.5):
    """Find peaks in the field evolution"""
    if len(field_values) == 0:
        return []
    
    # Find peaks in field oscillations
    max_field_magnitude = np.max(np.abs(field_values))
    threshold = threshold_fraction * max_field_magnitude
    
    peak_times = []
    for i in range(1, len(field_values)-1):
        # Look for local maxima in field magnitude
        if (np.abs(field_values[i]) > np.abs(field_values[i-1]) and 
            np.abs(field_values[i]) > np.abs(field_values[i+1]) and 
            np.abs(field_values[i]) > threshold):
            peak_times.append(times[i])
    
    return peak_times

def plot_field_analysis():
    """Create plots to visualize field evolution"""
    times, field_values = analyze_field_evolution()
    
    if len(times) == 0:
        print("No field data available for analysis")
        return
    
    # Convert field to GV/m for display
    field_GV_m = field_values / 1e9
    
    # Find field peaks
    peak_times = find_field_peaks(times, field_values)
    
    plt.figure(figsize=(12, 8))
    
    # Plot 1: Field evolution over time
    plt.subplot(2, 2, 1)
    plt.plot(times, field_GV_m, 'b-', linewidth=2)
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Mark field peaks
    for peak_time in peak_times[:5]:  # Show first 5 peaks
        plt.axvline(x=peak_time, color='r', linestyle=':', alpha=0.7, 
                   label=f'Field peak' if peak_time == peak_times[0] else "")
    
    plt.xlabel('Time Step')
    plt.ylabel('Electric Field (GV/m)')
    plt.title('Field Evolution at Monitor Point')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Plot 2: Phase analysis
    plt.subplot(2, 2, 2)
    # Calculate approximate phase from field oscillation
    if len(field_values) > 10:
        # Find zero crossings to estimate frequency
        zero_crossings = []
        for i in range(1, len(field_values)):
            if field_values[i-1] * field_values[i] < 0:
                zero_crossings.append(times[i])
        
        if len(zero_crossings) > 2:
            avg_period = np.mean(np.diff(zero_crossings)) * 2  # Full period
            phase = 2 * np.pi * times / avg_period
            
            plt.plot(phase % (2*np.pi), field_GV_m, 'g-', alpha=0.7)
            plt.xlabel('Phase (radians)')
            plt.ylabel('Electric Field (GV/m)')
            plt.title('Field vs Phase')
            plt.grid(True, alpha=0.3)
            
            # Mark peak phase regions
            peak_phases = []
            for peak_time in peak_times:
                if peak_time < len(times):
                    peak_phase = (2 * np.pi * peak_time / avg_period) % (2*np.pi)
                    peak_phases.append(peak_phase)
            
            for phase in peak_phases[:3]:
                plt.axvline(x=phase, color='r', linestyle=':', alpha=0.7)
    
    # Plot 3: Field amplitude analysis
    plt.subplot(2, 2, 3)
    if len(field_values) > 0:
        # Plot field amplitude over time
        plt.plot(times, np.abs(field_GV_m), 'g-', linewidth=2)
        plt.xlabel('Time Step')
        plt.ylabel('Field Amplitude (GV/m)')
        plt.title('Field Amplitude Evolution')
        plt.grid(True, alpha=0.3)
        
        # Plot 4: Field statistics
        plt.subplot(2, 2, 4)
        # Histogram of field values
        plt.hist(field_GV_m, bins=50, alpha=0.7, color='blue', edgecolor='black')
        plt.xlabel('Electric Field (GV/m)')
        plt.ylabel('Frequency')
        plt.title('Field Distribution')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('field_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print field statistics
    print("\n=== FIELD ANALYSIS ===")
    print(f"Field strength range: {np.min(field_GV_m):.2f} to {np.max(field_GV_m):.2f} GV/m")
    print(f"Field RMS: {np.sqrt(np.mean(field_GV_m**2)):.2f} GV/m")
    
    if peak_times:
        print(f"\nField peaks found (first 5):")
        for i, peak_time in enumerate(peak_times[:5]):
            field_at_time = field_GV_m[list(times).index(peak_time)] if peak_time in times else 0
            print(f"  Time step {peak_time}: Field = {field_at_time:.2f} GV/m")
    else:
        print("\nNo significant field peaks found with current criteria")

def estimate_field_power():
    """Estimate the power and energy in the electromagnetic field"""
    field_files = glob.glob('Ey_field_*.dat')
    if not field_files:
        print("Cannot analyze field power without field data")
        return
    
    # Load a representative field file
    latest_file = sorted(field_files)[-1]
    field_data = load_field_data(latest_file)
    if field_data is None:
        print("Could not load field data for power analysis")
        return
    
    # Calculate field energy density
    field_energy_density = 0.5 * 8.854e-12 * np.sum(field_data**2)  # J/m³
    total_volume = (len(field_data) * len(field_data[0])) * (40e-9)**2  # Approximate volume
    total_energy = field_energy_density * total_volume
    
    print(f"\n=== FIELD POWER ANALYSIS ===")
    print(f"Peak field strength: {np.max(np.abs(field_data))/1e9:.2f} GV/m")
    print(f"RMS field strength: {np.sqrt(np.mean(field_data**2))/1e9:.2f} GV/m")
    print(f"Total field energy: {total_energy*1e12:.2f} pJ")
    print(f"Field energy density: {field_energy_density:.2e} J/m³")

def main():
    print("FDTD Laser Simulation - Field Analysis")
    print("=" * 40)
    
    # Check if simulation data exists
    field_files = glob.glob('Ey_field_*.dat')
    
    if not field_files:
        print("No field files found. Run the simulation first.")
        return
    
    print(f"Found {len(field_files)} field files")
    
    # Perform analysis
    plot_field_analysis()
    estimate_field_power()
    
    print(f"\nResults saved to: field_analysis.png")

if __name__ == "__main__":
    main()
