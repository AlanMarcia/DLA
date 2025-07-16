#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

// --- Constants ---
const double PI = 3.14159265358979323846;
const double C0 = 299792458.0; // Speed of light in vacuum (m/s)
const double MU0 = 4.0 * PI * 1e-7; // Permeability of free space (H/m)
const double EPS0 = 1.0 / (MU0 * C0 * C0); // Permittivity of free space (F/m)

// --- Laser Source Parameters (defining wavelengths first) ---
const double lambda1 = 2e-6; 
const double lambda2 = 2e-6; // Wavelength of laser 2 (m) = 1500 nm
const double f1 = C0 / lambda1; // 
const double f2 = C0 / lambda2; // Frequency of laser 2 (Hz)
const double intensity1 = 1.0e9; // Intensity of laser 1 (V/m) - Target: 5 GV/m
const double intensity2 = 1.0e9; // Intensity of laser 2 (V/m) - Target: 5 GV/m

// --- Dielectric Bar Parameters ---
// Refractive indices at 2000 nm wavelength
const double n_air = 1.0; // Refractive index of air/vacuum
const double n_sio2 = 1.44; // Refractive index of SiO2 at 2000 nm (near-IR)
const double n_medium = n_air; // Medium where laser propagates (can be changed to other materials)

// Relative permittivity calculated from refractive index: ε_r = n²
const double eps_r_air = n_air * n_air; // Air/vacuum
const double eps_r_sio2 = n_sio2 * n_sio2; // SiO2 bars (~2.07 at 2000 nm)

// Bar geometry parameters (in SI units - meters)
const double bar_width_m = 500e-9; // Bar thickness = 500 nm
const double bar1_x_m = 5e-6; // Bar 1 center position = 5 μm
const double bar2_x_m = 7e-6; // Bar 2 center position = 7 μm
const double gap_between_bars_m = bar2_x_m - bar1_x_m; // Gap between bars = 2 μm

// --- Simulation Parameters ---
// Configurable grid resolution (change only this line to adjust fineness)
const double DX = 10e-9; // Grid spacing in meters (1 nm = ultra-fine, 2 nm = very fine, 5 nm = fine, 10 nm = standard)

// Speed of light in the propagation medium
const double c_medium = C0 / n_medium; // Speed of light in the medium
const double DT = DX / (c_medium * sqrt(2)); // Time step (s) - Courant stability condition for medium

// Wavelength in the medium (shorter than in vacuum)
const double lambda1_medium = lambda1 / n_medium; // Wavelength in the propagation medium
const double lambda2_medium = lambda2 / n_medium; // Wavelength in the propagation medium

// Configurable domain sizes (physical dimensions)
const double desired_domain_y_m = 12e-6; // Desired Y domain size in meters (12 μm)
const double margin_x = 1e-6; // 1 μm margin on each side in X direction

// Calculate grid size to fit structure with margins
const double tooth_height_m = 0.3 * lambda1_medium; // Tooth height based on wavelength in medium
const double structure_left = bar1_x_m - bar_width_m/2 - tooth_height_m; // Leftmost point
const double structure_right = bar2_x_m + bar_width_m/2 + tooth_height_m; // Rightmost point
const double total_domain_x = structure_right - structure_left + 2*margin_x; // Total domain width

// Grid sizes automatically calculated from physical dimensions and DX
const int SIZE_X = (int)(total_domain_x / DX); // Grid size in x direction (automatically scales with DX)
const int SIZE_Y = (int)(desired_domain_y_m / DX); // Grid size in y direction (automatically scales with DX)
const int MAX_TIME = 1000; // Total simulation time steps (dimensionless)

// Calculate new grid origin to center structure with margins
const double grid_origin_x = structure_left - margin_x; // New x=0 position in physical space

// Convert bar positions to grid points (relative to new origin)
const int bar_width = (int)(bar_width_m / DX);
const int bar1_x = (int)((bar1_x_m - grid_origin_x) / DX);
const int bar2_x = (int)((bar2_x_m - grid_origin_x) / DX);
const int gap_between_bars = bar2_x - bar1_x;

// --- Absorbing Boundary Conditions (PML) ---
const int PML_WIDTH = 10; // PML layer thickness (grid points)
const double PML_SIGMA_MAX = 0.8 * (3.0 + 1.0) / (377.0 * DX); // Maximum conductivity

// Teeth parameters (in terms of wavelength in medium)
const int num_teeth = 10; // Number of teeth per bar
const int tooth_width = (int)(0.3 * lambda1_medium / DX); // Width of each tooth = 0.3λ₁ in medium
const int tooth_height = (int)(0.3 * lambda1_medium / DX); // Height of teeth extending into gap
const int tooth_spacing = (int)(0.5 * lambda1_medium / DX); // Spacing between teeth centers = 0.5λ₁ in medium
const int first_tooth_y = (int)(1.5e-6 / DX); // Position of first tooth = 1 μm (SI units)

// Bar length (extends in y-direction) - in SI units
const double bar_start_y_m = 2e-6; // Start at 1 μm
const double bar_end_y_m = 12e-6; // End at 11 μm
const int bar_start_y = (int)(bar_start_y_m / DX); // Convert to grid points
const int bar_end_y = (int)(bar_end_y_m / DX); // Convert to grid points

// --- Source Positions ---
// External laser sources positioned outside the structure (in SI units)
const double left_source_x_m = 4.5e-6; // Left sources near left boundary
const double right_source_x_m = 7.5e-6; // Right sources near right boundary

// Laser intensity for external sources
const double external_source_intensity = 2.0e9; // 2 GV/m intensity for external sources

// Left side sources (positioned in the middle of spaces between teeth) - adjusted to new grid
const int left_source_x = (int)((left_source_x_m - grid_origin_x) / DX);
const int left_source1_y = first_tooth_y + (int)(0.5 * tooth_spacing); // Between teeth 1 and 2
const int left_source2_y = first_tooth_y + (int)(1.5 * tooth_spacing); // Between teeth 2 and 3
const int left_source3_y = first_tooth_y + (int)(2.5 * tooth_spacing); // Between teeth 3 and 4
const int left_source4_y = first_tooth_y + (int)(3.5 * tooth_spacing); // Between teeth 4 and 5
const int left_source5_y = first_tooth_y + (int)(4.5 * tooth_spacing); // Between teeth 5 and 6
const int left_source6_y = first_tooth_y + (int)(5.5 * tooth_spacing); // Between teeth 6 and 7
const int left_source7_y = first_tooth_y + (int)(6.5 * tooth_spacing); // Between teeth 7 and 8
const int left_source8_y = first_tooth_y + (int)(7.5 * tooth_spacing); // Between teeth 8 and 9
const int left_source9_y = first_tooth_y + (int)(8.5 * tooth_spacing); // Between teeth 9 and 10
const int left_source10_y = first_tooth_y + (int)(9.5 * tooth_spacing); // Between teeth 10 and 11

// Right side sources (positioned in the middle of spaces between teeth) - adjusted to new grid
const int right_source_x = (int)((right_source_x_m - grid_origin_x) / DX);
const int right_source1_y = first_tooth_y + (int)(0.5 * tooth_spacing); // Between teeth 1 and 2
const int right_source2_y = first_tooth_y + (int)(1.5 * tooth_spacing); // Between teeth 2 and 3
const int right_source3_y = first_tooth_y + (int)(2.5 * tooth_spacing); // Between teeth 3 and 4
const int right_source4_y = first_tooth_y + (int)(3.5 * tooth_spacing); // Between teeth 4 and 5
const int right_source5_y = first_tooth_y + (int)(4.5 * tooth_spacing); // Between teeth 5 and 6
const int right_source6_y = first_tooth_y + (int)(5.5 * tooth_spacing); // Between teeth 6 and 7
const int right_source7_y = first_tooth_y + (int)(6.5 * tooth_spacing); // Between teeth 7 and 8
const int right_source8_y = first_tooth_y + (int)(7.5 * tooth_spacing); // Between teeth 8 and 9
const int right_source9_y = first_tooth_y + (int)(8.5 * tooth_spacing); // Between teeth 9 and 10
const int right_source10_y = first_tooth_y + (int)(9.5 * tooth_spacing); // Between teeth 10 and 11

// Additional sources before the first tooth
const int left_source0_y = first_tooth_y - (int)(0.5 * tooth_spacing); // Before first tooth (left side)
const int right_source0_y = first_tooth_y - (int)(0.5 * tooth_spacing); // Before first tooth (right side)

// --- Main FDTD Class ---
class FDTDSimulator {
public:
    FDTDSimulator() :
        Ey(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        Hx(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        Hz(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        potential(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        eps(SIZE_X, std::vector<double>(SIZE_Y, EPS0)),
        sigma_x(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        sigma_y(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        time(0)
    {
        // Initialize PML absorbing boundaries
        initialize_pml();
        
        // Initialize the dielectric bars with teeth
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < SIZE_X; ++i) {
            for (int j = 0; j < SIZE_Y; ++j) {
                bool in_dielectric = false;
                
                // Initialize with medium permittivity (air by default)
                eps[i][j] = eps_r_air * EPS0;
                
                // Check if point is in Bar 1 (left bar)
                if (i >= (bar1_x - bar_width/2) && i <= (bar1_x + bar_width/2) &&
                    j >= bar_start_y && j <= bar_end_y) {
                    in_dielectric = true;
                }
                
                // Check if point is in Bar 2 (right bar)
                if (i >= (bar2_x - bar_width/2) && i <= (bar2_x + bar_width/2) &&
                    j >= bar_start_y && j <= bar_end_y) {
                    in_dielectric = true;
                }
                
                // Add teeth extending from Bar 1 towards the gap
                for (int t = 0; t < num_teeth; ++t) {
                    int tooth_center_y = first_tooth_y + t * tooth_spacing;
                    
                    // Tooth extending from Bar 1 (towards right)
                    if (i >= (bar1_x + bar_width/2) && i <= (bar1_x + bar_width/2 + tooth_height) &&
                        j >= (tooth_center_y - tooth_width/2) && j <= (tooth_center_y + tooth_width/2)) {
                        in_dielectric = true;
                    }
                    
                    // Tooth extending from Bar 2 (towards left)
                    if (i >= (bar2_x - bar_width/2 - tooth_height) && i <= (bar2_x - bar_width/2) &&
                        j >= (tooth_center_y - tooth_width/2) && j <= (tooth_center_y + tooth_width/2)) {
                        in_dielectric = true;
                    }
                }
                
                if (in_dielectric) {
                    eps[i][j] = eps_r_sio2 * EPS0; // SiO2 dielectric material
                }
            }
        }
    }

    void step() {
        update_H();
        update_E();
        apply_sources();
        //calculate_potential();
        time++;
    }

    void save_field_to_file(const std::string& filename) {
        std::ofstream file(filename);
        for (int i = 0; i < SIZE_X; ++i) {
            for (int j = 0; j < SIZE_Y; ++j) {
                file << Ey[i][j] << " ";
            }
            file << "\n";
        }
        file.close();
    }

    // void save_potential_to_file(const std::string& filename) {
    //     std::ofstream file(filename);
    //     for (int i = 0; i < SIZE_X; ++i) {
    //         for (int j = 0; j < SIZE_Y; ++j) {
    //             file << potential[i][j] << " ";
    //         }
    //         file << "\n";
    //     }
    //     file.close();
    // }

    void save_geometry_params(const std::string& filename) {
        std::ofstream file(filename);
        
        // Save wavelength-based parameters (fundamental design parameters)

        file << "lambda1 " << lambda1 << "\n"; // Wavelength 1 in meters
        file << "lambda2 " << lambda2 << "\n"; // Wavelength 2 in meters
        file << "lambda1_nm " << lambda1 * 1e9 << "\n"; // Wavelength 1 in nm
        file << "lambda2_nm " << lambda2 * 1e9 << "\n"; // Wavelength 2 in nm
        
        // Save grid parameters
        file << "SIZE_X " << SIZE_X << "\n";
        file << "SIZE_Y " << SIZE_Y << "\n";
        file << "DX " << DX << "\n";
        file << "DX_nm " << DX * 1e9 << "\n"; // Grid spacing in nanometers
        file << "DX_over_lambda1 " << DX / lambda1 << "\n"; // Grid spacing / λ₁
        file << "points_per_lambda1 " << lambda1 / DX << "\n"; // Points per λ₁
        file << "points_per_lambda2 " << lambda2 / DX << "\n"; // Points per λ₂
        
        // Physical domain size
        file << "domain_x_um " << SIZE_X * DX * 1e6 << "\n"; // Domain size in μm
        file << "domain_y_um " << SIZE_Y * DX * 1e6 << "\n"; // Domain size in μm
        file << "domain_x_lambda1 " << SIZE_X * DX / lambda1 << "\n"; // Domain size in λ₁
        file << "domain_y_lambda1 " << SIZE_Y * DX / lambda1 << "\n"; // Domain size in λ₁
        
        // Grid origin and margins (new compact grid with 1 μm margins)
        file << "grid_origin_x_m " << grid_origin_x << "\n"; // Grid origin x in meters
        file << "grid_origin_x_um " << grid_origin_x * 1e6 << "\n"; // Grid origin x in μm
        file << "margin_x_m " << margin_x << "\n"; // X-direction margin in meters
        file << "margin_x_um " << margin_x * 1e6 << "\n"; // X-direction margin in μm
        file << "structure_left_m " << structure_left << "\n"; // Structure left edge in meters
        file << "structure_left_um " << structure_left * 1e6 << "\n"; // Structure left edge in μm
        file << "structure_right_m " << structure_right << "\n"; // Structure right edge in meters
        file << "structure_right_um " << structure_right * 1e6 << "\n"; // Structure right edge in μm
        file << "total_domain_x_m " << total_domain_x << "\n"; // Total domain width in meters
        file << "total_domain_x_um " << total_domain_x * 1e6 << "\n"; // Total domain width in μm
        
        // Time parameters
        file << "DT " << DT << "\n"; // Time step in seconds
        file << "DT_fs " << DT * 1e15 << "\n"; // Time step in femtoseconds
        file << "MAX_TIME " << MAX_TIME << "\n";
        file << "total_time_fs " << MAX_TIME * DT * 1e15 << "\n"; // Total simulation time in fs
        file << "period1_fs " << (1.0 / f1) * 1e15 << "\n"; // Period of λ₁ in fs
        file << "period2_fs " << (1.0 / f2) * 1e15 << "\n"; // Period of λ₂ in fs
        file << "time_steps_per_period1 " << (1.0 / f1) / DT << "\n"; // Time steps per period λ₁
        
        // Save laser parameters
        file << "f1 " << f1 << "\n"; // Frequency in Hz
        file << "f2 " << f2 << "\n"; // Frequency in Hz
        file << "external_source_intensity " << external_source_intensity << "\n"; // External source intensity in V/m
        file << "external_source_intensity_GVm " << external_source_intensity * 1e-9 << "\n"; // External intensity in GV/m
        file << "total_sources " << 22 << "\n"; // Total number of laser sources (11 left + 11 right external sources)
        
        // Save dielectric parameters
        file << "eps_r_air " << eps_r_air << "\n"; // Relative permittivity of air (dimensionless)
        file << "eps_r_sio2 " << eps_r_sio2 << "\n"; // Relative permittivity of SiO2 (dimensionless)
        file << "n_air " << n_air << "\n"; // Refractive index of air
        file << "n_sio2 " << n_sio2 << "\n"; // Refractive index of SiO2
        file << "n_medium " << n_medium << "\n"; // Refractive index of propagation medium
        file << "num_teeth " << num_teeth << "\n"; // Number of teeth per bar
        
        // Bar dimensions (SI units and grid points)
        file << "bar_width_m " << bar_width_m << "\n"; // Bar width in meters
        file << "bar_width_nm " << bar_width_m * 1e9 << "\n"; // Bar width in nanometers
        file << "bar_width " << bar_width << "\n"; // Bar width in grid points
        file << "bar1_x_m " << bar1_x_m << "\n"; // Bar 1 position in meters
        file << "bar1_x_um " << bar1_x_m * 1e6 << "\n"; // Bar 1 position in micrometers
        file << "bar1_x " << bar1_x << "\n"; // Bar 1 position in grid points
        file << "bar2_x_m " << bar2_x_m << "\n"; // Bar 2 position in meters
        file << "bar2_x_um " << bar2_x_m * 1e6 << "\n"; // Bar 2 position in micrometers
        file << "bar2_x " << bar2_x << "\n"; // Bar 2 position in grid points
        file << "gap_between_bars_m " << gap_between_bars_m << "\n"; // Gap in meters
        file << "gap_between_bars_um " << gap_between_bars_m * 1e6 << "\n"; // Gap in micrometers
        file << "gap_between_bars " << gap_between_bars << "\n"; // Gap in grid points
        
        // Bar length (SI units and grid points)
        file << "bar_start_y_m " << bar_start_y_m << "\n"; // Bar start in meters
        file << "bar_start_y_um " << bar_start_y_m * 1e6 << "\n"; // Bar start in micrometers
        file << "bar_start_y " << bar_start_y << "\n"; // Bar start in grid points
        file << "bar_end_y_m " << bar_end_y_m << "\n"; // Bar end in meters
        file << "bar_end_y_um " << bar_end_y_m * 1e6 << "\n"; // Bar end in micrometers
        file << "bar_end_y " << bar_end_y << "\n"; // Bar end in grid points
        
        // Tooth dimensions (wavelength-based, grid points, and SI units)
        file << "tooth_width_lambda1 " << tooth_width * DX / lambda1 << "\n"; // Width in λ₁
        file << "tooth_width " << tooth_width << "\n"; // Width in grid points
        file << "tooth_width_nm " << tooth_width * DX * 1e9 << "\n"; // Width in nanometers
        file << "tooth_height_lambda1 " << tooth_height * DX / lambda1 << "\n"; // Height in λ₁
        file << "tooth_height " << tooth_height << "\n"; // Height in grid points
        file << "tooth_height_nm " << tooth_height * DX * 1e9 << "\n"; // Height in nanometers
        file << "tooth_spacing_lambda1 " << tooth_spacing * DX / lambda1 << "\n"; // Spacing in λ₁
        file << "tooth_spacing " << tooth_spacing << "\n"; // Spacing in grid points
        file << "tooth_spacing_nm " << tooth_spacing * DX * 1e9 << "\n"; // Spacing in nanometers
        
        // First tooth position (SI units and grid points)
        file << "first_tooth_y_m " << first_tooth_y * DX << "\n"; // Position in meters
        file << "first_tooth_y_um " << first_tooth_y * DX * 1e6 << "\n"; // Position in micrometers
        file << "first_tooth_y " << first_tooth_y << "\n"; // Position in grid points
        
        // External source positions (SI units and grid points)
        file << "left_source_x_m " << left_source_x_m << "\n"; // Left source x in meters
        file << "left_source_x_nm " << left_source_x_m * 1e9 << "\n"; // Left source x in nanometers
        file << "left_source_x " << left_source_x << "\n"; // Left source x in grid points
        file << "right_source_x_m " << right_source_x_m << "\n"; // Right source x in meters
        file << "right_source_x_um " << right_source_x_m * 1e6 << "\n"; // Right source x in micrometers
        file << "right_source_x " << right_source_x << "\n"; // Right source x in grid points
        
        // External source y-positions (before first tooth and between teeth)
        file << "left_source0_y " << left_source0_y << "\n"; // Left source 0 y in grid points (before first tooth)
        file << "left_source0_y_um " << left_source0_y * DX * 1e6 << "\n"; // Left source 0 y in micrometers
        file << "left_source1_y " << left_source1_y << "\n"; // Left source 1 y in grid points
        file << "left_source1_y_um " << left_source1_y * DX * 1e6 << "\n"; // Left source 1 y in micrometers
        file << "right_source0_y " << right_source0_y << "\n"; // Right source 0 y in grid points (before first tooth)
        file << "right_source0_y_um " << right_source0_y * DX * 1e6 << "\n"; // Right source 0 y in micrometers
        file << "right_source1_y " << right_source1_y << "\n"; // Right source 1 y in grid points
        file << "right_source1_y_um " << right_source1_y * DX * 1e6 << "\n"; // Right source 1 y in micrometers
        
        // Design ratios and scaling factors
        file << "tooth_fill_factor " << (double)tooth_width / (double)tooth_spacing << "\n"; // Tooth fill factor
        file << "domain_aspect_ratio " << (double)SIZE_X / (double)SIZE_Y << "\n"; // Domain aspect ratio
        
        // PML Absorbing Boundary Conditions
        file << "PML_WIDTH " << PML_WIDTH << "\n"; // PML layer thickness (grid points)
        file << "PML_WIDTH_nm " << PML_WIDTH * DX * 1e9 << "\n"; // PML thickness in nanometers
        file << "PML_WIDTH_lambda1 " << PML_WIDTH * DX / lambda1 << "\n"; // PML thickness in λ₁
        file << "PML_SIGMA_MAX " << PML_SIGMA_MAX << "\n"; // Maximum PML conductivity
        
        // Physical constants used
        file << "C0 " << C0 << "\n"; // Speed of light (m/s)
        file << "MU0 " << MU0 << "\n"; // Permeability of free space (H/m)
        file << "EPS0 " << EPS0 << "\n"; // Permittivity of free space (F/m)
        
        // Potential calculation parameters
        file << "potential_calculation " << "enabled" << "\n"; // Potential calculation status
        file << "potential_reference " << "boundary_zero" << "\n"; // Potential reference point
        file << "potential_method " << "field_integration_with_correction" << "\n"; // Calculation method
        file << "potential_solver " << "direct_integration" << "\n"; // Numerical solver
        file << "potential_smoothing " << "local_averaging" << "\n"; // Smoothing method
        
        file.close();
    }

private:
    std::vector<std::vector<double>> Ey;
    std::vector<std::vector<double>> Hx;
    std::vector<std::vector<double>> Hz;
    std::vector<std::vector<double>> potential; // Electric potential (V)
    std::vector<std::vector<double>> eps;
    std::vector<std::vector<double>> sigma_x; // PML conductivity in x direction
    std::vector<std::vector<double>> sigma_y; // PML conductivity in y direction
    int time;

    void initialize_pml() {
        // Initialize PML conductivity arrays
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < SIZE_X; ++i) {
            for (int j = 0; j < SIZE_Y; ++j) {
                double sigma_x_val = 0.0;
                double sigma_y_val = 0.0;
                
                // Left boundary
                if (i < PML_WIDTH) {
                    double dist = (PML_WIDTH - i) / (double)PML_WIDTH;
                    sigma_x_val = PML_SIGMA_MAX * pow(dist, 3.0);
                }
                // Right boundary  
                else if (i >= SIZE_X - PML_WIDTH) {
                    double dist = (i - (SIZE_X - PML_WIDTH - 1)) / (double)PML_WIDTH;
                    sigma_x_val = PML_SIGMA_MAX * pow(dist, 3.0);
                }
                
                // Bottom boundary
                if (j < PML_WIDTH) {
                    double dist = (PML_WIDTH - j) / (double)PML_WIDTH;
                    sigma_y_val = PML_SIGMA_MAX * pow(dist, 3.0);
                }
                // Top boundary
                else if (j >= SIZE_Y - PML_WIDTH) {
                    double dist = (j - (SIZE_Y - PML_WIDTH - 1)) / (double)PML_WIDTH;
                    sigma_y_val = PML_SIGMA_MAX * pow(dist, 3.0);
                }
                
                sigma_x[i][j] = sigma_x_val;
                sigma_y[i][j] = sigma_y_val;
            }
        }
    }

    void update_H() {
        // Update Hx with PML absorption (curl of Ey in z-direction)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < SIZE_X - 1; ++i) {
            for (int j = 0; j < SIZE_Y; ++j) {
                double curl_E = -(Ey[i+1][j] - Ey[i][j]);
                
                // PML absorption factor for Hx (x-direction absorption)
                double absorption_factor = 1.0 / (1.0 + sigma_x[i][j] * DT / (2.0 * MU0));
                double decay_factor = (1.0 - sigma_x[i][j] * DT / (2.0 * MU0)) * absorption_factor;
                
                Hx[i][j] = decay_factor * Hx[i][j] + absorption_factor * (DT / (MU0 * DX)) * curl_E;
            }
        }
        
        // Update Hz with PML absorption (curl of Ey in x-direction)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < SIZE_X; ++i) {
            for (int j = 0; j < SIZE_Y - 1; ++j) {
                double curl_E = (Ey[i][j+1] - Ey[i][j]);
                
                // PML absorption factor for Hz (y-direction absorption)
                double absorption_factor = 1.0 / (1.0 + sigma_y[i][j] * DT / (2.0 * MU0));
                double decay_factor = (1.0 - sigma_y[i][j] * DT / (2.0 * MU0)) * absorption_factor;
                
                Hz[i][j] = decay_factor * Hz[i][j] + absorption_factor * (DT / (MU0 * DX)) * curl_E;
            }
        }
    }

    void update_E() {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < SIZE_X - 1; ++i) {
            for (int j = 1; j < SIZE_Y - 1; ++j) {
                // Standard FDTD update for Ey field with PML absorption
                double curl_H = (Hz[i][j] - Hz[i][j-1]) - (Hx[i][j] - Hx[i-1][j]);
                
                // PML absorption factor
                double absorption_factor = 1.0 / (1.0 + (sigma_x[i][j] + sigma_y[i][j]) * DT / (2.0 * eps[i][j]));
                double decay_factor = (1.0 - (sigma_x[i][j] + sigma_y[i][j]) * DT / (2.0 * eps[i][j])) * absorption_factor;
                
                Ey[i][j] = decay_factor * Ey[i][j] + absorption_factor * (DT / (eps[i][j] * DX)) * curl_H;
            }
        }
    }

    void apply_sources() {
        // Left side external sources - source before first tooth
        Ey[left_source_x][left_source0_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        
        // Left side external sources (10 sources between teeth, frequency f1)
        Ey[left_source_x][left_source1_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source2_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source3_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source4_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source5_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source6_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source7_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source8_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source9_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[left_source_x][left_source10_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;

        // Right side external sources - source before first tooth
        Ey[right_source_x][right_source0_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        
        // Right side external sources (10 sources between teeth, frequency f1)
        Ey[right_source_x][right_source1_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source2_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source3_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source4_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source5_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source6_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source7_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source8_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source9_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
        Ey[right_source_x][right_source10_y] += sin(2.0 * PI * f1 * time * DT) * external_source_intensity;
    }

    // void calculate_potential() {
    //     // For electromagnetic waves, the scalar potential is not the primary quantity
    //     // Instead, we calculate a quasi-static potential from the instantaneous field
    //     // This gives a meaningful visualization of the potential landscape
        
    //     // Method 1: Direct integration of the electric field
    //     // V(x,y) = -∫ E⃗ · dl⃗ from reference point to (x,y)
    //     // Since we only have Ey, we integrate along field lines
        
    //     // Initialize potential with zero at left boundary
    //     #pragma omp parallel for
    //     for (int j = 0; j < SIZE_Y; ++j) {
    //         potential[0][j] = 0.0; // Reference potential at left boundary
    //     }
        
    //     // Integrate Ey field from left to right (this gives potential due to y-direction field)
    //     for (int i = 1; i < SIZE_X; ++i) {
    //         #pragma omp parallel for
    //         for (int j = 0; j < SIZE_Y; ++j) {
    //             // Simple integration: V(i) = V(i-1) - Ey * dy_effective
    //             // Since Ey creates potential differences in y-direction, we integrate along y
                
    //             // For visualization purposes, we create a smoothed potential based on local field
    //             double local_field_contribution = 0.0;
    //             int count = 0;
                
    //             // Average field in local neighborhood for smoother potential
    //             for (int di = -2; di <= 2; ++di) {
    //                 for (int dj = -2; dj <= 2; ++dj) {
    //                     int ni = i + di;
    //                     int nj = j + dj;
    //                     if (ni >= 0 && ni < SIZE_X && nj >= 0 && nj < SIZE_Y) {
    //                         local_field_contribution += Ey[ni][nj];
    //                         count++;
    //                     }
    //                 }
    //             }
                
    //             if (count > 0) {
    //                 local_field_contribution /= count;
    //             }
                
    //             // Integrate the field to get potential
    //             potential[i][j] = potential[i-1][j] - local_field_contribution * DX;
    //         }
    //     }
        
    //     // Method 2: Add contribution from field gradients (Poisson-like correction)
    //     // This helps capture potential variations due to charge concentrations
    //     std::vector<std::vector<double>> potential_correction(SIZE_X, std::vector<double>(SIZE_Y, 0.0));
        
    //     // Calculate correction based on field divergence (quasi-static approximation)
    //     #pragma omp parallel for collapse(2)
    //     for (int i = 1; i < SIZE_X - 1; ++i) {
    //         for (int j = 1; j < SIZE_Y - 1; ++j) {
    //             // Calculate field divergence ∇·E (indicates charge density)
    //             double dEy_dy = (Ey[i][j+1] - Ey[i][j-1]) / (2.0 * DX);
                
    //             // For TM mode (only Ey), we approximate Ex from wave equation
    //             // ∂Ex/∂t ≈ -c²∂Hy/∂z ≈ c²∂Ey/∂y (rough approximation)
    //             double Ex_approx = dEy_dy * 0.1; // Small coupling factor
    //             double dEx_dx = (i > 1 && i < SIZE_X - 2) ? 
    //                 ((Ey[i+1][j+1] - Ey[i+1][j-1]) - (Ey[i-1][j+1] - Ey[i-1][j-1])) / (4.0 * DX) * 0.1 : 0.0;
                
    //             // Divergence of E field
    //             double div_E = dEx_dx + dEy_dy;
                
    //             // Correction potential (small contribution)
    //             potential_correction[i][j] = -div_E * DX * DX * 0.01; // Small factor for stability
    //         }
    //     }
        
    //     // Apply correction with smoothing
    //     #pragma omp parallel for collapse(2)
    //     for (int i = 1; i < SIZE_X - 1; ++i) {
    //         for (int j = 1; j < SIZE_Y - 1; ++j) {
    //             potential[i][j] += potential_correction[i][j];
    //         }
    //     }
        
    //     // Smooth the potential to reduce numerical noise
    //     std::vector<std::vector<double>> potential_smooth = potential;
    //     #pragma omp parallel for collapse(2)
    //     for (int i = 1; i < SIZE_X - 1; ++i) {
    //         for (int j = 1; j < SIZE_Y - 1; ++j) {
    //             potential_smooth[i][j] = 0.6 * potential[i][j] + 
    //                                    0.1 * (potential[i-1][j] + potential[i+1][j] + 
    //                                           potential[i][j-1] + potential[i][j+1]);
    //         }
    //     }
        
    //     // Apply smoothed result
    //     potential = potential_smooth;
        
    //     // Ensure boundary conditions (grounded boundaries)
    //     #pragma omp parallel for
    //     for (int i = 0; i < SIZE_X; ++i) {
    //         potential[i][0] = 0.0;           // Bottom boundary
    //         potential[i][SIZE_Y-1] = 0.0;    // Top boundary
    //     }
    //     #pragma omp parallel for
    //     for (int j = 0; j < SIZE_Y; ++j) {
    //         potential[0][j] = 0.0;           // Left boundary
    //         potential[SIZE_X-1][j] = 0.0;    // Right boundary
    //     }
   // }
};


int main() {
    #ifdef _OPENMP
    std::cout << "OpenMP is enabled with " << omp_get_max_threads() << " threads" << std::endl;
    #else
    std::cout << "OpenMP is not available - running in serial mode" << std::endl;
    #endif
    
    // Print grid information
    std::cout << "=== GRID CONFIGURATION ===" << std::endl;
    std::cout << "Grid spacing (DX): " << DX * 1e9 << " nm" << std::endl;
    std::cout << "Domain X: " << SIZE_X << " points = " << SIZE_X * DX * 1e6 << " μm" << std::endl;
    std::cout << "Domain Y: " << SIZE_Y << " points = " << SIZE_Y * DX * 1e6 << " μm" << std::endl;
    std::cout << "Total grid points: " << (long long)SIZE_X * SIZE_Y << std::endl;
    std::cout << "Points per wavelength: " << lambda1 / DX << std::endl;
    std::cout << "==========================" << std::endl;
    
    FDTDSimulator sim;

    // Save geometry parameters at the beginning
    sim.save_geometry_params("geometry_params.txt");
    std::cout << "Geometry parameters saved to geometry_params.txt" << std::endl;

    for (int t = 0; t < MAX_TIME; ++t) {
        sim.step();

        // Save field data and print progress
        if (t % 50 == 0) {
            double current_time_s = t * DT; // Physical time in seconds
            double current_time_fs = current_time_s * 1e15; // Physical time in femtoseconds
            double current_time_periods = current_time_s * f1; // Time in periods of wavelength 1
            
            std::cout << "Time step: " << t 
                      << " | Physical time: " << current_time_fs << " fs"
                      << " (" << current_time_s << " s)"
                      << " | Periods: " << current_time_periods << std::endl;
            
            sim.save_field_to_file("Ey_field_" + std::to_string(t) + ".dat");
           // sim.save_potential_to_file("potential_" + std::to_string(t) + ".dat");
            
            // Print field statistics for monitoring
            if (t > 0) {
                std::cout << "  Field data saved." << std::endl;
            }
        }
    }

    std::cout << "Simulation finished." << std::endl;
    
    // Print final simulation summary
    double total_sim_time_s = MAX_TIME * DT;
    double total_sim_time_fs = total_sim_time_s * 1e15;
    double total_periods = total_sim_time_s * f1;
    
    std::cout << "\n=== SIMULATION SUMMARY ===" << std::endl;
    std::cout << "Total timesteps: " << MAX_TIME << std::endl;
    std::cout << "Total physical time: " << total_sim_time_fs << " fs (" << total_sim_time_s << " s)" << std::endl;
    std::cout << "Total periods simulated: " << total_periods << " periods of λ₁" << std::endl;
    std::cout << "Time step size: " << DT * 1e15 << " fs (" << DT << " s)" << std::endl;
    std::cout << "==========================" << std::endl;

    return 0;
}
