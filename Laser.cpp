#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

// --- Constants ---
const double PI = 3.14159265358979323846;
const double C0 = 299792458.0; // Speed of light in vacuum (m/s)
const double MU0 = 4.0 * PI * 1e-7; // Permeability of free space (H/m)
const double EPS0 = 1.0 / (MU0 * C0 * C0); // Permittivity of free space (F/m)

// --- Laser Source Parameters (defining wavelengths first) ---
const double lambda1 = 2e-6; // Wavelength of laser 1 (m) = 1500 nm
const double lambda2 = 2e-6; // Wavelength of laser 2 (m) = 1500 nm
const double f1 = C0 / lambda1; // Frequency of laser 1 (Hz)
const double f2 = C0 / lambda2; // Frequency of laser 2 (Hz)
const double intensity1 = 5.0e9; // Intensity of laser 1 (V/m) - Target: 5 GV/m
const double intensity2 = 5.0e9; // Intensity of laser 2 (V/m) - Target: 5 GV/m

// --- Simulation Parameters (in terms of wavelength) ---
const double DX = lambda1 / 50.0; // Spatial step = λ₁/30 ≈ 0.05 μm (excellent resolution)
const double DT = DX / (C0 * sqrt(2.0)); // Time step (s) - Courant stability condition

const int SIZE_X = (int)(6.0 * lambda1 / DX); // Grid size ≈ 6λ₁ in x direction (compact around structure)
const int SIZE_Y = (int)(6.0 * lambda1 / DX); // Grid size ≈ 6λ₁ in y direction (compact around structure)  
const int MAX_TIME = 1000; // Total simulation time steps (dimensionless)

// --- Absorbing Boundary Conditions (PML) ---
const int PML_WIDTH = 10; // PML layer thickness (grid points)
const double PML_SIGMA_MAX = 0.8 * (3.0 + 1.0) / (377.0 * DX); // Maximum conductivity

// --- Dielectric Pillar Parameters (in terms of wavelength) ---
// Relative permittivity of pillars of SiO2  
const double eps_r = 3.9; // Relative permittivity of the dielectric (dimensionless)
const int num_pillars_per_column = 5; // Number of pillars per column (dimensionless)
const int pillar_radius = (int)(0.2 * lambda1 / DX); // Pillar radius = 0.2λ₁
const int pillar_spacing_y = (int)(0.5 * lambda1 / DX); // Pillar spacing = 0.5λ₁
const int column1_x = (int)(3.0 * lambda1 / DX); // Column 1 at x = 3.1λ₁ (centered in domain)
const int column2_x = (int)(3.5 * lambda1 / DX); // Column 2 at x = 3.5λ₁ (centered in domain)
const int first_pillar_y = (int)(1.5 * lambda1 / DX); // First pillar at y = 1.5λ₁

// --- Source Positions (in terms of wavelength) ---
// Laser sources at the left boundary for x-direction propagation
const int source1_x = PML_WIDTH + 5; // Laser 1 near left boundary (after PML)
const int source1_y = (int)(2.5 * lambda1 / DX); // Laser 1 at y = 2.5λ₁ (centered on pillar array)
const int source2_x = PML_WIDTH + 5; // Laser 2 near left boundary (after PML)  
const int source2_y = (int)(3.5 * lambda1 / DX); // Laser 2 at y = 3.5λ₁ (offset from laser 1)

// --- Inter-pillar Source Parameters ---
const double inter_pillar_intensity = 2.0e9; // 2 GV/m intensity for inter-pillar sources

// Inter-pillar source positions in column 1 (exactly between pillars)
const int inter_source1_1_x = column1_x; // Centered on pillar column
const int inter_source1_1_y = first_pillar_y + (int)(0.5 * pillar_spacing_y); // Exactly between pillars 1 and 2
const int inter_source1_2_x = column1_x;
const int inter_source1_2_y = first_pillar_y + (int)(1.5 * pillar_spacing_y); // Exactly between pillars 2 and 3
const int inter_source1_3_x = column1_x;
const int inter_source1_3_y = first_pillar_y + (int)(2.5 * pillar_spacing_y); // Exactly between pillars 3 and 4
const int inter_source1_4_x = column1_x;
const int inter_source1_4_y = first_pillar_y + (int)(3.5 * pillar_spacing_y); // Exactly between pillars 4 and 5

// Inter-pillar source positions in column 2 (exactly between pillars)
const int inter_source2_1_x = column2_x; // Centered on pillar column
const int inter_source2_1_y = first_pillar_y + (int)(0.5 * pillar_spacing_y); // Exactly between pillars 1 and 2
const int inter_source2_2_x = column2_x;
const int inter_source2_2_y = first_pillar_y + (int)(1.5 * pillar_spacing_y); // Exactly between pillars 2 and 3
const int inter_source2_3_x = column2_x;
const int inter_source2_3_y = first_pillar_y + (int)(2.5 * pillar_spacing_y); // Exactly between pillars 3 and 4
const int inter_source2_4_x = column2_x;
const int inter_source2_4_y = first_pillar_y + (int)(3.5 * pillar_spacing_y); // Exactly between pillars 4 and 5

// --- Main FDTD Class ---
class FDTDSimulator {
public:
    FDTDSimulator() :
        Ey(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        Hx(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        Hz(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        eps(SIZE_X, std::vector<double>(SIZE_Y, EPS0)),
        sigma_x(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        sigma_y(SIZE_X, std::vector<double>(SIZE_Y, 0.0)),
        time(0)
    {
        // Initialize PML absorbing boundaries
        initialize_pml();
        
        // Initialize the dielectric pillars
        // We iterate through every grid point and check if it falls within any of the pillars.
        for (int i = 0; i < SIZE_X; ++i) {
            for (int j = 0; j < SIZE_Y; ++j) {
                bool in_pillar = false;
                // Check against all pillars in both columns
                for (int p = 0; p < num_pillars_per_column; ++p) {
                    int pillar_center_y = first_pillar_y + p * pillar_spacing_y;
                    
                    // Check if the point (i, j) is inside a pillar in the first column
                    if (sqrt(pow(i - column1_x, 2) + pow(j - pillar_center_y, 2)) < pillar_radius) {
                        in_pillar = true;
                        break; // Point is in a pillar, no need to check others
                    }
                    
                    // Check if the point (i, j) is inside a pillar in the second column
                    if (sqrt(pow(i - column2_x, 2) + pow(j - pillar_center_y, 2)) < pillar_radius) {
                        in_pillar = true;
                        break; // Point is in a pillar, no need to check others
                    }
                }
                
                if (in_pillar) {
                    eps[i][j] = eps_r * EPS0;
                }
            }
        }
    }

    void step() {
        update_H();
        update_E();
        apply_sources();
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
        file << "DX_um " << DX * 1e6 << "\n"; // Grid spacing in micrometers
        file << "DX_over_lambda1 " << DX / lambda1 << "\n"; // Grid spacing / λ₁
        file << "points_per_lambda1 " << lambda1 / DX << "\n"; // Points per λ₁
        file << "points_per_lambda2 " << lambda2 / DX << "\n"; // Points per λ₂
        
        // Physical domain size (in multiple units)
        file << "domain_x_um " << SIZE_X * DX * 1e6 << "\n"; // Domain size in μm
        file << "domain_y_um " << SIZE_Y * DX * 1e6 << "\n"; // Domain size in μm
        file << "domain_x_lambda1 " << SIZE_X * DX / lambda1 << "\n"; // Domain size in λ₁
        file << "domain_y_lambda1 " << SIZE_Y * DX / lambda1 << "\n"; // Domain size in λ₁
        
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
        file << "inter_pillar_intensity " << inter_pillar_intensity << "\n"; // Inter-pillar source intensity in V/m
        file << "inter_pillar_intensity_GVm " << inter_pillar_intensity * 1e-9 << "\n"; // Inter-pillar intensity in GV/m
        file << "total_sources " << 8 << "\n"; // Total number of laser sources (8 inter-pillar only)
        
        // Save dielectric parameters
        file << "eps_r " << eps_r << "\n"; // Relative permittivity (dimensionless)
        file << "num_pillars_per_column " << num_pillars_per_column << "\n";
        
        // Pillar dimensions (grid points, μm, and λ₁)
        file << "pillar_radius " << pillar_radius << "\n"; // Grid points
        file << "pillar_radius_um " << pillar_radius * DX * 1e6 << "\n"; // Radius in μm
        file << "pillar_radius_lambda1 " << pillar_radius * DX / lambda1 << "\n"; // Radius in λ₁
        file << "pillar_spacing_y " << pillar_spacing_y << "\n"; // Grid points
        file << "pillar_spacing_y_um " << pillar_spacing_y * DX * 1e6 << "\n"; // Spacing in μm
        file << "pillar_spacing_y_lambda1 " << pillar_spacing_y * DX / lambda1 << "\n"; // Spacing in λ₁
        
        // Column positions (grid points, μm, and λ₁)
        file << "column1_x " << column1_x << "\n"; // Grid position
        file << "column2_x " << column2_x << "\n"; // Grid position
        file << "column1_x_um " << column1_x * DX * 1e6 << "\n"; // Position in μm
        file << "column2_x_um " << column2_x * DX * 1e6 << "\n"; // Position in μm
        file << "column1_x_lambda1 " << column1_x * DX / lambda1 << "\n"; // Position in λ₁
        file << "column2_x_lambda1 " << column2_x * DX / lambda1 << "\n"; // Position in λ₁
        file << "column_separation_um " << (column2_x - column1_x) * DX * 1e6 << "\n"; // Separation in μm
        file << "column_separation_lambda1 " << (column2_x - column1_x) * DX / lambda1 << "\n"; // Separation in λ₁
        
        // First pillar position (grid points, μm, and λ₁)
        file << "first_pillar_y " << first_pillar_y << "\n"; // Grid position
        file << "first_pillar_y_um " << first_pillar_y * DX * 1e6 << "\n"; // Position in μm
        file << "first_pillar_y_lambda1 " << first_pillar_y * DX / lambda1 << "\n"; // Position in λ₁
        
        // Design ratios and scaling factors
        file << "pillar_diameter_over_spacing " << (2.0 * pillar_radius) / (double)pillar_spacing_y << "\n"; // Fill factor
        file << "domain_aspect_ratio " << (double)SIZE_X / (double)SIZE_Y << "\n"; // Domain aspect ratio
        
        // PML Absorbing Boundary Conditions
        file << "PML_WIDTH " << PML_WIDTH << "\n"; // PML layer thickness (grid points)
        file << "PML_WIDTH_um " << PML_WIDTH * DX * 1e6 << "\n"; // PML thickness in μm
        file << "PML_WIDTH_lambda1 " << PML_WIDTH * DX / lambda1 << "\n"; // PML thickness in λ₁
        file << "PML_SIGMA_MAX " << PML_SIGMA_MAX << "\n"; // Maximum PML conductivity
        
        // Physical constants used
        file << "C0 " << C0 << "\n"; // Speed of light (m/s)
        file << "MU0 " << MU0 << "\n"; // Permeability of free space (H/m)
        file << "EPS0 " << EPS0 << "\n"; // Permittivity of free space (F/m)
        
        file.close();
    }

private:
    std::vector<std::vector<double>> Ey;
    std::vector<std::vector<double>> Hx;
    std::vector<std::vector<double>> Hz;
    std::vector<std::vector<double>> eps;
    std::vector<std::vector<double>> sigma_x; // PML conductivity in x direction
    std::vector<std::vector<double>> sigma_y; // PML conductivity in y direction
    int time;

    void initialize_pml() {
        // Initialize PML conductivity arrays
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
        // Inter-pillar sources in column 1 (2 GV/m each, frequency f1)
        Ey[inter_source1_1_x][inter_source1_1_y] += sin(2.0 * PI * f1 * time * DT) * inter_pillar_intensity;
        Ey[inter_source1_2_x][inter_source1_2_y] += sin(2.0 * PI * f1 * time * DT) * inter_pillar_intensity;
        Ey[inter_source1_3_x][inter_source1_3_y] += sin(2.0 * PI * f1 * time * DT) * inter_pillar_intensity;
        Ey[inter_source1_4_x][inter_source1_4_y] += sin(2.0 * PI * f1 * time * DT) * inter_pillar_intensity;
        
        // Inter-pillar sources in column 2 (2 GV/m each, frequency f2)
        Ey[inter_source2_1_x][inter_source2_1_y] += sin(2.0 * PI * f2 * time * DT) * inter_pillar_intensity;
        Ey[inter_source2_2_x][inter_source2_2_y] += sin(2.0 * PI * f2 * time * DT) * inter_pillar_intensity;
        Ey[inter_source2_3_x][inter_source2_3_y] += sin(2.0 * PI * f2 * time * DT) * inter_pillar_intensity;
        Ey[inter_source2_4_x][inter_source2_4_y] += sin(2.0 * PI * f2 * time * DT) * inter_pillar_intensity;
    }
};


int main() {
    FDTDSimulator sim;

    // Save geometry parameters at the beginning
    sim.save_geometry_params("geometry_params.txt");
    std::cout << "Geometry parameters saved to geometry_params.txt" << std::endl;

    for (int t = 0; t < MAX_TIME; ++t) {
        sim.step();
        if (t % 50 == 0) {
            std::cout << "Time step: " << t << std::endl;
            sim.save_field_to_file("Ey_field_" + std::to_string(t) + ".dat");
        }
    }

    std::cout << "Simulation finished." << std::endl;

    return 0;
}
