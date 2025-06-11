# Viscoelastic J-Integral Calculator

A comprehensive C implementation for calculating the time-dependent J-integral in viscoelastic materials using contour integration methods.

## 🎯 Overview

This project implements a **time-dependent J-integral** for viscoelastic materials, providing a robust framework for fracture mechanics analysis in materials that exhibit significant time-dependent behavior such as asphalt concrete, polymers, and biological tissues.

### Key Features

- **🧮 Generalized Maxwell Viscoelastic Model**: Full Prony series implementation with recursive stress updating
- **🔄 Memory-Efficient Algorithm**: Avoids storing complete strain history using exponential decay recursion
- **🛣️ Path Independence Validation**: Numerical verification across multiple contour paths
- **⏱️ Time-Dependent Evolution**: Complete simulation framework for field evolution
- **📊 Multiple Test Scenarios**: Relaxation tests, load history application, and evolution analysis
- **🎯 High Performance**: Optimized C implementation for computational efficiency

## 📐 Theoretical Background

### J-Integral Formulation

The J-integral is defined as:

```
J = ∫_Γ [W δ₁ⱼ - σᵢⱼ (∂uᵢ/∂x₁)] nⱼ ds
```

Where:
- `W = ½ σᵢⱼ εᵢⱼ` is the strain energy density
- `σᵢⱼ` is the stress tensor
- `εᵢⱼ` is the strain tensor
- `uᵢ` are the displacement components
- `nⱼ` are the outward normal components
- `Γ` is the integration contour

### Viscoelastic Constitutive Model

The generalized Maxwell model with Prony series:

```
σᵢⱼ(t) = G∞ εᵢⱼ(t) + Σₖ Gₖ ∫₀ᵗ exp[-(t-τ)/τₖ] (dεᵢⱼ/dτ) dτ
```

### Recursive Stress Update

For computational efficiency, stress updates use:

```
σᵢⱼᵏ(t+Δt) = exp(-Δt/τₖ) σᵢⱼᵏ(t) + Gₖ[1-exp(-Δt/τₖ)] Δεᵢⱼ
```

## 🏗️ Project Structure

```
viscoelastic_j_integral/
├── viscoelastic_j_integral.h      # Main header with all declarations
├── memory_management.c            # Memory allocation and cleanup
├── contour_operations.c           # Contour initialization and geometry
├── field_initialization.c        # Analytical field initialization
├── viscoelastic_model.c          # Constitutive modeling and stress updates
├── j_integral_calculation.c      # J-integral computation and analysis
├── time_integration.c            # Time-stepping and simulation control
├── utilities.c                   # I/O, printing, and helper functions
├── main.c                        # Demonstration program
├── Makefile                      # Build system
└── README.md                     # This documentation
```

## 🚀 Quick Start

### Prerequisites

- GCC compiler (C99 standard or later)
- Make build system
- Math library (`libm`)

### Compilation

```bash
# Clone or download the project
cd viscoelastic_j_integral

# Build the project
make

# Or build with debug symbols
make debug

# Or build optimized release version
make release
```

### Running Demonstrations

```bash
# Run all demonstrations
./viscoelastic_j_integral 0

# Run specific demonstration
./viscoelastic_j_integral 1    # Basic J-integral calculation
./viscoelastic_j_integral 2    # Path independence validation
./viscoelastic_j_integral 3    # Viscoelastic evolution
./viscoelastic_j_integral 4    # Relaxation test
./viscoelastic_j_integral 5    # Load history application

# Or use Make targets
make demo1    # Basic J-integral
make demo2    # Path independence
make demo3    # Evolution simulation
make demo4    # Relaxation test
make demo5    # Load history
```

## 📋 Demonstration Scenarios

### 1. Basic J-Integral Calculation
- **Purpose**: Validate implementation against analytical solutions
- **Features**: Mode I crack tip fields, theoretical comparison
- **Output**: `demo1_contour.dat`

### 2. Path Independence Validation
- **Purpose**: Verify numerical path independence
- **Features**: Multiple concentric contours, statistical analysis
- **Output**: Console validation report

### 3. Viscoelastic Evolution Simulation
- **Purpose**: Time-dependent J-integral evolution
- **Features**: Asphalt concrete model, complete time history
- **Output**: `demo3_j_evolution.dat`, `demo3_report.txt`

### 4. Stress Relaxation Test
- **Purpose**: Validate relaxation behavior
- **Features**: Constant strain, stress decay analysis
- **Output**: `demo4_relaxation.dat`

### 5. Load History Application
- **Purpose**: Multi-step loading scenarios
- **Features**: Step loading, complex load histories
- **Output**: `demo5_load_history.dat`

## 📊 Output Files

### Data Files (.dat)
- **Format**: Space-separated columns with headers
- **Content**: Time series data for J-integral components
- **Usage**: Import into plotting software (Python, MATLAB, etc.)

### Report Files (.txt)
- **Content**: Simulation parameters, model details, summary statistics
- **Usage**: Documentation and result verification

## 🧪 Advanced Usage

### Custom Material Models

```c
// Create custom viscoelastic model
ViscoelasticModel *model = create_viscoelastic_model(3);
model->G_inf = 1.0e8;  // Long-term modulus
model->nu = 0.3;       // Poisson's ratio

// Define Prony series
double G_values[] = {2.0e8, 1.0e8, 0.5e8};
double tau_values[] = {1.0, 10.0, 100.0};
initialize_prony_series(model, G_values, tau_values);
```

### Custom Contours

```c
// Circular contour
Contour contour;
Point2D center = {0.0, 0.0};
initialize_circular_contour(&contour, center, 0.01, 32);

// Rectangular contour
Point2D bottom_left = {-0.02, -0.02};
Point2D top_right = {0.02, 0.02};
initialize_rectangular_contour(&contour, bottom_left, top_right, 64);
```

### Loading from External Data

```c
// Load field data from finite element results
load_fields_from_file(&contour, "fe_results.dat");

// Load contour geometry
load_contour_from_file(&contour, "contour_nodes.dat");
```

## 🔧 Build Options

```bash
# Standard build
make

# Debug build (with symbols and no optimization)
make debug

# Release build (optimized)
make release

# Memory leak checking
make memcheck

# Static analysis
make check

# Clean build artifacts
make clean

# Clean everything including outputs
make distclean

# Install system-wide
sudo make install
```

## 📈 Performance Considerations

- **Memory Usage**: O(N×K) where N = nodes, K = Prony terms
- **Computational Complexity**: O(N×K×T) for T time steps
- **Optimization**: Compiled with `-O2` for production builds
- **Parallel Processing**: Single-threaded (suitable for OpenMP extension)

## 🧮 Mathematical Validation

### Theoretical Benchmarks
- **Mode I Crack**: J = K₁²(1-ν²)/E for plane strain
- **Path Independence**: <5% variation across contours
- **Relaxation**: Exponential decay validation

### Numerical Accuracy
- **Integration**: Midpoint rule for improved accuracy
- **Time Stepping**: Adaptive time stepping available
- **Convergence**: Automatic validation checks

## 📚 References

1. **Rice, J.R.** (1968). A path independent integral and the approximate analysis of strain concentration by notches and cracks. *Journal of Applied Mechanics*, 35(2), 379-386.

2. **Schapery, R.A.** (1984). Correspondence principles and a generalized J integral for large deformation and fracture analysis of viscoelastic media. *International Journal of Fracture*, 25(3), 195-223.

3. **Christensen, R.M.** (2012). *Theory of Viscoelasticity: An Introduction*. Academic Press.

4. **Anderson, T.L.** (2017). *Fracture Mechanics: Fundamentals and Applications*. CRC Press.

## 🤝 Contributing

This implementation provides a solid foundation for research and development in viscoelastic fracture mechanics. Potential extensions include:

- **3D Implementation**: Extension to three-dimensional problems
- **Parallel Processing**: OpenMP or MPI parallelization
- **Advanced Materials**: Nonlinear viscoelasticity, damage coupling
- **GUI Interface**: Graphical user interface for parameter input
- **FE Integration**: Direct coupling with finite element codes

## 📄 License

This project is provided for educational and research purposes. Please cite appropriately in academic work.

## 🐛 Known Limitations

- **2D Only**: Current implementation limited to plane strain/stress
- **Small Strains**: Assumes small deformation theory
- **Simplified Loading**: Basic strain increment application
- **No Crack Growth**: Static crack configuration

## 💡 Tips for Users

1. **Start Simple**: Begin with Demo 1 to verify installation
2. **Check Path Independence**: Always validate with Demo 2
3. **Monitor Memory**: Use `make memcheck` for large simulations
4. **Visualize Results**: Import `.dat` files into plotting software
5. **Parameter Studies**: Modify material parameters for different materials

---

*For questions, suggestions, or collaboration opportunities, please refer to the project documentation or contact the development team.* 