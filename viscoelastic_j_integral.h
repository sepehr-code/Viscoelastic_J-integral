#ifndef VISCOELASTIC_J_INTEGRAL_H
#define VISCOELASTIC_J_INTEGRAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Mathematical constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Maximum number of Prony terms for viscoelastic model
#define MAX_PRONY_TERMS 10

// Data structure for 2D coordinates
typedef struct {
    double x;
    double y;
} Point2D;

// Data structure for stress tensor (2D)
typedef struct {
    double s11;  // σ₁₁
    double s12;  // σ₁₂
    double s21;  // σ₂₁
    double s22;  // σ₂₂
} StressTensor2D;

// Data structure for strain tensor (2D)
typedef struct {
    double e11;  // ε₁₁
    double e12;  // ε₁₂
    double e21;  // ε₂₁
    double e22;  // ε₂₂
} StrainTensor2D;

// Data structure for displacement field
typedef struct {
    double u1;   // u₁ displacement
    double u2;   // u₂ displacement
} Displacement2D;

// Data structure for displacement gradients
typedef struct {
    double du1_dx1;  // ∂u₁/∂x₁
    double du1_dx2;  // ∂u₁/∂x₂
    double du2_dx1;  // ∂u₂/∂x₁
    double du2_dx2;  // ∂u₂/∂x₂
} DisplacementGradient2D;

// Data structure for contour node
typedef struct {
    Point2D coord;                    // Node coordinates
    Displacement2D displacement;      // Displacements at node
    StressTensor2D stress;           // Stress tensor at node
    StrainTensor2D strain;           // Strain tensor at node
    DisplacementGradient2D du_dx;    // Displacement gradients
    Point2D normal;                  // Outward normal vector
    double ds;                       // Arc length element
} ContourNode;

// Data structure for contour
typedef struct {
    ContourNode *nodes;              // Array of contour nodes
    int num_nodes;                   // Number of nodes
    double total_length;             // Total contour length
} Contour;

// Data structure for Prony series term
typedef struct {
    double G_k;                      // Shear modulus for k-th term
    double tau_k;                    // Relaxation time for k-th term
    StressTensor2D *stress_history;  // Stress history for this term
} PronyTerm;

// Data structure for viscoelastic model
typedef struct {
    double G_inf;                    // Long-term shear modulus
    int num_prony_terms;             // Number of Prony series terms
    PronyTerm *prony_terms;          // Array of Prony terms
    double nu;                       // Poisson's ratio (assumed time-independent)
} ViscoelasticModel;

// Data structure for time-stepping parameters
typedef struct {
    double dt;                       // Time increment
    double current_time;             // Current simulation time
    double total_time;               // Total simulation time
    int time_step;                   // Current time step number
} TimeParameters;

// Data structure for J-integral results
typedef struct {
    double j_value;                  // J-integral value
    double time;                     // Time at which J was computed
    double strain_energy_integral;   // Strain energy contribution
    double stress_traction_integral; // Stress-traction contribution
} JIntegralResult;

// Function declarations

// Memory management
ContourNode* allocate_contour_nodes(int num_nodes);
void free_contour_nodes(ContourNode *nodes);
ViscoelasticModel* create_viscoelastic_model(int num_prony_terms);
void free_viscoelastic_model(ViscoelasticModel *model);
int allocate_stress_history(ViscoelasticModel *model, int num_nodes);
JIntegralResult* allocate_results_array(int num_time_steps);
void free_results_array(JIntegralResult *results);

// Contour initialization
void initialize_circular_contour(Contour *contour, Point2D center, double radius, int num_nodes);
void initialize_rectangular_contour(Contour *contour, Point2D bottom_left, Point2D top_right, int num_nodes);
void compute_contour_normals(Contour *contour);
void compute_arc_lengths(Contour *contour);
void create_concentric_contours(Contour **contours, int num_contours, Point2D center, 
                               double min_radius, double max_radius);

// Field initialization (analytical or from external data)
void initialize_displacement_field(Contour *contour, Point2D crack_tip);
void initialize_stress_field(Contour *contour, Point2D crack_tip);
void initialize_strain_field(Contour *contour, Point2D crack_tip);
void compute_displacement_gradients(Contour *contour);

// Viscoelastic constitutive modeling
void initialize_prony_series(ViscoelasticModel *model, double *G_values, double *tau_values);
void setup_asphalt_model(ViscoelasticModel *model);
void update_viscoelastic_stress(ViscoelasticModel *model, Contour *contour, 
                               StrainTensor2D *strain_increment, double dt);
void recursive_stress_update(PronyTerm *term, StrainTensor2D *strain_increment, 
                           StressTensor2D *total_stress, double dt, int node_index);
int validate_model_parameters(ViscoelasticModel *model);

// Strain energy calculations
double compute_strain_energy_density(StressTensor2D stress, StrainTensor2D strain);
double compute_total_strain_energy(Contour *contour);

// J-integral evaluation
double compute_j_integral(Contour *contour);
double compute_strain_energy_contribution(Contour *contour);
double compute_stress_traction_contribution(Contour *contour);
JIntegralResult evaluate_j_integral_at_time(Contour *contour, double time);
double compute_j_integral_midpoint(Contour *contour);
double compute_j_integral_rate(JIntegralResult *results, int num_results, int current_index);
void analyze_j_evolution(JIntegralResult *results, int num_results);

// Path independence verification
double compare_j_integral_paths(Contour *contour1, Contour *contour2);
void validate_path_independence(Contour **contours, int num_contours, double tolerance);

// Time integration
void time_step_update(ViscoelasticModel *model, Contour *contour, 
                     TimeParameters *time_params);
void simulate_viscoelastic_evolution(ViscoelasticModel *model, Contour *contour,
                                   TimeParameters *time_params, JIntegralResult **results);
void apply_load_history(ViscoelasticModel *model, Contour *contour, 
                       double *load_times, double *load_values, int num_loads,
                       TimeParameters *time_params, JIntegralResult **results);
void relaxation_test(ViscoelasticModel *model, Contour *contour,
                    double strain_level, TimeParameters *time_params,
                    JIntegralResult **results);
TimeParameters create_time_parameters(double dt, double total_time);
double adaptive_time_step(JIntegralResult *previous_results, int num_previous, 
                         double current_dt, double target_accuracy);

// Utility functions
void print_contour_info(Contour *contour);
void print_j_integral_result(JIntegralResult result);
void save_results_to_file(JIntegralResult *results, int num_steps, const char *filename);
void print_viscoelastic_model(ViscoelasticModel *model);
void save_contour_to_file(Contour *contour, const char *filename);
int load_contour_from_file(Contour *contour, const char *filename);
void create_simulation_report(ViscoelasticModel *model, Contour *contour,
                             JIntegralResult *results, int num_steps,
                             const char *filename);
void print_banner();

// Mathematical utilities
double kronecker_delta(int i, int j);
void matrix_multiply_2x2(double A[2][2], double B[2][2], double C[2][2]);
double tensor_double_contraction(StressTensor2D stress, StrainTensor2D strain);

#endif // VISCOELASTIC_J_INTEGRAL_H 