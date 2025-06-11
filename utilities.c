#include "viscoelastic_j_integral.h"

// Print contour information
void print_contour_info(Contour *contour) {
    if (contour == NULL) {
        printf("Contour: NULL\n");
        return;
    }
    
    printf("\n=== Contour Information ===\n");
    printf("Number of nodes: %d\n", contour->num_nodes);
    printf("Total length: %.6e m\n", contour->total_length);
    
    if (contour->nodes != NULL && contour->num_nodes > 0) {
        printf("First node: (%.3f, %.3f)\n", 
               contour->nodes[0].coord.x, contour->nodes[0].coord.y);
        printf("Last node: (%.3f, %.3f)\n", 
               contour->nodes[contour->num_nodes-1].coord.x, 
               contour->nodes[contour->num_nodes-1].coord.y);
        
        // Find bounding box
        double x_min = contour->nodes[0].coord.x;
        double x_max = contour->nodes[0].coord.x;
        double y_min = contour->nodes[0].coord.y;
        double y_max = contour->nodes[0].coord.y;
        
        for (int i = 1; i < contour->num_nodes; i++) {
            double x = contour->nodes[i].coord.x;
            double y = contour->nodes[i].coord.y;
            
            if (x < x_min) x_min = x;
            if (x > x_max) x_max = x;
            if (y < y_min) y_min = y;
            if (y > y_max) y_max = y;
        }
        
        printf("Bounding box: (%.3f, %.3f) to (%.3f, %.3f)\n", 
               x_min, y_min, x_max, y_max);
    }
    printf("===========================\n\n");
}

// Print J-integral result
void print_j_integral_result(JIntegralResult result) {
    printf("J-Integral Result at t = %.3e s:\n", result.time);
    printf("  Total J-integral: %.6e N/m\n", result.j_value);
    printf("  Strain energy:    %.6e N/m\n", result.strain_energy_integral);
    printf("  Stress-traction:  %.6e N/m\n", result.stress_traction_integral);
    printf("  Ratio SE/ST:      %.3f\n", 
           result.strain_energy_integral / (result.stress_traction_integral + 1e-12));
}

// Save results to file
void save_results_to_file(JIntegralResult *results, int num_steps, const char *filename) {
    if (results == NULL || num_steps <= 0 || filename == NULL) {
        fprintf(stderr, "Error: Invalid parameters for saving results\n");
        return;
    }
    
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", filename);
        return;
    }
    
    // Write header
    fprintf(file, "# Viscoelastic J-Integral Results\n");
    fprintf(file, "# Columns: time(s), J_total(N/m), J_strain_energy(N/m), J_stress_traction(N/m)\n");
    
    // Write data
    for (int i = 0; i < num_steps; i++) {
        fprintf(file, "%.6e %.6e %.6e %.6e\n",
                results[i].time,
                results[i].j_value,
                results[i].strain_energy_integral,
                results[i].stress_traction_integral);
    }
    
    fclose(file);
    printf("Results saved to %s (%d data points)\n", filename, num_steps);
}

// Print viscoelastic model parameters
void print_viscoelastic_model(ViscoelasticModel *model) {
    if (model == NULL) {
        printf("Viscoelastic Model: NULL\n");
        return;
    }
    
    printf("\n=== Viscoelastic Model ===\n");
    printf("Long-term modulus (G_inf): %.3e Pa\n", model->G_inf);
    printf("Poisson's ratio (nu): %.3f\n", model->nu);
    printf("Number of Prony terms: %d\n", model->num_prony_terms);
    
    printf("\nProny Series Terms:\n");
    printf("  Term    G_k (Pa)      tau_k (s)\n");
    printf("  ----    --------      ---------\n");
    
    for (int i = 0; i < model->num_prony_terms; i++) {
        printf("  %2d      %.3e     %.3e\n", 
               i, model->prony_terms[i].G_k, model->prony_terms[i].tau_k);
    }
    
    // Calculate total instantaneous modulus
    double G_0 = model->G_inf;
    for (int i = 0; i < model->num_prony_terms; i++) {
        G_0 += model->prony_terms[i].G_k;
    }
    
    printf("\nInstantaneous modulus (G_0): %.3e Pa\n", G_0);
    printf("Relaxation ratio (G_inf/G_0): %.3f\n", model->G_inf / G_0);
    printf("==========================\n\n");
}

// Mathematical utility: Kronecker delta
double kronecker_delta(int i, int j) {
    return (i == j) ? 1.0 : 0.0;
}

// Mathematical utility: 2x2 matrix multiplication
void matrix_multiply_2x2(double A[2][2], double B[2][2], double C[2][2]) {
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
}

// Mathematical utility: tensor double contraction
double tensor_double_contraction(StressTensor2D stress, StrainTensor2D strain) {
    return stress.s11 * strain.e11 + 
           stress.s22 * strain.e22 + 
           stress.s12 * strain.e12 + 
           stress.s21 * strain.e21;
}

// Save contour data to file
void save_contour_to_file(Contour *contour, const char *filename) {
    if (contour == NULL || filename == NULL) {
        fprintf(stderr, "Error: Invalid parameters for saving contour\n");
        return;
    }
    
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", filename);
        return;
    }
    
    // Write header
    fprintf(file, "# Contour Data\n");
    fprintf(file, "# Columns: x(m), y(m), u1(m), u2(m), s11(Pa), s22(Pa), s12(Pa), e11, e22, e12, nx, ny, ds(m)\n");
    
    // Write node data
    for (int i = 0; i < contour->num_nodes; i++) {
        ContourNode *node = &contour->nodes[i];
        fprintf(file, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                node->coord.x, node->coord.y,
                node->displacement.u1, node->displacement.u2,
                node->stress.s11, node->stress.s22, node->stress.s12,
                node->strain.e11, node->strain.e22, node->strain.e12,
                node->normal.x, node->normal.y, node->ds);
    }
    
    fclose(file);
    printf("Contour data saved to %s (%d nodes)\n", filename, contour->num_nodes);
}

// Load contour from file
int load_contour_from_file(Contour *contour, const char *filename) {
    if (contour == NULL || filename == NULL) {
        fprintf(stderr, "Error: Invalid parameters for loading contour\n");
        return -1;
    }
    
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for reading\n", filename);
        return -1;
    }
    
    // Count lines (excluding comments)
    char line[512];
    int node_count = 0;
    
    while (fgets(line, sizeof(line), file) != NULL) {
        if (line[0] != '#' && strlen(line) > 10) {
            node_count++;
        }
    }
    
    if (node_count == 0) {
        fprintf(stderr, "Error: No data found in file %s\n", filename);
        fclose(file);
        return -1;
    }
    
    // Allocate nodes
    contour->nodes = allocate_contour_nodes(node_count);
    if (contour->nodes == NULL) {
        fclose(file);
        return -1;
    }
    contour->num_nodes = node_count;
    
    // Rewind and read data
    rewind(file);
    int node_idx = 0;
    
    while (fgets(line, sizeof(line), file) != NULL && node_idx < node_count) {
        if (line[0] == '#') continue;
        
        ContourNode *node = &contour->nodes[node_idx];
        int items = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                          &node->coord.x, &node->coord.y,
                          &node->displacement.u1, &node->displacement.u2,
                          &node->stress.s11, &node->stress.s22, &node->stress.s12,
                          &node->strain.e11, &node->strain.e22, &node->strain.e12,
                          &node->normal.x, &node->normal.y, &node->ds);
        
        if (items >= 2) {
            // Set symmetry for stress and strain tensors
            node->stress.s21 = node->stress.s12;
            node->strain.e21 = node->strain.e12;
            node_idx++;
        }
    }
    
    fclose(file);
    
    // Compute total length
    compute_arc_lengths(contour);
    
    printf("Loaded contour from %s (%d nodes)\n", filename, node_idx);
    return 0;
}

// Create summary report
void create_simulation_report(ViscoelasticModel *model, Contour *contour,
                             JIntegralResult *results, int num_steps,
                             const char *filename) {
    if (model == NULL || contour == NULL || results == NULL || filename == NULL) {
        fprintf(stderr, "Error: Invalid parameters for creating report\n");
        return;
    }
    
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot create report file %s\n", filename);
        return;
    }
    
    fprintf(file, "# VISCOELASTIC J-INTEGRAL SIMULATION REPORT\n");
    fprintf(file, "# ==========================================\n\n");
    
    // Model information
    fprintf(file, "## VISCOELASTIC MODEL\n");
    fprintf(file, "Long-term modulus (G_inf): %.3e Pa\n", model->G_inf);
    fprintf(file, "Poisson's ratio: %.3f\n", model->nu);
    fprintf(file, "Number of Prony terms: %d\n\n", model->num_prony_terms);
    
    fprintf(file, "Prony Series:\n");
    for (int i = 0; i < model->num_prony_terms; i++) {
        fprintf(file, "  Term %d: G_k = %.3e Pa, tau_k = %.3e s\n",
                i, model->prony_terms[i].G_k, model->prony_terms[i].tau_k);
    }
    fprintf(file, "\n");
    
    // Contour information
    fprintf(file, "## CONTOUR INFORMATION\n");
    fprintf(file, "Number of nodes: %d\n", contour->num_nodes);
    fprintf(file, "Total length: %.6e m\n\n", contour->total_length);
    
    // Simulation results
    fprintf(file, "## SIMULATION RESULTS\n");
    fprintf(file, "Number of time steps: %d\n", num_steps);
    if (num_steps > 0) {
        fprintf(file, "Time range: %.3e to %.3e s\n", 
                results[0].time, results[num_steps-1].time);
        fprintf(file, "Initial J-integral: %.6e N/m\n", results[0].j_value);
        fprintf(file, "Final J-integral: %.6e N/m\n", results[num_steps-1].j_value);
        
        double change = results[num_steps-1].j_value - results[0].j_value;
        double percent_change = 100.0 * change / (fabs(results[0].j_value) + 1e-12);
        fprintf(file, "Total change: %.6e N/m (%.2f%%)\n", change, percent_change);
    }
    
    fclose(file);
    printf("Simulation report saved to %s\n", filename);
}

// Print program banner
void print_banner() {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════════╗\n");
    printf("║                   VISCOELASTIC J-INTEGRAL CALCULATOR                 ║\n");
    printf("║                                                                      ║\n");
    printf("║  A C implementation for time-dependent fracture mechanics analysis  ║\n");
    printf("║  Author: Advanced Fracture Mechanics Research Group                 ║\n");
    printf("║  Features:                                                           ║\n");
    printf("║    • Generalized Maxwell viscoelastic model                         ║\n");
    printf("║    • Recursive stress updating algorithm                            ║\n");
    printf("║    • Path-independent contour integral evaluation                   ║\n");
    printf("║    • Time-dependent field evolution simulation                      ║\n");
    printf("╚══════════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
} 