#include "viscoelastic_j_integral.h"

// Compute strain energy density at a node
double compute_strain_energy_density(StressTensor2D stress, StrainTensor2D strain) {
    // W = (1/2) * σᵢⱼ * εᵢⱼ
    double W = 0.5 * (stress.s11 * strain.e11 + 
                      stress.s22 * strain.e22 + 
                      2.0 * stress.s12 * strain.e12);
    return W;
}

// Compute total strain energy along contour
double compute_total_strain_energy(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for strain energy computation\n");
        return 0.0;
    }
    
    double total_energy = 0.0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        double W = compute_strain_energy_density(contour->nodes[i].stress, 
                                               contour->nodes[i].strain);
        total_energy += W * contour->nodes[i].ds;
    }
    
    return total_energy;
}

// Compute strain energy contribution to J-integral
double compute_strain_energy_contribution(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for strain energy contribution\n");
        return 0.0;
    }
    
    double J_strain_energy = 0.0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        ContourNode *node = &contour->nodes[i];
        
        // Strain energy density
        double W = compute_strain_energy_density(node->stress, node->strain);
        
        // J-integral strain energy term: W * δ₁ⱼ * nⱼ * ds
        // For 2D: W * (δ₁₁ * n₁ + δ₁₂ * n₂) = W * n₁
        double contribution = W * node->normal.x * node->ds;
        J_strain_energy += contribution;
    }
    
    return J_strain_energy;
}

// Compute stress-traction contribution to J-integral
double compute_stress_traction_contribution(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for stress-traction contribution\n");
        return 0.0;
    }
    
    double J_stress_traction = 0.0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        ContourNode *node = &contour->nodes[i];
        
        // Stress-displacement gradient terms: σᵢⱼ * (∂uᵢ/∂x₁) * nⱼ * ds
        double term1 = node->stress.s11 * node->du_dx.du1_dx1 * node->normal.x +
                       node->stress.s12 * node->du_dx.du1_dx1 * node->normal.y;
        
        double term2 = node->stress.s21 * node->du_dx.du2_dx1 * node->normal.x +
                       node->stress.s22 * node->du_dx.du2_dx1 * node->normal.y;
        
        double contribution = (term1 + term2) * node->ds;
        J_stress_traction -= contribution; // Note the negative sign in J-integral definition
    }
    
    return J_stress_traction;
}

// Compute complete J-integral
double compute_j_integral(Contour *contour) {
    if (contour == NULL) {
        fprintf(stderr, "Error: Invalid contour for J-integral computation\n");
        return 0.0;
    }
    
    double J_strain = compute_strain_energy_contribution(contour);
    double J_stress = compute_stress_traction_contribution(contour);
    
    double J_total = J_strain + J_stress;
    
    return J_total;
}

// Evaluate J-integral at specific time with detailed breakdown
JIntegralResult evaluate_j_integral_at_time(Contour *contour, double time) {
    JIntegralResult result;
    
    if (contour == NULL) {
        fprintf(stderr, "Error: Invalid contour for J-integral evaluation\n");
        memset(&result, 0, sizeof(JIntegralResult));
        return result;
    }
    
    result.time = time;
    result.strain_energy_integral = compute_strain_energy_contribution(contour);
    result.stress_traction_integral = compute_stress_traction_contribution(contour);
    result.j_value = result.strain_energy_integral + result.stress_traction_integral;
    
    return result;
}

// Compare J-integral values between different contours (path independence test)
double compare_j_integral_paths(Contour *contour1, Contour *contour2) {
    if (contour1 == NULL || contour2 == NULL) {
        fprintf(stderr, "Error: Invalid contours for comparison\n");
        return -1.0;
    }
    
    double J1 = compute_j_integral(contour1);
    double J2 = compute_j_integral(contour2);
    
    double relative_difference = fabs(J1 - J2) / (0.5 * (fabs(J1) + fabs(J2)) + 1e-12);
    
    printf("Path independence check: J1 = %.3e, J2 = %.3e, rel. diff = %.2e%%\n",
           J1, J2, relative_difference * 100.0);
    
    return relative_difference;
}

// Validate path independence for multiple contours
void validate_path_independence(Contour **contours, int num_contours, double tolerance) {
    if (contours == NULL || num_contours < 2) {
        fprintf(stderr, "Error: Need at least 2 contours for path independence validation\n");
        return;
    }
    
    printf("\n=== Path Independence Validation ===\n");
    
    double *J_values = (double*)malloc(num_contours * sizeof(double));
    if (J_values == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for J values\n");
        return;
    }
    
    // Compute J-integral for each contour
    for (int i = 0; i < num_contours; i++) {
        J_values[i] = compute_j_integral(contours[i]);
        printf("Contour %d: J = %.6e N/m\n", i, J_values[i]);
    }
    
    // Check pairwise differences
    double max_relative_diff = 0.0;
    for (int i = 0; i < num_contours; i++) {
        for (int j = i + 1; j < num_contours; j++) {
            double rel_diff = fabs(J_values[i] - J_values[j]) / 
                             (0.5 * (fabs(J_values[i]) + fabs(J_values[j])) + 1e-12);
            
            if (rel_diff > max_relative_diff) {
                max_relative_diff = rel_diff;
            }
            
            printf("Contours %d-%d: relative difference = %.2e%%\n", 
                   i, j, rel_diff * 100.0);
        }
    }
    
    printf("\nMaximum relative difference: %.2e%%\n", max_relative_diff * 100.0);
    
    if (max_relative_diff < tolerance) {
        printf("✓ Path independence VALIDATED (tolerance: %.2e%%)\n", tolerance * 100.0);
    } else {
        printf("✗ Path independence VIOLATED (exceeds tolerance: %.2e%%)\n", tolerance * 100.0);
    }
    
    free(J_values);
    printf("=====================================\n\n");
}

// Compute J-integral using midpoint rule for better accuracy
double compute_j_integral_midpoint(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL || contour->num_nodes < 2) {
        fprintf(stderr, "Error: Invalid contour for midpoint J-integral computation\n");
        return 0.0;
    }
    
    double J_total = 0.0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        int next = (i + 1) % contour->num_nodes;
        
        // Midpoint coordinates
        double x_mid = 0.5 * (contour->nodes[i].coord.x + contour->nodes[next].coord.x);
        double y_mid = 0.5 * (contour->nodes[i].coord.y + contour->nodes[next].coord.y);
        
        // Interpolated values at midpoint
        StressTensor2D stress_mid;
        StrainTensor2D strain_mid;
        DisplacementGradient2D du_dx_mid;
        Point2D normal_mid;
        
        stress_mid.s11 = 0.5 * (contour->nodes[i].stress.s11 + contour->nodes[next].stress.s11);
        stress_mid.s22 = 0.5 * (contour->nodes[i].stress.s22 + contour->nodes[next].stress.s22);
        stress_mid.s12 = 0.5 * (contour->nodes[i].stress.s12 + contour->nodes[next].stress.s12);
        stress_mid.s21 = 0.5 * (contour->nodes[i].stress.s21 + contour->nodes[next].stress.s21);
        
        strain_mid.e11 = 0.5 * (contour->nodes[i].strain.e11 + contour->nodes[next].strain.e11);
        strain_mid.e22 = 0.5 * (contour->nodes[i].strain.e22 + contour->nodes[next].strain.e22);
        strain_mid.e12 = 0.5 * (contour->nodes[i].strain.e12 + contour->nodes[next].strain.e12);
        strain_mid.e21 = 0.5 * (contour->nodes[i].strain.e21 + contour->nodes[next].strain.e21);
        
        du_dx_mid.du1_dx1 = 0.5 * (contour->nodes[i].du_dx.du1_dx1 + contour->nodes[next].du_dx.du1_dx1);
        du_dx_mid.du1_dx2 = 0.5 * (contour->nodes[i].du_dx.du1_dx2 + contour->nodes[next].du_dx.du1_dx2);
        du_dx_mid.du2_dx1 = 0.5 * (contour->nodes[i].du_dx.du2_dx1 + contour->nodes[next].du_dx.du2_dx1);
        du_dx_mid.du2_dx2 = 0.5 * (contour->nodes[i].du_dx.du2_dx2 + contour->nodes[next].du_dx.du2_dx2);
        
        normal_mid.x = 0.5 * (contour->nodes[i].normal.x + contour->nodes[next].normal.x);
        normal_mid.y = 0.5 * (contour->nodes[i].normal.y + contour->nodes[next].normal.y);
        
        // Normalize the midpoint normal
        double normal_mag = sqrt(normal_mid.x * normal_mid.x + normal_mid.y * normal_mid.y);
        if (normal_mag > 1e-12) {
            normal_mid.x /= normal_mag;
            normal_mid.y /= normal_mag;
        }
        
        // Segment length
        double dx = contour->nodes[next].coord.x - contour->nodes[i].coord.x;
        double dy = contour->nodes[next].coord.y - contour->nodes[i].coord.y;
        double ds = sqrt(dx*dx + dy*dy);
        
        // Strain energy contribution
        double W = compute_strain_energy_density(stress_mid, strain_mid);
        double J_strain = W * normal_mid.x * ds;
        
        // Stress-traction contribution
        double term1 = stress_mid.s11 * du_dx_mid.du1_dx1 * normal_mid.x +
                       stress_mid.s12 * du_dx_mid.du1_dx1 * normal_mid.y;
        
        double term2 = stress_mid.s21 * du_dx_mid.du2_dx1 * normal_mid.x +
                       stress_mid.s22 * du_dx_mid.du2_dx1 * normal_mid.y;
        
        double J_stress = -(term1 + term2) * ds;
        
        J_total += (J_strain + J_stress);
    }
    
    return J_total;
}

// Compute J-integral rate (dJ/dt) for time-dependent analysis
double compute_j_integral_rate(JIntegralResult *results, int num_results, int current_index) {
    if (results == NULL || num_results < 2 || current_index < 1 || current_index >= num_results) {
        return 0.0;
    }
    
    double dt = results[current_index].time - results[current_index-1].time;
    double dJ = results[current_index].j_value - results[current_index-1].j_value;
    
    if (fabs(dt) < 1e-12) {
        return 0.0;
    }
    
    return dJ / dt;
}

// Analyze J-integral evolution statistics
void analyze_j_evolution(JIntegralResult *results, int num_results) {
    if (results == NULL || num_results < 2) {
        fprintf(stderr, "Error: Insufficient data for J-integral evolution analysis\n");
        return;
    }
    
    printf("\n=== J-Integral Evolution Analysis ===\n");
    
    double J_initial = results[0].j_value;
    double J_final = results[num_results-1].j_value;
    double time_total = results[num_results-1].time - results[0].time;
    
    printf("Initial J-integral: %.6e N/m\n", J_initial);
    printf("Final J-integral: %.6e N/m\n", J_final);
    printf("Total change: %.6e N/m (%.2f%%)\n", 
           J_final - J_initial, 
           100.0 * (J_final - J_initial) / (fabs(J_initial) + 1e-12));
    printf("Time span: %.3e s\n", time_total);
    
    // Find maximum and minimum values
    double J_max = results[0].j_value;
    double J_min = results[0].j_value;
    int idx_max = 0, idx_min = 0;
    
    for (int i = 1; i < num_results; i++) {
        if (results[i].j_value > J_max) {
            J_max = results[i].j_value;
            idx_max = i;
        }
        if (results[i].j_value < J_min) {
            J_min = results[i].j_value;
            idx_min = i;
        }
    }
    
    printf("Maximum J: %.6e N/m at t = %.3e s\n", J_max, results[idx_max].time);
    printf("Minimum J: %.6e N/m at t = %.3e s\n", J_min, results[idx_min].time);
    
    // Compute average rate of change
    double avg_rate = (J_final - J_initial) / time_total;
    printf("Average rate: %.6e N/m/s\n", avg_rate);
    
    printf("=====================================\n\n");
} 