#include "viscoelastic_j_integral.h"

// Initialize Prony series parameters
void initialize_prony_series(ViscoelasticModel *model, double *G_values, double *tau_values) {
    if (model == NULL || G_values == NULL || tau_values == NULL) {
        fprintf(stderr, "Error: Invalid parameters for Prony series initialization\n");
        return;
    }
    
    for (int i = 0; i < model->num_prony_terms; i++) {
        model->prony_terms[i].G_k = G_values[i];
        model->prony_terms[i].tau_k = tau_values[i];
        
        printf("Prony term %d: G_k = %.2e Pa, tau_k = %.2e s\n", 
               i, G_values[i], tau_values[i]);
    }
    
    printf("Initialized Prony series with %d terms\n", model->num_prony_terms);
}

// Set up a typical viscoelastic model for asphalt concrete
void setup_asphalt_model(ViscoelasticModel *model) {
    if (model == NULL) {
        fprintf(stderr, "Error: Invalid model for asphalt setup\n");
        return;
    }
    
    // Typical asphalt concrete parameters
    model->G_inf = 1.0e8;  // Long-term shear modulus (Pa)
    model->nu = 0.35;      // Poisson's ratio
    
    // Prony series parameters for asphalt (example values)
    double G_values[] = {5.0e8, 3.0e8, 1.5e8, 0.8e8, 0.3e8};
    double tau_values[] = {1e-3, 1e-1, 1e1, 1e3, 1e5};  // Time constants (s)
    
    int num_terms = sizeof(G_values) / sizeof(G_values[0]);
    if (num_terms > model->num_prony_terms) {
        num_terms = model->num_prony_terms;
    }
    
    initialize_prony_series(model, G_values, tau_values);
    
    printf("Set up asphalt concrete model: G_inf = %.2e Pa, nu = %.3f\n",
           model->G_inf, model->nu);
}

// Update viscoelastic stress using recursive algorithm
void update_viscoelastic_stress(ViscoelasticModel *model, Contour *contour, 
                               StrainTensor2D *strain_increment, double dt) {
    if (model == NULL || contour == NULL || strain_increment == NULL) {
        fprintf(stderr, "Error: Invalid parameters for stress update\n");
        return;
    }
    
    for (int node = 0; node < contour->num_nodes; node++) {
        // Reset total stress for this node
        StressTensor2D total_stress = {0.0, 0.0, 0.0, 0.0};
        
        // Add long-term elastic contribution
        double G_inf = model->G_inf;
        total_stress.s11 += G_inf * strain_increment[node].e11;
        total_stress.s22 += G_inf * strain_increment[node].e22;
        total_stress.s12 += G_inf * strain_increment[node].e12;
        total_stress.s21 += G_inf * strain_increment[node].e21;
        
        // Add contributions from each Prony term
        for (int k = 0; k < model->num_prony_terms; k++) {
            recursive_stress_update(&model->prony_terms[k], 
                                   &strain_increment[node],
                                   &total_stress, dt, node);
        }
        
        // Add to existing stress (for incremental loading)
        contour->nodes[node].stress.s11 += total_stress.s11;
        contour->nodes[node].stress.s22 += total_stress.s22;
        contour->nodes[node].stress.s12 += total_stress.s12;
        contour->nodes[node].stress.s21 += total_stress.s21;
    }
}

// Recursive stress update for individual Prony term
void recursive_stress_update(PronyTerm *term, StrainTensor2D *strain_increment, 
                           StressTensor2D *total_stress, double dt, int node_index) {
    if (term == NULL || strain_increment == NULL || total_stress == NULL) {
        fprintf(stderr, "Error: Invalid parameters for recursive stress update\n");
        return;
    }
    
    if (term->stress_history == NULL) {
        fprintf(stderr, "Error: Stress history not allocated for Prony term\n");
        return;
    }
    
    double G_k = term->G_k;
    double tau_k = term->tau_k;
    
    // Exponential decay factor
    double decay_factor = exp(-dt / tau_k);
    
    // Memory coefficient for new strain increment
    double memory_coeff = G_k * (1.0 - decay_factor);
    
    // Update stress history for this term
    StressTensor2D *history = &term->stress_history[node_index];
    
    // Apply decay to existing stress history
    history->s11 *= decay_factor;
    history->s22 *= decay_factor;
    history->s12 *= decay_factor;
    history->s21 *= decay_factor;
    
    // Add contribution from new strain increment
    history->s11 += memory_coeff * strain_increment->e11;
    history->s22 += memory_coeff * strain_increment->e22;
    history->s12 += memory_coeff * strain_increment->e12;
    history->s21 += memory_coeff * strain_increment->e21;
    
    // Add this term's contribution to total stress
    total_stress->s11 += history->s11;
    total_stress->s22 += history->s22;
    total_stress->s12 += history->s12;
    total_stress->s21 += history->s21;
}

// Compute strain increment from displacement increment
void compute_strain_increment(Contour *contour, Contour *previous_contour, 
                             StrainTensor2D *strain_increment) {
    if (contour == NULL || previous_contour == NULL || strain_increment == NULL) {
        fprintf(stderr, "Error: Invalid parameters for strain increment computation\n");
        return;
    }
    
    for (int i = 0; i < contour->num_nodes; i++) {
        // Displacement increment
        double du1 = contour->nodes[i].displacement.u1 - previous_contour->nodes[i].displacement.u1;
        double du2 = contour->nodes[i].displacement.u2 - previous_contour->nodes[i].displacement.u2;
        
        // For simplicity, assume small strains and use current displacement gradients
        // In a full implementation, this would be computed more rigorously
        DisplacementGradient2D *du_dx = &contour->nodes[i].du_dx;
        
        // Strain increment (small strain assumption)
        strain_increment[i].e11 = du_dx->du1_dx1;
        strain_increment[i].e22 = du_dx->du2_dx2;
        strain_increment[i].e12 = 0.5 * (du_dx->du1_dx2 + du_dx->du2_dx1);
        strain_increment[i].e21 = strain_increment[i].e12; // Symmetry
        
        // Scale by some factor if needed (this is a simplified approach)
        double scale_factor = 0.01; // Adjust based on loading rate
        strain_increment[i].e11 *= scale_factor;
        strain_increment[i].e22 *= scale_factor;
        strain_increment[i].e12 *= scale_factor;
        strain_increment[i].e21 *= scale_factor;
    }
}

// Apply relaxation modulus function
double relaxation_modulus(ViscoelasticModel *model, double time) {
    if (model == NULL || time < 0.0) {
        return 0.0;
    }
    
    double G_t = model->G_inf;
    
    for (int k = 0; k < model->num_prony_terms; k++) {
        double G_k = model->prony_terms[k].G_k;
        double tau_k = model->prony_terms[k].tau_k;
        
        G_t += G_k * exp(-time / tau_k);
    }
    
    return G_t;
}

// Apply creep compliance function
double creep_compliance(ViscoelasticModel *model, double time) {
    if (model == NULL || time < 0.0) {
        return 0.0;
    }
    
    // For a generalized Maxwell model, the creep compliance is more complex
    // This is a simplified approximation
    double G_total = model->G_inf;
    for (int k = 0; k < model->num_prony_terms; k++) {
        G_total += model->prony_terms[k].G_k;
    }
    
    // Initial compliance
    double J_0 = 1.0 / G_total;
    
    // Add time-dependent terms (simplified)
    double J_t = J_0;
    for (int k = 0; k < model->num_prony_terms; k++) {
        double G_k = model->prony_terms[k].G_k;
        double tau_k = model->prony_terms[k].tau_k;
        
        J_t += (G_k / (G_total * G_total)) * (1.0 - exp(-time / tau_k));
    }
    
    return J_t;
}

// Reset stress history (for new loading scenarios)
void reset_stress_history(ViscoelasticModel *model, int num_nodes) {
    if (model == NULL || num_nodes <= 0) {
        return;
    }
    
    for (int k = 0; k < model->num_prony_terms; k++) {
        if (model->prony_terms[k].stress_history != NULL) {
            memset(model->prony_terms[k].stress_history, 0, 
                   num_nodes * sizeof(StressTensor2D));
        }
    }
    
    printf("Reset stress history for %d nodes\n", num_nodes);
}

// Validate viscoelastic model parameters
int validate_model_parameters(ViscoelasticModel *model) {
    if (model == NULL) {
        fprintf(stderr, "Error: Null model pointer\n");
        return -1;
    }
    
    if (model->G_inf < 0.0) {
        fprintf(stderr, "Error: Negative long-term modulus\n");
        return -1;
    }
    
    if (model->nu < 0.0 || model->nu >= 0.5) {
        fprintf(stderr, "Error: Invalid Poisson's ratio (%.3f)\n", model->nu);
        return -1;
    }
    
    for (int k = 0; k < model->num_prony_terms; k++) {
        if (model->prony_terms[k].G_k < 0.0) {
            fprintf(stderr, "Error: Negative modulus in Prony term %d\n", k);
            return -1;
        }
        if (model->prony_terms[k].tau_k <= 0.0) {
            fprintf(stderr, "Error: Non-positive relaxation time in Prony term %d\n", k);
            return -1;
        }
    }
    
    printf("Model parameters validation: PASSED\n");
    return 0;
} 