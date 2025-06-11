#include "viscoelastic_j_integral.h"

// Perform single time step update
void time_step_update(ViscoelasticModel *model, Contour *contour, 
                     TimeParameters *time_params) {
    if (model == NULL || contour == NULL || time_params == NULL) {
        fprintf(stderr, "Error: Invalid parameters for time step update\n");
        return;
    }
    
    // Allocate strain increment array
    StrainTensor2D *strain_increment = (StrainTensor2D*)calloc(contour->num_nodes, 
                                                               sizeof(StrainTensor2D));
    if (strain_increment == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for strain increment\n");
        return;
    }
    
    // For this simplified implementation, we'll apply a small strain increment
    // In a full implementation, this would come from solving the equilibrium equations
    double load_factor = 1.0 + 0.01 * sin(2.0 * M_PI * time_params->current_time / 10.0);
    
    for (int i = 0; i < contour->num_nodes; i++) {
        // Apply small incremental strain based on current loading
        strain_increment[i].e11 = 1e-6 * load_factor * time_params->dt;
        strain_increment[i].e22 = -0.3 * strain_increment[i].e11; // Poisson effect
        strain_increment[i].e12 = 0.5e-6 * load_factor * time_params->dt;
        strain_increment[i].e21 = strain_increment[i].e12;
    }
    
    // Update viscoelastic stress
    update_viscoelastic_stress(model, contour, strain_increment, time_params->dt);
    
    // Update total strain (accumulate increments)
    for (int i = 0; i < contour->num_nodes; i++) {
        contour->nodes[i].strain.e11 += strain_increment[i].e11;
        contour->nodes[i].strain.e22 += strain_increment[i].e22;
        contour->nodes[i].strain.e12 += strain_increment[i].e12;
        contour->nodes[i].strain.e21 += strain_increment[i].e21;
    }
    
    // Update displacement field (simplified - based on strain)
    // In a full implementation, this would be solved from equilibrium
    for (int i = 0; i < contour->num_nodes; i++) {
        double x = contour->nodes[i].coord.x;
        double y = contour->nodes[i].coord.y;
        
        // Simple displacement update based on accumulated strain
        contour->nodes[i].displacement.u1 += strain_increment[i].e11 * x * 0.1;
        contour->nodes[i].displacement.u2 += strain_increment[i].e22 * y * 0.1;
    }
    
    // Recompute displacement gradients
    compute_displacement_gradients(contour);
    
    // Update time parameters
    time_params->current_time += time_params->dt;
    time_params->time_step++;
    
    free(strain_increment);
}

// Simulate complete viscoelastic evolution
void simulate_viscoelastic_evolution(ViscoelasticModel *model, Contour *contour,
                                   TimeParameters *time_params, JIntegralResult **results) {
    if (model == NULL || contour == NULL || time_params == NULL || results == NULL) {
        fprintf(stderr, "Error: Invalid parameters for viscoelastic simulation\n");
        return;
    }
    
    // Calculate number of time steps
    int num_steps = (int)(time_params->total_time / time_params->dt) + 1;
    
    // Allocate results array
    *results = allocate_results_array(num_steps);
    if (*results == NULL) {
        return;
    }
    
    // Allocate stress history for all Prony terms
    if (allocate_stress_history(model, contour->num_nodes) != 0) {
        free_results_array(*results);
        *results = NULL;
        return;
    }
    
    printf("\n=== Starting Viscoelastic Evolution Simulation ===\n");
    printf("Total time: %.3e s\n", time_params->total_time);
    printf("Time step: %.3e s\n", time_params->dt);
    printf("Number of steps: %d\n", num_steps);
    printf("Contour nodes: %d\n", contour->num_nodes);
    printf("Prony terms: %d\n", model->num_prony_terms);
    printf("==================================================\n\n");
    
    // Initialize time parameters
    time_params->current_time = 0.0;
    time_params->time_step = 0;
    
    // Compute initial J-integral
    (*results)[0] = evaluate_j_integral_at_time(contour, time_params->current_time);
    
    printf("Step %4d: t = %.3e s, J = %.6e N/m\n", 
           0, time_params->current_time, (*results)[0].j_value);
    
    // Time-stepping loop
    for (int step = 1; step < num_steps; step++) {
        // Update fields for current time step
        time_step_update(model, contour, time_params);
        
        // Evaluate J-integral at current time
        (*results)[step] = evaluate_j_integral_at_time(contour, time_params->current_time);
        
        // Print progress every 10% of simulation
        if (step % (num_steps / 10) == 0 || step == num_steps - 1) {
            double progress = 100.0 * step / (num_steps - 1);
            printf("Step %4d: t = %.3e s, J = %.6e N/m (%.1f%%)\n", 
                   step, time_params->current_time, (*results)[step].j_value, progress);
        }
        
        // Check for convergence issues
        if (!isfinite((*results)[step].j_value)) {
            fprintf(stderr, "Error: Non-finite J-integral at step %d\n", step);
            break;
        }
    }
    
    printf("\n=== Simulation Complete ===\n");
    printf("Final time: %.3e s\n", time_params->current_time);
    printf("Final J-integral: %.6e N/m\n", (*results)[num_steps-1].j_value);
    
    // Analyze evolution
    analyze_j_evolution(*results, num_steps);
}

// Apply load history to simulation
void apply_load_history(ViscoelasticModel *model, Contour *contour, 
                       double *load_times, double *load_values, int num_loads,
                       TimeParameters *time_params, JIntegralResult **results) {
    if (model == NULL || contour == NULL || load_times == NULL || 
        load_values == NULL || num_loads <= 0) {
        fprintf(stderr, "Error: Invalid parameters for load history application\n");
        return;
    }
    
    printf("\n=== Applying Load History ===\n");
    printf("Number of load steps: %d\n", num_loads);
    
    // Calculate total number of time steps
    int total_steps = 0;
    for (int i = 0; i < num_loads; i++) {
        int steps = (int)(load_times[i] / time_params->dt);
        total_steps += steps;
        printf("Load %d: magnitude = %.3e, duration = %.3e s (%d steps)\n",
               i, load_values[i], load_times[i], steps);
    }
    
    // Allocate results array
    *results = allocate_results_array(total_steps);
    if (*results == NULL) {
        return;
    }
    
    // Allocate stress history
    if (allocate_stress_history(model, contour->num_nodes) != 0) {
        free_results_array(*results);
        *results = NULL;
        return;
    }
    
    printf("Total simulation steps: %d\n", total_steps);
    printf("=============================\n\n");
    
    int result_index = 0;
    time_params->current_time = 0.0;
    time_params->time_step = 0;
    
    // Apply each load step
    for (int load_step = 0; load_step < num_loads; load_step++) {
        double load_magnitude = load_values[load_step];
        double load_duration = load_times[load_step];
        int steps_in_load = (int)(load_duration / time_params->dt);
        
        printf("Applying load step %d (magnitude: %.3e)\n", load_step, load_magnitude);
        
        for (int step = 0; step < steps_in_load && result_index < total_steps; step++) {
            // Apply load-dependent strain increment
            StrainTensor2D *strain_increment = (StrainTensor2D*)calloc(contour->num_nodes, 
                                                                       sizeof(StrainTensor2D));
            if (strain_increment == NULL) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                continue;
            }
            
            for (int i = 0; i < contour->num_nodes; i++) {
                // Load-proportional strain increment
                strain_increment[i].e11 = load_magnitude * 1e-6 * time_params->dt;
                strain_increment[i].e22 = -0.3 * strain_increment[i].e11;
                strain_increment[i].e12 = 0.5 * strain_increment[i].e11;
                strain_increment[i].e21 = strain_increment[i].e12;
                
                // Update total strain
                contour->nodes[i].strain.e11 += strain_increment[i].e11;
                contour->nodes[i].strain.e22 += strain_increment[i].e22;
                contour->nodes[i].strain.e12 += strain_increment[i].e12;
                contour->nodes[i].strain.e21 += strain_increment[i].e21;
            }
            
            // Update viscoelastic stress
            update_viscoelastic_stress(model, contour, strain_increment, time_params->dt);
            
            // Update time
            time_params->current_time += time_params->dt;
            time_params->time_step++;
            
            // Evaluate J-integral
            (*results)[result_index] = evaluate_j_integral_at_time(contour, time_params->current_time);
            
            if (step % (steps_in_load / 5) == 0) {
                printf("  Step %d/%d: t = %.3e s, J = %.6e N/m\n",
                       step, steps_in_load, time_params->current_time, 
                       (*results)[result_index].j_value);
            }
            
            result_index++;
            free(strain_increment);
        }
    }
    
    printf("\nLoad history application complete\n");
    analyze_j_evolution(*results, result_index);
}

// Perform relaxation test (constant strain)
void relaxation_test(ViscoelasticModel *model, Contour *contour,
                    double strain_level, TimeParameters *time_params,
                    JIntegralResult **results) {
    if (model == NULL || contour == NULL || time_params == NULL) {
        fprintf(stderr, "Error: Invalid parameters for relaxation test\n");
        return;
    }
    
    printf("\n=== Relaxation Test ===\n");
    printf("Applied strain level: %.3e\n", strain_level);
    printf("Test duration: %.3e s\n", time_params->total_time);
    printf("=======================\n\n");
    
    // Apply initial strain instantly
    for (int i = 0; i < contour->num_nodes; i++) {
        contour->nodes[i].strain.e11 = strain_level;
        contour->nodes[i].strain.e22 = -0.3 * strain_level; // Poisson effect
        contour->nodes[i].strain.e12 = 0.0;
        contour->nodes[i].strain.e21 = 0.0;
        
        // Apply initial elastic stress
        double E = 2.0 * model->G_inf * (1.0 + model->nu);
        contour->nodes[i].stress.s11 = E * strain_level / (1.0 - model->nu * model->nu);
        contour->nodes[i].stress.s22 = E * (-0.3 * strain_level) / (1.0 - model->nu * model->nu);
        contour->nodes[i].stress.s12 = 0.0;
        contour->nodes[i].stress.s21 = 0.0;
    }
    
    // Calculate number of time steps
    int num_steps = (int)(time_params->total_time / time_params->dt) + 1;
    
    // Allocate results and stress history
    *results = allocate_results_array(num_steps);
    if (*results == NULL) {
        return;
    }
    
    if (allocate_stress_history(model, contour->num_nodes) != 0) {
        free_results_array(*results);
        *results = NULL;
        return;
    }
    
    // Initialize time parameters
    time_params->current_time = 0.0;
    time_params->time_step = 0;
    
    // Initial J-integral
    (*results)[0] = evaluate_j_integral_at_time(contour, 0.0);
    printf("Initial J-integral: %.6e N/m\n", (*results)[0].j_value);
    
    // Relaxation simulation (no new strain, just stress relaxation)
    for (int step = 1; step < num_steps; step++) {
        // Zero strain increment for relaxation
        StrainTensor2D *zero_strain = (StrainTensor2D*)calloc(contour->num_nodes, 
                                                               sizeof(StrainTensor2D));
        if (zero_strain != NULL) {
            // Update stress with zero strain increment (pure relaxation)
            update_viscoelastic_stress(model, contour, zero_strain, time_params->dt);
            free(zero_strain);
        }
        
        // Update time
        time_params->current_time += time_params->dt;
        time_params->time_step++;
        
        // Evaluate J-integral
        (*results)[step] = evaluate_j_integral_at_time(contour, time_params->current_time);
        
        // Print progress
        if (step % (num_steps / 10) == 0 || step == num_steps - 1) {
            double progress = 100.0 * step / (num_steps - 1);
            printf("Step %4d: t = %.3e s, J = %.6e N/m (%.1f%%)\n",
                   step, time_params->current_time, (*results)[step].j_value, progress);
        }
    }
    
    printf("\nRelaxation test complete\n");
    analyze_j_evolution(*results, num_steps);
}

// Create time parameters structure
TimeParameters create_time_parameters(double dt, double total_time) {
    TimeParameters params;
    params.dt = dt;
    params.total_time = total_time;
    params.current_time = 0.0;
    params.time_step = 0;
    
    return params;
}

// Adaptive time stepping (basic implementation)
double adaptive_time_step(JIntegralResult *previous_results, int num_previous, 
                         double current_dt, double target_accuracy) {
    if (previous_results == NULL || num_previous < 2 || target_accuracy <= 0.0) {
        return current_dt;
    }
    
    // Estimate error based on J-integral rate of change
    double dJ_dt = compute_j_integral_rate(previous_results, num_previous, num_previous - 1);
    double estimated_error = fabs(dJ_dt * current_dt);
    
    // Adjust time step based on error estimate
    double safety_factor = 0.8;
    double new_dt = current_dt;
    
    if (estimated_error > target_accuracy) {
        // Reduce time step
        new_dt = current_dt * safety_factor * sqrt(target_accuracy / estimated_error);
    } else if (estimated_error < 0.1 * target_accuracy) {
        // Increase time step
        new_dt = current_dt * 1.2;
    }
    
    // Limit time step changes
    if (new_dt > 2.0 * current_dt) new_dt = 2.0 * current_dt;
    if (new_dt < 0.1 * current_dt) new_dt = 0.1 * current_dt;
    
    return new_dt;
} 