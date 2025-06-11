#include "viscoelastic_j_integral.h"

// Function prototypes for demonstration scenarios
void demo_basic_j_integral();
void demo_path_independence();
void demo_viscoelastic_evolution();
void demo_relaxation_test();
void demo_load_history();

int main(int argc, char *argv[]) {
    print_banner();
    
    printf("🔬 DEMONSTRATION SCENARIOS:\n");
    printf("1. Basic J-integral calculation with analytical fields\n");
    printf("2. Path independence validation\n");
    printf("3. Viscoelastic evolution simulation\n");
    printf("4. Stress relaxation test\n");
    printf("5. Load history application\n");
    printf("\n");
    
    // Check for command line arguments
    int demo_choice = 0;
    if (argc > 1) {
        demo_choice = atoi(argv[1]);
    }
    
    if (demo_choice == 0) {
        printf("Enter demo number (1-5) or 0 for all: ");
        scanf("%d", &demo_choice);
    }
    
    printf("\n============================================================\n");
    
    switch (demo_choice) {
        case 1:
            demo_basic_j_integral();
            break;
        case 2:
            demo_path_independence();
            break;
        case 3:
            demo_viscoelastic_evolution();
            break;
        case 4:
            demo_relaxation_test();
            break;
        case 5:
            demo_load_history();
            break;
        case 0:
        default:
            printf("Running all demonstrations...\n\n");
            demo_basic_j_integral();
            demo_path_independence();
            demo_viscoelastic_evolution();
            demo_relaxation_test();
            demo_load_history();
            break;
    }
    
    printf("\n🎉 All demonstrations completed successfully!\n");
    printf("Check the generated output files for detailed results.\n\n");
    
    return 0;
}

// Demo 1: Basic J-integral calculation
void demo_basic_j_integral() {
    printf("\n🔍 DEMO 1: Basic J-Integral Calculation\n");
    printf("=======================================\n");
    
    // Create contour around crack tip
    Contour contour;
    Point2D crack_tip = {0.0, 0.0};
    Point2D center = {0.0, 0.0};
    double radius = 0.01; // 1 cm radius
    int num_nodes = 32;
    
    initialize_circular_contour(&contour, center, radius, num_nodes);
    
    // Initialize fields with analytical crack tip solution
    initialize_displacement_field(&contour, crack_tip);
    initialize_stress_field(&contour, crack_tip);
    initialize_strain_field(&contour, crack_tip);
    compute_displacement_gradients(&contour);
    
    // Print contour information
    print_contour_info(&contour);
    
    // Compute J-integral
    double J_value = compute_j_integral(&contour);
    double J_strain = compute_strain_energy_contribution(&contour);
    double J_stress = compute_stress_traction_contribution(&contour);
    
    printf("J-Integral Results:\n");
    printf("  Strain energy contribution: %.6e N/m\n", J_strain);
    printf("  Stress-traction contribution: %.6e N/m\n", J_stress);
    printf("  Total J-integral: %.6e N/m\n", J_value);
    
    // Compare with theoretical value (K_I = 1000 Pa√m, E = 2e9 Pa, nu = 0.3)
    double K_I = 1000.0;
    double E = 2.0e9;
    double nu = 0.3;
    double J_theoretical = K_I * K_I * (1.0 - nu*nu) / E;
    
    printf("  Theoretical J-integral: %.6e N/m\n", J_theoretical);
    printf("  Relative error: %.2f%%\n", 
           100.0 * fabs(J_value - J_theoretical) / J_theoretical);
    
    // Save contour data
    save_contour_to_file(&contour, "demo1_contour.dat");
    
    // Clean up
    free_contour_nodes(contour.nodes);
    
    printf("✓ Demo 1 completed\n");
}

// Demo 2: Path independence validation
void demo_path_independence() {
    printf("\n🔍 DEMO 2: Path Independence Validation\n");
    printf("======================================\n");
    
    Point2D crack_tip = {0.0, 0.0};
    Point2D center = {0.0, 0.0};
    int num_contours = 5;
    double min_radius = 0.005; // 5 mm
    double max_radius = 0.025; // 25 mm
    
    // Create multiple concentric contours
    Contour **contours = (Contour**)malloc(num_contours * sizeof(Contour*));
    if (contours == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return;
    }
    
    create_concentric_contours(contours, num_contours, center, min_radius, max_radius);
    
    // Initialize fields for all contours
    for (int i = 0; i < num_contours; i++) {
        initialize_displacement_field(contours[i], crack_tip);
        initialize_stress_field(contours[i], crack_tip);
        initialize_strain_field(contours[i], crack_tip);
        compute_displacement_gradients(contours[i]);
    }
    
    // Validate path independence
    double tolerance = 0.05; // 5% tolerance
    validate_path_independence(contours, num_contours, tolerance);
    
    // Clean up
    for (int i = 0; i < num_contours; i++) {
        free_contour_nodes(contours[i]->nodes);
        free(contours[i]);
    }
    free(contours);
    
    printf("✓ Demo 2 completed\n");
}

// Demo 3: Viscoelastic evolution simulation
void demo_viscoelastic_evolution() {
    printf("\n🔍 DEMO 3: Viscoelastic Evolution Simulation\n");
    printf("===========================================\n");
    
    // Create viscoelastic model (asphalt concrete)
    ViscoelasticModel *model = create_viscoelastic_model(5);
    if (model == NULL) {
        fprintf(stderr, "Error: Failed to create viscoelastic model\n");
        return;
    }
    
    setup_asphalt_model(model);
    print_viscoelastic_model(model);
    
    // Validate model parameters
    if (validate_model_parameters(model) != 0) {
        free_viscoelastic_model(model);
        return;
    }
    
    // Create contour
    Contour contour;
    Point2D crack_tip = {0.0, 0.0};
    Point2D center = {0.0, 0.0};
    initialize_circular_contour(&contour, center, 0.01, 32);
    
    // Initialize fields
    initialize_displacement_field(&contour, crack_tip);
    initialize_stress_field(&contour, crack_tip);
    initialize_strain_field(&contour, crack_tip);
    compute_displacement_gradients(&contour);
    
    // Set up time parameters
    double dt = 0.1;          // 0.1 second time steps
    double total_time = 100.0; // 100 second simulation
    TimeParameters time_params = create_time_parameters(dt, total_time);
    
    // Run simulation
    JIntegralResult *results = NULL;
    simulate_viscoelastic_evolution(model, &contour, &time_params, &results);
    
    if (results != NULL) {
        int num_steps = (int)(total_time / dt) + 1;
        
        // Save results
        save_results_to_file(results, num_steps, "demo3_j_evolution.dat");
        create_simulation_report(model, &contour, results, num_steps, "demo3_report.txt");
        
        // Print key results
        printf("\nKey Results:\n");
        printf("Initial J-integral: %.6e N/m\n", results[0].j_value);
        printf("Final J-integral: %.6e N/m\n", results[num_steps-1].j_value);
        
        free_results_array(results);
    }
    
    // Clean up
    free_contour_nodes(contour.nodes);
    free_viscoelastic_model(model);
    
    printf("✓ Demo 3 completed\n");
}

// Demo 4: Stress relaxation test
void demo_relaxation_test() {
    printf("\n🔍 DEMO 4: Stress Relaxation Test\n");
    printf("================================\n");
    
    // Create viscoelastic model
    ViscoelasticModel *model = create_viscoelastic_model(3);
    if (model == NULL) {
        fprintf(stderr, "Error: Failed to create viscoelastic model\n");
        return;
    }
    
    // Set up simplified model for relaxation test
    model->G_inf = 1.0e8;  // 100 MPa
    model->nu = 0.3;
    
    double G_values[] = {2.0e8, 1.0e8, 0.5e8};
    double tau_values[] = {1.0, 10.0, 100.0};
    initialize_prony_series(model, G_values, tau_values);
    
    print_viscoelastic_model(model);
    
    // Create contour
    Contour contour;
    Point2D center = {0.0, 0.0};
    initialize_circular_contour(&contour, center, 0.01, 24);
    
    // Set up relaxation test parameters
    double strain_level = 0.001; // 0.1% strain
    double dt = 1.0;              // 1 second time steps
    double total_time = 1000.0;   // 1000 second relaxation
    TimeParameters time_params = create_time_parameters(dt, total_time);
    
    // Run relaxation test
    JIntegralResult *results = NULL;
    relaxation_test(model, &contour, strain_level, &time_params, &results);
    
    if (results != NULL) {
        int num_steps = (int)(total_time / dt) + 1;
        
        // Save results
        save_results_to_file(results, num_steps, "demo4_relaxation.dat");
        
        // Print relaxation statistics
        printf("\nRelaxation Statistics:\n");
        printf("Applied strain: %.4f\n", strain_level);
        printf("Initial J-integral: %.6e N/m\n", results[0].j_value);
        printf("Final J-integral: %.6e N/m\n", results[num_steps-1].j_value);
        
        double relaxation_ratio = results[num_steps-1].j_value / results[0].j_value;
        printf("Relaxation ratio: %.3f\n", relaxation_ratio);
        
        free_results_array(results);
    }
    
    // Clean up
    free_contour_nodes(contour.nodes);
    free_viscoelastic_model(model);
    
    printf("✓ Demo 4 completed\n");
}

// Demo 5: Load history application
void demo_load_history() {
    printf("\n🔍 DEMO 5: Load History Application\n");
    printf("==================================\n");
    
    // Create viscoelastic model
    ViscoelasticModel *model = create_viscoelastic_model(4);
    if (model == NULL) {
        fprintf(stderr, "Error: Failed to create viscoelastic model\n");
        return;
    }
    
    setup_asphalt_model(model);
    
    // Create contour
    Contour contour;
    Point2D center = {0.0, 0.0};
    initialize_circular_contour(&contour, center, 0.01, 32);
    
    // Initialize fields
    Point2D crack_tip = {0.0, 0.0};
    initialize_displacement_field(&contour, crack_tip);
    initialize_stress_field(&contour, crack_tip);
    initialize_strain_field(&contour, crack_tip);
    compute_displacement_gradients(&contour);
    
    // Define load history (step loading)
    int num_loads = 4;
    double load_times[] = {10.0, 20.0, 30.0, 40.0};    // Duration of each load step
    double load_values[] = {1.0, 2.0, 1.5, 0.5};       // Load magnitude for each step
    
    printf("Load History:\n");
    for (int i = 0; i < num_loads; i++) {
        printf("  Step %d: Load = %.1f, Duration = %.1f s\n", 
               i+1, load_values[i], load_times[i]);
    }
    printf("\n");
    
    // Set up time parameters
    double dt = 0.5; // 0.5 second time steps
    TimeParameters time_params = create_time_parameters(dt, 0.0); // total_time calculated in function
    
    // Apply load history
    JIntegralResult *results = NULL;
    apply_load_history(model, &contour, load_times, load_values, num_loads, 
                      &time_params, &results);
    
    if (results != NULL) {
        // Calculate total steps
        int total_steps = 0;
        for (int i = 0; i < num_loads; i++) {
            total_steps += (int)(load_times[i] / dt);
        }
        
        // Save results
        save_results_to_file(results, total_steps, "demo5_load_history.dat");
        
        free_results_array(results);
    }
    
    // Clean up
    free_contour_nodes(contour.nodes);
    free_viscoelastic_model(model);
    
    printf("✓ Demo 5 completed\n");
} 