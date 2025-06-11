#include "viscoelastic_j_integral.h"

// Initialize displacement field using analytical crack tip solution
void initialize_displacement_field(Contour *contour, Point2D crack_tip) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for displacement field initialization\n");
        return;
    }
    
    // Parameters for analytical crack tip field
    double K_I = 1000.0;  // Mode I stress intensity factor (Pa·√m)
    double K_II = 0.0;    // Mode II stress intensity factor (for now, pure Mode I)
    double E = 2.0e9;     // Young's modulus (Pa)
    double nu = 0.3;      // Poisson's ratio
    double mu = E / (2.0 * (1.0 + nu));  // Shear modulus
    double kappa = 3.0 - 4.0 * nu;       // Kolosov constant (plane strain)
    
    for (int i = 0; i < contour->num_nodes; i++) {
        // Relative position from crack tip
        double x = contour->nodes[i].coord.x - crack_tip.x;
        double y = contour->nodes[i].coord.y - crack_tip.y;
        
        double r = sqrt(x*x + y*y);
        double theta = atan2(y, x);
        
        if (r < 1e-12) {
            // At crack tip, displacements are zero
            contour->nodes[i].displacement.u1 = 0.0;
            contour->nodes[i].displacement.u2 = 0.0;
            continue;
        }
        
        double sqrt_r = sqrt(r);
        double theta_half = theta / 2.0;
        double cos_theta_half = cos(theta_half);
        double sin_theta_half = sin(theta_half);
        double cos_3theta_half = cos(3.0 * theta_half);
        double sin_3theta_half = sin(3.0 * theta_half);
        
        // Mode I displacement field
        double u1_I = (K_I / (2.0 * mu)) * sqrt(r / (2.0 * M_PI)) * 
                      cos_theta_half * (kappa - cos(theta));
        double u2_I = (K_I / (2.0 * mu)) * sqrt(r / (2.0 * M_PI)) * 
                      sin_theta_half * (kappa - cos(theta));
        
        // Mode II displacement field
        double u1_II = (K_II / (2.0 * mu)) * sqrt(r / (2.0 * M_PI)) * 
                       sin_theta_half * (kappa + cos(theta) + 2.0);
        double u2_II = -(K_II / (2.0 * mu)) * sqrt(r / (2.0 * M_PI)) * 
                       cos_theta_half * (kappa - cos(theta) - 2.0);
        
        // Total displacement
        contour->nodes[i].displacement.u1 = u1_I + u1_II;
        contour->nodes[i].displacement.u2 = u2_I + u2_II;
    }
    
    printf("Initialized displacement field with K_I = %.1f Pa√m\n", K_I);
}

// Initialize stress field using analytical crack tip solution
void initialize_stress_field(Contour *contour, Point2D crack_tip) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for stress field initialization\n");
        return;
    }
    
    // Parameters for analytical crack tip field
    double K_I = 1000.0;  // Mode I stress intensity factor
    double K_II = 0.0;    // Mode II stress intensity factor
    
    for (int i = 0; i < contour->num_nodes; i++) {
        // Relative position from crack tip
        double x = contour->nodes[i].coord.x - crack_tip.x;
        double y = contour->nodes[i].coord.y - crack_tip.y;
        
        double r = sqrt(x*x + y*y);
        double theta = atan2(y, x);
        
        if (r < 1e-12) {
            // At crack tip, stresses are singular (set to large finite value)
            contour->nodes[i].stress.s11 = 1e10;
            contour->nodes[i].stress.s22 = 1e10;
            contour->nodes[i].stress.s12 = 0.0;
            contour->nodes[i].stress.s21 = 0.0;
            continue;
        }
        
        double sqrt_r = sqrt(r);
        double theta_half = theta / 2.0;
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        double cos_theta_half = cos(theta_half);
        double sin_theta_half = sin(theta_half);
        double cos_3theta_half = cos(3.0 * theta_half);
        double sin_3theta_half = sin(3.0 * theta_half);
        
        double coeff = 1.0 / sqrt(2.0 * M_PI * r);
        
        // Mode I stress field
        double s11_I = K_I * coeff * cos_theta_half * (1.0 - sin_theta_half * sin_3theta_half);
        double s22_I = K_I * coeff * cos_theta_half * (1.0 + sin_theta_half * sin_3theta_half);
        double s12_I = K_I * coeff * sin_theta_half * cos_theta_half * cos_3theta_half;
        
        // Mode II stress field
        double s11_II = -K_II * coeff * sin_theta_half * (2.0 + cos_theta_half * cos_3theta_half);
        double s22_II = K_II * coeff * sin_theta_half * cos_theta_half * cos_3theta_half;
        double s12_II = K_II * coeff * cos_theta_half * (1.0 - sin_theta_half * sin_3theta_half);
        
        // Total stress
        contour->nodes[i].stress.s11 = s11_I + s11_II;
        contour->nodes[i].stress.s22 = s22_I + s22_II;
        contour->nodes[i].stress.s12 = s12_I + s12_II;
        contour->nodes[i].stress.s21 = contour->nodes[i].stress.s12; // Symmetry
    }
    
    printf("Initialized stress field with analytical crack tip solution\n");
}

// Initialize strain field from stress using elastic relations
void initialize_strain_field(Contour *contour, Point2D crack_tip) {
    if (contour == NULL || contour->nodes == NULL) {
        fprintf(stderr, "Error: Invalid contour for strain field initialization\n");
        return;
    }
    
    double E = 2.0e9;     // Young's modulus
    double nu = 0.3;      // Poisson's ratio
    
    // Elastic compliance matrix components
    double C11 = 1.0 / E;
    double C12 = -nu / E;
    double C66 = 2.0 * (1.0 + nu) / E;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        StressTensor2D *stress = &contour->nodes[i].stress;
        StrainTensor2D *strain = &contour->nodes[i].strain;
        
        // Plane strain relations
        strain->e11 = C11 * stress->s11 + C12 * stress->s22;
        strain->e22 = C12 * stress->s11 + C11 * stress->s22;
        strain->e12 = C66 * stress->s12;
        strain->e21 = strain->e12; // Symmetry
    }
    
    printf("Initialized strain field from stress using elastic relations\n");
}

// Compute displacement gradients using finite differences
void compute_displacement_gradients(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL || contour->num_nodes < 3) {
        fprintf(stderr, "Error: Invalid contour for displacement gradient computation\n");
        return;
    }
    
    for (int i = 0; i < contour->num_nodes; i++) {
        int prev = (i - 1 + contour->num_nodes) % contour->num_nodes;
        int next = (i + 1) % contour->num_nodes;
        
        // Coordinate differences
        double dx_prev = contour->nodes[i].coord.x - contour->nodes[prev].coord.x;
        double dy_prev = contour->nodes[i].coord.y - contour->nodes[prev].coord.y;
        double dx_next = contour->nodes[next].coord.x - contour->nodes[i].coord.x;
        double dy_next = contour->nodes[next].coord.y - contour->nodes[i].coord.y;
        
        // Displacement differences
        double du1_prev = contour->nodes[i].displacement.u1 - contour->nodes[prev].displacement.u1;
        double du2_prev = contour->nodes[i].displacement.u2 - contour->nodes[prev].displacement.u2;
        double du1_next = contour->nodes[next].displacement.u1 - contour->nodes[i].displacement.u1;
        double du2_next = contour->nodes[next].displacement.u2 - contour->nodes[i].displacement.u2;
        
        // Central difference approximation
        double h_prev = sqrt(dx_prev*dx_prev + dy_prev*dy_prev);
        double h_next = sqrt(dx_next*dx_next + dy_next*dy_next);
        
        if (h_prev > 1e-12 && h_next > 1e-12) {
            // Average gradients from both sides
            double du1_dx1_prev = du1_prev * dx_prev / (h_prev * h_prev);
            double du1_dx2_prev = du1_prev * dy_prev / (h_prev * h_prev);
            double du2_dx1_prev = du2_prev * dx_prev / (h_prev * h_prev);
            double du2_dx2_prev = du2_prev * dy_prev / (h_prev * h_prev);
            
            double du1_dx1_next = du1_next * dx_next / (h_next * h_next);
            double du1_dx2_next = du1_next * dy_next / (h_next * h_next);
            double du2_dx1_next = du2_next * dx_next / (h_next * h_next);
            double du2_dx2_next = du2_next * dy_next / (h_next * h_next);
            
            contour->nodes[i].du_dx.du1_dx1 = 0.5 * (du1_dx1_prev + du1_dx1_next);
            contour->nodes[i].du_dx.du1_dx2 = 0.5 * (du1_dx2_prev + du1_dx2_next);
            contour->nodes[i].du_dx.du2_dx1 = 0.5 * (du2_dx1_prev + du2_dx1_next);
            contour->nodes[i].du_dx.du2_dx2 = 0.5 * (du2_dx2_prev + du2_dx2_next);
        } else {
            // Set to zero if unable to compute
            memset(&contour->nodes[i].du_dx, 0, sizeof(DisplacementGradient2D));
        }
    }
    
    printf("Computed displacement gradients using finite differences\n");
}

// Initialize fields from external data file (e.g., FE results)
int load_fields_from_file(Contour *contour, const char *filename) {
    if (contour == NULL || filename == NULL) {
        fprintf(stderr, "Error: Invalid parameters for field loading\n");
        return -1;
    }
    
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return -1;
    }
    
    char line[256];
    int node_count = 0;
    
    // Skip header if present
    if (fgets(line, sizeof(line), file) != NULL) {
        if (strncmp(line, "#", 1) == 0 || strncmp(line, "x", 1) == 0) {
            // Header detected, continue
        } else {
            // No header, rewind
            rewind(file);
        }
    }
    
    while (fgets(line, sizeof(line), file) != NULL && node_count < contour->num_nodes) {
        double x, y, u1, u2, s11, s22, s12, e11, e22, e12;
        
        int items = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                          &x, &y, &u1, &u2, &s11, &s22, &s12, &e11, &e22, &e12);
        
        if (items >= 4) {
            // Update node coordinates if provided
            contour->nodes[node_count].coord.x = x;
            contour->nodes[node_count].coord.y = y;
            contour->nodes[node_count].displacement.u1 = u1;
            contour->nodes[node_count].displacement.u2 = u2;
            
            if (items >= 7) {
                contour->nodes[node_count].stress.s11 = s11;
                contour->nodes[node_count].stress.s22 = s22;
                contour->nodes[node_count].stress.s12 = s12;
                contour->nodes[node_count].stress.s21 = s12; // Assume symmetry
            }
            
            if (items >= 10) {
                contour->nodes[node_count].strain.e11 = e11;
                contour->nodes[node_count].strain.e22 = e22;
                contour->nodes[node_count].strain.e12 = e12;
                contour->nodes[node_count].strain.e21 = e12; // Assume symmetry
            }
            
            node_count++;
        }
    }
    
    fclose(file);
    
    if (node_count != contour->num_nodes) {
        fprintf(stderr, "Warning: Loaded %d nodes, expected %d\n", node_count, contour->num_nodes);
    }
    
    printf("Loaded field data from %s for %d nodes\n", filename, node_count);
    return 0;
} 