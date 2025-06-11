#include "viscoelastic_j_integral.h"

// Initialize circular contour around crack tip
void initialize_circular_contour(Contour *contour, Point2D center, double radius, int num_nodes) {
    if (contour == NULL || num_nodes <= 0 || radius <= 0.0) {
        fprintf(stderr, "Error: Invalid parameters for circular contour initialization\n");
        return;
    }
    
    contour->nodes = allocate_contour_nodes(num_nodes);
    if (contour->nodes == NULL) {
        return;
    }
    
    contour->num_nodes = num_nodes;
    
    // Generate nodes along circular contour
    double theta_increment = 2.0 * M_PI / num_nodes;
    
    for (int i = 0; i < num_nodes; i++) {
        double theta = i * theta_increment;
        
        // Node coordinates
        contour->nodes[i].coord.x = center.x + radius * cos(theta);
        contour->nodes[i].coord.y = center.y + radius * sin(theta);
        
        // Outward normal vector (for circular contour, same as radial direction)
        contour->nodes[i].normal.x = cos(theta);
        contour->nodes[i].normal.y = sin(theta);
        
        // Arc length element (constant for circular contour)
        contour->nodes[i].ds = radius * theta_increment;
        
        // Initialize other fields to zero
        contour->nodes[i].displacement.u1 = 0.0;
        contour->nodes[i].displacement.u2 = 0.0;
        
        memset(&contour->nodes[i].stress, 0, sizeof(StressTensor2D));
        memset(&contour->nodes[i].strain, 0, sizeof(StrainTensor2D));
        memset(&contour->nodes[i].du_dx, 0, sizeof(DisplacementGradient2D));
    }
    
    contour->total_length = 2.0 * M_PI * radius;
    
    printf("Initialized circular contour: center=(%.3f, %.3f), radius=%.3f, %d nodes\n",
           center.x, center.y, radius, num_nodes);
}

// Initialize rectangular contour
void initialize_rectangular_contour(Contour *contour, Point2D bottom_left, Point2D top_right, int num_nodes) {
    if (contour == NULL || num_nodes < 4) {
        fprintf(stderr, "Error: Invalid parameters for rectangular contour (need at least 4 nodes)\n");
        return;
    }
    
    contour->nodes = allocate_contour_nodes(num_nodes);
    if (contour->nodes == NULL) {
        return;
    }
    
    contour->num_nodes = num_nodes;
    
    double width = top_right.x - bottom_left.x;
    double height = top_right.y - bottom_left.y;
    double perimeter = 2.0 * (width + height);
    
    // Distribute nodes along perimeter
    int nodes_per_side = num_nodes / 4;
    int extra_nodes = num_nodes % 4;
    
    int node_idx = 0;
    
    // Bottom edge (left to right)
    int bottom_nodes = nodes_per_side + (extra_nodes > 0 ? 1 : 0);
    for (int i = 0; i < bottom_nodes && node_idx < num_nodes; i++, node_idx++) {
        double t = (double)i / (bottom_nodes - 1);
        contour->nodes[node_idx].coord.x = bottom_left.x + t * width;
        contour->nodes[node_idx].coord.y = bottom_left.y;
        contour->nodes[node_idx].normal.x = 0.0;
        contour->nodes[node_idx].normal.y = -1.0; // Outward normal (downward)
        contour->nodes[node_idx].ds = width / (bottom_nodes - 1);
    }
    
    // Right edge (bottom to top)
    int right_nodes = nodes_per_side + (extra_nodes > 1 ? 1 : 0);
    for (int i = 1; i < right_nodes && node_idx < num_nodes; i++, node_idx++) {
        double t = (double)i / (right_nodes - 1);
        contour->nodes[node_idx].coord.x = top_right.x;
        contour->nodes[node_idx].coord.y = bottom_left.y + t * height;
        contour->nodes[node_idx].normal.x = 1.0; // Outward normal (rightward)
        contour->nodes[node_idx].normal.y = 0.0;
        contour->nodes[node_idx].ds = height / (right_nodes - 1);
    }
    
    // Top edge (right to left)
    int top_nodes = nodes_per_side + (extra_nodes > 2 ? 1 : 0);
    for (int i = 1; i < top_nodes && node_idx < num_nodes; i++, node_idx++) {
        double t = (double)i / (top_nodes - 1);
        contour->nodes[node_idx].coord.x = top_right.x - t * width;
        contour->nodes[node_idx].coord.y = top_right.y;
        contour->nodes[node_idx].normal.x = 0.0;
        contour->nodes[node_idx].normal.y = 1.0; // Outward normal (upward)
        contour->nodes[node_idx].ds = width / (top_nodes - 1);
    }
    
    // Left edge (top to bottom)
    for (int i = 1; i < nodes_per_side && node_idx < num_nodes; i++, node_idx++) {
        double t = (double)i / (nodes_per_side - 1);
        contour->nodes[node_idx].coord.x = bottom_left.x;
        contour->nodes[node_idx].coord.y = top_right.y - t * height;
        contour->nodes[node_idx].normal.x = -1.0; // Outward normal (leftward)
        contour->nodes[node_idx].normal.y = 0.0;
        contour->nodes[node_idx].ds = height / (nodes_per_side - 1);
    }
    
    // Initialize other fields to zero
    for (int i = 0; i < num_nodes; i++) {
        contour->nodes[i].displacement.u1 = 0.0;
        contour->nodes[i].displacement.u2 = 0.0;
        memset(&contour->nodes[i].stress, 0, sizeof(StressTensor2D));
        memset(&contour->nodes[i].strain, 0, sizeof(StrainTensor2D));
        memset(&contour->nodes[i].du_dx, 0, sizeof(DisplacementGradient2D));
    }
    
    contour->total_length = perimeter;
    
    printf("Initialized rectangular contour: (%.3f,%.3f) to (%.3f,%.3f), %d nodes\n",
           bottom_left.x, bottom_left.y, top_right.x, top_right.y, num_nodes);
}

// Compute outward normal vectors for general contour
void compute_contour_normals(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL || contour->num_nodes < 3) {
        fprintf(stderr, "Error: Invalid contour for normal computation\n");
        return;
    }
    
    for (int i = 0; i < contour->num_nodes; i++) {
        int prev = (i - 1 + contour->num_nodes) % contour->num_nodes;
        int next = (i + 1) % contour->num_nodes;
        
        // Tangent vector (average of adjacent segments)
        double dx = contour->nodes[next].coord.x - contour->nodes[prev].coord.x;
        double dy = contour->nodes[next].coord.y - contour->nodes[prev].coord.y;
        
        // Normalize tangent
        double magnitude = sqrt(dx*dx + dy*dy);
        if (magnitude > 1e-12) {
            dx /= magnitude;
            dy /= magnitude;
        }
        
        // Normal vector (rotate tangent by 90 degrees, pointing outward)
        contour->nodes[i].normal.x = dy;   // Outward normal
        contour->nodes[i].normal.y = -dx;
    }
}

// Compute arc length elements for contour
void compute_arc_lengths(Contour *contour) {
    if (contour == NULL || contour->nodes == NULL || contour->num_nodes < 2) {
        fprintf(stderr, "Error: Invalid contour for arc length computation\n");
        return;
    }
    
    contour->total_length = 0.0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        int next = (i + 1) % contour->num_nodes;
        
        double dx = contour->nodes[next].coord.x - contour->nodes[i].coord.x;
        double dy = contour->nodes[next].coord.y - contour->nodes[i].coord.y;
        
        contour->nodes[i].ds = sqrt(dx*dx + dy*dy);
        contour->total_length += contour->nodes[i].ds;
    }
}

// Utility function to check if point is inside contour (for validation)
int point_inside_contour(Contour *contour, Point2D point) {
    if (contour == NULL || contour->nodes == NULL) {
        return 0;
    }
    
    int crossings = 0;
    
    for (int i = 0; i < contour->num_nodes; i++) {
        int j = (i + 1) % contour->num_nodes;
        
        Point2D p1 = contour->nodes[i].coord;
        Point2D p2 = contour->nodes[j].coord;
        
        // Ray casting algorithm
        if (((p1.y > point.y) != (p2.y > point.y)) &&
            (point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x)) {
            crossings++;
        }
    }
    
    return (crossings % 2 == 1);
}

// Create multiple concentric contours for path independence testing
void create_concentric_contours(Contour **contours, int num_contours, Point2D center, 
                               double min_radius, double max_radius) {
    if (contours == NULL || num_contours <= 0) {
        fprintf(stderr, "Error: Invalid parameters for concentric contours\n");
        return;
    }
    
    for (int i = 0; i < num_contours; i++) {
        double radius = min_radius + (max_radius - min_radius) * i / (num_contours - 1);
        int num_nodes = 32; // Fixed number of nodes per contour
        
        contours[i] = (Contour*)malloc(sizeof(Contour));
        if (contours[i] == NULL) {
            fprintf(stderr, "Error: Failed to allocate contour %d\n", i);
            continue;
        }
        
        initialize_circular_contour(contours[i], center, radius, num_nodes);
    }
    
    printf("Created %d concentric contours with radii from %.3f to %.3f\n",
           num_contours, min_radius, max_radius);
} 