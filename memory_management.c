#include "viscoelastic_j_integral.h"

// Allocate memory for contour nodes
ContourNode* allocate_contour_nodes(int num_nodes) {
    if (num_nodes <= 0) {
        fprintf(stderr, "Error: Invalid number of nodes (%d)\n", num_nodes);
        return NULL;
    }
    
    ContourNode *nodes = (ContourNode*)calloc(num_nodes, sizeof(ContourNode));
    if (nodes == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for %d contour nodes\n", num_nodes);
        return NULL;
    }
    
    return nodes;
}

// Free memory allocated for contour nodes
void free_contour_nodes(ContourNode *nodes) {
    if (nodes != NULL) {
        free(nodes);
    }
}

// Create and initialize viscoelastic model
ViscoelasticModel* create_viscoelastic_model(int num_prony_terms) {
    if (num_prony_terms <= 0 || num_prony_terms > MAX_PRONY_TERMS) {
        fprintf(stderr, "Error: Invalid number of Prony terms (%d). Must be between 1 and %d\n", 
                num_prony_terms, MAX_PRONY_TERMS);
        return NULL;
    }
    
    ViscoelasticModel *model = (ViscoelasticModel*)malloc(sizeof(ViscoelasticModel));
    if (model == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for viscoelastic model\n");
        return NULL;
    }
    
    model->num_prony_terms = num_prony_terms;
    model->G_inf = 0.0;
    model->nu = 0.3; // Default Poisson's ratio
    
    // Allocate memory for Prony terms
    model->prony_terms = (PronyTerm*)calloc(num_prony_terms, sizeof(PronyTerm));
    if (model->prony_terms == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for Prony terms\n");
        free(model);
        return NULL;
    }
    
    // Initialize each Prony term
    for (int i = 0; i < num_prony_terms; i++) {
        model->prony_terms[i].G_k = 0.0;
        model->prony_terms[i].tau_k = 1.0;
        model->prony_terms[i].stress_history = NULL;
    }
    
    return model;
}

// Free memory allocated for viscoelastic model
void free_viscoelastic_model(ViscoelasticModel *model) {
    if (model == NULL) return;
    
    if (model->prony_terms != NULL) {
        // Free stress history for each Prony term
        for (int i = 0; i < model->num_prony_terms; i++) {
            if (model->prony_terms[i].stress_history != NULL) {
                free(model->prony_terms[i].stress_history);
            }
        }
        free(model->prony_terms);
    }
    
    free(model);
}

// Allocate stress history for Prony terms (called after contour is initialized)
int allocate_stress_history(ViscoelasticModel *model, int num_nodes) {
    if (model == NULL || num_nodes <= 0) {
        fprintf(stderr, "Error: Invalid parameters for stress history allocation\n");
        return -1;
    }
    
    for (int i = 0; i < model->num_prony_terms; i++) {
        model->prony_terms[i].stress_history = (StressTensor2D*)calloc(num_nodes, sizeof(StressTensor2D));
        if (model->prony_terms[i].stress_history == NULL) {
            fprintf(stderr, "Error: Failed to allocate stress history for Prony term %d\n", i);
            
            // Clean up previously allocated memory
            for (int j = 0; j < i; j++) {
                free(model->prony_terms[j].stress_history);
                model->prony_terms[j].stress_history = NULL;
            }
            return -1;
        }
    }
    
    return 0;
}

// Allocate results array for time simulation
JIntegralResult* allocate_results_array(int num_time_steps) {
    if (num_time_steps <= 0) {
        fprintf(stderr, "Error: Invalid number of time steps (%d)\n", num_time_steps);
        return NULL;
    }
    
    JIntegralResult *results = (JIntegralResult*)calloc(num_time_steps, sizeof(JIntegralResult));
    if (results == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for %d results\n", num_time_steps);
        return NULL;
    }
    
    return results;
}

// Free results array
void free_results_array(JIntegralResult *results) {
    if (results != NULL) {
        free(results);
    }
} 