#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;


typedef struct particle_s {
    double x, y, z;  
    double vx, vy, vz;  
    double mass;  
} particle_t;


typedef struct octree_node_s {
    double mass;  
    double cx, cy, cz;  
    struct octree_node_s* children[8];  
    particle_t* particle;  
    bounds_t bounds;
} octree_node_t;


typedef struct bounds_s {
    double minX, maxX;
    double minY, maxY;
    double minZ, maxZ;
} bounds_t;

/*
 * Définit les limites d'un nœud de l'octree.
 *
 * @param node: Pointeur vers le nœud de l'octree à configurer.
 * @param bounds: Les limites à définir pour ce nœud.
 */
void set_bounds(octree_node_t* node, bounds_t bounds) {
    node->bounds = bounds;
}

/*
 * Crée un nouveau nœud de l'octree.
 *
 * @return Pointeur vers le nouveau nœud de l'octree.
 */
octree_node_t* create_octree_node() {
    octree_node_t* node = (octree_node_t*)malloc(sizeof(octree_node_t));
    if (!node) {
        
        return NULL;
    }

    
    node->mass = 0.0;
    node->cx = node->cy = node->cz = 0.0;
    for (int i = 0; i < 8; ++i) {
        node->children[i] = NULL;
    }
    node->particle = NULL;

    return node;
}

/*
 * Détermine dans quel octant d'un nœud donné une particule se trouve.
 *
 * @param node: Le nœud pour lequel déterminer l'octant.
 * @param particle: La particule pour laquelle déterminer l'octant.
 * @return L'indice de l'octant dans lequel la particule se trouve.
 */
int get_octant(octree_node_t* node, particle_t* particle) {
    int octant = 0;
    if (particle->x > node->cx) octant |= 4;
    if (particle->y > node->cy) octant |= 2;
    if (particle->z > node->cz) octant |= 1;
    return octant;
}

/*
 * Ajoute une particule à l'octree.
 *
 * @param node: Nœud de l'octree auquel ajouter la particule.
 * @param particle: Particule à ajouter à l'octree.
 * @param boundarySize: Taille des limites du nœud actuel.
 */
void add_particle_to_octree(octree_node_t* node, particle_t* particle, double boundarySize) {
    if (node->particle == NULL && node->mass == 0) {
        
        node->particle = particle;
    } else {
        
    if (node->particle != NULL) {
        
        divide_and_distribute(node);
        
        int octant = get_octant(node, particle);
        add_particle_to_octree(node->children[octant], particle, ...);
    }

        
        int octant = get_octant(node, particle);
        if (node->children[octant] == NULL) {
            node->children[octant] = create_octree_node();
            
            
        }
        add_particle_to_octree(node->children[octant], particle, boundarySize / 2);
    }

    
    
}

/*
 * Divise un nœud de l'octree et distribue les particules aux sous-nœuds.
 *
 * @param node: Nœud de l'octree à diviser.
 */
void divide_and_distribute(octree_node_t* node) {
    
    if (node->particle != NULL && node->children[0] == NULL) {
        
        double centerX = (node->bounds.minX + node->bounds.maxX) / 2;
        double centerY = (node->bounds.minY + node->bounds.maxY) / 2;
        double centerZ = (node->bounds.minZ + node->bounds.maxZ) / 2;

        
        for (int i = 0; i < 8; i++) {
            node->children[i] = create_octree_node();

            
            bounds_t childBounds;
            childBounds.minX = (i & 4) ? centerX : node->bounds.minX;
            childBounds.maxX = (i & 4) ? node->bounds.maxX : centerX;
            childBounds.minY = (i & 2) ? centerY : node->bounds.minY;
            childBounds.maxY = (i & 2) ? node->bounds.maxY : centerY;
            childBounds.minZ = (i & 1) ? centerZ : node->bounds.minZ;
            childBounds.maxZ = (i & 1) ? node->bounds.maxZ : centerZ;

            set_bounds(node->children[i], childBounds);
        }

        
        int octant = get_octant(node, node->particle);
        add_particle_to_octree(node->children[octant], node->particle, /*...*/);

        
        node->particle = NULL;
    }
}

/*
 * Met à jour la masse et le centre de masse d'un nœud de l'octree.
 *
 * @param node: Nœud de l'octree pour lequel mettre à jour la masse et le centre de masse.
 * @param particle: Particule contribuant à la masse et au centre de masse du nœud.
 */
void update_mass_and_center(octree_node_t* node, particle_t* particle) {
    if (node->mass == 0.0) {
        
        node->mass = particle->mass;
        node->cx = particle->x;
        node->cy = particle->y;
        node->cz = particle->z;
    } else {
        
        double totalMass = node->mass + particle->mass;
        node->cx = (node->cx * node->mass + particle->x * particle->mass) / totalMass;
        node->cy = (node->cy * node->mass + particle->y * particle->mass) / totalMass;
        node->cz = (node->cz * node->mass + particle->z * particle->mass) / totalMass;
        node->mass = totalMass;
    }
}


void add_particle_to_octree(octree_node_t* node, particle_t* particle, double boundarySize) {
    if (node->particle == NULL && node->mass == 0) {
        node->particle = particle;
    } else {
        if (node->particle != NULL) {
            
            divide_and_distribute(node);
        }

        int octant = get_octant(node, particle);
        if (node->children[octant] == NULL) {
            node->children[octant] = create_octree_node();
            
            
        }
        add_particle_to_octree(node->children[octant], particle, boundarySize / 2);
    }

    
    update_mass_and_center(node, particle);
}


/*
 * Construit l'octree à partir d'un ensemble de particules.
 *
 * @param node: Racine de l'octree à construire ou à mettre à jour.
 * @param particles: Tableau de particules à ajouter à l'octree.
 * @param n: Nombre de particules dans le tableau.
 * @param boundarySize: Taille des limites de l'octree.
 */
void build_octree(octree_node_t* node, particle_t* particles, int n, double boundarySize) {
    for (int i = 0; i < n; ++i) {
        add_particle_to_octree(node, &particles[i], boundarySize);
    }
}

/*
 * Calcule récursivement le centre de masse pour chaque nœud de l'octree.
 *
 * @param node: Nœud de l'octree pour lequel calculer le centre de masse.
 */
void compute_center_of_mass(octree_node_t* node) {
    if (node == NULL || (node->particle != NULL && node->mass != 0)) {
        
        return;
    }

    double totalMass = 0;
    double cx = 0, cy = 0, cz = 0;

    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != NULL) {
            compute_center_of_mass(node->children[i]);
            totalMass += node->children[i]->mass;
            cx += node->children[i]->mass * node->children[i]->cx;
            cy += node->children[i]->mass * node->children[i]->cy;
            cz += node->children[i]->mass * node->children[i]->cz;
        }
    }

    if (totalMass > 0) {
        node->mass = totalMass;
        node->cx = cx / totalMass;
        node->cy = cy / totalMass;
        node->cz = cz / totalMass;
    }
}

/*
 * Calcule la force exercée sur une particule en utilisant l'octree.
 *
 * @param particle: Particule pour laquelle calculer la force.
 * @param node: Nœud de l'octree à utiliser pour le calcul de la force.
 * @param theta: Critère d'ouverture pour l'approximation de Barnes-Hut.
 */
void compute_force(particle_t *particle, octree_node_t *node, double theta) {
    if (node == NULL || (node->mass == 0 && node->particle == NULL)) {
        
        return;
    }

    if (node->particle == particle) {
        
        return;
    }

    double dx = node->cx - particle->x;
    double dy = node->cy - particle->y;
    double dz = node->cz - particle->z;
    double distance = sqrt(dx * dx + dy * dy + dz * dz);

    if (distance == 0) {
        
        return;
    }

    double d = (node->bounds.maxX - node->bounds.minX); 
    if (node->particle != NULL || (d / distance) < theta) {
        
        double F = (G * particle->mass * node->mass) / (distance * distance * distance);
        particle->fx += F * dx;
        particle->fy += F * dy;
        particle->fz += F * dz;
    } else {
        
        for (int i = 0; i < 8; i++) {
            compute_force(particle, node->children[i], theta);
        }
    }
}

/*
 * Met à jour la position et la vitesse d'une particule en fonction des forces appliquées.
 *
 * @param particle: Particule à mettre à jour.
 * @param dt: Pas de temps pour l'intégration.
 */
void update_particle(particle_t *particle, double dt) {
    
    particle->vx += particle->fx * dt / particle->mass;
    particle->vy += particle->fy * dt / particle->mass;
    particle->vz += particle->fz * dt / particle->mass;

    
    particle->x += particle->vx * dt;
    particle->y += particle->vy * dt;
    particle->z += particle->vz * dt;

    
    particle->fx = 0;
    particle->fy = 0;
    particle->fz = 0;
}



int main() {
    const int num_particles = 1000; 
    const double theta = 0.5; 
    const double dt = 0.01; 
    const int num_steps = 1000; 

    
    particle_t* particles = malloc(num_particles * sizeof(particle_t));
    

    
    octree_node_t* root = create_octree_node();
    
    bounds_t root_bounds = {/* ... */};
    set_bounds(root, root_bounds);

    
    for (int step = 0; step < num_steps; ++step) {
        
        build_octree(root, particles, num_particles, root_bounds);

        
        compute_center_of_mass(root);

        
        for (int i = 0; i < num_particles; ++i) {
            compute_force(&particles[i], root, theta);
        }

        
        for (int i = 0; i < num_particles; ++i) {
            update_particle(&particles[i], dt);
        }

        
      //bench
        
        
    }

    
    free(particles);
    

    return 0;
}
