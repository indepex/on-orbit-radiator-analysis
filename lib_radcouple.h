#include <stdint.h>
#include <string.h>
#include "lib_vec.h"

#define WORDLENGTH 64

#define HIT_NONE     0 // Ray escapes without hitting another surface
#define HIT_ABSORB   1 // Ray hits an absorbing surface
#define HIT_REFLECT  2 // Ray hits a reflecting surface

// Surface reflection flags (bitmask); zero means absorption
#define REF_FRONT 1  // Reflect from the front face
#define REF_BACK  2  // Reflect from the back face

#define RECT   0
#define CYLIND 1
#define CUBE   2
#define TRI    3
#define DISC   4

#define SHELL -1


// Ray path: vertices, termination status, intersected surfaces, and outgoing direction ID
typedef struct {
  int32_t n;     // Number of stored points along the ray
  int32_t f;     // Ray termination status (HIT_*)
  VEC * points;  // Ray path points
  int32_t * id;  // Surface ID associated with each path point
  int final_dir; // Index of the final outgoing direction
} RAY;

typedef struct { char name[WORDLENGTH]; int type; int reflect; double r; double h;
                 VEC p0; VEC p1; VEC p2; VEC p3; VEC p4;} SURFACE;

// Reusable workspace for ray tracing
typedef struct {
  double *dist;
  VEC *Qq;
  VEC *N;
  int32_t nalloc;
} TRACE_WORK;


void init_trace_work(TRACE_WORK *w, int32_t n);
void free_trace_work(TRACE_WORK *w);
void free_vec_array(VEC *);

// Accumulate cosine-weighted emission by outgoing direction
void count_out_dirs(
  VEC *StartPoints,      // Array of ray origins
  int32_t num_start,     // Number of ray origins
  int32_t start_surf_id, // Surface ID of ray origins
  VEC *Directions,       // Array of initial ray directions
  int32_t num_dir,       // Number of initial directions
  double *dir_emission,     // Emission escaping to infinity, binned by direction
  double *total_emission,   // Total emitted weight, including absorbed rays
  double *absorbed_emission // Weight absorbed by surfaces
);

void trace_one_ray(
  VEC start, int32_t start_surface,       // Ray origin and originating surface ID
  VEC initial_dir,                        // Initial ray direction
  int32_t add_shell, double shell_radius, // Append intersection with outer sphere; sphere radius
  RAY *ray,                               // Output ray path
  VEC *final,                             // Final outgoing direction
  RAY *workray,                           // Reusable ray-path buffer
  TRACE_WORK *tw                          // Ray-tracing workspace
);
void free_rays(RAY *rays, int32_t n);
int32_t trace_step(
  int32_t current_surface, // ID of the current (originating) surface
  int32_t *next_surface, // ID of the next intersected surface
  VEC  P, // Ray origin
  VEC  U, // Incident ray direction
  VEC *Q, // Intersection point
  VEC *V, // Outgoing (reflected) direction
  TRACE_WORK *tw // Ray-tracing workspace
);


int32_t get_fibonacci_id(VEC V);
void set_surfaces(SURFACE *, int, int);
void print_surfaces(void);

VEC * fibonacci_sphere(int32_t);
VEC * fibonacci_sphere_half(int32_t, int32_t);
#define GOLDEN_ANGLE (M_PI * (3.0 - sqrt(5.0)))  // Approximately 137.508 degrees
