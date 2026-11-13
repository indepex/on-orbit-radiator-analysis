#define EPS_NORM 1e-10

// Data types ----------------------------------------------------------------------------------
// Three-dimensional vector
typedef struct { double x, y, z; } VEC;


// lib_vec.c ------------------------------------------------------------------------------------
// Geometry utilities for ray tracing
// Ray-surface intersection routines
int32_t cross_sphere(VEC, double, VEC, VEC,  VEC *,  VEC *);

// Intersection routines
// Intersect ray (P, U) with a cylinder (center C, axis D, height h, radius r).
// Return the intersection Q, surface normal N, and ray parameter k.
int32_t cylinder_cross(VEC C,VEC D,double r,double h,VEC P,VEC U, VEC *N, VEC *Q, double *k);
// Intersect a ray (P, U) with quadrilateral ABCD; return Q, N, and distance dist
int32_t rect_cross(VEC A, VEC B, VEC C, VEC D, VEC P, VEC U, VEC *N, VEC *Q, double *dist);
// Intersect a ray (P, U) with triangle ABC; return Q, N, and distance dist
int32_t tri_cross(VEC A, VEC B, VEC C, VEC P, VEC U, VEC *N, VEC *Q, double *dist);
// Intersect ray (P, U) with a disk (center C, normal D, radius r).
// Return the intersection Q, surface normal N, and ray parameter k.
int32_t disc_cross(VEC C,VEC D,double r,VEC P,VEC U, VEC *N, VEC *Q, double *k);
int32_t inside_triangle(VEC P, VEC A, VEC B, VEC C);

// Specular reflection
VEC reflection(const VEC, const VEC);
// Basic vector operations
VEC linear1(double, VEC);
VEC linear2(double, VEC, double, VEC);
VEC linear3(double, VEC, double, VEC, double, VEC);
double dot_product(VEC, VEC);
VEC cross_product(VEC, VEC);
VEC normalize(VEC X);
void printvec(char *, VEC);
