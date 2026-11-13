#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include"lib_radcouple.h"

#define MAXPOINTS 1000  // Maximum number of ray path points (including reflections).

// Store surface geometry and Fibonacci direction samples globally.
SURFACE * SUR = NULL;
int nsurf = 0;
VEC * fibo_dirs = NULL;
int ndirs = 0;

//-------------------------------------------------------------------------------------------------
// Allocate per-ray tracing workspace.
void init_trace_work(TRACE_WORK *w, int32_t n)
{
  w->nalloc = n;
  w->dist = malloc(n * sizeof(double));
  w->Qq   = malloc(n * sizeof(VEC));
  w->N    = malloc(n * sizeof(VEC));

  if(w->dist == NULL || w->Qq == NULL || w->N == NULL){
    fprintf(stderr, "TRACE_WORK allocation error\n");
    exit(1);
  }
}
void free_trace_work(TRACE_WORK *w)
{
  free(w->dist);
  free(w->Qq);
  free(w->N);

  w->dist = NULL;
  w->Qq = NULL;
  w->N = NULL;
  w->nalloc = 0;
}

void free_vec_array(VEC *p)
{
  free(p);
}
//-------------------------------------------------------------------------------------------------
// Accumulate cosine-weighted emission by outgoing direction.
void count_out_dirs(
  VEC *StartPoints,      // Array of ray origins
  int32_t num_start,     // Number of ray origins
  int32_t start_surf_id, // Index of the emitting surface
  VEC *Directions,       // Array of initial ray directions
  int32_t num_dir,       // Number of initial directions
  double *dir_emission,     // Cosine-weighted emission escaping in each direction
  double *total_emission,   // Total emitted weight, including absorbed rays
  double *absorbed_emission // Cosine-weighted emission absorbed by surfaces
){
int32_t total, absorbed;

  total = num_start * num_dir;
  absorbed = 0;

#pragma omp parallel
  {
    double *local_emission = calloc(ndirs, sizeof(double));
    double local_tot_emission = 0.;
    double local_absorbed = 0.;
    RAY work;
    TRACE_WORK tw;

    work.points = malloc(MAXPOINTS * sizeof(VEC));
    work.id = malloc(MAXPOINTS * sizeof(int32_t));
    work.n = 0;
    work.f = HIT_NONE;
    work.final_dir = -1;
    init_trace_work(&tw, nsurf);

#pragma omp for schedule(dynamic)
    for(int32_t k = 0; k < total; ++k){
      int32_t i = k / num_dir;   // Index of the ray origin
      int32_t j = k % num_dir;   // Index of the initial direction
      // Trace rays individually.
      RAY ray;
      VEC final;
      int32_t out_id;

      // Initialize the ray record.
      ray.n = 0; ray.f = HIT_NONE; ray.points = NULL; ray.id = NULL; ray.final_dir = -1;

      // Trace one ray without adding a distant spherical-shell intersection.
      trace_one_ray(StartPoints[i],start_surf_id, Directions[j],0,3000.0, &ray, &final, &work,&tw);

      local_tot_emission += Directions[j].z ;  // Count the initial cos(theta) weight regardless of the tracing outcome.
      if(ray.f == HIT_ABSORB){ // Ray absorbed by a surface
        local_absorbed += Directions[j].z ;  // Add its cos(theta) weight to the absorbed total.
      }else if(ray.f == HIT_NONE){  // Ray escapes to infinity
        out_id = get_fibonacci_id(final); // Find the closest Fibonacci direction sample.
        if(out_id >= 0){
          local_emission[out_id] += Directions[j].z ;  // Accumulate the escaping cos(theta) weight.
        }
      }

      free(ray.points);
      free(ray.id);
    }
    free(work.points);
    free(work.id);

#pragma omp critical
    {
      for(int32_t i = 0; i < ndirs; ++i){ dir_emission[i] += local_emission[i]; }
      absorbed += local_absorbed;
      *total_emission += local_tot_emission;
    }
    free(local_emission);
  }

  *absorbed_emission = absorbed;
}
//-------------------------------------------------------------------------------------------------
// Generate and trace num_start * num_dir rays.
RAY * gen_rays(
  VEC *StartPoints,      // Array of ray origins
  int32_t num_start,     // Number of ray origins
  int32_t start_surf_id, // Index of the emitting surface
  VEC *Directions,       // Array of initial ray directions
  int32_t num_dir,       // Number of initial directions
  int32_t out_id_flag,   // Whether to identify the outgoing Fibonacci hemisphere direction
  int32_t add_shell, double shell_radius // Whether to append a distant shell hit, and the shell radius
){
RAY * rays;
int32_t i,j,k;
VEC final;
RAY work;
TRACE_WORK tw;

  // Allocate ray records and reusable tracing workspace.
  rays = (RAY *)malloc(num_start*num_dir*sizeof(RAY));
  work.points = malloc(MAXPOINTS * sizeof(VEC));
  work.id = malloc(MAXPOINTS * sizeof(int32_t));
  work.n = 0; work.f = HIT_NONE;
  work.final_dir = -1;
  init_trace_work(&tw, nsurf);
  // Loop over ray origins and initial directions.
  for(i=0;i<num_start;++i){
  for(j=0;j<num_dir;++j){
    k = j + i*num_dir;
    // Initialize the output record defensively.
    rays[k].n = 0; rays[k].f = HIT_NONE; rays[k].points = NULL; rays[k].id = NULL;
    // Trace the ray.
    trace_one_ray(
      StartPoints[i], start_surf_id, Directions[j], // Ray origin, emitting surface, and initial direction
      add_shell, shell_radius , // Optional distant-shell intersection and radius
      &rays[k],                 // Output ray record
      &final,                   // Final outgoing direction
      &work,
      &tw
    );
    if(out_id_flag && rays[k].f == HIT_NONE){ rays[k].final_dir = get_fibonacci_id(final); }
    else{rays[k].final_dir = -1;}
  }}
  free(work.points);
  free(work.id);
  free_trace_work(&tw);
  return(rays);
}
void free_rays(RAY *rays, int32_t n)
{
  if(rays == NULL) return;

  for(int32_t i = 0; i < n; ++i){
    free(rays[i].points);
    free(rays[i].id);
  }

  free(rays);
}
// Trace a single ray through successive surface interactions. ----------------------------------------------------------------------
void trace_one_ray(
  VEC start, int32_t start_surface,       // Ray origin and emitting surface index
  VEC initial_dir,                        // Initial ray direction
  int32_t add_shell, double shell_radius, // Optional distant-shell intersection and radius
  RAY *ray,                               // Output ray record
  VEC *final,                             // Final outgoing direction
  RAY *workray,                           // Reusable ray-path workspace
  TRACE_WORK *tw                          // Intersection workspace
){
int32_t n, ret,next_surf;
VEC U,Q,V,N,Org;

  Org.x = Org.y = Org.z = 0;

  // Initialize the path length.
  n = 0;
  // Record the initial point and its surface index.
  workray->points[n] = start; workray->id[n] = start_surface; n++;
  U = initial_dir;
  // Continue until the ray is absorbed or escapes.
  while(1){
    // printf("trace_one_ray loop\n");
    ret = trace_step(workray->id[n-1], &next_surf, workray->points[n-1], U, &Q, &V, tw );
    // printf("return from trace_step is %d REFLECT:2 ABSORB:1 NONE:0 from %d to %d\n",ret,
    //          workray->id[n-1], next_surf);
    if(ret == HIT_NONE) { // No further intersection: the ray escapes.
      *final = U ; // The current direction is the final outgoing direction.
      break; // Exit without appending another path point.
    }else if(ret == HIT_ABSORB){// Ray reaches a non-reflecting, absorbing surface.
      final->x = final->y = final->z = 0.; // Use a zero vector for the final direction of an absorbed ray.
      workray->points[n] = Q;workray->id[n] = next_surf; n++; break; // Record the terminal intersection and stop.
    }else { // A reflected ray continues toward another surface.
      if(n >= MAXPOINTS-1) { // Stop if the ray path would exceed the allocated capacity.
        fprintf(stderr,"Too many reflection \n"); exit(1); }
      // Append the reflection point and continue tracing.
      workray->points[n] = Q; workray->id[n] = next_surf; U = V ; n++;
    }
  }
  // printf("trace_one_ray retrun code %d\n",ret);
  // printvec("trace_one_ray",*final);
  // Ray tracing is complete.
  // Optionally intersect an escaping ray with a distant spherical shell.
  if(ret == HIT_NONE && add_shell){
    if(cross_sphere(Org, shell_radius, workray->points[n-1], U, &Q, &N) != 1){ // The ray does not intersect the requested spherical shell.
      workray->n = n-1;
      fprintf(stderr,"Error not in the sphere shell\n");exit(1);
    }
    // Append the shell intersection without changing HIT_NONE.
    workray->points[n] = Q; workray->id[n] = SHELL; n++;
  }
  // Copy the n recorded path points into the output ray.
  ray->n = n; ray->f = ret;
  ray->points = (VEC *)malloc(n*sizeof(VEC));
  ray->id     = (int32_t *)malloc(n*sizeof(int32_t));
  memcpy(ray->points,workray->points,n * sizeof(VEC));
  memcpy(ray->id,workray->id,n * sizeof(int32_t));

  return;
}
//-------------------------------------------------------------------------------------------------
// Advance a ray by one surface interaction. -------------------------------------------------------------------
int32_t trace_step(
  int32_t current_surface, // Index of the current surface
  int32_t *next_surface, // Index of the next intersected surface
  VEC  P, // Current ray origin
  VEC  U, // Current ray direction
  VEC *Q, // Next intersection point
  VEC *V,  // Direction after the interaction
  TRACE_WORK *tw // Reusable intersection workspace
){
/*
 * Find the nearest surface intersected by the ray from P along U.
 * Return codes:
 *   HIT_NONE:   no surface intersection
 *   HIT_ABSORB: intersection with an absorbing surface
 *   HIT_REFLECT: intersection with a reflecting surface
 */
int32_t k, min_k, r;
double min_dist, un;
double *dist = tw->dist;
VEC *Qq = tw->Qq;
VEC *N = tw->N;

  // printf("Entering trace_step  nsurf=%d  current_surface,=%d\n",nsurf,current_surface);
  min_dist = DBL_MAX; min_k = -1;
  // Test each surface and calculate the intersection distance and normal.
  for(k=1;k<nsurf;++k){  // Surface 0 is a dummy entry.
    if(k==current_surface) { continue; } // Skip the surface from which the ray departs.
    // printf("trace_step surface %d\n",k);
    if(SUR[k].type == RECT){
      r = rect_cross(SUR[k].p0,SUR[k].p1,SUR[k].p2,SUR[k].p3,P,U,&(N[k]),&(Qq[k]),&(dist[k]));
    }else if(SUR[k].type == CYLIND){
      r = cylinder_cross(SUR[k].p0,SUR[k].p1,SUR[k].r,SUR[k].h,P,U,&(N[k]),&(Qq[k]),&(dist[k]));
    }else if(SUR[k].type == TRI){
      r = tri_cross(SUR[k].p0,SUR[k].p1,SUR[k].p2,P,U,&(N[k]),&(Qq[k]),&(dist[k]));
    }else if(SUR[k].type == DISC){
      r = disc_cross(SUR[k].p0,SUR[k].p1,SUR[k].r,P,U,&(N[k]),&(Qq[k]),&(dist[k]));
    }else{
      fprintf(stderr,"Surface type error \n") ; exit(1);
    }
    if(r == 1){ // A valid intersection has a positive distance.
      // printf("Type %d dist %f\n",SUR[k].type ,  dist[k]);
      if(min_dist > dist[k]) { min_dist = dist[k] ; min_k = k; } // Keep the closest intersection found so far.
    }
  }
  // min_k identifies the nearest intersected surface.
   // printf("Surface %d is the nearest\n",min_k);
  if(min_k<0) { V->x = V->y = V->z = 0.; r = HIT_NONE;  // No surface is intersected.
  }else{
    k = min_k ; // Index of the intersected surface
    *Q = Qq[k];*next_surface = k; // Return the intersection point and surface index.
    // Determine whether the surface reflects this ray.
    un = dot_product(U,N[k]);
    // printf("Reflect ? suface.reflect %d , (U,N) %f  at %d %s\n",SUR[k].reflect,un,k,SUR[k].name);
    if((SUR[k].reflect & REF_FRONT) && (un <= 0.0)){ // Reflect when incident on the front face.
      r=HIT_REFLECT; // printf("REFLECT\n");
    }else if((SUR[k].reflect & REF_BACK) && (un >= 0.0)){ // Reflect when incident on the back face.
      r=HIT_REFLECT; // printf("REFLECT\n");
    }else { // Absorb when neither reflection condition is met.
      r = HIT_ABSORB; // printf("ABSORB\n");
    }
    if(r == HIT_REFLECT) { // Reflected ray
      *V = reflection(U,N[k]); // Calculate the reflected direction.
    }else{
      V->x = V->y = V->z = 0. ; // An absorbed ray has no outgoing direction.
    }
  }
  return(r);
}
//-------------------------------------------------------------------------------------------------
int32_t get_fibonacci_id(VEC F){
int32_t i0,i;
double c,bestc,maxc;
double rf;
int32_t best;

  if(ndirs == 0) { //fprintf(stderr,"No fibonacci data\n");
    return(-1);} // No direction samples have been registered.

  if(fabs(dot_product(F,F))<EPS_NORM) return(-1);
  F = normalize(F);

  if(F.z < 0.0){
    fprintf(stderr,"WARNING :: Output Ray goes under hemisphere\n"); return(-1);
  }

  // For Fibonacci hemisphere samples, z is approximately 1 - i/N; start near F.z.
  i0 = (int32_t)(ndirs*(1.0 - F.z));
  if(i0<0) i0=0;
  if(i0>=ndirs) i0=ndirs-1;
  // Initial best candidate
  best = i0;
  bestc = dot_product(F,fibo_dirs[i0]);
  // Magnitude of the projection onto the XY plane
  rf = sqrt(fmax(0.0,1.0 - F.z*F.z));

  // Search toward increasing indices.
  for(i=i0+1;i<ndirs;++i){
    maxc = F.z*fibo_dirs[i].z + rf*sqrt(fmax(0.0,1.0 - fibo_dirs[i].z*fibo_dirs[i].z));
    if(maxc < bestc) break;
    c = dot_product(F,fibo_dirs[i]);
    if(c > bestc){ bestc = c; best = i; }
  }

  // Search toward decreasing indices.
  for(i=i0-1;i>=0;--i){
    maxc = F.z*fibo_dirs[i].z + rf*sqrt(fmax(0.0,1.0 - fibo_dirs[i].z*fibo_dirs[i].z));
    if(maxc < bestc) break;
    c = dot_product(F,fibo_dirs[i]);
    if(c > bestc){ bestc = c; best = i; }
  }

  return(best);
}
//-------------------------------------------------------------------------------------------------
// Register the surface geometry.
void set_surfaces(SURFACE *src, int n, int debug)
{
  nsurf = n;
  free(SUR); // Release the previously registered surfaces, if any.
  SUR = malloc(n * sizeof(SURFACE));

  if(SUR == NULL){ fprintf(stderr, "surface allocation error\n"); exit(1); }

  memcpy(SUR, src, n * sizeof(SURFACE));

  if(debug)  print_surfaces();
}

void print_surfaces(void)
{
  printf("Number of surfaces = %d\n", nsurf);

  for(int i = 0; i < nsurf; ++i){

    printf("--------------------------------------------------\n");

    printf("Surface %d\n", i);

    printf("  name    = %s\n", SUR[i].name);
    printf("  type    = %d\n", SUR[i].type);
    printf("  reflect = %d\n", SUR[i].reflect);

    printf("  r       = %f\n", SUR[i].r);
    printf("  h       = %f\n", SUR[i].h);

    printvec("p0", SUR[i].p0);
    printvec("p1", SUR[i].p1);
    printvec("p2", SUR[i].p2);
    printvec("p3", SUR[i].p3);
    printvec("p4", SUR[i].p4);
  }

  printf("--------------------------------------------------\n");
}
//-------------------------------------------------------------------------------------------------
// Generate approximately uniform directions on a sphere using a Fibonacci lattice.
// Allocate and return an array of N three-dimensional unit vectors.
VEC * fibonacci_sphere(int32_t N){
// int32_t fibonacci_sphere(int32_t N, VEC **Pnt){
int32_t i;
double w,theta,r;
VEC *P;

  // Allocate the vector array.
  if((P = (VEC *)malloc(N*sizeof(VEC)))==NULL){
    fprintf(stderr,"struct P allocation error\n");
    return(NULL);
  }

  for (i = 0; i < N; i++) {
    w = 1.0 - (2.0 * i / (double)N);  // Uniform spacing in the z coordinate
    P[i].z = w;
    theta = i * GOLDEN_ANGLE;  // Distribute azimuths using the golden angle.
    r = sqrt(1.0 - w * w);  // Radius of the horizontal circle at z = w
    P[i].x = r * cos(theta);
    P[i].y = r * sin(theta);
  }
  return P;
}
//-------------------------------------------------------------------------------------------------
// Generate approximately uniform directions on the upper hemisphere using a Fibonacci lattice.
// Allocate and return an array of N three-dimensional unit vectors.
VEC * fibonacci_sphere_half(int32_t N, int32_t keep_flag){
int32_t i;
double w,theta,r;
VEC *P;

  // Allocate the vector array.
  if((P = (VEC *)malloc(N*sizeof(VEC)))==NULL){
    fprintf(stderr,"struct P allocation error\n");
    // return(-1);
    return(NULL);
  }

  for (i = 0; i < N; i++) {
    w = 1.0 - (1.0 * i / (double)N);  // Uniform spacing in the z coordinate
    P[i].z = w;
    theta = i * GOLDEN_ANGLE;  // Distribute azimuths using the golden angle.
    r = sqrt(1.0 - w * w);  // Radius of the horizontal circle at z = w
    P[i].x = r * cos(theta);
    P[i].y = r * sin(theta);
  }
  if(keep_flag) {   // Retain the generated direction array for subsequent lookups.
     fibo_dirs = P ;
     ndirs = N;
  }

  return P;
}
