#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"lib_vec.h"

// Ray origin U, direction V, intersection P, and surface normal N.
//-----------------------------------------------------------------------------------------------
// Ray-sphere intersection.
int32_t cross_sphere(VEC O, double r, VEC U, VEC V,  VEC *P,  VEC *N){
// Sphere: center O, radius r.
// Ray: origin U, unit direction V.
// Intersection: P = U + k V.
// Since |P-O|^2 = r^2, we have |k V + U - O|^2 = r^2.
// Thus k^2 + 2 (U-O,V) k + |U-O|^2 - r^2 = 0, assuming |V|^2 = 1.
VEC W;
double b, c, d, k1, k2, k, kk;

  W = linear2(1, U, -1, O);
  b = 2*dot_product(W,V); c = dot_product(W,W) - r*r;
  d = b*b-4*c;
  if(d < 0){ //printf("d=%f\n",d);
    return(0);} // No intersection.
  k1 = (-b + sqrt(d) ) / 2;
  k2 = (-b - sqrt(d) ) / 2;
  if(k2>k1) { kk = k1; k1 = k2; k2 = kk;}  // Ensure k1 >= k2.
  if(k2>=0){ k = k2;      // Use the nearest intersection along the ray.
  }else if(k1>=0) {k=k1;  // Otherwise, use the farther intersection if it is forward.
  }else {return(0);} // Both intersections lie behind the ray origin.
  *P = linear2(1., U, k, V);
  *N = normalize(linear2(1., *P, -1, O));
  return(1);

}
//-----------------------------------------------------------------------------------------------
// Ray-cylinder intersection.
// Find intersection Q, normal N, and ray parameter k for a cylinder centered at C, with axis D, radius r, and height h.
int32_t cylinder_cross(VEC C,VEC D,double r,double h,VEC P,VEC U, VEC *N, VEC *Q, double *k){
// Cylinder: center C, axis Z = normalize(D), radius r, axial extent [-h/2, h/2].
// Ray: origin P, direction U.
// Intersection: Q = P + k U.
// At axial coordinate z, the vector from the cylinder axis to Q is
// (P + k U) - (C + z Z), which is perpendicular to Z.
// Therefore (P,Z) + k (U,Z) - (C,Z) = z.
// The intersection also lies at radial distance r from the axis:
//      r^2 = |P + k U - C - z Z|^2
// Substitute z from the axial projection into the radial constraint:
//      r^2 = |P + k U - C - ((P,Z) + k (U,Z) - (C,Z)) Z|^2
// Define W = P - C:
//      r^2 = | W + k U - ((W,Z) + k (U,Z)) Z|^2
//          = | W - (W,Z) Z + k (U - (U,Z) Z) |^2
// Define S = W - (W,Z) Z and T = U - (U,Z) Z:
//      r^2 = | S + k T |^2
//     (T,T)k^2 + 2(S,T)k + (S,S)-r*r = 0
// Solve the resulting quadratic equation for k.
VEC Z,W, S, T;
double wz,uz,a,b,c,d,k1,k2,kk,z1,z2,z;
  // printf("Entering cylinder_cross h=%f r=%f\n",h,r);
  Z = normalize(D);
  W = linear2(1, P, -1, C); wz = dot_product(W, Z); uz = dot_product(U, Z);
  S = linear2(1, W, -wz, Z); T = linear2(1, U, -uz, Z);
  a = dot_product(T,T); b = 2*dot_product(S,T); c = dot_product(S,S) - r*r;
  d = b*b-4*a*c;
  if(d < 0){ //printf("d<0\n");
    return(0);} // No real intersection.
  if( fabs(a) < EPS_NORM ) { // Degenerate quadratic coefficient (ray nearly parallel to the axis).
    if( fabs(b)  < EPS_NORM ) { // printf("b=0");
       return(0);} // No intersection.
    *k = -c/b;
    z = wz + (*k)*uz;  // printf("Z = %f\n",z);
    if(!(-h/2 <= z && z <= h/2)) { // printf("z %f out of range %f\n",z,h);a
       return(0);} // Intersection lies outside the finite cylinder.
  }else{
    k1 = (-b + sqrt(d) ) / 2 / a;
    k2 = (-b - sqrt(d) ) / 2 / a;
    if(k2>k1) { kk = k1; k1 = k2; k2 = kk;}  // Ensure k1 >= k2.
    z1 = wz + k1*uz; z2 = wz + k2*uz;//  printf("Z1 = %f  Z2 = %f\n",z1,z2);
    if(k2>=0 && z2>=-h/2 && z2<=h/2 ){ *k = k2; z=z2;  // Prefer the nearer forward intersection within the cylinder height.
    }else if(k1>=0 && z1>=-h/2 && z1<=h/2 ) {*k=k1;z=z1;  // Otherwise, accept the farther forward intersection.
    }else { // printf("z1 %f z2 %f out of range %f\n",z1,z2,h);
      return(0);} // Neither intersection lies on the finite cylinder.
  }
  *Q = linear2(1., P, *k, U);
  *N = normalize(linear3(1., *Q, -1, C, -z, Z));
  // printvec("Cylinder Q",*Q);
  // printvec("Cylinder N",*N);
  return(1); // Intersection found.
}
//-------------------------------------------------------------------------------------------------
int32_t rect_cross(VEC A, VEC B, VEC C, VEC D, VEC P, VEC U, VEC *N, VEC *Q, double *dist){
// Intersect a ray from P in direction U with quadrilateral A-B-C-D; return Q, N, and distance.
VEC X, Y;
double un,pn,an;

  X = linear2(1, B, -1, A);
  Y = linear2(1, D, -1, A);  // Orthogonalization is not needed here.
  *N = normalize(cross_product(X,Y));
  // Q = P + dist U = A + x X + y Y -> (P,N) + dist (U,N) = (A,N)
  un = dot_product(U, (*N)); pn = dot_product(P, (*N)); an = dot_product(A, (*N));
  if(fabs(un) < EPS_NORM ) {return(0);} // Ray is parallel to the plane.
  *dist = (an - pn)/un ;
  if(*dist < 0.0){ return(0);} // Intersection lies behind the ray origin.
  *Q = linear2( 1, P, *dist, U);

  // Test whether Q lies inside either of the two triangles.
  if( inside_triangle(*Q,A,B,C) || inside_triangle(*Q,A,C,D) ){
    return(1); // Intersection found.
  }else{ return(0);} // No intersection within the quadrilateral.
}
//-------------------------------------------------------------------------------------------------
int32_t tri_cross(VEC A, VEC B, VEC C, VEC P, VEC U, VEC *N, VEC *Q, double *dist){
// Intersect a ray from P in direction U with triangle A-B-C; return Q, N, and distance.
VEC X, Y;
double un,pn,an;

  X = linear2(1, B, -1, A);
  Y = linear2(1, C, -1, A);  // Orthogonalization is not needed here.
  *N = normalize(cross_product(X,Y));
  // Q = P + dist U = A + x X + y Y -> (P,N) + dist (U,N) = (A,N)
  un = dot_product(U, (*N)); pn = dot_product(P, (*N)); an = dot_product(A, (*N));
  if(fabs(un) < EPS_NORM ) {return(0);} // Ray is parallel to the plane.
  *dist = (an - pn)/un ;
  if(*dist < 0.0){ return(0);} // Intersection lies behind the ray origin.
  *Q = linear2( 1, P, *dist, U);

  // Test whether Q lies inside the triangle.
  if( inside_triangle(*Q,A,B,C) ){
    return(1); // Intersection found.
  }else{ return(0);} // No intersection within the triangle.
}
//-------------------------------------------------------------------------------------------------
// Intersect a ray from P in direction U with a disk centered at C, normal D, and radius r.
int32_t disc_cross(VEC C,VEC D,double r,VEC P,VEC U, VEC *N, VEC *Q, double *dist){
// The disk plane satisfies (P + dist U - C, D) = 0, so dist (U,D) = (C-P,D).
VEC X,R;
double ud, w;

  X = linear2(1, C, -1, P);
  ud = dot_product(U, D); w = dot_product(X, D);
  if(fabs(ud) < EPS_NORM ) {return(0);} // Ray is parallel to the disk plane.
  *dist = w/ud ;
  *N = normalize(D);
  *Q = linear2(1, P, *dist, U);
  R = linear2(1, *Q, -1, C);
  // Check radial distance and positive ray parameter.
  if(dot_product(R, R)<r*r && *dist > 0){ return(1); }else{ return(0);}
}
//-------------------------------------------------------------------------------------------------
// Determine whether P lies in triangle ABC using barycentric coordinates.
int inside_triangle(VEC P, VEC A, VEC B, VEC C)
{
VEC v0,v1,v2;
double d00,d01,d02,d11,d12,invDenom,u,v;

  v0 = linear2(1.0, C, -1.0, A); v1 = linear2(1.0, B, -1.0, A); v2 = linear2(1.0, P, -1.0, A);

  d00 = dot_product(v0, v0); d01 = dot_product(v0, v1); d02 = dot_product(v0, v2);
  d11 = dot_product(v1, v1); d12 = dot_product(v1, v2);

  invDenom = 1.0 / (d00 * d11 - d01 * d01);
  u = (d11 * d02 - d01 * d12) * invDenom;
  v = (d00 * d12 - d01 * d02) * invDenom;

  return (u >= -EPS_NORM) && (v >= -EPS_NORM) && (u + v <= 1.0 + EPS_NORM);
}
//-------------------------------------------------------------------------------------------------
// Compute the reflection of incident vector D about unit normal N_norm; |D| need not be 1.
VEC reflection(const VEC D, const VEC N_norm){
// The component (D,N_norm) N_norm parallel to the normal changes sign.
// The tangential component D - (D,N_norm) N_norm is unchanged.
// Thus D_reflected = D - 2 (D,N_norm) N_norm.
  double dcos = dot_product(D, N_norm);
  return( linear2( -2*dcos, N_norm, 1., D) );
}
//-------------------------------------------------------------------------------------------------
// Linear combinations of vectors: a A + b B (and extensions).
VEC linear1(double a, VEC A){
VEC W;
  W.x = a*A.x; W.y = a*A.y; W.z = a*A.z;
  return(W);
}
VEC linear2(double a, VEC A, double b, VEC B){
VEC W;
  W.x = a*A.x + b*B.x; W.y = a*A.y + b*B.y; W.z = a*A.z + b*B.z;
  return(W);
}
VEC linear3(double a, VEC A, double b, VEC B, double c, VEC C){
VEC W;
  W.x = a*A.x + b*B.x + c*C.x; W.y = a*A.y + b*B.y + c*C.y; W.z = a*A.z + b*B.z + c*C.z;
  return(W);
}
//-------------------------------------------------------------------------------------------------
// Dot product.
double dot_product(const VEC A, const VEC B){

  return( A.x * B.x + A.y * B.y + A.z * B.z );
}
//-------------------------------------------------------------------------------------------------
// Cross product.
VEC cross_product(const VEC A, const VEC B){
VEC C;

  C.x = A.y * B.z - A.z * B.y;
  C.y = A.z * B.x - A.x * B.z;
  C.z = A.x * B.y - A.y * B.x;

  return(C);
}
//-------------------------------------------------------------------------------------------------
// Normalize a vector.
VEC normalize(VEC A){
VEC W;
double n;

  W.x = W.y = W.z = 0.0 ;
  n = sqrt(dot_product(A, A));
  if(n>EPS_NORM){
    W.x = A.x/n; W.y = A.y/n; W.z = A.z/n;
  }else{
    fprintf(stderr,"Normalize error\n");
  }
  return(W);
}
//-------------------------------------------------------------------------------------------------
// Print a vector.
void printvec(char *s, VEC V){
   printf("%-5s %6.4f %6.4f %6.4f\n", s, V.x,  V.y,  V.z);
}
