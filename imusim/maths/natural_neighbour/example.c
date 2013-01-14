#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "delaunay.h"
#include "natural.h"


int main(int argc, char **argv)
{
  // the number of points in our point set.
  int n = 3;

  // the points.
  double x[] = {0,1,2};
  double y[] = {3,4,5};
  double z[] = {6,7,8};
  
  double u[] = {0,1,2};
  double v[] = {3,4,5};
  double w[] = {6,7,8};
  
  // initialise the points.
  vertex* ps = initPoints(x,y,z,  u,v,w,  n);
  
  // create and build the mesh.
  mesh *exampleMesh = newMesh(); 
  buildMesh(ps ,n, exampleMesh);
  
  // perform interpolation.
  double u_out, v_out, w_out;  
  interpolate3_3(1,4,7, &u_out, &v_out, &w_out, exampleMesh);
  
  printf("Interpolated value at (1,2,3) is (%lf, %lf %lf).\n", 
                                                          u_out, v_out, w_out);
  
  // free the mesh.  
  freeMesh(exampleMesh);
  free(ps);

}
