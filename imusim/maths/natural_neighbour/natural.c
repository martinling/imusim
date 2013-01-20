/*******************************************************************************
*
*  interpolate.c - By Ross Hemsley Aug. 2009 - rh7223@bris.ac.uk.
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "delaunay.h"
#include "assert.h"
#include "utils.h"
#include "natural.h"

/******************************************************************************/
/* Set this to run tests, and include debugging information.                  */
// #define DEBUG 0

/* set this to disable all asserts.                                           */
// #define NDEBUG 0 

//#define _TEST_
/******************************************************************************/

#define SQR(x)  (x)*(x)
#define ABS(x)  (x) >0 ? x : -(x)
#define MAX(x,y)  x<y ? y : x
#define MIN(x,y)  x>y ? y : x

/******************************************************************************/

vertex *loadPoints(char *filename, int *n)
{
  int i;
  FILE *f = fopen(filename, "r");
  
  if (!f)
  {
    fprintf(stderr, "Could not open file: '%s'.\n", filename);
    exit(1);
  }
  
  // Get the number of points.
  fscanf(f, "%d", n);

  // Allocate enough memory for all of our points.
  // and also the interpolant.
  vertex *ps = malloc(sizeof(vertex) **n);

  for (i=0; i<*n; i++)
  {
    // Initialise the voronoiVolume to an impossible value, so that
    // we can tell whether or not it has been calculated.
    ps[i].voronoiVolume = -1;
    ps[i].index = i;
    // Load the known point of this vector field at this vertex.
    fscanf(f, "%lf %lf %lf   %lf %lf %lf", &ps[i].X, &ps[i].Y, &ps[i].Z, 
                                           &ps[i].U, &ps[i].V, &ps[i].W );
  }
  fclose(f);
  
  return ps;  
}

/******************************************************************************/

vertex *initPoints(double *x, double *y, double *z, 
                   double *u, double *v, double *w, int n)
{
  vertex* ps = malloc(sizeof(vertex) *n);

  int i;
  for (i=0; i<n; i++)
  {
    ps[i].X = x[i];
    ps[i].Y = y[i];
    ps[i].Z = z[i];
    
    ps[i].U = u[i];
    ps[i].V = v[i];
    ps[i].W = w[i]; 
    
    ps[i].voronoiVolume = -1;   
  }
  
  return ps;
}

/******************************************************************************/

void writePointsToFile(vertex *ps, int n)
{
  FILE *f = fopen("./points.mat", "wt");
  if (!f)
  {
    fprintf(stderr, "Could not open points file for writing.\n");
    exit(1);
  }
  
  int i;
  
  for (i=0; i<n; i++)
    fprintf(f, "%lf %lf %lf %lf %lf %lf\n", ps[i].X, ps[i].Y, ps[i].Z, 
                                            ps[i].U, ps[i].V, ps[i].W);  
  fclose(f);
}

/******************************************************************************/

void lastNaturalNeighbours(vertex *v, mesh *m, arrayList *neighbours, 
                                               arrayList *neighbourSimplicies)
{
  int i, j;
  for (i=0; i<arrayListSize(m->updates); i++)
  {
    simplex *this = getFromArrayList(m->updates,i); 
    for (j=0; j<4; j++)
    {     
      if (this->p[j] != v && (! arrayListContains(neighbours, this->p[j])) )
      {
        if ((! pointOnSimplex(this->p[j], m->super)))
        {
          addToArrayList(neighbours, this->p[j]);      
          addToArrayList(neighbourSimplicies, this);
        }
      }      
    }
  }
}

/******************************************************************************/

// This function will interpolate the value of a new vertex in a given 
// vector field.

void interpolate3_3( double  x, double  y, double  z, 
                     double *u, double *v, double *w, mesh *m )
{
  int i;
  
  // Set up a temporary vertex to add to this mesh.
  vertex p;
  p.X             =  x;
  p.Y             =  y;
  p.Z             =  z;
  p.index         = -1;
  p.voronoiVolume = -1;
  
  // The verticies which form the natural neighbours of this point.
  arrayList *neighbours;
  // The The list of neighbouring simplicies, attached to each one
  // of the given neighbours. This means makes neighbour lookup much faster
  // later on.
  arrayList *neighbourSimplicies;
  
  // The array of neighbour volumes for each voronoi cell of each
  // natural neighbour.
  double *neighbourVolumes;
  
  // The volume of the point to be interpolated.
  double pointVolume;
  
  // The interpolated value.
  double value[3] = {0,0,0};
  
  // The sum of the weighing function: may sometimes be less than 1, when
  // we have points on the super simplex.
  double sum, weight;
  
  // Add the point to the Delaunay Mesh - storing the original state.
  addPoint(&p, m);    

  // Find the natural neighbours of the inserted point, and also keep 
  // a list of an arbitrary neighbouring simplex, this will give us faster
  // neighbour lookup later.
  neighbours          = newArrayList();  
  neighbourSimplicies = newArrayList();  
  lastNaturalNeighbours(&p, m, neighbours, neighbourSimplicies);

  // Calculate the volumes of the Voronoi Cells of the natural neighbours.
  neighbourVolumes = malloc(arrayListSize(neighbours) * sizeof(double));

  // Calculate the 'before' volumes of each voronoi cell.
  for (i=0; i<arrayListSize(neighbours); i++)
  {
    vertex  *thisVertex  = getFromArrayList(neighbours, i);
    simplex *thisSimplex = getFromArrayList(neighbourSimplicies,i);  
    voronoiCell *vc      = getVoronoiCell(thisVertex, thisSimplex, m);    
    neighbourVolumes[i]  = voronoiCellVolume(vc, thisVertex);  
    freeVoronoiCell(vc,m); 
  }

  // Calculate the volume of the new point's Voronoi Cell.
  // We just need any neighbour simplex to use as an entry point into the
  // mesh.
  simplex *s             = getFromArrayList(neighbourSimplicies,0);
  voronoiCell *pointCell = getVoronoiCell(&p, s, m);
  pointVolume            = voronoiCellVolume(pointCell, &p);
  freeVoronoiCell(pointCell,m);
         
  // Remove the last point.
  removePoint(m);

  // Calculate the 'stolen' volume of each neighbouring Voronoi Cell,
  // by calculating the original volumes, and subtracting the volumes
  // given when the point was added.
  for (i=0; i<arrayListSize(neighbours); i++)
  {
    vertex *thisVertex   = getFromArrayList(neighbours, i);  
    
    // All verticies have -1 here to start with, so we can tell if 
    // we have already calculated this value, and use it again here.
    if (thisVertex->voronoiVolume < 0)
    {
      simplex *s           = findAnyNeighbour(thisVertex, m->conflicts);
      voronoiCell *vc      = getVoronoiCell(thisVertex, s, m);
      thisVertex->voronoiVolume = voronoiCellVolume(vc, thisVertex);
      freeVoronoiCell(vc,m);
    }
    neighbourVolumes[i]  = thisVertex->voronoiVolume-neighbourVolumes[i];
  }
   
  // Weight the data values of each natural neighbour using the volume
  // ratios.
  sum   = 0;

  for (i=0; i<arrayListSize(neighbours); i++)
  {
    vertex *thisVertex = getFromArrayList(neighbours, i);
    assert (neighbourVolumes[i]>= -0.001);
    
    // Get the weight of this vertex.
    weight = neighbourVolumes[i]/pointVolume;
    
    // Add this componenet to the result.
    sum      += weight;
    value[0] += weight * thisVertex->U;   
    value[1] += weight * thisVertex->V;   
    value[2] += weight * thisVertex->W;           
  }
  
  // Normalise the output.
  vertexByScalar(value, (double)1/(double)sum, value);

  // If the sum is 0 or less, we will get meaningless output. 
  // If it is slightly greater than 1, this could be due to rounding errors.
  // We tolerate up to 0.1 here.  
  if (sum <= 0 || sum > 1.1)
  {
    fprintf(stderr, "Error: sum value: %lf, expected range (0,1].\n",sum);
    fprintf(stderr, "There could be a degenerecy in the mesh, either retry "
                    "(since input is randomised this may resolve the problem), "
                    "or try adding a random peterbation to every point.\n");
    exit(1);
  }

  // Put the dead simplicies in the memory pool.
  for (i=0; i<arrayListSize(m->updates); i++)
    push(m->deadSimplicies, getFromArrayList(m->updates, i));

  // Free all the memory that we allocated whilst interpolating this point.
  emptyArrayList(m->conflicts);
  emptyArrayList(m->updates);
  
  // Free memory associated with adding this point.
  freeArrayList(neighbours,          NULL);
  freeArrayList(neighbourSimplicies, NULL); 
  free(neighbourVolumes);
  
  // set the output.
  *u = value[0];
  *v = value[1];
  *w = value[2];

}

/******************************************************************************/

/* Unit testing. */
#ifdef _TEST_

  #include <sys/time.h>

  /* The number of points to create in our test data. */
  #define NUM_TEST_POINTS 1e4
  
  /* How detailed the interpolated output sohuld be: the cube of this value
     is the number of points we create.                                     */
  #define INTERPOLATE_DETAIL 20
  
  /* Do we print the output to file? */
  #define OUTPUT_TO_FILE

/******************************************************************************/

double getTime()
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec + tv.tv_usec/1.0e6;
}

/******************************************************************************/

int main(int argc, char **argv)
{  
  int i; 
  int n = NUM_TEST_POINTS;
  srand ( time(NULL) );
  
  // Create a random pointset for testing.
  vertex *ps = malloc(sizeof(vertex)*NUM_TEST_POINTS);

  for (i=0; i<n; i++)
  {
    ps[i].X = 100*(double)rand() / ((double)RAND_MAX + 1);
    ps[i].Y = 100*(double)rand() / ((double)RAND_MAX + 1);
    ps[i].Z = 100*(double)rand() / ((double)RAND_MAX + 1);

    // We chose the value at every point of the vector field, to be the same
    // as the coordinate, this will enable us to have a good idea of
    // the error margins.
    ps[i].U =  ps[i].X;
    ps[i].V =  ps[i].Y;
    ps[i].W =  ps[i].Z;
    ps[i].index = i;
    ps[i].voronoiVolume = -1;
  }

  mesh *m = newMesh();
  buildMesh(ps, n, m);

  // Display some information about the mesh.
  printf("Number of Verticies: %d.\n", n);
  printf("Number of Simplicies: %d.\n", getNumSimplicies(m));
  printf("Co-planar degenerecies fixed: %d.\n", numPlanarDegenerecies(m));
  printf("Co-spherical degenerecies fixed: %d.\n", numSphericalDegenerecies(m));

  // Write output to files for plotting.
  #ifdef OUTPUT_TO_FILE
  writeTetsToFile(m);  
  writePointsToFile(ps, n);
  #endif

  vertex min, max, range;
  getRange(ps, n, &min, &max, &range, 0);

  // We will store the component-wise sum over all errors, and max error
  // so that we can get an idea for the performance of our interpolator.
  double error_sum = 0;

  double x, y, z, d1, d2, d3;
  int detail = INTERPOLATE_DETAIL;  

  d1 = range.X / detail;
  d2 = range.Y / detail;
  d3 = range.Z / detail;
  
  // This will tell us how many steps are in our main loop:
  // so that we can show an indication as to how far we have gone through
  // the calculation.
  int to_do = detail*detail*detail;
  
  #ifdef OUTPUT_TO_FILE
  // Interpolate a set of points.
  FILE *f = fopen("./out.mat","wt");
  if (!f)
  {
    fprintf(stderr, "Could not open point file for writing.\n");
    exit(1);
  }
  #endif

  int j,k,done=1;
  double t1 = getTime();
  for (i=0; i<detail; i++)
  {
    for (j=0; j<detail; j++)
    {
      for (k=0; k<detail; k++, done++)
      {
        double u, v, w;
        
        x = min.X + i*d1;
        y = min.Y + j*d2;
        z = min.Z + k*d3;
        
        // Do the interpolation.        
        interpolate3_3(x, y, z, &u, &v, &w, m);
        
        // Add this error to the sum.
        error_sum += SQR(x-u) + SQR(y-v) + SQR(z-w);

        // Show status.
        #ifdef OUTPUT_TO_FILE
        fprintf(f,"%lf %lf %lf %lf %lf %lf\n", x, y, z, u, v, w);  
        #endif
        printf("Interpolating: %d/%d (%d%% complete).\n%c[1A", done,to_do,
                                         (int)(done/(double)to_do *100), 27);   
      }
    }
  }
  double t2 = getTime();
  
  #ifdef OUTPUT_TO_FILE
  fclose(f);
  #endif
  
  printf("\n\n");

  int pps = (double)to_do/(double)(t2-t1);
  
  printf("Interpolated %d points in %lf seconds (%d points/second).\n\n", 
                                         to_do, t2-t1,pps);
    
  printf("Sum of Squared errors over vector norms: %lf.\n", error_sum);

  #if DEBUG >= 0
  printf("\nSimplex pool size at end: %d\n", stackSize(m->deadSimplicies));
  printf("Voronoi pool size at end: %d\n", stackSize(m->deadVoronoiCells));
  #endif
  
  // Free up used memory.
  freeMesh(m);
  free(ps);
 
  return 0;
}

/******************************************************************************/
#endif

