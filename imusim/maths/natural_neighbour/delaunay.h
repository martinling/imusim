/******************************************************************************/
/*

  delaunay.h - By Ross Hemsley Aug. 2009 - rh7223@bris.ac.uk.
  
  This module will compute the delaunay triangulation of a set of uniformly
  distributed points in R^3. We will use the iterative edge flipping
  algorithm to add points one at a time.
  
  To store the triangulation, we start by just storing simplicies with pointers
  to their respective coordinates. 
  
  To simplify insertion, we first create a super-simplex with contains our
  entire dataset, this means that we don't have to specify any special
  conditions for the creation of our original simplicies.

  To make our algorithm robust we use Jonathan Shewchuk's Arbitrary Precision 
  Floating-poing Arithmetic predicates[1]  


  [1]  Routines for Arbitrary Precision Floating-point Arithmetic               
       and Fast Robust Geometric Predicates. May 18, 1996.
       Jonathan Richard Shewchuk.      
       
*/
/******************************************************************************/
#ifndef delaunay_h
#define delaunay_h
#include "utils.h"
/******************************************************************************/

/* These macros make code more readable. They allow us to access  
   the indexed elements of verticies directly.                    */
#define X v[0]
#define Y v[1]
#define Z v[2]
#define U data[0]
#define V data[1]
#define W data[2]

/******************************************************************************/

typedef struct
{
  // This is the location of this point.
  double v[3];
  
  // This is the point index in the point list.
  // We only store this so that it is convenient for using the
  // tet-mesh function in Matlab.
  // We can remove it for when the code is actually used.
  int    index;

  // These are the values for our vector field at this location.
  double data[3];
  
  // We use this for caching the voronoi volume of this point.
  // it will provide a good speed-up!
  double voronoiVolume;
  
} vertex;

/*******************************************************************************
* This is how we represent a Voronoi Cell in memory.                          
*******************************************************************************/
typedef struct 
{
  // The number of points on this cell, and the amount
  // of memory allocated for points on this cell.
  int n, nallocated;
  // The array of points on this cell.
  double **points;

  // This defines the cell, it contains a list of faces, each one is
  // consistantly oriented relative to itself (so traversing the pionts gives 
  // the convex hull of the face). Each face is seperated by a NULL pointer.
  // No gaurentee is made about the consistancy of orientations between 
  // different faces.
  arrayList *verticies;
  
} voronoiCell;

/******************************************************************************/
/* This is how we store an individual simplex: 4 pointers to the coordinates. */
/* We should try storing this without pointers probably.                      */
/******************************************************************************/
typedef struct _simplex
{
  // The verticies of this simplex.
  vertex  *p[4];
  // The neighbouring simlpicies of this simplex.
  // These are ordered in accordance with our 'get face' routine: 
  // so that the i'th face is shared with the i'th neighbour.
  struct _simplex *s[4];
  // This is the node in our auxillary list structure that holds this simplex.
  // It's far from an elegant solution: but it is fast and space-efficient.
  listNode *node;
} simplex;

/******************************************************************************/
/* We want to efficiently change back the neighbour pointers when             */
/* we remove a point.                                                         */
/******************************************************************************/
typedef struct
{
  stack  *ptrs;
  stack  *old;
} neighbourUpdate;

/******************************************************************************/
// We will keep all details of the mesh in this structure: thus hiding
// the complexities of memory pooling from the user. We also want to store
// the lists of most recently destroyed and allocated simplicies/neighbour
// updates, so that we speed up the 'remove last point' operation,
// which is _crucial_ to fast natural neighbour interpolation.
/******************************************************************************/

typedef struct
{
  // a linked list of all the simplicies.
  linkedList    *tets;
  
  // The simplex which contains all of the points.
  // its verticies contain no data values.
  simplex *super;
  vertex   superVerticies[4];
  
  // Memory pool.
  stack   *deadSimplicies;
  stack   *deadVoronoiCells;

  // We modify these when a point is inserted/removed.
  arrayList       *conflicts;
  arrayList       *updates;
  neighbourUpdate *neighbourUpdates;
  
  // Keep count of the number of degenerecies we find in the mesh, 
  // so that we can spot errors, and be aware of particularly degenerate data.
  int coplanar_degenerecies;
  int cospherical_degenerecies;

} mesh;

/******************************************************************************/
mesh*            newMesh();
//------------------------------------------------------------------------------
void             freeMesh(mesh *m);
//------------------------------------------------------------------------------
void             removePoint(mesh *m);
//------------------------------------------------------------------------------
vertex*          loadPoints(char *filename, int *n);
//------------------------------------------------------------------------------
void             getRange(vertex *ps, int n, vertex *min,
                                             vertex *max, vertex *range, int r);
//------------------------------------------------------------------------------
void             initSuperSimplex(vertex *ps, int n, mesh *m);
//------------------------------------------------------------------------------
void             writePointsToFile(vertex *ps, int n);
//------------------------------------------------------------------------------
void             writeTetsToFile(mesh *m);
//------------------------------------------------------------------------------
int              simplexContainsPoint(simplex *s, vertex *p);
//------------------------------------------------------------------------------
void             getFaceVerticies(simplex *s, int i, vertex **p1, vertex **p2, 
                                                     vertex **p3, vertex **p4 );
//------------------------------------------------------------------------------
int              vercmp(vertex *v1, vertex *v2);
//------------------------------------------------------------------------------
void             faceTest(mesh *m);
//------------------------------------------------------------------------------
void             orientationTest(linkedList *tets);
//------------------------------------------------------------------------------
void             allTests(linkedList *tets);
//------------------------------------------------------------------------------
void             addSimplexToMesh(mesh *m, simplex *s);
//------------------------------------------------------------------------------
void             removeSimplexFromMesh(mesh *m, simplex *s);
//------------------------------------------------------------------------------
simplex*         findContainingSimplex(mesh *m, vertex *p);
//------------------------------------------------------------------------------
int              isDelaunay(simplex *s, vertex *p);
//------------------------------------------------------------------------------
simplex**        swapSimplexNeighbour(simplex *s, simplex *old, simplex *new);
//------------------------------------------------------------------------------
simplex*         findNeighbour(simplex *s, vertex *p);
//------------------------------------------------------------------------------
void             setOrientationBits(simplex *s);
//------------------------------------------------------------------------------
void             buildMesh(vertex* ps, int n, mesh *m);
//------------------------------------------------------------------------------
void             addPointToMesh(vertex *p, linkedList *tets);
//------------------------------------------------------------------------------
int              pointOnSimplex(vertex *p, simplex *s);
//------------------------------------------------------------------------------
void             printEdge(vertex *v1, vertex* v2, FILE *stream);
//------------------------------------------------------------------------------
int              isConvex(vertex *v1, vertex *v2, vertex *v3, 
                                      vertex *t,  vertex *b);
//------------------------------------------------------------------------------
void             setNeighbourIndex(simplex *s, int i, int newIndex);
//------------------------------------------------------------------------------
int              getNeighbourIndex(simplex *s, int i);
//------------------------------------------------------------------------------
arrayList*       findNeighbours(vertex *v, simplex *s);
//------------------------------------------------------------------------------
simplex*         newSimplex(mesh *m);
//------------------------------------------------------------------------------
void             addPoint(vertex *p, mesh *m);
//------------------------------------------------------------------------------
int              delaunayTest(mesh *m, vertex *ps, int n);
//------------------------------------------------------------------------------
void             circumCenter(simplex *s, double *out);
//------------------------------------------------------------------------------
voronoiCell*     getVoronoiCell(vertex *point, simplex *s0, mesh *m);
//------------------------------------------------------------------------------
void             setNeighbours(arrayList *newTets);
//------------------------------------------------------------------------------
int              shareThreePoints(simplex *s0, int i, simplex *s1);
//------------------------------------------------------------------------------
void             vertexAdd(double *a, double *b, double *out);
//------------------------------------------------------------------------------
void             vertexByScalar(double *a, double b, double *out);
//------------------------------------------------------------------------------
void             vertexSub(double *a, double *b, double *out);
//------------------------------------------------------------------------------
void             crossProduct(double *b, double *c, double *out);
//------------------------------------------------------------------------------
double           squaredDistance(double *a);
//------------------------------------------------------------------------------
double           scalarProduct(double *a, double *b);
//------------------------------------------------------------------------------
double           volumeOfTetrahedron(double *a,double *b, double *c, double *d);
//------------------------------------------------------------------------------
double           voronoiCellVolume(voronoiCell *vc, vertex *p);
//------------------------------------------------------------------------------
void             removeExternalSimplicies(mesh *m);
//------------------------------------------------------------------------------
arrayList*       naturalNeighbours(vertex *v, mesh *m);
//------------------------------------------------------------------------------
void             writeVoronoiCellToFile(FILE* f, voronoiCell *vc);
//------------------------------------------------------------------------------
void             freeVoronoiCell(voronoiCell *vc, mesh *m);
//------------------------------------------------------------------------------
neighbourUpdate* initNeighbourUpdates();
//------------------------------------------------------------------------------
void             resetNeighbourUpdates(neighbourUpdate *nu);
//------------------------------------------------------------------------------
void             undoNeighbourUpdates(neighbourUpdate *nu);
//------------------------------------------------------------------------------
void             pushNeighbourUpdate(neighbourUpdate *nu, simplex **ptr,
                                                          simplex  *old);
//------------------------------------------------------------------------------
void             freeNeighbourUpdates(neighbourUpdate *nu);
//------------------------------------------------------------------------------
simplex*         findAnyNeighbour(vertex *v, arrayList *tets);
//------------------------------------------------------------------------------
int              getNumSimplicies(mesh *m);
//------------------------------------------------------------------------------
void             randomPerturbation(vertex *v, int attempt);
//------------------------------------------------------------------------------
int              numSphericalDegenerecies(mesh *m);
//------------------------------------------------------------------------------
int              numPlanarDegenerecies(mesh *m);
/******************************************************************************/
#endif

