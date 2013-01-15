/*******************************************************************************
*
*        delaunay.c - By Ross Hemsley Aug. 2009 - rh7223@bris.ac.uk.
*
* This file implements Delaunay meshing in 3D, using the edge flipping
* algorithm. To stop degenerecies arising from floating point errors, we use
* the geometical predicates provided in predicates.c - giving adaptive 
* floating point arithmetic. We also remove degenerecies present in data 
* caused by points which are coplanar, or cospherical. These points are removed
* by gradually adding random peterbations until the degenerecies are removed.
*
* This file has unit testing, which can be done by defining _TEST_ as shown
* seen below. The file can then be compiled by running:
*
*   >gcc -O3 delaunay.c utils.c
*
* The executible created can be run to create a set of random points, which are
* then meshed and checked for Delaunayness.
* 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "./utils.h"
#include "delaunay.h"
#include "predicates.c"
#include "assert.h"
#include <time.h>

/******************************************************************************/

/* Set this to be lower than the average distance between points. It is the
   amount that we will shift points by when we detect degenerecies. We
   gradually increase the value until the degenercy is removed                */
  #define PERTURBATION_VALUE  0.0001

/* Do we show status of meshing? */
  #define VERBOSE

/* This dictates the amount of debugging information to print. Starting
   at level 0. If not defined, no debugging information is printed.           */
//  #define DEBUG 0
                                                                            
/* This allows us to turn off all error checking. */
//#define NDEBUG

/* Turn on Unit testing. */
// #define _TEST_

/******************************************************************************/
/* - DO NOT EDIT - */
#ifdef _TEST_
  #undef NDEBUG
#endif

#ifdef DEBUG
int SIMPLEX_MALLOC = 0;
int VORONOI_MALLOC = 0;
int VERTEX_MALLOC  = 0;
int COPLANAR_DEGENERECIES = 0;
int COSPHERICAL_DEGENERECIES = 0;
#endif
/******************************************************************************/
#define MAX(x,y)  x<y ? y : x
#define MIN(x,y)  x>y ? y : x
#define SWAP(x,y)                                                              \
{                                                                              \
  double tmp;                                                                  \
  tmp = x;                                                                     \
  x   = y;                                                                     \
  y   = tmp;                                                                   \
}
/******************************************************************************/

simplex *newSimplex(mesh *m)
{
  simplex *s = pop(m->deadSimplicies);
 
  // Obviously, we aren't going to re-use the super simplex..
  if (s==m->super) s=0;
 
  if (!s)
  {
    s = malloc(sizeof(simplex));
    #ifdef DEBUG
    SIMPLEX_MALLOC ++;
    #endif
  }
  s->s[0] = 0;
  s->s[1] = 0;
  s->s[2] = 0;
  s->s[3] = 0;
  
  return s;
}

/******************************************************************************/
// This will take a list of points, and a mesh struct, and create a 
// Delaunay Tetrahedralisation.

void buildMesh(vertex* ps, int n, mesh *m)
{
  // Seed the random function, we will use the random function to remove
  // any degenerecies as we find them.
  srand ( time(NULL) );
    
  // We have no degenerecies to start with.
  m->coplanar_degenerecies  = 0;
  m->cospherical_degenerecies = 0;

  // This simplex will contain our entire point-set.
  initSuperSimplex(ps, n, m);
  addSimplexToMesh(m, m->super);
  
  int i,j; 
  // Add each point to the mesh 1-by-1 using the Edge Flipping technique.
  for (i=0; i<n; i++)
  {
    addPoint(&ps[i], m);
    
    // Push conflicts to the memory pool.
    for (j=0; j<arrayListSize(m->conflicts); j++)
      push(m->deadSimplicies, getFromArrayList(m->conflicts, j));    
    
    // Reset the conflict and update lists.
    emptyArrayList(m->conflicts);
    emptyArrayList(m->updates);
    
    // Clear out the old neighobur update structs. (We don't use them here).
    resetNeighbourUpdates(m->neighbourUpdates);
    
    #ifdef VERBOSE
    // Show status of meshing.
    printf("Meshing: %d%%.\n%c[1A", (int)((i+1)/(double)n *100),27);   
    #endif 
  }
}

/******************************************************************************/
// This will allow us to remove all the simplicies which are connected
// to the super simplex.
void removeExternalSimplicies(mesh *m)
{
  listNode *iter = topOfLinkedList(m->tets);
  simplex *s;
 
  // Remove all simplicies which connect to the super simplex
  while ((s = nextElement(m->tets, &iter)))
  {
    if (  simplexContainsPoint(s, m->super->p[0]) || 
          simplexContainsPoint(s, m->super->p[1]) || 
          simplexContainsPoint(s, m->super->p[2]) ||
          simplexContainsPoint(s, m->super->p[3])     )      
    {     
      swapSimplexNeighbour(s->s[0], s, NULL);
      swapSimplexNeighbour(s->s[1], s, NULL);
      swapSimplexNeighbour(s->s[2], s, NULL);
      swapSimplexNeighbour(s->s[3], s, NULL);
      
      removeSimplexFromMesh(m, s);
    }
  }
}

/******************************************************************************/
// return the value that we modified.

simplex** swapSimplexNeighbour(simplex *s, simplex *old, simplex *new)
{
  // If this neighbour is on the exterior, we don't need to do anything.
  if (!s) return NULL;
 
  int i,found=0;
  
  // We are going to go through each of the elements children to see which one
  // points to the old simplex. When we find that value, we are going to swap 
  // it for the new simplex value.
  for (i=0;i<4;i++)
  {
    if (s->s[i] == old) 
    {
      found=1;
      break;    
    }
  }

  s->s[i] = new;

  assert(found);
  return &s->s[i];
}

/******************************************************************************/
// we are going to go through every face of every simplex to see if the
// orientation is consistent.

void orientationTest(linkedList *tets)
{
  #if DEBUG >= 1
  printf("Running orientation test: ---------------------------------------\n");
  #endif
  
  int i;  
  listNode *iter = topOfLinkedList(tets);
  simplex  *s;
  
  while((s = nextElement(tets, &iter)))
  {
    vertex *p1, *p2, *p3, *p4;
    
    #if DEBUG >= 1
    printf("Checking orientation of %p\n", s);
    #endif
    
    // Go through every face of this simplex
    for (i=0;i<4;i++)
    {
      getFaceVerticies(s, i, &p1, &p2, &p3, &p4);
      double o =  orient3dfast(p1->v, p2->v, p3->v, p4->v);
      assert (o>0);
    }
  }
}

/******************************************************************************/

int delaunayTest(mesh *m, vertex *ps, int n)
{
 
 listNode *iter = topOfLinkedList(m->tets);
 simplex  *s;
 
 int isDel=0;
 int notDel=0;
  
 while ((s = nextElement(m->tets, &iter)))
 {
    #if DEBUG >= 2
    printf("Checking Delaunayness of %p.\n",s);
    #endif    
    
    // we want to see if this simplex is delaunay  
    int i, succes=1;
    for (i=0; i<n; i++)
    {
      // if this point is not on the simplex, then it should not be within 
      // the circumsphere of this given simplex.
      if (! pointOnSimplex(&ps[i], s))
      {
        #if DEBUG >= 0
        double orientation = orient3dfast(s->p[0]->v, 
                                          s->p[1]->v, 
                                          s->p[2]->v, 
                                          s->p[3]->v);

        assert(orientation != 0); 
        assert(orientation  > 0);
        #endif
        
        double inSph = inspherefast(s->p[0]->v,
                                    s->p[1]->v, 
                                    s->p[2]->v, 
                                    s->p[3]->v,
                                    ps[i].v);
        if (inSph >= 0)
        {
          notDel++;
          succes = 0;
          break;
        }
      }
    
    } 
    if (succes) isDel ++;  
  }
 
  #if DEBUG >=2
  printf("Non-Delaunay Simplicies: %d.\n", notDel);
  printf("There are %f%% non-Delaunay Simplicies.\n", 
                                         (notDel/(double)(isDel + notDel))*100);
 #endif
  return notDel == 0;                                        
}

/******************************************************************************/
// This function is purely to test whether the set of neighbours of each 
// simplex is correct - If it is not reliable, then the program behaviour will
// be undeterministic: potentially giving a very difficult bug. 
// We only need to run this test when the code is modified.

void faceTest(mesh *m)
{
  int j;
  
  #if DEBUG >= 1
  printf("Running face test: ----------------------------------------------\n");
  #endif
  
  // Set our iterator to point to the top of the tet list.
  listNode *iter = topOfLinkedList(m->tets);
  // The pointre to the current simplex that we are considering.
  simplex  *s;
  
  // Go through every simplex in the list.
  while ((s = nextElement(m->tets, &iter)))
  {
    #if DEBUG >= 1
    printf("# Checking simplex: %p.\n", s);
    #endif
    
    // Go through each neighbour of the simplex (this is equivilent
    // to going through every face of the simplex).
    for (j=0;j<4;j++)
    {
      vertex  *p1, *p2, *p3, *p4, *t1, *t2, *t3, *t4;; 
      
      // Get the verticies of the face we are considering.      
      getFaceVerticies(s, j, &p1, &p2, &p3, &p4);
            
      // This is the neighbour that should share the given verticies.
      simplex *neighbour = s->s[j];    
             
      #if DEBUG >=1
      printf("  Neighbour: s->[%d]: %p\n", j, neighbour);
      #endif
      
      // This could be an outer-face: in which case, there is no neighbour here.
      if (neighbour != NULL)
      {          
        int x,found=0;
        //assert(!s->s[j]->dead);
        
        // Go through each neighbour and see if it points to us. 
        // if it does (which it should) check the points match.
        for (x=0; x<4; x++)
        {
          #if DEBUG >=1
          printf("  Checking: %p\n", neighbour->s[x]);
          #endif
          if (neighbour && neighbour->s[x] && neighbour->s[x] == s)
          {
            found = 1;
            
            // Get the verticies of the face that we share with the current
            // simplex.
            getFaceVerticies(neighbour, x, &t1, &t2, &t3, &t4);
            
            // We want to check that these two simplicies share their first 
            // three verticies. 
            getFaceVerticies(neighbour, x, &t1, &t2, &t3, &t4);
             
            assert (vercmp(t1,p1) || vercmp(t2, p1) || vercmp(t3,p1));
            assert (vercmp(t1,p2) || vercmp(t2, p2) || vercmp(t3,p2));
            assert (vercmp(t1,p3) || vercmp(t2, p3) || vercmp(t3,p3));                          
          }
        }
        // We have a pointer to a neighbour which does not point back to us.
        assert(found);
      }
    }
  }
}

/******************************************************************************/

int vercmp(vertex *v1, vertex *v2)
{
  int i;
  for (i=0; i<3; i++)
    if ( v1->v[i] != v2->v[i] ) return 0; 
  return 1;
}

/******************************************************************************/
// This is a slightly optimised method to find the containing simplex
// of a point. We go through each simplex, check to see which faces, if any
// face the point we are looking for. The first one we find that does, we
// follow that neighbour. If all the faces are oriented so that the point is
// not in front of them, then we know that we have found the containing simplex.
// It is likely to be provably O(n^1/2).


simplex* findContainingSimplex(mesh *m, vertex *p)
{
  // This will arbitrarily get the first simplex to consider.
  // ideally we want to start from the middle, but chosing a random 
  // simplex will give us good performance in general.
  
  listNode *iter = topOfLinkedList(m->tets);
  simplex  *s    = nextElement(m->tets,&iter); 
  vertex *v1, *v2, *v3, *v4;
  
  int i;
  for (i=0; i<4; i++)
  {
    // get the orientation of this face.
    getFaceVerticies(s, i, &v1, &v2, &v3, &v4);
    
    if ((orient3dfast(v1->v, v2->v, v3->v, p->v) < 0) && s->s[i])
    {
      // Go to the next simplex, and start the loop again.
      s = s->s[i];
      i = -1;
    }
  }
    
  // All the orientation tests passed: the point lies within/on the simplex.
  return s;
}

/******************************************************************************/
// Return, as 3 arrays of double, the verticies of the face i of this simplex.
// This function aims to help us ensure consistant orientation.
// The last value is that of the remaining vertex which is left over.

void getFaceVerticies(simplex *s, int i, vertex **p1, vertex **p2, 
                                         vertex **p3, vertex **p4  )
{
  switch (i)
  {
    case 0:
      *p1 = s->p[0];
      *p2 = s->p[1];
      *p3 = s->p[2];      
      *p4 = s->p[3];
      break;
    case 1:
      *p1 = s->p[3];
      *p2 = s->p[1];
      *p3 = s->p[0];      
      *p4 = s->p[2];   
      break;
    case 2:
      *p1 = s->p[0];
      *p2 = s->p[2];
      *p3 = s->p[3];      
      *p4 = s->p[1];  
      break;
    case 3:  
      *p1 = s->p[3];
      *p2 = s->p[2];
      *p3 = s->p[1];      
      *p4 = s->p[0];  
      break;
  } 
}

/******************************************************************************/
// This routine will tell us whether or not a simplex contains a given point.
// To perform this test robustly, we will use the code provided by
// Jonathan Richard Shewchuk[3]. This code allows us to scale the precision 
// of our calculations so that we can be sure of valid results whilst retaining
// good performance in the general case.

int simplexContainsPoint(simplex *s, vertex *p)
{
  // To perform this test, we check the orientation of our point against 
  // the plane defined by each triangular face of our given simplex.
  // if the sign is always negative then the point lies within the simplex.
  
  int i;
  
  // The points on this face.
  vertex *p1, *p2, *p3, *p4;  
  
  for (i=0; i<4; i++)
  {
    // Get the face values for this simplex.
    getFaceVerticies(s, i, &p1, &p2, &p3, &p4);
    if (orient3dfast(p1->v, p2->v, p3->v, p->v) < 0) return 0;
  }
  
  return 1;
}

/******************************************************************************/
// Write out all the tets in the list, except for those ones connected to 
// the points on S0: which we can use as the super simplex.

void writeTetsToFile(mesh *m)
{
  FILE *f = fopen("./tets.mat", "wt");
  if (!f)
  {
    fprintf(stderr, "Could not open tet. file for writing.\n");
    exit(1);
  } 
  
  simplex  *s;
  listNode *iter = topOfLinkedList(m->tets);
  
 int i;
  while ((s = nextElement(m->tets, &iter)))
  {
    int super =0;
    for (i=0; i<4; i++)
      if (pointOnSimplex(s->p[i],m->super)) super =1;

    if (!super) 
      fprintf(f,"%d %d %d %d\n", s->p[0]->index, s->p[1]->index, 
                                 s->p[2]->index, s->p[3]->index);
  }
  fclose(f);  
}

/******************************************************************************/
// Add gradually larger random perturbations to this point, until we can
// get a sphere which is not degenerate.

void randomPerturbation(vertex *v, int attempt)
{
  int i;
  for (i=0;i<3;i++)
  {
    // Get a [0,1] distributed random variable.
    double rand01 = (double)rand()/((double)RAND_MAX + 1);
    // add a random perturbation to each component of this vertex.
    double p = (rand01-0.5) * PERTURBATION_VALUE * (attempt+1);
    v->v[i] +=  p;
  }
}

/******************************************************************************/
// This routine will return 0 if the simplex is no longer Delaunay with 
// the addition of this new point, 1 if this simplex is still Delaunay
// with the addition of this new point, and -1 if this simplex is 
// degenerate with the addition of this new point (i.e. if the simplex is
// co-spherical.)

int isDelaunay(simplex *s, vertex *p)
{ 
  // If the orientation is incorrect, then the output will be indeterministic.
 // #if DEBUG >= 0
  double orientation = orient3dfast(s->p[0]->v, 
                                    s->p[1]->v, 
                                    s->p[2]->v, 
                                    s->p[3]->v);


  if (orientation <= 0)
  {
    printf("orientation error: %p, %lf\n",s,orientation);

    exit(1);
  }
//  assert(orientation != 0);
//  assert(orientation >  0);

  //#endif
  double inSph = inspherefast(  s->p[0]->v,
                                s->p[1]->v, 
                                s->p[2]->v, 
                                s->p[3]->v, p); 

            
  // We have a degenerecy.
  if (inSph == 0) return -1;
                
  return inSph < 0;

}

/******************************************************************************/
// We assume that the list is correct on starting (i.e. contains no
// non-conflicting simplicies).

void updateConflictingSimplicies(vertex *p, mesh *m)
{
  int i;
  // Get at least one simplex which contains this point.
  simplex *s0 = findContainingSimplex(m, p);
  simplex *current;
  
  // Go through each simplex, if it contains neighbours which are
  // not already present, which are not in the list already, 
  // and which are not delaunay, we add them to the list of conflicts
  stack *toCheck = newStack();
  push(toCheck, s0);

  while (!isEmpty(toCheck))
  {
    // pop the next one to check from the stack.
    current = pop(toCheck);
    
    int isDel = isDelaunay(current,p); 
    
    // Check to see whether or not we have a degenerecy
    if (isDel == -1) 
    {     
      m->cospherical_degenerecies ++;
      int i=0;
      while( isDel == -1 )
      {
        randomPerturbation(p,i);
        isDel = isDelaunay(current,p);
        //printf("Degenerecy removing for %p, attempt: %d\n",current,i);
        i++;
      }   
      
      // Start this function again now that we have moved the point.
      freeStack(toCheck,NULL);
      emptyArrayList(m->conflicts);
      updateConflictingSimplicies(p,m);
      return;
    }
    
    if ((!isDel) && (!arrayListContains(m->conflicts, current)))
    {
      // add this simplex, and check its neighbours.
      addToArrayList(m->conflicts, current);
      for (i=0; i<4;i++)
        if (current->s[i])
          push(toCheck, current->s[i]);
    }
  }
  freeStack(toCheck,NULL);
}

/******************************************************************************/
// Add a point by using the edge flipping algorithm.

void addPoint(vertex *p, mesh *m)
{
  // If the list arguments are NULL, then we create local lists.
  // Otherwise, we return the list of updates we did. 
  // This is so that we can easily perform point removal.

  // This will set a list of conflicting non-Delaunay simplicies in the mesh
  // structure.
  updateConflictingSimplicies(p,m);
  
  // We now have a list of simplicies which contain the point p within
  // their circum-sphere.
  // We now want to create a new tetrahedralisation in the polytope formed
  // by removing the old simplicies.
  // We know which faces we should connect to our point, by deleting every
  // face which is shared by another conflicting simplex.
  
  int i,j;
  for (j=0; j< arrayListSize(m->conflicts); j++)
  {
     simplex *s = getFromArrayList(m->conflicts,j);
     
    // Now go through each face, if it is not shared by any other face
    // on the stack, we will create a new simplex which is joined to 
    // our point.
    for (i=0; i<4; i++)
    {
      vertex *v1, *v2, *v3, *v4;
      getFaceVerticies(s, i, &v1, &v2, &v3, &v4);
      
      // Now, check to see whether or not this face is shared with any 
      // other simplicies in the list.
      if (! arrayListContains(m->conflicts, s->s[i]))
      {
        // We will create a new simplex connecting this face to our point. 
        simplex *new = newSimplex(m);
        new->p[0] = v1;
        new->p[1] = v2;
        new->p[2] = v3;
        new->p[3] =  p;
        
        int attempt = 0;
        // Detect degenerecies resulting from coplanar points.
        double o = orient3dfast(v1->v, v2->v, v3->v, p->v);
        if (o<=0)
        {
          m->coplanar_degenerecies ++;
          while (o<=0)
          {
            randomPerturbation(p, attempt);
            o = orient3dfast(v1->v, v2->v, v3->v, p->v);
            attempt ++;
          }
          // We are going to have to start adding this point again.
          // That means removing all changes we have done so far.
          undoNeighbourUpdates(m->neighbourUpdates);
          int k;
          for (k=0; k<arrayListSize(m->updates); k++)
          {
            removeSimplexFromMesh(m, getFromArrayList(m->updates,k));
            push(m->deadSimplicies, getFromArrayList(m->updates, k));    
          }
          emptyArrayList(m->updates);
          emptyArrayList(m->conflicts);
          // Start adding this point again, now that we have
          // (hopefully) removed the coplanar dependencies.
          addPoint(p,m);
          return;
        }
        
        // We know that every successful face will be pointing
        // outwards from the point. We can therefore directly set the neighbour
        // to be the same as the one that was with this face before.
        new->s[0] = s->s[i];
         
        // update, storing each neighbour pointer change we make.
        simplex **update = swapSimplexNeighbour(s->s[i], s, new);
        pushNeighbourUpdate(m->neighbourUpdates, update, s);

        // This is a list of all the new tets created whilst adding
        // this point.
        addToArrayList(m->updates, new);
        addSimplexToMesh(m, new);        
      }      
    }    
  }

  // Connect up the internal neighbours of all our new simplicies.
  setNeighbours(m->updates);  
  
  // Remove the conflicting simplicies.  
  for (i=0; i<arrayListSize(m->conflicts); i++)
  {
    simplex *s = getFromArrayList(m->conflicts, i);
    removeSimplexFromMesh(m,s);
  }
}

/******************************************************************************/
// Slightly quick and dirty way to connect up all the neighbours of the 
// new simplicies.

void setNeighbours(arrayList *newTets)
{
  simplex *s, *s2;
  vertex *v1, *v2, *v3, *t1, *t2, *t3, *tmp;
  // Go through each new simplex.
  int j;
  
  for (j=0; j<arrayListSize(newTets); j++)
  {
    s = getFromArrayList(newTets,j);

    // These are the verticies on the 2-simplex pointing outwards from 
    // the current point.   
    v1 = s->p[0];
    v2 = s->p[1];
    v3 = s->p[2];

    // We need to find neighbours for the edges (v1,v2) (v2,v3) (v3,v1)
    // We will do this by going through every other simplex in the list, 
    // and checking to see if its outward pointing face shares any of these
    // pairs. If it does, then we connect up the neighbours.
    
    int k;
    for (k=0; k<arrayListSize(newTets); k++)
    {
      s2 = getFromArrayList(newTets,k); 
      if (s == s2) continue;
      int i;
      // NOTE: we don't consider the outside face.
      // We want to know which side the neighbours are on.
      for (i=1; i<4; i++)
      {
        getFaceVerticies(s2, i, &t1, &t2, &t3, &tmp);
        // We now want to see if any of the edges (v1,v2) (v2,v3) (v3,v1) are 
        // on this triangle:
        if (      (v1 == t1 || v1 == t2 || v1 == t3) && 
                  (v2 == t1 || v2 == t2 || v2 == t3)    )
          s->s[1] = s2;        
        else if ( (v2 == t1 || v2 == t2 || v2 == t3) && 
                  (v3 == t1 || v3 == t2 || v3 == t3)    )
          s->s[3] = s2;
        else if ( (v3 == t1 || v3 == t2 || v3 == t3) && 
                  (v1 == t1 || v1 == t2 || v1 == t3)    )
          s->s[2] = s2;           
      }
    }
  }
}

/******************************************************************************/

int shareThreePoints(simplex *s0, int i, simplex *s1)
{
  vertex *v1, *v2, *v3, *v4;
  
  getFaceVerticies(s0, i, &v1, &v2, &v3, &v4);

  return (pointOnSimplex(v1,s1) && pointOnSimplex(v2,s1) &&
          pointOnSimplex(v3,s1) );
}

/******************************************************************************/
// Print an edge of a simplex to an output stream.

void printEdge(vertex *v1, vertex* v2, FILE *stream)
{
  fprintf(stream, "%lf %lf %lf   %lf %lf %lf\n",
                                   v1->X, v2->X, v1->Y, v2->Y, v1->Z, v2->Z);
}

/******************************************************************************/
// Does this simplex have the point p? - notice that we are comparing pointersa
// and not coordinates: so that duplicate coordinates will evaluate to not 
// equal.

int pointOnSimplex(vertex *p, simplex *s)
{
  if (!s) return 0;
  
  if (p == s->p[0] || p == s->p[1] || p == s->p[2] || p == s->p[3])   
    return 1;

  return 0;
}

/******************************************************************************/
// This routine tell us the neighbour of a simplex which is _not_ connected
// to the given point. 

simplex *findNeighbour(simplex *s, vertex *p)
{
  vertex *t1, *t2, *t3, *t4;
  int i,found=0;
  for (i=0; i<4; i++)
  {
    getFaceVerticies(s, i, &t1, &t2, &t3, &t4);
    if (t4 == p) 
    {
      found = 1;
      break;
    }
  }
  
  // If this fails then we couldn't find this point on this simplex.
  assert(found);
  
  return s->s[i];
}

/******************************************************************************/
// Check to see if the two simplicies sharing face v1, v2, v3 with top point
// t, and bottom point b are convex: i.e. can we draw a line between t and b
// which passes through the 2-simplex v1,v2,v3.

int isConvex(vertex *v1, vertex *v2, vertex *v3, vertex *t, vertex *b)
{
  int i=0;
  if (orient3dfast(v3->v, t->v, v1->v, b->v) < 0) i++; 
  if (orient3dfast(v1->v, t->v, v2->v, b->v) < 0) i++;  
  if (orient3dfast(v2->v, t->v, v3->v, b->v) < 0) i++;
  
  return (i==0);
}

/******************************************************************************/
// This will return an arrayList of verticies which are the Natural
// Neighbours of a point. This is currently not used, and is slow.

arrayList *naturalNeighbours(vertex *v, mesh *m)
{
  simplex *s;
  // User is responsible for freeing this structure.
  arrayList *l = newArrayList();
  listNode  *iter = topOfLinkedList(m->tets);

  int i;
  while ((s = nextElement(m->tets, &iter)))
    if (pointOnSimplex(v, s))
      for (i=0;i<4;i++)
        if ( (s->p[i] != v) && (! pointOnSimplex(s->p[i], m->super) )
                            && (! arrayListContains(l, s->p[i])) )
          addToArrayList(l, s->p[i]);

  return l;
}

/******************************************************************************/
// Given a point and a list of simplicies, we want to find any valid 
// neighbour of this point.

simplex *findAnyNeighbour(vertex *v, arrayList *tets)
{
  int i;
  
  for (i=0; i<arrayListSize(tets); i++)
  {
    simplex *s = getFromArrayList(tets,i);
    if (pointOnSimplex(v, s)) return s; 
  }
  return NULL;
}

/******************************************************************************/
// This function will find the neighbours of a given point.
// given a simplex and at least one neighbour.
// This is much more efficient than the previous Natural Neighobour method, 
// because we take a local simplex and then check the neighbourhood for 
// matching simplicies.

arrayList *findNeighbours(vertex *v, simplex *s)
{
  int i;
  arrayList *l   = newArrayList(); 
  stack *toCheck = newStack();

  simplex *current;
  push(toCheck, s);

  while (!isEmpty(toCheck))
  {
    // pop the next one to check from the stack.
    current = pop(toCheck);
    
    // We try to chose the things most likely to fail first, to take 
    // advantage of lazy evaluation.
    if ( pointOnSimplex(v, current) && (! arrayListContains(l, current)) )
    {
      // add this simplex, and check its neighbours.
      addToArrayList(l, current);
      for (i=0; i<4;i++)
        if (current->s[i])
          push(toCheck, current->s[i]);
    }   
  }
  freeStack(toCheck,NULL);

  return l;   
}

/******************************************************************************/
// Given a simplex, we want to find the correctly oriented verticies which are
// not connected

void getRemainingFace(simplex *s, vertex *p, vertex **v1, 
                                             vertex **v2, 
                                             vertex **v3  )
{
  int i,found=0;
  vertex *tmp;
  for (i=0; i<4; i++)
  {
    getFaceVerticies(s, i, v1, v2, v3, &tmp);
    if (tmp == p) 
    {
      found = 1;
      break; 
    }
  }
  // Make sure that we found the point.
  assert(found);
}


/******************************************************************************/

int isNeighbour(simplex *s0, simplex *s1)
{
  int i;
  for (i=0;i<4;i++)
    if (s0->s[i] == s1) return 1;
    
  return 0;
}

/******************************************************************************/

voronoiCell *newVoronoiCell(mesh *m, int n)
{

  voronoiCell *vc ;
  vc = pop(m->deadVoronoiCells);

  if (!vc)
  {
    vc = malloc(sizeof(voronoiCell));
    vc->verticies   = newArrayList();
    vc->nallocated  = 0;
    vc->points      = 0;
    #ifdef DEBUG
    VORONOI_MALLOC ++;  
    #endif
  } else {
    emptyArrayList(vc->verticies);
  }

  // Allocate memory for the point list.
  // We do a realloc, because we want to expand the array to the required size,
  // and then not have to do any more alloc's later. - This is basically
  // a memory pooling technique.
  if (n > vc->nallocated)
  {
    vc->points = realloc(vc->points, sizeof(double)*n);
    
    int i;
    for (i=vc->nallocated; i<n; i++)
    {
      #ifdef DEBUG
      VERTEX_MALLOC++;
      #endif
      vc->points[i] = malloc(sizeof(double)*3);
    }
    vc->nallocated = n;
  }
  vc->n = n;
  return vc;
}

/******************************************************************************/

void addVertexToVoronoiCell(voronoiCell *vc, double *v)
{
  addToArrayList(vc->verticies, v); 
}

/******************************************************************************/
// We use a NULL pointer as a seperator between different faces.

void startNewVoronoiFace(voronoiCell *vc)
{
  addToArrayList(vc->verticies, NULL);
}

/******************************************************************************/
// Given a list of conflicts from the last insert, a list of updates
// from the last insert, and the mesh. We can 'roll-back' the mesh to its
// previous state.

void removePoint(mesh *m)
{
  int i;
  simplex  *s;  

  for (i=0; i< arrayListSize(m->conflicts); i++)
  {
    s = getFromArrayList(m->conflicts,i);
    addSimplexToMesh(m,s);
  }
   
  undoNeighbourUpdates(m->neighbourUpdates);
  
  for (i=0; i<arrayListSize(m->updates); i++)
  {
    s = getFromArrayList(m->updates, i);
    removeSimplexFromMesh(m,s);
  }
}

/******************************************************************************/
// This will take a voronoi cell and calculate the volume.
// the point p is the point which the voronoi cell is defined about.

double voronoiCellVolume(voronoiCell *vc, vertex *p)
{
  int i,j;
  double volume = 0;
  
  for (i=0; i<arrayListSize(vc->verticies); i++)
  {
    double   *thisV;
    double   *firstV;
    double   *lastV = NULL;
    
    // Calculate the center point of this face.
    double center[3] = {0,0,0};
    
    // Find the center point of this vertex.
    for (j=i; j<arrayListSize(vc->verticies); j++)
    {     
      thisV = getFromArrayList(vc->verticies, j);
           
      // We have reached the next face.
      if (!thisV) break;
      vertexAdd(thisV, center, center);          
    }

    // This will give us the center point of the face.
    vertexByScalar(center, 1/(double)(j-i), center);
       
    // First vertex on the face.
    firstV = getFromArrayList(vc->verticies, i);
    lastV  = NULL;
    
    for (j=i; j<arrayListSize(vc->verticies); j++)
    {
      // Get the current vertex from the face.
      thisV = getFromArrayList(vc->verticies,j);
      
      // We've reached the end of this face.
      if (thisV == NULL)
      {
        i=j;
        break;
      }
      // If we have two points to join up, add the volume.
      if (lastV)
        volume += volumeOfTetrahedron(thisV, lastV, p->v, center);   
      else 
        firstV = thisV;
      lastV = thisV;
    }
    // Add the first segment.
    volume += volumeOfTetrahedron(lastV, firstV, p->v, center);      
  }

  assert(volume>0);
  
  return volume;

}

/******************************************************************************/
// a different function for getting Voronoi Cells - actually slower than
// original, so has been removed. In theory does less computation, but has
// bigger overheads for dealing with memory etc.

voronoiCell* getVoronoiCell2(vertex *point, simplex *s0, mesh *m)
{
  arrayList *neighbours = findNeighbours(point, s0);
  int n                 = arrayListSize(neighbours);
  
  int i,j,k;
  
  // Alloc the memory for our new cell.
  voronoiCell *vc = newVoronoiCell(m,n);
  
  // Set all the points to be used in this Voronoi cell. 
  // We do this by going through all n neighbouring simplicies, and calculating
  // their circum centers.
  for (i=0;i<n;i++)
    circumCenter((simplex*)getFromArrayList(neighbours,i), vc->points[i]);
  
  /* The following two lists must be updated atomically */
  
  // This is the list of incident edges.
  arrayList *incidentEdges = newArrayList();
  
  // This is a list of lists of the simplicies attached to those incident edges.
  arrayList *incidentSimplexLists = newArrayList();
  
  // This will extract every edge from the neighoburing simplicies.
  for (i=0; i<arrayListSize(neighbours); i++)
  {
    // get this simplex.
    int sIndex = i;
    simplex *s = getFromArrayList(neighbours,i);

    // Go through every point on this simplex. We ignore one of these: which is
    // the point we are interpolating.
    for (j=0;j<4;j++)
    {
      // The vertex we are considering on the simplex.
      vertex *thisVertex = s->p[j];
      
      // If this is the point we are interpolating, ignore it.
      if (thisVertex == point) continue;
      
      int index = arrayListGetIndex(incidentEdges, thisVertex);
      // This edge is already in the list of incident edges.
      if (index != -1)
      {
        // We the simplex incident to this edge to index'th list 
        // contained within the list incidentSimplicies 
        arrayList *thisList = getFromArrayList(incidentSimplexLists, index);
        int* thisInt = malloc(sizeof(int));
        *thisInt = sIndex;
        addToArrayList(thisList, thisInt);     
      // This edge has not yet been added.
      } else {

        // This edge has not been seen yet, create a new entry in our 
        // edge list, and create a list to contain the simplicies which 
        // are incident to it.
        arrayList *thisList = newArrayList();
        int* thisInt = malloc(sizeof(int));
        *thisInt = sIndex;
        addToArrayList(thisList, thisInt);
        // Note atomic adds: these must always be coherent!
        addToArrayList(incidentEdges, thisVertex);
        addToArrayList(incidentSimplexLists, thisList);
      }
    }
  }
  /* We now have a list of edges which are incident to the point
     being interpolated. (these are defined by one vertex only, as we 
     know that all edges are connected to one common point). For each of these
     entries in the edgeList, there is a list in the list incidentSimplexlists
     which contains a list of all simplicies which are incident to that edge. */

    for (i=0; i<arrayListSize(incidentEdges); i++)
    {
      // the current edge we are considering (defined by the edge which is not
      // the point being interpolated).
      arrayList *incidentSimplicies = getFromArrayList(incidentSimplexLists,i);
      
      // We now want to get the list of simplicies which are incidient to this
      // edge, with a consistant orientation, so that we can calculate
      // the volume.
      // To do this we are going to fetch one simplex at a time, see if it is
      // a neighbour to the simplex we last used. if it is, we make sure it is
      // not a neighbour that we already visited (because there are two
      // valid neighbours for each corresponding direction).
      
      int      currentIndex   = *(int*)getFromArrayList(incidentSimplicies, 0);
      simplex *currentSimplex = getFromArrayList(neighbours, currentIndex);
      simplex *firstSimplex   = currentSimplex;
      simplex *lastConsidered = NULL;
      
      do
      {
        for (k=0; k<arrayListSize(incidentSimplicies); k++)
        {
          int      thisIndex   = *(int*)getFromArrayList(incidentSimplicies,k);
          simplex *thisSimplex = getFromArrayList(neighbours, thisIndex);
          
          if (thisSimplex != lastConsidered && 
                                        isNeighbour(thisSimplex,currentSimplex))
          {

            addVertexToVoronoiCell(vc, vc->points[thisIndex]);
            
            lastConsidered = currentSimplex;  
            currentSimplex = thisSimplex;
            break;
          }
        }        
      } while(currentSimplex != firstSimplex);
     
      startNewVoronoiFace(vc);
    }
    
    
  // Free all of the lists.
  for (i=0; i<arrayListSize(incidentEdges); i++)
  {
    arrayList *thisList = getFromArrayList(incidentSimplexLists, i);
    freeArrayList(thisList,free);
  }
  freeArrayList(incidentEdges, NULL);
  freeArrayList(incidentSimplexLists, NULL);  
  
  return vc;
}

/******************************************************************************/
// This will give us the volume of the voronoi cell about the point p.
// We pass a point, at least one simplex containing that point, and the mesh.

voronoiCell* getVoronoiCell(vertex *point, simplex *s0, mesh *m)
{
  simplex  *s;
  // Find the Natural Neighbour verticies of this point.
  arrayList *neighbours = findNeighbours(point, s0);
  int n = arrayListSize(neighbours);
  
  // If no neighbours were found, it could be because we are trying to 
  // get a cell outside of the points
  if (n==0)
  {
    if (!simplexContainsPoint(m->super, point))
      fprintf(stderr,"Error: point outside of delaunay triangulation. - "
                     "try extending the super-simplex and re-starting.\n");
    else
     fprintf(stderr, "Error: No neighbours found for point! - mesh appears "
                     "to be degenerate.\n");
    exit(1);
  }
  
  // Create a new voronoi cell.
  voronoiCell *vc = newVoronoiCell(m,n);

  // This is the list of simplicies incident to each edge.
  simplex  *simps[n];

  // An edge will always start at our given point, therefore we only need
  // to store the second point to fully define each edge.
  // There are going to be three edges for each simplex.
  vertex      *edges[3*n];
  // Since edges must contain duplicate edges, we want to know 
  // when a particular edge has been visited.
  int i, j = 0, done[3*n];

  for (i=0; i<arrayListSize(neighbours); i++)
  {
    s = getFromArrayList(neighbours, i);
    vertex *v1, *v2, *v3;
    getRemainingFace(s, point, &v1, &v2, &v3);
    
    // Add this simplex to the list note we add three points for each.
    simps[i]        = s;
    edges[3*i]      = v1;   
    edges[3*i+1]    = v2;
    edges[3*i+2]    = v3;
    
    done[3*i]       = 0;
    done[3*i+1]     = 0;
    done[3*i+2]     = 0;
 
    // Calculate the circumcenter of this simplex.
    circumCenter(s, vc->points[i]);
  }

  // For every edge that is in the list, we are going to get the first simplex
  // which is incident to it, and then draw a line from it to the next
  // neighbour that it is incident to it. This next neighbour will be chosen 
  // because it is NOT the last one we considered.  
  // We are effectively rotating around each edge which spans from the point
  // to one of its natural neighbours, and storing the circum-centers of 
  // the simplicies that we found.
  for (i=0; i<3*n; i++)
  {   
    // We don't want to recompute edges.
    if (done[i]) continue;
    
    // This is the current simplex.
    // We are going to find a neighbour for it, which shares an edge, 
    // and is NOT equal to lastConsidered.
    int first   = i;
    int current = i;
    int lastConsidered = -1;
    // Create this voronoi face.

    int match;
    do {
      match=0;
      for (j=0; j < 3*n; j++)
      {
        if (done[j]) continue;
       // if (done[j]) continue;
        // Is this edge shared?
        // Is this simplex a neighbour of the current simplex?
        // Are we making progress: is this a new neighbour?

        if ((edges[i] == edges[j]) && j != lastConsidered
                                   && isNeighbour(simps[current/3], simps[j/3]))
        {
          done[j] = 1;  
          match   = 1;      
          // Add this vertex to this face of this cell.
          addVertexToVoronoiCell(vc, vc->points[j/3]);
          lastConsidered = current;
          current = j;
          break;
        }      
      }  
    } while (match && (current != first));
    
    startNewVoronoiFace(vc);
  }
   
  freeArrayList(neighbours, NULL);

  return vc;
}

/******************************************************************************/

void freeVoronoiCell(voronoiCell *vc, mesh *m)
{
  // We just push the cell to the memory pool.
  // We can free the memory pools manually, or let the program do it 
  // automatically at the end.
  push(m->deadVoronoiCells, vc);
}

/******************************************************************************/
// This writes a voronoi cell to a file. This uses (-1,-1,-1) as a special
// value to seperate faces, so: things will obviously break if there is a 
// circum center on this value... We only really use this function for testing
// so we allow this bug to become a "feature"!

void writeVoronoiCellToFile(FILE* f, voronoiCell *vc)
{
  int i;
  for (i=0; i<arrayListSize(vc->verticies); i++)
  {
    vertex *v = getFromArrayList(vc->verticies, i);
    
    if (!v)
      fprintf(f,"-1 -1 -1\n");
    
    else
      fprintf(f, "%lf %lf %lf\n", v->X, v->Y, v->Z);
  } 
}

/******************************************************************************/
// We should make sure that we only use these two functions for interacting
// with the global simplex list, otherwise the program will behave
// indeterministically.

void addSimplexToMesh(mesh *m, simplex *s)
{
  s->node = addToLinkedList(m->tets, s);  
}

/******************************************************************************/

void removeSimplexFromMesh(mesh *m, simplex *s)
{
  // The simplex has a special pointer which gives its location in the mesh
  // linked list. This allows us to easily remove it from the list.
  removeFromLinkedList(m->tets, s->node);
}

/******************************************************************************/
// This will create a 'super simplex' that contains all of our data to form a
// starting point for our triangulation.

void initSuperSimplex(vertex *ps, int n, mesh *m)
{
  int i;
  m->super = newSimplex(m);
  
  // Get the range of our data set.
  vertex min,max,range;
  getRange(ps, n, &min, &max, &range,1);
  
  // Make the super simplex bigger!
 // vertexByScalar(range.v, 4, range.v);
  
  // We will go clockwise around the base, and then do the top.
  m->superVerticies[0].X =  min.X + range.X/2;
  m->superVerticies[0].Y =  max.Y + 3*range.Y; 
  m->superVerticies[0].Z =  min.Z - range.Z;


  m->superVerticies[1].X =  max.X + 2*range.X;
  m->superVerticies[1].Y =  min.Y - 2*range.Y;
  m->superVerticies[1].Z =  min.Z - range.Z;
  
  m->superVerticies[2].X =  min.X - 2*range.X;
  m->superVerticies[2].Y =  min.Y - 2*range.Y;
  m->superVerticies[2].Z =  min.Z - range.Z;

  m->superVerticies[3].X = min.X + range.X/2;
  m->superVerticies[3].Y = min.Y + range.Y/2;
  m->superVerticies[3].Z = max.Z + 2*range.Z;
  
  // The super-simplex doesn't have any neighbours.
  for (i=0;i<4;i++)
  {
    m->superVerticies[i].index = 0;
    m->super->p[i] = &m->superVerticies[i];
    m->super->s[i] = NULL;
  }
}

/******************************************************************************/
// We are using two stacks, instead of a struct, because it gives us good
// memory advantages. - We will always want about the same number of 
// neighbour updates: using two array-based stacks means that we can always have
// that memory allocated: we should not have to do any memory reallocation
// for the neighbour updating.
// We use a function, so that the process becomes atomic: we don't want 
// to end up with the two stacks being incoherent!

void pushNeighbourUpdate(neighbourUpdate *nu, simplex **ptr, simplex *old)
{
  push(nu->ptrs, ptr);
  push(nu->old,  old);
}

/******************************************************************************/

void freeNeighbourUpdates(neighbourUpdate *nu)
{
  freeStack(nu->ptrs, free);
  freeStack(nu->old,  free);
  free(nu);
}

/******************************************************************************/
// We will go through, and use our neighbour update list to change the 
// neighbour values back to their originals.

void undoNeighbourUpdates(neighbourUpdate *nu)
{
  simplex **thisPtr;
  simplex  *thisSimplex;
  
  // We use isEmpty, because the pointers might sometimes be NULL.
  while (!isEmpty(nu->ptrs))
  {
    thisPtr     = pop(nu->ptrs);
    thisSimplex = pop(nu->old);  
    
    if (thisPtr)
      *thisPtr = thisSimplex;
  }
}

/******************************************************************************/

void resetNeighbourUpdates(neighbourUpdate *nu)
{
  // This will empty the stacks, without freeing any memory. This is the key
  // to this 'memory-saving hack'.
  emptyStack(nu->ptrs);
  emptyStack(nu->old);
}

/******************************************************************************/

neighbourUpdate *initNeighbourUpdates()
{
  neighbourUpdate *nu = malloc(sizeof(neighbourUpdate));
  nu->ptrs = newStack();
  nu->old  = newStack();
  return nu;
}

/******************************************************************************/
// Allocate all the strucutres required to maintain a mesh in memory.

mesh *newMesh()
{
  // Create the struct to hold all of the data strucutres.
  mesh *m             = malloc(sizeof(mesh));
  // Pointer to the super simplex.
  m->super            = NULL;
  // A linked list of simplicies: We can actually remove this without losing
  // any functionality (but it is useful for testing etc.).
  m->tets             = newLinkedList();
  // Instead of freeing old simplicies/voronoi cells, we put them on a stack
  // and reuse them as necesary.
  m->deadSimplicies   = newStack();
  m->deadVoronoiCells = newStack();
  // This is an array of currently conflicting simplicies.
  m->conflicts        = newArrayList();
  // This is an array of the most recently added simplicies.
  m->updates          = newArrayList();
  // This is an array describing the most recent neighbour updates performed.
  m->neighbourUpdates = initNeighbourUpdates();

  return m;
}

/******************************************************************************/

void freeMesh(mesh *m)
{
  #ifdef DEBUG
  printf("Mallocs for vertex: %d.\n", VERTEX_MALLOC);
  printf("Mallocs for simplex: %d.\n", SIMPLEX_MALLOC);
  printf("Mallocs for voronoi: %d.\n", VORONOI_MALLOC);
  #endif
  
  free(m->super);
  freeStack(m->deadSimplicies,   free);
  
  while(!isEmpty(m->deadVoronoiCells))
  {
    voronoiCell *vc = pop(m->deadVoronoiCells);
    int i;
    for (i=0;i<vc->nallocated; i++)
      free(vc->points[i]);
    free(vc->points);
    freeArrayList(vc->verticies, NULL);
    free(vc);
  }
  
  freeStack(m->deadVoronoiCells, NULL);
  freeLinkedList(m->tets,        free);
  freeArrayList(m->conflicts, free);
  freeArrayList(m->updates, NULL);
  freeNeighbourUpdates(m->neighbourUpdates);
  free(m); 
}

/******************************************************************************/
// This will give us the volume of the arbitrary tetrahedron formed by 
// v1, v2, v3, v4
// All arguments are arrays of length three of doubles.

double volumeOfTetrahedron(double *a, double *b, double *c, double *d)
{
  double a_d[3], b_d[3], c_d[3], cross[3];
  
  vertexSub(a,d, a_d);
  vertexSub(b,d, b_d);
  vertexSub(c,d, c_d);
  
  crossProduct(b_d, c_d, cross);  
  double v = scalarProduct(a_d, cross)/(double)6;
   
  return (v >= 0) ? v : -v;
}

/******************************************************************************/

double squaredDistance(double *a)
{
  return scalarProduct(a,a);
}

/******************************************************************************/
// Take the cross product of two verticies and put it in the vertex 'out'.
void crossProduct(double *b, double *c, double *out)
{
  out[0] = b[1] * c[2] - b[2] * c[1];
  out[1] = b[2] * c[0] - b[0] * c[2];
  out[2] = b[0] * c[1] - b[1] * c[0];  
}

/******************************************************************************/

double scalarProduct(double *a, double *b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/******************************************************************************/

void vertexSub(double *a, double *b, double *out)
{
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

/******************************************************************************/

void vertexAdd(double *a, double *b, double *out)
{
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
}

/******************************************************************************/
// Note that this modifies the actual value of the given vertex.

void vertexByScalar(double *a, double b, double *out)
{
  out[0] = a[0] * b;
  out[1] = a[1] * b;
  out[2] = a[2] * b;
}

/******************************************************************************/
// This function will compute the circumcenter of a given simplex.
// -it returns the radius.-

void circumCenter(simplex *s, double *out)
{
  vertex *a, *b, *c, *d;
  getFaceVerticies(s, 0, &a, &b, &c, &d);
 
  double b_a[3]   , c_a[3]   , d_a[3], 
         cross1[3], cross2[3], cross3[3], 
         mult1[3] , mult2[3] , mult3[3], 
         sum[3];
  double denominator;
  
  // Calculate diferences between points.
  vertexSub(b->v, a->v, b_a);
  vertexSub(c->v, a->v, c_a);
  vertexSub(d->v, a->v, d_a);
  
  // Calculate first cross product.
  crossProduct(b_a, c_a, cross1);
  
  // Calculate second cross product.
  crossProduct(d_a, b_a, cross2);
  
  // Calculate third cross product.
  crossProduct(c_a, d_a, cross3);

  vertexByScalar(cross1, squaredDistance(d_a), mult1);
  vertexByScalar(cross2, squaredDistance(c_a), mult2);
  vertexByScalar(cross3, squaredDistance(b_a), mult3);
  
  // Add up the sum of the numerator.
  vertexAdd(mult1, mult2, sum);
  vertexAdd(mult3, sum  , sum);

  // Calculate the denominator.
  denominator = 2*scalarProduct(b_a, cross3);
  
  // Do the division, and output to out.
  vertexByScalar(sum, 1/(double)(denominator), out);
  
  vertexAdd(out, a->v, out);
  
  // Calculate the radius of this sphere. - We don't actually need this.
  // But if we need it for debugging, we can add it back in.
  // return sqrt((double)squaredDistance(sum))/(double)denominator;
}

/******************************************************************************/

int getNumSimplicies(mesh *m)
{
  return linkedListSize(m->tets);
}

/******************************************************************************/

int numSphericalDegenerecies(mesh *m)
{
  return m->cospherical_degenerecies;
}

/******************************************************************************/

int numPlanarDegenerecies(mesh *m)
{
  return m->coplanar_degenerecies;
}

/******************************************************************************/

void getRange(vertex *ps, int n, vertex *min, vertex *max, vertex *range, int r)
{
  int i;
  
  *min = ps[0];
  *max = ps[0];
  
  for (i=0; i<n; i++)
  {
    if (0)
    {
      ps[i].X +=  ((double)rand() / ((double)RAND_MAX + 1) -0.5);
      ps[i].Y +=  ((double)rand() / ((double)RAND_MAX + 1) -0.5);
      ps[i].Z +=  ((double)rand() / ((double)RAND_MAX + 1) -0.5);
    }
    
    max->X = MAX(max->X, ps[i].X);
    max->Y = MAX(max->Y, ps[i].Y);
    max->Z = MAX(max->Z, ps[i].Z);
    
    min->X = MIN(min->X, ps[i].X);
    min->Y = MIN(min->Y, ps[i].Y);
    min->Z = MIN(min->Z, ps[i].Z);   
  }
  
  for (i=0;i<3;i++)
    range->v[i] = max->v[i] - min->v[i];
}

/*******************************************************************************
* Due to the complexity of the Delaunay Meshing, we provide more extensive     *
* debugging tests. We perform checks to ensure that the mesh is Delaunay, to   *
* check that all of the simplicies are consistantly oriented, and to check     *
* that all neighbour relations are correctly assigned.                         *
*******************************************************************************/

/* Turn on the unit testing for this file.                            */
/* We can then compile this file and run it to perform self-testing.  */
#ifdef _TEST_

  #include <sys/time.h>
  /* Most tests rely on asserts: so we make sure these are turned on. */
  #undef NDEBUG
  /* Set this to give more or less debugging information. */                                                                            
  /* The number of points to use in testing. */
  #define NUM_TEST_POINTS 1e4
  
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
  srand ( time(NULL) );
  
  // Create a random pointset for testing.
  vertex *ps = malloc(sizeof(vertex)*NUM_TEST_POINTS);

  for (i=0; i<NUM_TEST_POINTS; i++)
  {
    ps[i].X = (double)rand() / ((double)RAND_MAX + 1);
    ps[i].Y = (double)rand() / ((double)RAND_MAX + 1);
    ps[i].Z = (double)rand() / ((double)RAND_MAX + 1);
    ps[i].U =  ps[i].X;
    ps[i].V =  ps[i].Y;
    ps[i].W =  ps[i].Z;
    ps[i].index = i;
    ps[i].voronoiVolume = -1;
  }

  mesh *delaunayMesh = newMesh();
  
  // Build the mesh, timing how long it takes.
  double t1 = getTime();
  buildMesh(ps, NUM_TEST_POINTS, delaunayMesh);
  double t2 = getTime();
  
  int n = NUM_TEST_POINTS;
  printf("\nMeshed %d points using %d simplicies in %lf seconds.\n", n,
                                     getNumSimplicies(delaunayMesh), t2-t1);
  printf("Co-planar degenerecies fixed: %d.\n",
                                     numPlanarDegenerecies(delaunayMesh));
  printf("Co-spherical degenerecies fixed: %d.\n", 
                                     numSphericalDegenerecies(delaunayMesh));
  
  printf("Now testing mesh...\n");
  
  orientationTest(delaunayMesh->tets);
  delaunayTest(delaunayMesh, ps, NUM_TEST_POINTS);
  faceTest(delaunayMesh);
  
  freeMesh(delaunayMesh);
  printf("Testing Complete.\n");
  return 0;
}

/******************************************************************************/
#endif

