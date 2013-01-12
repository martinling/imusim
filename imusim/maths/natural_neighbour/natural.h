/*******************************************************************************
*
*  interpolate.h - By Ross Hemsley Aug. 2009 - rh7223@bris.ac.uk.
*  
*  This unit will perform Natural-Neighbour interpolation. To do this, we first
*  need a Mesh of Delaunay Tetrahedrons for our input points. Each point to be
*  interpolated is then inserted into the mesh (remembering the steps that were
*  taken to insert it) and then the volume of the modified Voronoi cells 
*  (easily computed from the Delaunay Mesh) are used to weight the neighbouring
*  points. We can then revert the Delaunay mesh back to the original mesh by 
*  reversing the flips required to insert the point.
*
*******************************************************************************/

#ifndef natural_h
#define natural_h

/******************************************************************************/
vertex *loadPoints(char *filename, int *n);
//------------------------------------------------------------------------------
vertex *initPoints(double *x, double *y, double *z, 
                   double *u, double *v, double *w, int n);
//------------------------------------------------------------------------------
void    writePointsToFile(vertex *ps, int n);
//------------------------------------------------------------------------------
void    lastNaturalNeighbours(vertex *v, mesh *m, arrayList *neighbours, 
                                                arrayList *neighbourSimplicies);
//------------------------------------------------------------------------------                                               
void     interpolate3_3(double  x, double  y, double  z, 
                        double *u, double *v, double *w, mesh *m);
                        
/******************************************************************************/
#endif

