#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>

/* Set up variable type for coordinates, depending on chosen precision */
#ifdef DOUBLEPRECISION
typedef double MyFloat;
#else
typedef float MyFloat;
#endif

extern long NumPart;
extern int  DesDensNgb;
extern float* Hsml;
extern float* Quantity;
extern float* Vel;
extern float* Mass;
extern float Hmax;
extern float *Value, *ValueQuantity;
extern float *Cube, *CubeQuantity;

/* Variables that are only used by hsml-finding part */
extern int *Head, *Next, *Tail, *Len;
extern int *GroupLen, *GroupOffset;
extern FILE *Logfile;
extern double BoxSize, BoxHalf;
extern float Softening;

extern struct particle_data {
  MyFloat Pos[3];	  /*!< Particle position */
}
  *P;                     /*!< points to particles on this processor */

extern struct particle_vel_data {
  float Vel[3];
}
  *PVel;

extern struct r2data {
  float r2;
  int   index;
} 
  *R2list;

extern int    AnzNodes;
extern int    MaxNodes;

extern int *Nextnode;
extern int *Father;

extern struct NODE
{
  float len;			/*!< sidelength of treenode */
  float center[3];		/*!< geometrical center of node */
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      float s[3];               /*!< center of mass of node */
      int mass;            /*!< mass of node */
      int cost;            /*!< counts the number of interactions in which this node is used */
      int sibling;         /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;        /*!< this gives the next node in case the current node needs to be opened */
      int father;          /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */


#endif     /* ALLVARS_H */








