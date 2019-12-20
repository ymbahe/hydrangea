#include "allvars.h"

long NumPart;


int *Head, *Next, *Tail, *Len;
int *GroupLen, *GroupOffset;
FILE *Logfile;
double BoxSize, BoxHalf;
int DesDensNgb;

float *Hsml;
float *Quantity;
float *Vel;
float *Mass;
float Softening;
float *Value, *ValueQuantity;
float *Cube, *CubeQuantity;
float Hmax;

struct particle_data *P;     /*!< points to particles on this processor */
struct particle_vel_data *PVel;

struct r2data *R2list;

int    AnzNodes;
int    MaxNodes;

int *Nextnode;
int *Father;

struct NODE *Nodes_base;                    /*!< points to the actual memory allocted for the nodes */ 
struct NODE *Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */

