/* C version of katamaran search */

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
  long val;
  long index;
} vi;

int compare(const void* pa, const void* pb)
{
  vi va = *((vi*)pa);
  vi vb = *((vi*)pb);
  
  if (va.val == vb.val) return 0;
  else if (va.val < vb.val) return -1;
  else return 1;
}


int ckat(int argc, void *argv[])
{

  long* vA = (long *) argv[0];
  long* vB = (long *) argv[1];

  long n_a = *(long *) argv[2];
  long n_b = *(long *) argv[3];

  long* locs_a = (long *) argv[4];

  vi* varr_a = malloc(n_a*sizeof(vi));
  vi* varr_b = malloc(n_b*sizeof(vi));

  long ii;

  for (ii = 0; ii < n_a; ii++) {
    varr_a[ii].index = ii;
    varr_a[ii].val = vA[ii];
  }

  for (ii = 0; ii < n_b; ii++) {
    varr_b[ii].index = ii;
    varr_b[ii].val = vB[ii];
  }


  qsort(varr_a, n_a, sizeof(vi), compare);
  qsort(varr_b, n_b, sizeof(vi), compare);
  

  /* Start actual katamaran loop */

  long ind_a = 0;
  long ind_b = 0;
    
  long val_a = varr_a[ind_a].val;
  long val_b = varr_b[ind_b].val;

  while(1)
    {
            
      if (val_a > val_b)
	{
	  ind_b++;
	  if (ind_b >= n_b)
	    break;
	  val_b = varr_b[ind_b].val;
	  continue;
	}

      if (val_a < val_b)
	{
	  ind_a++;
	  if (ind_a >= n_a)
	    break;
	  val_a = varr_a[ind_a].val;
	  continue;
	}

      if (val_a == val_b)
	{
	  locs_a[varr_a[ind_a].index] = varr_b[ind_b].index;
	  ind_a++;
	  ind_b++;
	  if (ind_a >= n_a || ind_b >= n_b)
	    break;
	  val_a = varr_a[ind_a].val;
	  val_b = varr_b[ind_b].val;
	  continue;
	}
    }

  return 0;
}


  

