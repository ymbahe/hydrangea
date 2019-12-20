/* C version of katamaran search */

#include <stdio.h>

int ckat(int argc, void *argv[])
{

  long* vA = (long *) argv[0];
  long* vB = (long *) argv[1];

  long n_a = *(long *) argv[2];
  long n_b = *(long *) argv[3];

  long* locs_a = (long *) argv[4];

  /* Start actual katamaran loop */

  long ind_a = 0;
  long ind_b = 0;
    
  long val_a = vA[ind_a];
  long val_b = vB[ind_b];

  while(1)
    {
            
      if (val_a > val_b)
	{
	  ind_b++;
	  if (ind_b >= n_b)
	    break;
	  val_b = vB[ind_b];
	  continue;
	}

      if (val_a < val_b)
	{
	  ind_a++;
	  if (ind_a >= n_a)
	    break;
	  val_a = vA[ind_a];
	  continue;
	}

      if (val_a == val_b)
	{
	  locs_a[ind_a] = ind_b;
	  ind_a++;
	  ind_b++;
	  if (ind_a >= n_a || ind_b >= n_b)
	    break;
	  val_a = vA[ind_a];
	  val_b = vB[ind_b];
	  continue;
	}
    }

  return 0;
}


  

