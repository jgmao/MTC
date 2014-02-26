//  outputfunction for images
//  range 0	127 + (127 / (maxin - minin)
//  range 1	255 / (maxin - minin)
//  range 3	just clipped between 0 and 255
//
//================================================== 
 #include <stdio.h>
 #include <Imlib.h>
 #include <math.h>
 #include <stdlib.h>

namespace metric{
void
output (double *imin, int W, int H, ImlibData * idout,
	ImlibImage * imout, char *nameout, int range)
{
  int i1, j1;
  double minin, maxin;
  minin = maxin = 0.0;
  for (j1 = 0; j1 < H; j1++)
    {
      for (i1 = 0; i1 < W; i1++)
	{
	  int index = i1 + W * j1;
	  if (imin[index] > maxin) {maxin = imin[index];}
	  if (imin[index] < minin) {minin = imin[index];}
	}
    }
  //printf ("Nameout %s \tmaxim = %lf \t minim = %lf\n", nameout, maxin, minin);
  //printf ("Red pixels represent nan-values\n"); 
  //============================================================================
  
  if (range == 0)
    {
      for (j1 = 0; j1 < H; j1++)
	{
	  for (i1 = 0; i1 < W; i1++)
	    {
	      int index  = 3 * (i1 + W * j1);
	      int index2 = i1 + W * j1;
	      
	      imout->rgb_data[index] = imout->rgb_data[index + 1] =
		imout->rgb_data[index + 2] =
		127 + (127 / (maxin - minin)) * imin[index2];
		
	      if (isnan(imin[index2])) 
	        {
		  imout->rgb_data[index]     = 255;
		  imout->rgb_data[index + 1] = 0;
		  imout->rgb_data[index + 2] = 0;
		}
	    }
	}
    }
    
  if (range == 1)
    {
      for (j1 = 0; j1 < H; j1++)
	{
	  for (i1 = 0; i1 < W; i1++)
	    {
	      int index  = 3 * (i1 + W * j1);
	      int index2 = i1 + W * j1;
	      imout->rgb_data[index] = imout->rgb_data[index + 1] =
		imout->rgb_data[index + 2] =
		(255 / (maxin - minin)) * (imin[index2] - minin);
		
	      if (isnan(imin[index2])) 
	        {
		  imout->rgb_data[index]     = 255;
		  imout->rgb_data[index + 1] = 0;
		  imout->rgb_data[index + 2] = 0;
		}
		
	    }
	}
    }
    
  if (range == 2)
    {
      for (j1 = 0; j1 < H; j1++)
	{
	  for (i1 = 0; i1 < W; i1++)
	    {
	      int index  = 3 * (i1 + W * j1);
	      int index2 = i1 + W * j1;
	      if ((imin[index2] - minin) >= 0.0)
		{
		  imout->rgb_data[index] = imout->rgb_data[index + 1] =
		    imout->rgb_data[index + 2] = (int) (imin[index2] - minin);
		}
	      if ((imin[index2] - minin) < 0.0)
		{
		  imout->rgb_data[index] = 255;
		  imout->rgb_data[index + 1] = imout->rgb_data[index + 2] = 0;
		}
	      if ((imin[index2] - minin) > 255.0)
		{
		  imout->rgb_data[index] = 255;
		  imout->rgb_data[index + 1] =
		    imout->rgb_data[index + 2] = 255;
		}    
	      if (isnan(imin[index2])) 
	        {
		  imout->rgb_data[index]     = 255;
		  imout->rgb_data[index + 1] = 0;
		  imout->rgb_data[index + 2] = 0;
		}
	      
	    }
	}
    }
    
  if (range == 3)
    {
      for (j1 = 0; j1 < H; j1++)
	{
	  for (i1 = 0; i1 < W; i1++)
	    {
	      int index  = 3 * (i1 + W * j1);
	      int index2 = i1 + W * j1;
	      if ((imin[index2] >= 0.0) && (imin[index2] <= 255.0))
		{
		  imout->rgb_data[index] = imout->rgb_data[index + 1] 
		   = imout->rgb_data[index + 2] = (int) (imin[index2] + 0.5);
		}
	      if (imin[index2] < 0.0)
		{
		  imout->rgb_data[index] = imout->rgb_data[index + 1] 
		   = imout->rgb_data[index + 2] = 0;
		}
	      if (imin[index2] > 255.0)
		{
		  imout->rgb_data[index] = imout->rgb_data[index + 1] 
		   = imout->rgb_data[index + 2] = 255;
		}    
	      if (isnan(imin[index2])) 
	        {
		  imout->rgb_data[index]     = 255;
		  imout->rgb_data[index + 1] = 0;
		  imout->rgb_data[index + 2] = 0;
		}
	      
	    }
	}
    }
  Imlib_save_image_to_ppm (idout, imout, nameout);
}

}
