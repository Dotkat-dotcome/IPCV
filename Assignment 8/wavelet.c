#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                     2-D HAAR WAVELET DECOMPOSITION                       */
/*                                                                          */
/*(Copyright by Pascal Peter, Laurent Hoeltgen and Joachim Weickert, 8/2014)*/
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Remarks: 
  - Note that the program probably only works for square images where
  each dimension is a power of 2. 
  - We use the non-standard decomposition scheme.
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n1)         /* size */

     /* allocates memory for a vector of size n1 */


{
*vector = (float *) malloc (n1 * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough memory available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2)         /* size in direction 2 */

     /* allocates memory for matrix of size n1 * n2 */


{
long i;

*matrix = (float **) malloc (n1 * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (float *) malloc (n2 * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough memory available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n1)         /* size */

     /* disallocates memory for a vector of size n1 */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  n1,         /* size in direction 1 */
      long  n2)         /* size in direction 2 */

     /* disallocates memory for matrix of size n1 * n2 */

{
long i;

for (i=0; i<n1; i++)
    free(matrix[i]);

free(matrix);

return;
}

/*--------------------------------------------------------------------------*/

void read_string

     (char *v)         /* string to be read */

/*
 reads a long value v
*/

{
fgets (v, 80, stdin);
if (v[strlen(v)-1] == '\n')
   v[strlen(v)-1] = 0;
return;
}

/*--------------------------------------------------------------------------*/

void read_long

     (long *v)         /* value to be read */

/*
 reads a long value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%ld", &*v);
return;
}

/*--------------------------------------------------------------------------*/

void read_float

     (float *v)         /* value to be read */

/*
 reads a float value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%f", &*v);
return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

     (const char  *file_name,    /* name of pgm file */ 
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      float       ***u)          /* image, output */   

/* 
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u; 
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
FILE   *inimage;    /* input file */
char   row[80];     /* for reading data */
long   i, j;        /* loop variables */

/* open file */
inimage = fopen (file_name, "rb");
if (NULL == inimage) 
   {
   printf ("could not open file '%s' for reading, aborting.\n", file_name);
   exit (1);
   }

/* read header */
fgets (row, 80, inimage);          /* skip format definition */
fgets (row, 80, inimage);        
while (row[0]=='#')                /* skip comments */
      fgets (row, 80, inimage);
sscanf (row, "%ld %ld", nx, ny);   /* read image size */
fgets (row, 80, inimage);          /* read maximum grey value */

/* allocate memory */
alloc_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=0; j<(*ny); j++) 
 for (i=0; i<(*nx); i++) 
     (*u)[i][j] = (float) getc(inimage);

/* close file */
fclose(inimage);

return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/* 
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",float_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start(arguments,lineformat);

/* convert format string and arguments to plain text line string */
vsprintf(line,lineformat,arguments);

/* add line to total commentary string */
strncat(comment,line,80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   sprintf(comment,"%s\n",comment);

/* close argument list */
va_end(arguments);

return;

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_pgm

     (float  **u,          /* image, unchanged */ 
      long   nx,           /* image size in x direction */
      long   ny,           /* image size in y direction */
      char   *file_name,   /* name of pgm file */
      char   *comments)    /* comment string (set 0 for no comments) */

/* 
  writes a greyscale image into a pgm P5 file;
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
float          aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage) 
   {
   printf("Could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fprintf (outimage, comments);             /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=0; j<ny; j++)
 for (i=0; i<nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

static int int_log2 

     (long val)           /* the value to be taken log2 */

/* Computes the base to logarithm of an integer input value. */

{
long ret = -1;
while (val > 0) 
      {
      val >>= 1;
      ret++;
      }
return ret;
}

/*--------------------------------------------------------------------------*/

void trans_Haar_level 

     (float   *a,         /* input vector */
      long    n)          /* length of the vector */
      
/*
  Performs a 1D Haar Wavelet transform on the vector a. n is the length of the
  vector. At return, the first half of a contains the approximation
  coefficients and the second half the detail coefficients. 
*/
      
{
float  h0, h1;       /* time saving variables */
float  g0, g1;       /* time saving variables */
long   half;         /* time saving variables */
long   i, j;         /* loop variables */

h0 = 1.0 / sqrt (2.0);
h1 = 1.0 / sqrt (2.0);
g0  = -h1;
g1  =  h0;

if (n >= 2)
   {
   half = n >> 1;         
   float *tmp = calloc (n, sizeof(float));
   for (i=0; i<half; i++) 
       {
       j = 2 * i;
       tmp[i]      = a[j] * h0 + a[j+1] * h1;
       tmp[half+i] = a[j] * g0 + a[j+1] * g1;
       }
   for (i=0; i<n; i++)
       a[i] = tmp[i];
   free(tmp);
   }
return;
}

/*--------------------------------------------------------------------------*/

void inv_trans_Haar_level 

     (float   *a,         /* input vector */
      long    n)          /* length of the vector */ 
      
/*
  Performs a 1D inverse Haar Wavelet transform on the vector a. n is the length
  of the vector. The first half of a must contain the approximation
  coefficients and the second half the detail coefficients
*/
      
{
float  h0, h1;       /* time saving variables */
float  g0, g1;       /* time saving variables */
long   half;         /* time saving variables */
long   i, j;         /* loop variables */

h0 = 1.0 / sqrt (2.0);
h1 = 1.0 / sqrt (2.0);
g0  = -h1;
g1  =  h0;

if (n >= 2) 
   {
   half = n >> 1;
   float *tmp = calloc (n, sizeof(float));
   for (i=0; i<half; i++)
       {
       j = 2 * i;
       tmp[j]     = a[i] * h0  + a[half+i] * g0;
       tmp[j+1]   = a[i] * h1  + a[half+i] * g1;
       }
   for (i=0; i<n; i++) 
       a[i] = tmp[i];
   free(tmp);
   }
return;
}

/*--------------------------------------------------------------------------*/

void transpose_data 

     (float   **in,       /* input image */
      long    maxw,       /* maximum width */
      long    maxh)       /* maximum height */
      
/*
  Transposes the data inside the matrix in from 0 to maxw and from 0 to maxh.
  All the remaining coefficients are untouched. 
*/
      
{
long   i, j;        /* loop variables */
float  dummy = 0;   /* swapping variable */
  
for (j=0; j<maxw; j++)
 for (i=j; i<maxh; i++) 
     {
     dummy    = in[i][j];
     in[i][j] = in[j][i];
     in[j][i] = dummy;
     }
return;
}

/*--------------------------------------------------------------------------*/

void decompose

     (long    nx,                  /* image size in x direction */
      long    ny,                  /* image size in y direction */
      long    levels,              /* levels */
      float   **u,                 /* input image */
      float   **coefficients)      /* Haar coefficients */

/* Perform a 2-D Haar Wavelet Decomposition */

{
long   i=0, k=0;      /* loop variables */
float  *tmp = NULL;   /* temporary array for 1-D wavelet transform */
  
for (i=0; i<nx; i++)
 for (k=0; k<ny; k++)
     coefficients[i][k] = u[i][k];

for (k=0; k<=levels; k++) 
    {
    /* Perform the Haar wavelet transform column-wise */
    for (i=0; i<nx/pow(2,k); i++) 
        {
        /* Copy the data to be transformed */
        tmp = calloc (nx / pow (2,k), sizeof (float));
        memcpy (tmp, coefficients[i], (nx / pow (2,k)) * sizeof (float));
    
        /* Apply Haar transform. */
        trans_Haar_level (tmp, nx / pow (2, k));
    
        /* Copy the coefficients back  */
        memcpy (coefficients[i], tmp, (nx / pow (2, k)) * sizeof (float));
        free (tmp);
        }

    /* Perform the Haar wavelet transform row-wise */
    transpose_data (coefficients, nx / pow (2, k), ny / pow (2, k));
    for (i=0; i<ny/pow(2,k); i++) 
        {
        /* Copy the data to be transformed */
        tmp = calloc (ny / pow (2, k), sizeof (float));
        memcpy (tmp, coefficients[i], (ny / pow (2, k)) * sizeof (float));
      
        /* Apply Haar transform. */
        trans_Haar_level (tmp, ny / pow (2, k));
      
        /* Copy the coefficients back. */
        memcpy (coefficients[i], tmp, (ny / pow (2, k)) * sizeof (float));
        free (tmp);
        }
        
    transpose_data (coefficients, ny / pow (2, k), nx / pow (2, k));
    }

return;

} /* decompose */

/*--------------------------------------------------------------------------*/

void reconstruct

     (long    nx,                  /* image size in x direction */
      long    ny,                  /* image size in y direction */
      long    levels,              /* levels */
      float   **u,                 /* input image */
      float   **coefficients)      /* Haar coefficients */

/* Reconstruct image from Haar wavelet coefficients */
 
{
long   i, k;     /* loop variables */
float  *tmp;     /* temporary array for 1-D wavelet transform */
   
for (i=0; i<nx; i++)
 for (k=0; k<ny; k++)
     u[i][k] = coefficients[i][k];
  
for (k=levels; k>=0; k--)
    {
    /* Perform the inverse Haar Wavelet transform row-wise */
    transpose_data (u, nx / pow (2, k), ny / pow (2, k));
    for (i=0; i<ny/pow(2,k); i++) 
        {
        /* Copy the data to be transformed from in into data. */
        tmp = calloc (ny / pow (2, k), sizeof (float));
        memcpy (tmp, u[i], (ny / pow (2, k)) * sizeof (float));
      
        /* Apply inverse Haar transform. */
        inv_trans_Haar_level (tmp, ny / pow (2, k));
      
        /* Copy the coefficients back into in. */
        memcpy (u[i], tmp, (ny / pow (2, k)) * sizeof (float));
        free (tmp);
        }

    /* Perform the inverse Haar Wavelet transform column-wise. */
    transpose_data (u, ny / pow (2, k), nx / pow (2, k));
    for (i=0; i<nx/pow(2,k); i++)
        {
        /* Copy the data to be transformed from in into data. */
        tmp = calloc (nx / pow (2, k), sizeof (float));
        memcpy (tmp, u[i], (nx / pow (2, k)) * sizeof (float));
      
        /* Apply inverse Haar transform. */
        inv_trans_Haar_level (tmp, nx / pow (2, k));
      
        /* Copy the coefficients back into in. */
        memcpy (u[i], tmp, (nx / pow (2, k)) * sizeof (float));
        free (tmp);
        }
    }
  
return;

} /* reconstruct */

/*--------------------------------------------------------------------------*/

void hard_shrinkage

     (long    nx,              /* image size in x direction */
      long    ny,              /* image size in Y direction */
      float   **coefficients,  /* Haar coefficients */
      float   threshold)       /* threshold of shrinkage */
      
/* applies hard shrinkage to Haar wavelet coefficients */      
      
{
long  i, j;         /* loop variables */
long  count = 0;    /* number of zero coefficients */
  
for (j=0; j<nx; j++)
 for (i=0; i<ny; i++) 
     {
     /* SUPPLEMENT CODE HERE */
     //printf("%f",coefficients[i][j]);
     if (coefficients[i][j]<=threshold && coefficients[i][j]>=(-threshold))
     {
         coefficients[i][j] = 0.0;
         count++;
     }
     }
  
printf ("Set %d out of %d coefficients to 0.\n\n", count, nx*ny);

return;
}

/*--------------------------------------------------------------------------*/

void soft_shrinkage

     (long    nx,              /* image size in x direction */
      long    ny,              /* image size in Y direction */
      float   **coefficients,  /* Haar coefficients */
      float   threshold)       /* threshold of shrinkage */ 
      
/* applies soft shrinkage to Haar wavelet coefficients */ 
      
{
long  i, j;         /* loop variables */
long  count = 0;    /* number of zero coefficients */
  
for (j=0; j<nx; j++)
 for (i=0; i<ny; i++) 
     {
     /* SUPPLEMENT CODE HERE */
     if (coefficients[i][j]<=threshold && coefficients[i][j]>=(-threshold))
     {
         coefficients[i][j] = 0.0;
         count++;
     }
     else
     {
         if (coefficients[i][j]>0)
         {
             coefficients[i][j] = coefficients[i][j] - threshold;
         }
         if (coefficients[i][j]<0)
         {
             coefficients[i][j] = coefficients[i][j] + threshold;
         }
         else
         {
             coefficients[i][j] = coefficients[i][j];
         }
     }
     }
  
printf ("Set %d out of %d coefficients to 0.\n", count, nx*ny);
printf ("Rescaled remaining coefficients according to soft shrinkage\n\n");

return;
}
 
/*--------------------------------------------------------------------------*/

void garrote_shrinkage

     (long    nx,              /* image size in x direction */
      long    ny,              /* image size in Y direction */
      float   **coefficients,  /* Haar coefficients */
      float   threshold)       /* threshold of shrinkage */ 
      
/* applies garrote shrinkage to Haar wavelet coefficients */ 
      
{
long  i, j;         /* loop variables */
long  count = 0;    /* number of zero coefficients */
  
if (threshold < 0.000001f)
   {
   printf ("Threshold too low for Garrote-shrinkage, aborting\n");
   return;
   }
  
for (j=0; j<nx; j++)
 for (i=0; i<ny; i++) 
     {
     /* SUPPLEMENT CODE HERE */
     if (coefficients[i][j]<=threshold && coefficients[i][j]>=(-threshold))
     {
         coefficients[i][j] = 0.0;
         count++;
     }
     else
     {
         coefficients[i][j] = coefficients[i][j] - ((threshold*threshold)/coefficients[i][j]);
     }
     }
  
printf ("Set %d out of %d coefficients to 0.\n", count, nx*ny);
printf ("Rescaled remaining coefficients according to Garrote shrinkage\n\n");

return;
}

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **u;                  /* image */
float  **coefficients;       /* wavelet coefficients */
long   nx, ny;               /* image size in x, y direction */ 
long   goal;                 /* choice variable for shrinkage */
float  threshold;            /* threshold for shrinkage */
long   levels;               /* levels */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("WAVELET SHRINKAGE\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert,           \n");
printf ("    Pascal Peter and Laurent Hoeltgen             \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("Non-standard Haar Wavelet decomposition.\n");
printf ("Only works for square images with sizes that are a power of 2.\n\n");

printf ("input image (pgm):                ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("\n");
printf ("Shrinkage types:                \n");
printf ("  (0) hard shrinkage            \n");
printf ("  (1) soft shrinkage            \n");
printf ("  (2) Garrote shrinkage         \n");
printf ("Your choice:                      ");
read_long (&goal);

printf ("\n");
printf ("threshold (>0) (float):           ");
read_float (&threshold);

printf ("output image (pgm):               ");
read_string (out);
printf ("\n");


/* ---- allocate memory for coefficients ---- */

alloc_matrix (&coefficients, nx, ny);


/* ---- process image with your favourite filter ---- */

if (int_log2(nx) < int_log2(ny)) 
   levels = int_log2(nx);
else 
   levels = int_log2(ny);

/* compute wavelet coefficients */
decompose (nx, ny, levels, u, coefficients);

/* perform shrinkage */
if (goal == 0) 
   hard_shrinkage (nx, ny, coefficients, threshold);
if (goal == 1) 
   soft_shrinkage (nx, ny, coefficients, threshold);
if (goal == 2) 
   garrote_shrinkage (nx, ny, coefficients, threshold);
  
/* compute reconstruction */
reconstruct (nx, ny, levels, u, coefficients); 


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
if (goal == 0) comment_line (comments, "# hard shrinkage\n");
if (goal == 1) comment_line (comments, "# soft shrinkage\n");
if (goal == 2) comment_line (comments, "# Garrote shrinkage\n");
comment_line (comments, "# threshold: %8.4f\n", threshold);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u,            nx, ny);
disalloc_matrix (coefficients, nx, ny);

return(0);
}
