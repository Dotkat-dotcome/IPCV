#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                          POINT TRANSFORMATIONS                           */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 8/2014)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Point Transformations:
 - affine rescaling
 - gamma correction
 - histogram equalization
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (double **vector,   /* vector */
      long   n1)         /* size */

     /* allocates memory for a vector of size n1 */


{
*vector = (double *) malloc (n1 * sizeof(double));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough memory available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (double ***matrix,  /* matrix */
      long   n1,         /* size in direction 1 */
      long   n2)         /* size in direction 2 */

     /* allocates memory for matrix of size n1 * n2 */


{
long i;

*matrix = (double **) malloc (n1 * sizeof(double *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
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

     (double *vector,    /* vector */
      long   n1)         /* size */

     /* disallocates memory for a vector of size n1 */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (double **matrix,   /* matrix */
      long   n1,         /* size in direction 1 */
      long   n2)         /* size in direction 2 */

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

void read_double

     (double *v)         /* value to be read */

/*
 reads a double value v
*/

{
char   row[80];    /* string for reading data */

fgets (row, 80, stdin);
if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%lf", &*v);
return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

     (const char  *file_name,    /* name of pgm file */ 
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      double      ***u)          /* image, output */   

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
for (j=1; j<=(*ny); j++) 
 for (i=1; i<=(*nx); i++) 
     (*u)[i][j] = (double) getc(inimage);

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

     (double **u,          /* image, unchanged */ 
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
double         aux;        /* auxiliary variable */
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
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
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

void rescale 

     (double  **u,        /* input image, range [0,255] */
      long    nx,         /* size in x direction */
      long    ny,         /* size in y direction */
      double  a,          /* smallest transformed grey level */
      double  b,          /* largest transformed grey level */
      double  *g)         /* transformed grey levels */

/* 
 affine rescaling of the grey values of u such that 
 min(u) -> a, and max(u) -> b. 
*/

{
long    i, j, k;       /* loop variables */
double  min, max;      /* extrema of u */
double  factor;        /* time saver */

/* determine extrema of u */

/*
 INSERT CODE
*/
min = 255;
max = 0;
for (i=1;i<=nx;i++)
{
  for (j=1;j<=ny;j++)
  {
    if (u[i][j] < min)
    {
      min = u[i][j];
    }
    if (u[i][j] > max)
    {
      max = u[i][j];
    }
  }
}
/* rescale */

/*
 INSERT CODE
*/
for (k=0;k<=255;k++)
{
  g[k] = a + (b-a)*((k-min)/(max-min));
}
return;

} /* rescale */

/*--------------------------------------------------------------------------*/

void gamma_correct

     (double  gamma,      /* gamma correction factor */
      double  *g)         /* transformed grey levels */

/* 
 applies gamma correction to the 256 grey levels that may appear
 in byte wise coded images
*/

{
long  k;   /* loop variable */

/*
 INSERT CODE
*/
for (k=0;k<=255;k++)
{
  g[k] = 255.0*pow((k/255.0),(1/gamma));
}
return;
}

/*--------------------------------------------------------------------------*/

void hist_equal

     (double  **u,        /* input image, range [0,255] */
      long    nx,         /* size in x direction */
      long    ny,         /* size in y direction */
      double  *g)         /* transformed grey levels */

/* performs histogram equalization on u. */

{
long    i, j, k;    /* loop variables */
long    r;          /* current summation index r */
long    k_r;        /* current summation index k_r */
long    n;          /* pixel number */
double  sum;        /* for summing up */
double  psum, qsum; /* sums in equalisation algorithm */
double  hist[256];  /* histogramme */  

/* initialise histogram with zeros */

for (k = 0; k <= 255; k++)
  hist[k]=0;

/* create histogram of u with bin width 1 */

/*
 INSERT CODE
*/
for (i=1;i<=nx;i++)
{
  for (j=1;j<=ny;j++)
  {
    hist[(int)u[i][j]] = hist[(int)u[i][j]] + 1;
  }
}
/* equalization */
r = 0;
k_r = 0;
psum = 0.0;
qsum = 0.0;
/*
 INSERT CODE
*/
for (k_r=0;k_r<=255;k_r++)
{
  qsum = 0.0;
  psum = psum + hist[k_r];
  for (r=0;r<=255;r++)
  {
    qsum = qsum + ((nx*ny)/256.0);
    if (qsum >= psum)
    {
      g[k_r] = r;
      break;
    }
  }
}

return;

} /* hist_equal */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
double **u;                  /* image */
double *g;                   /* grey level mapping */
long   nx, ny;               /* image size in x, y direction */ 
long   i, j;                 /* loop variables */ 
long   transform;            /* type of point transformation */
double a, b;                 /* rescaling bounds */
double gamma;                /* gamma correction factor */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("POINT TRANSFORMATIONS\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("\n");
printf ("available point transformations:\n");
printf ("(0) affine rescaling\n");
printf ("(1) gamma correction\n");
printf ("(2) histogram equalization\n\n");
printf ("your choice:                      ");
read_long (&transform);

if (transform == 0)
   {
   printf ("smallest grey value:              ");
   read_double (&a);
   printf ("largest  grey value:              ");
   read_double (&b);
   }
if (transform == 1)
   {
   printf ("gamma correction factor:          ");
   read_double (&gamma);
   }

printf ("output image (pgm):               ");
read_string (out);
printf ("\n");


/* ---- greyscale transformation ---- */

/* allocate storage for greyscale transformation vector */
alloc_vector (&g, 256);

/* calculate greyscale transformation vector */
if (transform == 0) 
   rescale (u, nx, ny, a, b, g);
if (transform == 1) 
   gamma_correct (gamma, g);
if (transform == 2) 
   hist_equal (u, nx, ny, g);

/* apply greyscale transformation to the image */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = g[(long)(u[i][j])];


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
if (transform == 0)
   {
   comments[0]='\0';
   comment_line (comments, "# affine rescaling\n");
   comment_line (comments, "# a: %8.4lf\n", a);
   comment_line (comments, "# b: %8.4lf\n", b);
   }
if (transform == 1)
   {
   comments[0]='\0';
   comment_line (comments, "# gamma correction\n");
   comment_line (comments, "# gamma: %8.4lf\n", gamma);
   }
if (transform == 2)
   {
   comments[0]='\0';
   comment_line (comments, "# histogram equalization\n");
   }

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_vector (g, 256);
disalloc_matrix (u, nx+2, ny+2);

return(0);
}
