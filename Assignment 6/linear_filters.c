#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define PRECISION 3.f
#define BOUNDARY 1
#define GRID_SIZE 1.f


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*          GAUSSIAN-BASED HIGHPASS, LOWPASS AND BANDPASS FILTERS           */
/*                                                                          */
/*         (Copyright by Joachim Weickert and Pascal Peter, 8/2014)         */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Gaussian-based highpass, lowpass and bandpass filters
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
for (j=1; j<=(*ny); j++) 
 for (i=1; i<=(*nx); i++) 
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

     (float   **u,        /* input image */
      long    nx,         /* size in x direction */
      long    ny,         /* size in y direction */
      float   a,          /* smallest transformed grey level */
      float   b)          /* largest transformed grey level */

/* 
 affine rescaling of the grey values of u such that 
 min(u) -> a, and max(u) -> b. 
*/

{
long   i, j;       /* loop variables */
float  min, max;   /* extrema of u */
float  m,n;        /* time saver */
float  g[256];     /* mapping of grey values */

/* determine extrema of u */
min = max = u[1][1];
for (i=1; i<=nx; i++)
    {
    for (j=1; j<=ny; j++)
        {
        if (u[i][j] < min) min = u[i][j];
        if (u[i][j] > max) max = u[i][j];
        }
    }

/* rescale */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = a + (u[i][j] - min) / (max - min) * (b - a);

return;

} /* rescale */

/*--------------------------------------------------------------------------*/

void gauss_conv 

     (float    sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      float    **f)       /* input: original image ;  output: smoothed */


/*  Gaussian convolution. */

{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */
      

/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hx) + 1;

if ((bc != 0) && (length > nx))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
              * exp (- (i * i * hx * hx) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[nx+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[nx+length-1+p] = help[nx+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[nx+length-p];
           help[nx+length-1+p] = help[length+p-1];
           }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);


/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hy) + 1;
if ((bc != 0) && (length > ny))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) 
              * exp (- (j * j * hy * hy) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[ny+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[ny+length-1+p] = help[ny+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[ny+length-p];
           help[ny+length-1+p] = help[length+p-1];
           } 
 
    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        f[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);

return;

} /* gauss_conv */

/*--------------------------------------------------------------------------*/

void lowpass

     (float   sigma,      /* standard deviation of Gaussian */
      long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      float   **f)        /* input: original; output: lowpass filtered */

/* lowpass filter */

{  
/* Apply Gaussian for lowpass filtering */
gauss_conv (sigma, nx, ny, GRID_SIZE, GRID_SIZE, PRECISION, BOUNDARY, f);

return;
}

/*--------------------------------------------------------------------------*/

void highpass

     (float   sigma,      /* standard deviation of Gaussian */
      long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      float   **f)        /* input: original; output: highpass filtered */
      
/* highpass filter */

{  
long   i,j;       /* loop variables */
float  **gauss;   /* Gaussian smoothed image */

/* allocate memory */
alloc_matrix (&gauss, nx+2, ny+2);

/* copy original image to temporary array */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     gauss[i][j] = f[i][j];

/* apply Gaussian to temporary array */
gauss_conv(sigma, nx, ny, GRID_SIZE, GRID_SIZE, PRECISION, BOUNDARY, gauss);

/* compute highpass filter */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++) 
     {
     /* SUPPLEMENT YOUR CODE HERE */
        f[i][j] = f[i][j] - gauss[i][j];
     }
  
/* free memory */
disalloc_matrix (gauss, nx+2, ny+2);

return;
}

/*--------------------------------------------------------------------------*/

void bandpass

     (float   sigma1,     /* standard deviation of first Gaussian */
      float   sigma2,     /* standard deviation of second Gaussian */
      long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      float   **f)        /* input: original; output: bandpass filtered */

/* bandpass filter: sigma1 > sigma2 */

{
long   i,j;        /* loop variables */
float  **gauss1;   /* Gaussian smoothed image (stddev sigma1) */
float  **gauss2;   /* Gaussian smoothed image (stddev sigma2) */

/* allocate memory */
alloc_matrix (&gauss1, nx+2, ny+2);
alloc_matrix (&gauss2, nx+2, ny+2);

/* copy original image to temporary arrays */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++) 
     {
     gauss1[i][j] = f[i][j];
     gauss2[i][j] = f[i][j];
     }

/* apply Gaussian to temporary arrays */
gauss_conv (sigma1, nx, ny, GRID_SIZE, GRID_SIZE, PRECISION, BOUNDARY, gauss1);
gauss_conv (sigma2, nx, ny, GRID_SIZE, GRID_SIZE, PRECISION, BOUNDARY, gauss2);

/* compute bandpass filter */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++) 
     {
     /* SUPPLEMENT YOUR CODE HERE */
        f[i][j] = gauss1[i][j] - gauss2[i][j];
     }

/* free memory */
disalloc_matrix (gauss1, nx+2, ny+2);
disalloc_matrix (gauss2, nx+2, ny+2);

return;
}

/*--------------------------------------------------------------------------*/

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *std)        /* standard deviation, output */

/*
 computes minimum, maximum, mean, and standard deviation of an image u
*/

{
long    i, j;       /* loop variables */
double  help1;      /* auxiliary variable */
float   help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help1 = help1 + (double)u[i][j];
     }
*mean = (float)help1 / (nx * ny);

*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help2  = u[i][j] - *mean;
     *std = *std + help2 * help2;
     }
*std = sqrt(*std / (nx * ny));

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **u;                  /* image */
long   nx, ny;               /* image size in x, y direction */ 
float  sigma1;               /* standard deviation for first Gaussian */
float  sigma2;               /* standard deviation for second Gaussian */
long   filter;               /* variable for filter choice */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("GAUSSIAN-BASED HIGHPASS, LOWPASS AND BANDPASS FILTERS\n\n");
printf ("*****************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert               \n");
printf ("    and Pascal Peter                                 \n");
printf ("    Dept. of Mathematics and Computer Science        \n");
printf ("    Saarland University, Saarbruecken, Germany       \n\n");
printf ("    All rights reserved. Unauthorized usage,         \n");
printf ("    copying, hiring, and selling prohibited.         \n\n");
printf ("    Send bug reports to                              \n");
printf ("    weickert@mia.uni-saarland.de                     \n\n");
printf ("*****************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                                           ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("\navailable linear filters:\n");
printf ("(0) lowpass                \n");
printf ("(1) highpass               \n");
printf ("(2) bandpass               \n");
printf ("your choice:                                                 ");
read_long (&filter);

printf ("\n");
printf ("standard deviation sigma1 for Gaussian (float):              ");
read_float (&sigma1);

if (filter == 2) 
   {
   printf ("standard deviation sigma2 for second Gaussian (<sigma1):     ");
   read_float (&sigma2);
   }

printf ("output image (pgm):                                          ");
read_string (out);
printf ("\n");


/* ---- process image with linear filter ---- */

if (filter == 0) 
   {
   printf ("Applying lowpass filter...\n\n");
   lowpass (sigma1, nx, ny, u);
   }
if (filter == 1) 
   {
   printf ("Applying highpass filter...\n\n");
   highpass (sigma1, nx, ny, u);
   } 
if (filter == 2) 
   {
   printf ("Applying bandpass filter...\n\n");
   bandpass (sigma1, sigma2, nx, ny, u);
   }

/* perform affine rescaling */
if (filter >= 1)
   rescale(u, nx, ny, 0, 255);


/* ---- analyse filtered image ---- */

analyse (u, nx, ny, &min, &max, &mean, &std);
printf ("filtered image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
if (filter == 0)
   {
   comment_line (comments, "# lowpass filter\n");
   comment_line (comments, "# sigma1: %8.4f\n", sigma1);
   }
if (filter == 1)
   {
   comment_line (comments, "# highpass filter\n");
   comment_line (comments, "# sigma1: %8.4f\n", sigma1);
   }
if (filter == 2)
   {
   comment_line (comments, "# bandpass filter\n");
   comment_line (comments, "# sigma1: %8.4f\n", sigma1);
   comment_line (comments, "# sigma2: %8.4f\n", sigma2);
   }

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);

return(0);
}
