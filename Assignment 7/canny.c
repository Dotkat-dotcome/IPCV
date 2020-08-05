#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define EPSILON 10E-7


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                 CANNY EDGE DETECTOR FOR GREY VALUE IMAGES                */
/*                                                                          */
/*         (Copyright by Markus Mainberger and Pascal Peter, 8/2014)        */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Canny edge detector for grey value images
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

void dummies
 
     (float **u,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = u[i][1];
    u[i][ny+1] = u[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = u[1][j];
    u[nx+1][j] = u[nx][j];
    }
return;
}  

/*--------------------------------------------------------------------------*/

void zero_boundaries
 
     (float **u,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* sets boundary pixels to zero */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = 0.0f;
    u[i][ny+1] = 0.0f;
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = 0.0f;
    u[nx+1][j] = 0.0f;
    }
return;
}  

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

void presmooth

     (float   **f,        /* original image */
      float   **u,        /* smoothed image */
      long    nx,         /* size in x-direction */
      long    ny,         /* size in y-direction */
      float   sigma)      /* std. dev. of Gaussian */
      
/* smoothes an image f (f has to have a border of 1 pixel) */

{
long  i, j;    /* loop variables */

/* copy image */
for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     u[i][j] = f[i][j];
   
if (sigma > 0.0)
   gauss_conv (sigma, nx, ny, 1.0, 1.0, 5.0, 1, u);

return;
}

/*--------------------------------------------------------------------------*/

void apply_nms

     (float  x1,           /* left neighbour */
      float  x2,           /* central value */
      float  *x3,          /* new central value */
      float  x4)           /* right neighbour */

/* sets a point to 0 when one of its neighbours (x1,x4) is larger */

{
if (x1 > x2 || x4 > x2)
   *x3 = 0.0f;

return;
}

/*--------------------------------------------------------------------------*/

float get_direction

     (float x,            /* first component of vector */
      float y)            /* second component of vector */

/* 
  computes the edge direction in radian for a vector (x,y)^T
  used for grey value images 
*/

{
float  direction;     /* direction */

if (fabsf (x) < EPSILON)
   if (fabsf (y) < EPSILON)
      direction = 0;
   else
      direction = 1.5708;
else
   direction = atanf (y / x);

return direction;
}

/*--------------------------------------------------------------------------*/

void nonmaxima_suppression

     (float   **u,        /* image */
      float   **dx,       /* x-derivatives */
      float   **dy,       /* y-derivatives */
      long    nx,         /* x-dimension of u */
      long    ny)         /* y-dimension of u */

/* applies nonmaxima suppression */

{
long   i,j;       /* loop variables */
float  **tmpu;    /* work copy of image u */

alloc_matrix (&tmpu, nx+2, ny+2);

/* copy image u to tmpu */
for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     tmpu[i][j] = u[i][j];

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (tmpu[i][j] !=0.0f)
        {
        float direction;  /* directional information */

        /* get edge direction in pixel (i,j) */
        direction = get_direction (dx[i][j], dy[i][j]);

        /* apply nonmaxima suppression considering neighbours
           in edge direction */
        if (direction >= -0.3927 && direction <= 0.3927)
           apply_nms (tmpu[i-1][j], tmpu[i][j], &u[i][j], tmpu[i+1][j]);
        else if (direction > 0.3927 && direction <= 1.1781)
           apply_nms (tmpu[i-1][j-1], tmpu[i][j], &u[i][j], tmpu[i+1][j+1]);
        else if ((direction > 1.1781  && direction <= 1.5708) ||
                 (direction < -1.1781 && direction >= -1.5708))
           apply_nms (tmpu[i][j-1], tmpu[i][j], &u[i][j], tmpu[i][j+1]);
        else if (direction < -0.3927 && direction >= -1.1781)
           apply_nms (tmpu[i+1][j-1], tmpu[i][j], &u[i][j], tmpu[i-1][j+1]);
        else
           {
           printf ("Error in nonmaxima_suppression(): ");
           printf ("\nUndefined direction\n");
           exit (-1);
           }
        }

/* disallocate memory */
disalloc_matrix (tmpu, nx+2, ny+2);

return;

} /* nonmaxima_suppression */

/*--------------------------------------------------------------------------*/

void trace_edge

     (long    i,          /* x component of current pixel */
      long    j,          /* y component of current pixel */
      float   **u,        /* image */
      float   T1,         /* lower threshold */
      float   T2)         /* higher threshold */

/*
  Consider the neighbours of pixel (i,j). If a neighbour is under the higher
  threshold T2 and over the lower threshold T1 add this pixel to the 
  edge pixels and repeat all steps for that pixel.
*/

{
u[i][j] = 255;

if (u[i+1][j+1] <= T2 && u[i+1][j+1] > T1)
   trace_edge(i+1, j+1, u, T1, T2);

if (u[i+1][j] <= T2 && u[i+1][j] > T1)
   trace_edge(i+1, j, u, T1, T2);

if (u[i+1][j-1] <= T2 && u[i+1][j-1] > T1)
   trace_edge(i+1, j-1, u, T1, T2);

if (u[i-1][j+1] <= T2 && u[i-1][j+1] > T1)
   trace_edge(i-1, j+1, u, T1, T2);

if (u[i-1][j] <= T2 && u[i-1][j] > T1)
   trace_edge(i-1, j, u, T1, T2);

if (u[i-1][j-1] <= T2 && u[i-1][j-1] > T1)
   trace_edge(i-1, j-1, u, T1, T2);

if (u[i][j+1] <= T2 && u[i][j+1] > T1)
   trace_edge(i, j+1, u, T1, T2);

if (u[i][j-1] <= T2 && u[i][j-1] > T1)
   trace_edge(i, j-1, u, T1, T2);

return;

} /* trace_edge */

/*--------------------------------------------------------------------------*/

void get_derivatives

     (float   **u,        /* image */
      float   **dx,       /* derivatives in x-direction */
      float   **dy,       /* derivatives in y-direction */
      long    nx,         /* x-dimension of u */
      long    ny)         /* y-dimension of u */

/*
  compute the derivatives in x and y direction;
  the grid size h is assumed to be 1
*/

{
long  i,j;     /* loop variables */

/* compute the derivatives of u using the sobel operator */
float sobelX;
float sobelY;


for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        sobelX = u[i+1][j+1] + u[i+1][j-1] + 2 * u[i+1][j] - u[i-1][j-1] - 2 * u[i-1][j] - u[i-1][j-1];
        sobelX = sobelX / 8.0f;

        sobelY = u[i-1][j+1] + 2 * u[i][j+1] + u[i+1][j+1] - u[i-1][j-1] - 2 * u[i][j-1] - u[i+1][j-1];
        sobelY = sobelY / 8.0f;

        dx[i][j] = sobelX;
        dy[i][j] = sobelY;
    }
}
/*
 INSERT CODE
*/

return;
}

/*--------------------------------------------------------------------------*/

void hysteresis_thresholding

     (float   **u,        /* image */
      float   T1,         /* lower threshold */
      float   T2,         /* higher threshold */
      long    nx,         /* x-dimension of u */
      long    ny)         /* y-dimension of u */

/* applies hysteresis thresholding and creates the binary edge image */

{
long  i,j;    /* loop variables */

/* use pixels that are larger than T2 as seed points to find
   additional edge pixels that are larger than T1 */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] >= T2) 
        trace_edge (i, j, u, T1, T2);
     }

/* create a binary image */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < 255) 
        u[i][j] = 0.0f;
     }
     
return;
}

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **u;                  /* binary edge image */
float  **f;                  /* original image */
float  **dx;                 /* derivative in x direction */
float  **dy;                 /* derivative in y direction */
long   nx, ny;               /* image size in x, y direction */ 
long   i, j;                 /* loop variables */
float  sigma;                /* Gaussian standard deviation */
float  T1;                   /* lower threshold */
float  T2;                   /* higher threshold */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("CANNY EDGE DETECTOR FOR GREY VALUE IMAGES\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by                             \n");
printf ("    Markus Mainberger and Pascal Peter            \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                      ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &f);


/* ---- read parameters ---- */

printf ("std. dev. of Gaussian sigma (float):    ");
read_float (&sigma);

printf ("lower threshold T1 (float):             ");
read_float (&T1);

printf ("higher threshold T2 (float):            ");
read_float (&T2);

printf ("output image (pgm):                     ");
read_string (out);
printf ("\n");


/* ---- initialization of binary edge image f ---- */

alloc_matrix (&u, nx+2, ny+2);

for (i=0; i<nx+2; i++)
 for (j=0; j<ny+2; j++)
     u[i][j] = 0.0f;


/* ---- apply canny edge detection ---- */

alloc_matrix (&dx, nx+2, ny+2);
alloc_matrix (&dy, nx+2, ny+2);

/* smooth the image with a Gaussian of standard deviation sigma */
presmooth (f, u, nx, ny, sigma);

/* mirror image boundaries */
dummies (u, nx, ny);

/* get the image derivatives */
get_derivatives(u, dx, dy, nx, ny);

/* set boundaries to zero */
zero_boundaries (u, nx, ny);

/* compute the gradient magnitude and threshold
   with T1 to get edge candidates */
for (i=1; i<nx+1; ++i)
 for (j=1; j<ny+1; ++j)
     {
     u[i][j]  = 0.0f;
     u[i][j] += dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j];
     u[i][j]  = sqrtf (u[i][j]);
     u[i][j]  = (u[i][j] >= T1) ? u[i][j] : 0.0f;
     }

/* apply nonmaxima suppression */
nonmaxima_suppression (u, dx, dy, nx, ny);

/* apply hysteresis thresholding */
hysteresis_thresholding (u, T1, T2, nx, ny);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# canny edge detector for grey value images\n");
comment_line (comments, "# sigma: %8.4f\n", sigma);
comment_line (comments, "# T1:    %8.4f\n", T1);
comment_line (comments, "# T2:    %8.4f\n", T2);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u,  nx+2, ny+2);
disalloc_matrix (f,  nx+2, ny+2);
disalloc_matrix (dx, nx+2, ny+2);
disalloc_matrix (dy, nx+2, ny+2);

return(0);
}
