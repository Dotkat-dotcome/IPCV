#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                    OPTIC FLOW USING LUCAS AND KANADE                     */
/*                                                                          */
/*  (Copyright by and Stephan Didas, 1/2007 and Joachim Weickert, 8/2014)   */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - Optic flow calculation with method by Lucas and Kanade
*/


/*----------------------------------------------------------------------------*/

float maxf(float a,float b)
{
    return a>b ? a : b;
}

/*---------------------------------------------------------------------------*/

float minf(float a,float b)
{
    return a<b ? a : b;
}

/*---------------------------------------------------------------------------*/

int maxi(int a,int b)
{
  return (a>b) ? a : b;
}

/*---------------------------------------------------------------------------*/

int mini(int a,int b)
{
  return (a<b) ? a : b;
}

/*---------------------------------------------------------------------------*/

int clamped (int min, int val, int max)
{
    assert(min<=max);
    return mini(maxi(min,val),max);
}

/*----------------------------------------------------------------------------*/

int byte_range(int a)
/* restricts number to unsigned char range */
{
    return clamped(0,a,255);
}

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

void alloc_cubix

     (float ****cubix,  /* cubix */
      long  n1,         /* size in direction 1 */
      long  n2,         /* size in direction 2 */
      long  n3)         /* size in direction 3 */

     /* allocates memory for cubix of size n1 * n2 * n3 */


{
long i, j;

*cubix = (float ***) malloc (n1 * sizeof(float **));
if (*cubix == NULL) 
   {
   printf("alloc_cubix: not enough memory available\n");
   exit(1);
   }
for (i=0; i<n1; i++) 
    {
    (*cubix)[i] = (float **) malloc (n2 * sizeof(float *));
    if ((*cubix)[i] == NULL) 
       {
       printf("alloc_cubix: not enough memory available\n");
       exit(1);
       }
    for (j=0; j<n2; j++) 
        {
        (*cubix)[i][j] = (float *) malloc (n3 * sizeof(float));
        if ((*cubix)[i][j] == NULL) 
           {
           printf("alloc_cubix: not enough memory available\n");
           exit(1);
           }
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

void disalloc_cubix

     (float ***cubix,   /* cubix */
      long  n1,         /* size in direction 1 */
      long  n2,         /* size in direction 2 */
      long  n3)         /* size in direction 3 */

     /* disallocates memory for cubix of size n1 * n2 * n3 */

{
long i, j;

for (i=0; i<n1; i++) 
 for (j=0; j<n2; j++) 
     free(cubix[i][j]);

for (i=0; i<n1; i++)
    free(cubix[i]);

free(cubix);

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
      char   *comments,    /* comment string (set 0 for no comments) */
      long   flag,         /* is rescaling required? */
      float  max)          /* max grey value */

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
     if (flag == 0)
        aux = u[i][j] + 0.499999;    /* for correct rounding */
     else
        aux = ((255.0f * u[i][j]) / max) + 0.499999;
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

/*----------------------------------------------------------------------------*/

void write_ppm

     (float ***u,          /* colour image, unchanged */
      int   nx,            /* size in x direction */
      int   ny,            /* size in y direction */
      char  *file_name,    /* name of ppm file */
      char  *comments)     /* comment string (set 0 for no comments) */

     /* writes an image into a pgm P5 (greyscale) or ppm P6 (colour) file */

{
FILE           *outimage;  /* output file */
int            i, j, m;    /* loop variables */
float          aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage) 
   {
   printf("Could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

fprintf (outimage, "P6\n");                /* colour format */

if (comments != 0)
   fprintf (outimage, comments);           /* comments */

fprintf (outimage, "%d %d\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");               /* maximal value */

/* write image data */
for (j = 1; j <= ny; j++) 
 for (i = 1; i <= nx; i++) 
  for (m = 0; m < 3; m++) 
      {
      aux = u[i][j][m] + 0.499999;    /* for correct rounding */
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
  
} /* write_ppm */

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

/*----------------------------------------------------------------------------*/

void vector_to_RGB

     (float x,   /* x-component */
      float y,   /* y-component */
      int   *R,  /* red component */
      int   *G,  /* green component */
      int   *B)  /* blue component */

/* 
  Computes the color representation of a vector. 
*/

{               
float Pi;          /* pi */
float amp;         /* amplitude (magnitude) */
float phi;         /* phase (angle) */
float alpha, beta; /* weights for linear interpolation */

/* set pi */
Pi = 2.0 * acos(0.0);

/* determine amplitude and phase (cut amp at 1) */
amp = sqrt (x * x + y * y);
if (amp > 1) amp = 1;
if (x == 0.0)
  if (y >= 0.0) phi = 0.5 * Pi;
  else phi = 1.5 * Pi;
else if (x > 0.0)
  if (y >= 0.0) phi = atan (y/x);
  else phi = 2.0 * Pi + atan (y/x);
else phi = Pi + atan (y/x);

phi = phi / 2.0;

// interpolation between red (0) and blue (0.25 * Pi)
if ((phi >= 0.0) && (phi < 0.125 * Pi)) 
   {
   beta  = phi / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *G = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
   *B = (int)floor(amp * (alpha *   0.0 + beta * 255.0));
   }
if ((phi >= 0.125 * Pi) && (phi < 0.25 * Pi)) 
   {
   beta  = (phi-0.125 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta *  64.0));
   *G = (int)floor(amp * (alpha *   0.0 + beta *  64.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   }

// interpolation between blue (0.25 * Pi) and green (0.5 * Pi)
if ((phi >= 0.25 * Pi) && (phi < 0.375 * Pi)) 
   {
   beta  = (phi - 0.25 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha *  64.0 + beta *   0.0));
   *G = (int)floor(amp * (alpha *  64.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   }
if ((phi >= 0.375 * Pi) && (phi < 0.5 * Pi)) 
   {
   beta  = (phi - 0.375 * Pi) / (0.125 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
   }
   
// interpolation between green (0.5 * Pi) and yellow (0.75 * Pi)
if ((phi >= 0.5 * Pi) && (phi < 0.75 * Pi)) 
   {
   beta  = (phi - 0.5 * Pi) / (0.25 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 0.0   + beta * 255.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *B = (int)floor(amp * (alpha * 0.0   + beta * 0.0));
   }

// interpolation between yellow (0.75 * Pi) and red (Pi)
if ((phi >= 0.75 * Pi) && (phi <= Pi)) 
   {
   beta  = (phi - 0.75 * Pi) / (0.25 * Pi);
   alpha = 1.0 - beta;
   *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
   *G = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
   *B = (int)floor(amp * (alpha * 0.0   + beta *   0.0));
   }

/* check RGB range */
*R = byte_range(*R);
*G = byte_range(*G);
*B = byte_range(*B);

return;

} /* vector_to_RGB */

/*----------------------------------------------------------------------------*/

void flow_to_color

     (float **u,            /* flow field, first channel (input) */
      float **v,            /* flow field, second channel (input) */
      float ***color_img,   /* color representation (output) */
      float max_disp,       /* maximal disparity (set to -1 if not used) */
      int nx,               /* size in x direction */
      int ny)               /* size in y direction */

/* 
  Computes a color representation of a flow field. 
*/

{
int i,j;
int R,G,B;
float maximum_length = 0;

for (i=1; i<=nx; i++) 
 for (j=1; j<=ny; j++) 
     maximum_length = maxf(maximum_length, u[i][j] * u[i][j]
                                         + v[i][j] * v[i][j]);

if(max_disp==-1.0f)
   maximum_length = 1.0f / sqrt(maximum_length);
else 
   maximum_length = 1.0f / max_disp;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++) 
     {
     if(u[i][j]!=100.0f && v[i][j]!=100.0f) 
        vector_to_RGB(maximum_length*u[i][j],
                      maximum_length*v[i][j],&R,&G,&B);
     else
        R=G=B=128.0f;

     color_img[i][j][0]=(float)R;
     color_img[i][j][1]=(float)G;
     color_img[i][j][2]=(float)B;
     }

return;

} /* flow_to_color */

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

/*
 Gaussian convolution.
*/

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

/* ---------------------------------------------------------------------- */

void create_eq_systems

     (float    **f1,      /* frame 1, input */
      float    **f2,      /* frame 2, input */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    ht,        /* distance between two frames */
      float    rho,       /* integration scale */
      float    **dxx,     /* element of structure tensor, output */
      float    **dxy,     /* element of structure tensor, output */
      float    **dyy,     /* element of structure tensor, output */
      float    **dxz,     /* element of the right-hand side, output */
      float    **dyz)     /* element of the right-hand side, output */

/*
 Calculates the entries of the structure tensor and the right-hand side
 for the linear systems of equations.
*/

{
long    i, j;                 /* loop variables */
float   df_dx, df_dy, df_dz;  /* derivatives of f */
float   w1, w2, w3, w4, w5;   /* time savers */


/* ---- calculate gradient and its tensor product ---- */

dummies(f1, nx, ny); 
dummies(f2, nx, ny); 

for(i=1; i<=nx; i++) 
 for(j=1; j<= ny; j++) 
    {
    /* compute the spatial derivatives using Sobel operators */
    w1 = 1.0 / (8.0 * hx);
    w2 = 1.0 / (4.0 * hx);
 
    df_dx = w1 * (f1[i+1][j+1] - f1[i-1][j+1]
                + f1[i+1][j-1] - f1[i-1][j-1])
          + w2 * (f1[i+1][j]   - f1[i-1][j]);

    w3 = 1.0 / (8.0 * hy);
    w4 = 1.0 / (4.0 * hy);
 
    df_dy = w3 * (f1[i+1][j+1] - f1[i+1][j-1]
                + f1[i-1][j+1] - f1[i-1][j-1])
          + w4 * (f1[i][j+1]   - f1[i][j-1]);

    /* time derivative with forward difference, ht = 1 */
    /* 
      SUPPLEMENT CODE HERE
    */
    w5 = 1.0;
    df_dz = w5 * (f2[i][j] - f1[i][j]);

    /* calculate matrix entries and right-hand side */
    dxx[i][j] = df_dx * df_dx;
    dxy[i][j] = df_dx * df_dy;
    dyy[i][j] = df_dy * df_dy;
    /* 
      SUPPLEMENT CODE HERE
    */
    dxz[i][j] = df_dx * df_dz;
    dyz[i][j] = df_dy * df_dz;
    }


/* ---- smoothing at integration scale, Dirichlet b.c. ---- */

if (rho > 0.0)
   {
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxx);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxy);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dyy);
   /* 
     SUPPLEMENT CODE HERE
   */
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxz);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dyz);
   }

return;

} /* create_eq_systems */

/* ---------------------------------------------------------------------- */

void lucas_kanade 

     (long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      float    **dxx,       /* component of system matrix, input */
      float    **dxy,       /* component of system matrix, input */
      float    **dyy,       /* component of system matrix, input */
      float    **dxz,       /* component of right-hand side, input */
      float    **dyz,       /* component of right-hand side, input */
      float    eps,         /* threshold */
      float    **c,         /* confidence of the flow, output */
      float    **u,         /* x component of optic flow, output */
      float    **v)         /* v component of optic flow, output */

/* 
   Performs optic flow estimation with Lucas-Kanade method.
*/

{
long  i, j;      /* loop variables */
float trace;     /* trace of the matrix at pixel [i][j] */
float det;       /* determinant of the matrix at pixel [i][j] */

for(i=1; i<=nx; i++) 
 for(j=1; j<=ny; j++) 
    {
    /* 
      SUPPLEMENT CODE HERE
    */
    trace = dxx[i][j] + dyy[i][j];
    det = (dxx[i][j] * dyy[i][j]) - (dxy[i][j] * dxy[i][j]);
    if (trace <= eps && det <= eps) 
    {
        u[i][j] = 0;
        v[i][j] = 0;
        c[i][j] = 0;
    }
    if (det <= eps && trace > eps)
    {
        u[i][j] = -dxz[i][j]/(dxx[i][j]+dyy[i][j]);
        v[i][j] = -dyz[i][j]/(dxx[i][j]+dyy[i][j]);
        c[i][j] = 128;
    }
    if (det > eps && trace > eps)
    {
        u[i][j] = (1/det)*(-dyy[i][j]*dxz[i][j] + dxy[i][j]*dyz[i][j]);
        v[i][j] = (1/det)*(dxy[i][j]*dxz[i][j] - dxx[i][j]*dyz[i][j]);
        c[i][j] = 255.0;
    }
    }

return;

} /* lucas_kanade */

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
char   out1[80];             /* for reading data */
char   out2[80];             /* for reading data */
float  **f1, **f2;           /* images */
float  **dxx, **dxy, **dyy;  /* components of matrix */
float  **dxz, **dyz;         /* components of right-hand side */
float  **u, **v;             /* optic flow components */
float  **w, **c;             /* optic flow magnitude and confidence */
float  ***colour;            /* colour representation of optic flow */
long   i, j;                 /* loop variables */
long   nx, ny;               /* image size in x, y direction */ 
float  rho;                  /* integration scale */
float  eps;                  /* threshold for Lucas-Kanade method */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("OPTIC FLOW COMPUTATION WITH THE METHOD OF LUCAS AND KANADE\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert            \n");
printf ("    and 2007 by Stephan Didas                     \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image 1 (pgm):                      ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &f1);

printf ("input image 2 (pgm):                      ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &f2);


/* ---- read parameters ---- */

printf ("integration scale rho (>0) (float):       ");
read_float (&rho);

printf ("threshold epsilon (>0) (float):           ");
read_float (&eps);

printf ("output image (flow classification) (pgm): ");
read_string (out1);

printf ("output image (colour coded flow) (ppm):   ");
read_string (out2);
printf ("\n");


/* ---- allocate storage ---- */

alloc_matrix (&dxx, nx+2, ny+2);
alloc_matrix (&dxy, nx+2, ny+2);
alloc_matrix (&dyy, nx+2, ny+2);
alloc_matrix (&dxz, nx+2, ny+2);
alloc_matrix (&dyz, nx+2, ny+2);
alloc_matrix (&u, nx+2, ny+2);
alloc_matrix (&v, nx+2, ny+2);
alloc_matrix (&w, nx+2, ny+2);
alloc_matrix (&c, nx+2, ny+2);
alloc_cubix (&colour, nx+2, ny+2, 3);


/* ---- process images ---- */

create_eq_systems(f1, f2, nx, ny, 1.0, 1.0, 1.0, rho, dxx, dxy, dyy, 
                  dxz, dyz);
lucas_kanade(nx, ny, dxx, dxy, dyy, dxz, dyz, eps, c, u, v);


/* ---- calculate flow magnitude ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     w[i][j] = sqrt (u[i][j] * u[i][j] + v[i][j] * v[i][j]);
    
    
/* ---- analyse filtered image ---- */

analyse (w, nx, ny, &min, &max, &mean, &std);
printf ("flow magnitude image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write flow confidence image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# optic flow confidence, Lucas and Kanade method\n");
comment_line (comments, "# rho: %8.4f\n", rho);
comment_line (comments, "# threshold eps: %8.4f\n", eps);

/* write image */
write_pgm (c, nx, ny, out1, comments, 0, max);
printf ("output image %s successfully written\n\n", out1);

/* ---- write colour coded flow image (ppm format P6) ---- */

flow_to_color(u, v, colour, 4, nx, ny);

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# colour coded optic flow, Lucas and Kanade method\n");
comment_line (comments, "# rho: %8.4f\n", rho);
comment_line (comments, "# threshold eps: %8.4f\n", eps);

/* write image */
write_ppm (colour, nx, ny, out2, comments);
printf ("output image %s successfully written\n\n", out2);

/* ---- free memory  ---- */

disalloc_matrix (f1, nx+2, ny+2);
disalloc_matrix (f2, nx+2, ny+2);
disalloc_matrix (u,  nx+2, ny+2);
disalloc_matrix (v,  nx+2, ny+2);
disalloc_matrix (w,  nx+2, ny+2);
disalloc_matrix (dxx, nx+2, ny+2);
disalloc_matrix (dxy, nx+2, ny+2);
disalloc_matrix (dyy,  nx+2, ny+2);
disalloc_matrix (dxz,  nx+2, ny+2);
disalloc_matrix (dyz,  nx+2, ny+2);
disalloc_matrix (c,  nx+2, ny+2);
disalloc_cubix (colour, nx+2, ny+2, 3);

return(0);
}
