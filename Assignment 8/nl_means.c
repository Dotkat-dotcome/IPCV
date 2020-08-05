#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                    NL-MEANS AND GAUSSIAN CONVOLUTION                     */
/*                                                                          */
/* (Copyright by Sven Grewenig, Pascal Peter and Joachim Weickert, 8/2014)  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 consists of:
 - Gaussian convolution
 - NL-Means filtering
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

     (double  **u,          /* image, unchanged */ 
      long    nx,           /* image size in x direction */
      long    ny,           /* image size in y direction */
      char    *file_name,   /* name of pgm file */
      char    *comments)    /* comment string (set 0 for no comments) */

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

void gauss_conv 

     (double   sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      double   hx,        /* pixel size in x direction */
      double   hy,        /* pixel size in y direction */
      double   precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      double   **f)       /* input: original image ;  output: smoothed */

/*  Gaussian convolution. */

{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
double  sum;                  /* for summing up */
double  *conv;                /* convolution vector */
double  *help;                /* row or column with dummy boundaries */
      

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

void smooth

     (double  **f,        /* original image */
      double  **u,        /* smoothed image */
      long    nx,         /* size in x-direction */
      long    ny,         /* size in y-direction */
      double  sigma)      /* std. dev. of Gaussian */
      
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

void dummies
 
     (double   **u,       /* image matrix */
      long     nx,        /* size in x direction */
      long     ny,        /* size in y direction */
      long     p)         /* dummy-size */     

/* creates dummy boundaries of size p by mirroring */

{
long  i, j;       /* loop variables */
long  imax,jmax;  /* maximal indices */

imax = nx + 2 * p - 1;
jmax = ny + 2 * p - 1;

for (j=0; j<=jmax; j++)
    {
    for (i=0; i<p; i++)
        {
        if (j < p)
           u[i][j] = u[2 * p - 1 - i][2 * p - 1 - j];
        else if ((j >= p) && (j <= ny+p-1))
           u[i][j] = u[2 * p - 1 - i][j];
        else
           u[i][j] = u[2 * p - 1 - i][jmax + ny - j];
        }
    for (i=imax-p+1; i<=imax; i++)
        {
        if (j < p)
           u[i][j] = u[imax + nx - i][2 * p - 1 - j];
        else if ((j >= p) && (j <= ny+p-1))
           u[i][j] = u[imax + nx - i][j];
        else
           u[i][j] = u[imax + nx - i][jmax + ny - j];
        }
    }

for (i=p; i<=nx+p-1; i++)
    {
    for (j=0; j<p; j++)
        u[i][j] = u[i][2 * p - 1 - j];
    for (j=ny+p; j<=jmax; j++)
        u[i][j] = u[i][jmax + ny - j];
    } 

return;

} /* dummies */

/*--------------------------------------------------------------------------*/

void nl_means 

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      long     m,         /* patch size */
      double   sigma,     /* std. dev. */
      long     t,         /* search window size */
      double   **u)       /* input: original image; output: processed */

/*
  performs NL-Means filtering
*/

{
long    i, j, k, l, p, r;          /* loop variables */
long    imax, jmax;                /* extended image boundaries */
double  weight, weightsum, temp;   /* time saving variables */
double  **f;                       /* work copy of image u */
long    rmin, rmax, smin, smax;    /* boundaries of search window */

/* create reflecting dummy boundaries */
dummies (u, nx, ny, m);

imax = nx + m - 1;
jmax = ny + m - 1;

alloc_matrix (&f, nx+2*m, ny+2*m);

/* copy u into f  */
for (i=0; i<nx+2*m; i++)
 for (j=0; j<ny+2*m; j++)
     f[i][j]=u[i][j];


/* ---- NL-Means Filter ---- */

/* iterate over all pixels u[i][j] */
for (i=m; i<=imax; i++) 
    {
    for (j=m; j<=jmax; j++) 
        {
        u[i][j]   = 0.0;
        weightsum = 0.0;

        /* compute boundaries of the search window */
        if (i-t < m)     rmin = m;
        else             rmin = i-t;
        if (i+t > imax)  rmax = imax;
        else             rmax = i+t;
        if (j-t < m)     smin = m;
        else             smin = j-t;
        if (j+t > jmax)  smax = jmax;
        else             smax = j+t;

        /* for each pixel in the search window... */
        for (k=rmin; k<=rmax; k++) 
            {
            for (l=smin; l<=smax; l++) 
                {
                temp=0.0;
                
                /* compute patch distance (squared Euclidean norm) */
                for (p=-m; p<=m; p++) 
                 for (r=-m; r<=m; r++) 
                     temp+=(f[i-p][j-r]-f[k-p][l-r])*(f[i-p][j-r]-f[k-p][l-r]);
            
                /* compute the corresponding Gaussian weight. note that temp
                   already contains the squared Euclidean norm */
                weight = exp (-temp / (2 * sigma * sigma));
                     
                /* averaging */
                u[i][j]   += weight * f[k][l];
                weightsum += weight;
                }  
             }

         /* normalise result */
         if (weightsum > 0)
            u[i][j]=u[i][j]/weightsum;
         else
            u[i][j]=f[i][j];
         }
    }

disalloc_matrix (f, nx+2*m, ny+2*m);

return;

} /* nl_means */

/*--------------------------------------------------------------------------*/

void analyse

     (double  **u,         /* image, unchanged */
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
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **u;                  /* filtered image */
double  **f;                  /* original image */
double  **tmp;                /* work copy of f */
double  **method_noise;       /* shifted difference between f and u */
long    i, j, k, l;           /* loop variable */ 
long    nx, ny;               /* image size in x, y direction */ 
long    m;                    /* patch radius */
long    t;                    /* search window size */
double  sigma1;               /* std. deviation sigma for nl-means */
double  sigma2;               /* std. deviation sigma for Gaussian */
long    imax, jmax;           /* help indices */
long    goal;                 /* choice variable for different algorithms */
float   max, min;             /* largest, smallest grey value */
float   mean;                 /* average grey value */
float   std;                  /* standard deviation */
char    comments[1600];       /* string for comments */
long    length;               /* string length */

printf ("\n");
printf ("NL-MEANS AND GAUSSIAN CONVOLUTION\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert,           \n");
printf ("    Sven Grewenig and Pascal Peter                \n");
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

printf ("\n");
printf ("available algorithms:                 \n");
printf ("  (0) denoising with NL-means         \n");
printf ("  (1) denoising with Gauss-convolution\n");
printf ("your choice:                            ");
read_long (&goal);
printf ("\n");

/* ---- read parameters ---- */

if (goal == 0)
   {
   printf ("patch size m (integer):                 ");
   read_long (&m);

   printf ("search window size n (integer):         ");
   read_long (&t);

   printf ("std. dev. sigma (float):                ");
   read_double (&sigma1);
   }
else
   {
   printf ("std. dev. of Gaussian sigma (float):    ");
   read_double (&sigma2);
   m = 1;
   }

printf ("output image (pgm):                     ");
read_string (out);
printf ("\n");


/* ---- allocate storage for u and method_noise ---- */

alloc_matrix (&u,            nx+2*m, ny+2*m);
alloc_matrix (&method_noise, nx+2*m, ny+2*m);
alloc_matrix (&tmp,          nx+2*m, ny+2*m);


/* ---- initialisations ---- */

imax = nx + m - 1;
jmax = ny + m - 1;

for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     tmp[i][j] = u[i][j] = f[k][l];
  
  
/* ---- process image with your favourite filter ---- */

if (goal == 0) 
   {
   printf ("Applying NL-means filter with patch size %d, \n", m);
   printf ("window size %d and std. dev. %8.4lf\n\n", t, sigma1);
   nl_means (nx, ny, m, sigma1, t, u);
   } 
else 
   {
   printf ("Applying Gauss-Convolution with std. dev. %8.4lf\n\n", sigma2);
   smooth (tmp, u, nx, ny, sigma2);
   }

/* copy result to f */
for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     f[k][l] = u[i][j];


/* ---- analyse filtered image ---- */

analyse (f, nx, ny, &min, &max, &mean, &std);
printf ("filtered image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- compute method noise ---- */

for (j=m; j<=jmax; j++)
 for (i=m; i<=imax; i++)
     method_noise[i][j] = 127.5 + tmp[i][j] - u[i][j];


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
if (goal == 0)
   {
   comment_line (comments, "# NL-means\n");
   comment_line (comments, "# patch size m:         %8ld\n", m);
   comment_line (comments, "# search window size n: %8.4f\n", t);
   comment_line (comments, "# std. dev. sigma:      %8.4f\n", sigma1);
   }
else
   {
   comment_line (comments, "# Gaussian convolution\n");
   comment_line (comments, "# std. dev. sigma:      %8.4f\n", sigma2);
   }

/* write image */
write_pgm (f, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- write method noise image (pgm format P5) ---- */

length = strlen(out);
if (length < 5) 
   {
   printf("Error: output file name does not match the format filename.pgm\n");
   }
out[length-4] = '\0';
sprintf(out,"%s_method_noise.pgm",out);

/* copy result to f */
for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     f[k][l] = method_noise[i][j];

/* generate comment string */
comments[0]='\0';
if (goal == 0)
   {
   comment_line (comments, "# method noise of NL-means\n");
   comment_line (comments, "# patch size m:         %8ld\n", m);
   comment_line (comments, "# search window size n: %8.4f\n", t);
   comment_line (comments, "# std. dev. sigma:      %8.4f\n", sigma1);
   }
else
   {
   comment_line (comments, "# method noise of Gaussian convolution\n");
   comment_line (comments, "# std. dev. sigma:      %8.4f\n", sigma2);
   }

/* write image */
write_pgm (f, nx, ny, out, comments);
printf ("method noise %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (f,            nx+2,    ny+2);
disalloc_matrix (tmp,          nx+2*m, ny+2*m);
disalloc_matrix (u,            nx+2*m, ny+2*m);
disalloc_matrix (method_noise, nx+2*m, ny+2*m);

return(0);
}
