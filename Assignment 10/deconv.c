#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                   DECONVOLUTION VIA THE FOURIER DOMAIN                   */
/*                                                                          */
/*         (Copyright by Martin Welk and Joachim Weickert, 8/2014)          */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - DFT or FFT, whatever is possible 
 - reflecting extension at the boundaries 
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

void filter

     (float    **ur,        /* real part of Fourier coeffs, changed */
      float    **ui,        /* imag. part of Fourier coeffs, changed */
      float    **hr,        /* real part of Fourier kernel, unchanged */
      float    **hi,        /* imag. part of Fourier kernel, unchanged */
      float    param,       /* filter parameter */
      long     nx,          /* pixel number in x direction */
      long     ny)          /* pixel number in y direction */

/*
 allows to filter the Fourier coefficients
*/

{
long   i, j;              /* loop variables */
float  N;                 /* denominator */
float  vr,vi;             /* auxiliary variables for cplx arithm. */


/* ---- compute filtered coefficients ---- */

/* SUPPLEMENT CODE HERE */
for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            {
                vr = ((ur[i][j]*hr[i][j]) + (ui[i][j]*hi[i][j]))/((hr[i][j]*hr[i][j]) + (hi[i][j]*hi[i][j]) + param);
                vi = ((ui[i][j]*hr[i][j]) - (ur[i][j]*hi[i][j]))/((hr[i][j]*hr[i][j]) + (hi[i][j]*hi[i][j]) + param);
                ur[i][j] = vr;
                ui[i][j] = vi;
            }
        }
    }

return;
}

/*--------------------------------------------------------------------------*/

long mylog2 

     (long n)               /* should be positive */

/*
 returns ld(n) if n is a power of 2; 
 returns -1    in all other cases. 
*/

{
long  ld;     /* ld(n) */
long  m;      /* auxiliary variable */

if (n <= 0)
   ld = -1;
else if (n == 1)
   ld = 0;
else
   {
   ld = 0;
   m  = 1; 
   do {
      m = 2 * m; 
      ld = ld + 1;
      }
   while (m < n);
   if (m > n) ld = -1;
   }

return (ld);
}

/*--------------------------------------------------------------------------*/

void FFT 

     (float    *vr,         /* real part of signal / Fourier coeff. */
      float    *vi,         /* imaginary part of signal / Fourier coeff. */
      long     n)           /* signal length, has to be power of 2 */ 

/*
 Fast Fourier Transform of a (complex) 1-D signal. 
 Based on the description in the Bronstein book.
 The signal length has to be a power of 2.
*/

{
const    float fa = sqrt(0.5);    /* frequently used renormalization factor */
long     p, q, r, j, k;           /* variables for indices and sizes */
long     nh, qq, qh;              /* variables for indices and sizes */
long     jq, jqh, jv, jn;         /* variables for indices and sizes */
long     logn;                    /* ld(n) */
long     m;                       /* auxiliary variable */
float    rh, ih, ch, chih;        /* for intermediate results */
float    *scrr, *scri, *exh;      /* auxiliary vectors */
float*   srr;                     /* point at source arrays, real part */ 
float*   sri;                     /* point at source arrays, imag. part */ 
float*   der;                     /* point at dest. arrays, real part */ 
float*   dei;                     /* point at dest. arrays, imag. part */
float*   swpp;                    /* used for pointer swapping */


/* ---- memory allocations ---- */

alloc_vector (&scrr, n);
alloc_vector (&scri, n);
alloc_vector (&exh,  n);


/* ---- initializations ----*/

/* init pointers */
srr = vr; 
sri = vi; 
der = scrr; 
dei = scri; 

/* logn := ld(n) */
m = n;
logn = -1;
while (m >= 1)
      {
      m >>= 1;              /* m = m / 2 */
      logn++;               
      }

/* init sizes */
nh = n>>1;                  /* n / 2 */ 
qh = nh>>1;                 /* n / 4 */ 

/* trigonometric values */
exh[0]  = 1.0;   exh[nh]    = 0.0;    /* exp (2 pi i 0.0 ) */
exh[qh] = 0.0;   exh[nh+qh] = 1.0;    /* exp (2 pi i 0.25) */
ch = -1.0;                            /* cos pi */
/* further trigonometric values will be computed by recursion */

/* other initializations */
qq = n; 
q  = nh; 
qh = q>>1; 
p  = 1;


/* ---- loop through all levels ----*/

for (r=logn; r>=1; r--) 
    {
    /* iterate through levels */
    if (r < logn) 
       ch = sqrt (0.5 * (1.0 + ch));    /* recursion for cosines */
    for (j=0; j<p; j++) 
        {         
        /* iterate through columns */
        jq = j * qq; 
        jqh = jq >> 1;
        if ((j&1==1) && (r<logn-1)) 
           {       
           /* recursion for exp(i*angle) */
           chih = 0.5 / ch;                    /* cosine inverse half */
           jv = jqh - q;                       
           jn = jqh + q; 
           if (jn == nh) 
              {                   
              /* use half-angle formula for exp */
              exh[jqh]    = (exh[jv] - 1.0) * chih;
              exh[jqh+nh] = exh[jv+nh] * chih;
              }
           else 
              {
              exh[jqh]    = (exh[jv]    + exh[jn]   ) * chih;
              exh[jqh+nh] = (exh[jv+nh] + exh[jn+nh]) * chih;
              }
           } /* if */
        for (k=0; k<q; k++) 
            {               
            /* iterate through rows */
            rh =  exh[jqh]    * srr[jq+k+q] + exh[jqh+nh] * sri[jq+k+q];
            ih = -exh[jqh+nh] * srr[jq+k+q] + exh[jqh]    * sri[jq+k+q];
            der[jqh+k]    = fa * (srr[jq+k] + rh);
            dei[jqh+k]    = fa * (sri[jq+k] + ih);
            der[jqh+nh+k] = fa * (srr[jq+k] - rh);
            dei[jqh+nh+k] = fa * (sri[jq+k] - ih);
            }
        } /* for j */
        
        /* swap array pointers */
        swpp = srr;   srr = der;   der = swpp;        
        swpp = sri;   sri = dei;   dei = swpp;

        /* adjust sizes */ 
        qq = q; 
        q = qh; 
        qh >>= 1; 
        p <<= 1;          
    }  /* for r */

/* copy f^ back, if ld(n) is odd */
if (logn&1 == 1)                             
   for (j=0; j<n; j++) 
       {                 
       der[j] = srr[j]; 
       dei[j] = sri[j]; 
       }


/* ---- disallocate memory ----*/

disalloc_vector (scrr, n);
disalloc_vector (scri, n);
disalloc_vector (exh,  n);

return;

} /* FFT */
  
/*--------------------------------------------------------------------------*/

void DFT  

     (float    *vr,         /* real part of signal / Fourier coeff. */
      float    *vi,         /* imaginary part of signal / Fourier coeff. */
      long     n)           /* signal length (>0) */ 


/* 
 Discrete Fourier transform of a (complex) 1-D signal. Quadratic complexity.
 Does not require powers of 2 as signal length.
*/

{
long    i, j;              /* loop variables */
float   help1, help2;      /* time savers */
float   help3, c, s;       /* time savers */
float   *fr, *fi;          /* auxiliary vectors (real / imaginary part) */
     
 
/* ---- allocate storage ---- */

alloc_vector (&fr, n);
alloc_vector (&fi, n);


/* ---- copy (vr,vi) into (fr,fi) ---- */

for (i=0; i<=n-1; i++)
    {
    fr[i] = vr[i];
    fi[i] = vi[i];
    }


/* ---- time savers ---- */

help1 = - 2.0 * 3.1415927 / (float)n;
help2 = 1.0 / sqrt ((float)n);

 
/* ---- perform DFT ---- */

for (i=0; i<=n-1; i++)
    {
    vr[i] = 0.0;
    vi[i] = 0.0;
    for (j=0; j<=n-1; j++)
        {
        help3 = help1 * i * j;
        c     = cos (help3);
        s     = sin (help3);
        vr[i] = vr[i] + fr[j] * c - fi[j] * s;
        vi[i] = vi[i] + fi[j] * c + fr[j] * s;
        }
    vr[i] = vr[i] * help2;
    vi[i] = vi[i] * help2;
    }


/* ---- disallocate storage ---- */

disalloc_vector (fr, n);
disalloc_vector (fi, n);
return;

} /* DFT */

/* ---------------------------------------------------------------------- */

void FT2D  

     (float    **ur,        /* real part of image / Fourier coeff. */
      float    **ui,        /* imaginary part of image / Fourier coeff. */
      long     nx,          /* pixel number in x direction */ 
      long     ny)          /* pixel number in y direction */ 


/* 
 Two-dimensional discrete Fourier transform of a (complex) image.
 This algorithm exploits the separability of the Fourier transform. 
 Uses FFT when the pixel numbers are powers of 2, DFT otherwise.
*/


{
long   i, j;              /* loop variables */
long   n;                 /* max (nx, ny) */
long   logn;              /* ld(n) */
float  *vr, *vi;          /* real / imaginary signal or Fourier data */


/* ---- allocate auxiliary vectors vr, vi ---- */

if (nx > ny) 
   n = nx; 
else 
   n = ny;

alloc_vector (&vr, n);
alloc_vector (&vi, n);


/* ---- transform along x direction ---- */

logn = mylog2 (nx);
for (j=0; j<=ny-1; j++)
    {
    /* write in 1-D vector; exclude boundary pixels of ur and ui */
    for (i=0; i<=nx-1; i++)
        {
        vr[i] = ur[i+1][j+1];
        vi[i] = ui[i+1][j+1];
        }

    /* apply Fourier transform */
    if (logn == -1)
       /* nx is not a power of 2; use DFT */
       DFT (vr, vi, nx); 
    else
       /* nx is a power of 2; use FFT */
       FFT (vr, vi, nx); 

    /* write back in 2-D image; include boundary pixels */
    for (i=0; i<=nx-1; i++)
        {
        ur[i+1][j+1] = vr[i];
        ui[i+1][j+1] = vi[i];
        }
    }


/* ---- transform along y direction ---- */

logn = mylog2 (ny);
for (i=0; i<=nx-1; i++)
    {
    /* write in 1-D vector; exclude boundary pixels of ur and ui */
    for (j=0; j<=ny-1; j++)
        {
        vr[j] = ur[i+1][j+1];
        vi[j] = ui[i+1][j+1];
        }

    /* apply Fourier transform */
    if (logn == -1) 
       /* ny is not a power of 2; use DFT */
       DFT (vr, vi, ny);
    else
       /* ny is a power of 2; use FFT */
       FFT (vr, vi, ny); 

    /* write back in 2-D image; include boundary pixels */
    for (j=0; j<=ny-1; j++)
        {
        ur[i+1][j+1] = vr[j];
        ui[i+1][j+1] = vi[j];
        }
    }


/* ---- disallocate storage ---- */

disalloc_vector (vr, n);
disalloc_vector (vi, n);

return;

} /* FT2D */

/* ---------------------------------------------------------------------- */

void periodic_shift  

     (float    **u,         /* image, changed */
      long     nx,          /* pixel number in x direction */ 
      long     ny,          /* pixel number in y direction */
      long     xshift,      /* shift in x direction */ 
      long     yshift)      /* shift in y direction */ 

/*
 shifts an image u by the translation vector (xshift,yshift) 
 with 0 <= xshift <= nx-1 and 0 <= yshift <= ny-1.
*/

{
long    i, j;         /* loop variables */
float   **f;          /* auxiliary image */

/* allocate storage */
alloc_matrix (&f, nx+2, ny+2);

/* shift in x direction */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (i-xshift >= 1)
        f[i][j] = u[i-xshift][j];
     else 
        f[i][j] = u[i+nx-xshift][j];

/* shift in y direction */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (j-yshift >= 1)
        u[i][j] = f[i][j-yshift];
     else 
        u[i][j] = f[i][j+ny-yshift];

/* disallocate storage */
disalloc_matrix (f, nx+2, ny+2);

return;
}

/* ---------------------------------------------------------------------- */

void FT_deconv  

     (float    **ur,        /* input image, changed */
      float    **hr,        /* kernel image, changed */
      float    param,       /* parameter for filtering */
      long     nx,          /* pixel number in x direction, unchanged */ 
      long     ny)          /* pixel number in y direction, unchanged */ 
      

/*
 computes the Gaussian derivative of order (xord, yord) and writes
 it back into **ur
*/


{
long   i, j;           /* loop variables */
float  **ui;           /* imaginary image or Fourier data (image)  */
float  **hi;           /* imaginary image or Fourier data (kernel) */
float  hweight;        /* weight of kernel, for normalisation */


/* ---- allocate storage ---- */

alloc_matrix (&ui, nx+2, ny+2);
alloc_matrix (&hi, nx+2, ny+2);


/* ---- initialise imaginary part with 0 ---- */

for (j=0; j<=ny+1; j++)
 for (i=0; i<=nx+1; i++) {
     ui[i][j] = 0.0;
     hi[i][j] = 0.0;
 }


/* ---- compute discrete Fourier transformation ---- */

printf ("computing Fourier transformation of image...\n");
FT2D (ur, ui, nx, ny);

printf ("computing Fourier transformation of kernel...\n");
FT2D (hr, hi, nx, ny);

/* normalise kernel */
hweight=hr[1][1];
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++) {
     hr[i][j] = hr[i][j] / hweight;
     hi[i][j] = hi[i][j] / hweight;
 }


/* ---- shift lowest frequency in the centre ----*/

/* periodic_shift (ur, nx, ny, nx/2, ny/2); */
/* periodic_shift (ui, nx, ny, nx/2, ny/2); */


/* ---- filter the Fourier coefficients ---- */

printf ("filtering the Fourier coefficients...\n");
filter (ur, ui, hr, hi, param, nx, ny);

/* ---- shift lowest frequency back to the corners ----*/

/* periodic_shift (ur, nx, ny, nx-nx/2, ny-ny/2); */
/* periodic_shift (ui, nx, ny, nx-nx/2, ny-ny/2); */


/* ---- compute discrete Fourier backtransformation ---- */

printf ("computing Fourier backtransformation...\n\n");

/* backtransformation = DFT of complex conjugated Fourier coefficients */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     ui[i][j] = - ui[i][j];
FT2D (ur, ui, nx, ny);


/* ---- free storage ---- */

disalloc_matrix (ui, nx+2, ny+2);
disalloc_matrix (hi, nx+2, ny+2);

return;

} /* FT_deconv */

/* ---------------------------------------------------------------------- */

void rescale
 
     (float   **u,         /* image, changed */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   u1,          /* grey value of u */
      float   u2,          /* grey value of u */
      float   v1,          /* corresponding grey value of v */
      float   v2)          /* corresponding grey value of v */

/*
 rescales image u such that the grey value u1 is mapped into v1,
 and u2 into v2.
*/

{
long    i, j;       /* loop variables */
float   m, b;       /* parameters of the transformation */

m = (v2 - v1) / (u2 - u1);
b = v1 - m * u1;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = m * u[i][j] + b;

return;

} /* rescale */

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
float  **ur;                 /* image data */
float  **hr;                 /* kernel data */
float  **hr_tmp;             /* helping data to hr */
long   i, j;                 /* loop variables */
long   nx, ny;               /* image size in x, y direction */ 
long   hnx, hny;             /* kernel image size in x, y direction */
float  param;                /* filter parameter */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("DECONVOLUTION IN THE FOURIER DOMAIN\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by                             \n");
printf ("    Joachim Weickert and Martin Welk              \n");
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
read_pgm_and_allocate_memory (in, &nx, &ny, &ur);

printf ("input kernel (pgm):               ");
read_string (in);
read_pgm_and_allocate_memory (in, &hnx, &hny, &hr_tmp);


/* ---- initialise hr ---- */

alloc_matrix (&hr, nx+2, ny+2);

for (j=0; j<=ny+1; j++)
 for (i=0; i<=nx+1; i++) 
  hr[i][j]=0.0;
 
for (j=1; j<hny; j++)
 for (i=1; i<hnx; i++)
     if (i<=nx && j<=ny)
        hr[i][j] = hr_tmp[i][j];


/* ---- read parameters ---- */

printf ("filter parameter (float):         ");
read_float (&param);

printf ("output image (pgm):               ");
read_string (out);
printf ("\n");


/* ---- shift kernel to the corner ----*/

periodic_shift (hr, nx, ny, nx-hnx/2, ny-hny/2);


/* ---- compute filter in the Fourier domain ---- */

FT_deconv (ur, hr, param, nx, ny);


/* ---- analyse filtered image ---- */

analyse (ur, nx, ny, &min, &max, &mean, &std);
printf ("filtered image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# deconvolution via the Fourier domain\n");
comment_line (comments, "# param: %8.4f\n", param);

/* write image */
write_pgm (ur, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (ur,     nx+2,  ny+2);
disalloc_matrix (hr_tmp, hnx+2, hny+2);
disalloc_matrix (hr,     nx+2,  ny+2);

return(0);
}
