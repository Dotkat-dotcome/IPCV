#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*            MORPHOLOGY WITH A SQUARE-SHAPED STRUCTURING ELEMENT           */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 8/2014)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Osher-Sethian scheme (explicit, one-sided differences) 
 Satisfies extremum principle for time steps ht <= 0.5. 
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

void dilation 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Dilation with a square of size (2m + 1) * (2m + 1) as structuring element.
*/

{
long   i, j, k;    /* loop variables */
long   n;          /* max(nx,ny) */
float  max;        /* time saver */
float  *f;         /* auxiliary vector */


/* ---- allocate storage ---- */

if (nx > ny) n = nx; else n = ny;
alloc_vector (&f, n+2*m+1);


/* ---- dilation in x direction ---- */

for (j=1; j<=ny; j++)
    {
    /* copy row in vector with reflected boundary layer */
    for (i=1; i<=m; i++)
        f[i] = u[m+1-i][j];
    for (i=1; i<=nx; i++)
        f[m+i] = u[i][j];
    for (i=1; i<=m; i++)
        f[m+nx+i] = u[nx-i+1][j];
 
    /* perform dilation for each pixel */
    for (i=1; i<=nx; i++)  
        {
        max = f[i];
        for (k=i+1; k<=i+2*m; k++)
            if (f[k] > max) max = f[k];
        u[i][j] = max; 
        }
    }


/* ---- dilation in y direction ---- */

for (i=1; i<=nx; i++)
    {
    /* copy column in vector with refelcted boundary layer */
    for (j=1; j<=m; j++)
        f[j] = u[i][m+1-j];
    for (j=1; j<=ny; j++)
        f[m+j] = u[i][j];
    for (j=1; j<=m; j++)
        f[m+ny+j] = u[i][ny-j+1];

    /* perform dilation for each pixel */
    for (j=1; j<=ny; j++) 
        {
        max = f[j];
        for (k=j+1; k<=j+2*m; k++)
            if (f[k] > max) max = f[k];
        u[i][j] = max;
        }
    }


/* ---- disallocate storage for f ---- */

disalloc_vector (f, n+2*m+1);
return;

} /* dilation */

/*--------------------------------------------------------------------------*/

void erosion 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Erosion with a square of size (2m + 1) * (2m + 1) as structuring element.
*/

{
long   i, j, k;    /* loop variables */
long   n;          /* max(nx,ny) */
float  min;        /* time saver */
float  *f;         /* auxiliary vector */


/* ---- allocate storage ---- */

if (nx > ny) n = nx; else n = ny;
alloc_vector (&f, n+2*m+1);


/* ---- erosion in x direction ---- */

for (j=1; j<=ny; j++)
    {
    /* copy row in vector with reflected boundary layer */
    for (i=1; i<=m; i++)
        f[i] = u[m+1-i][j];
    for (i=1; i<=nx; i++)
        f[m+i] = u[i][j];
    for (i=1; i<=m; i++)
        f[m+nx+i] = u[nx-i+1][j];
 
    /* perform erosion for each pixel */
    for (i=1; i<=nx; i++)  
        {
        min = f[i];
        for (k=i+1; k<=i+2*m; k++)
            if (f[k] < min) min = f[k];
        u[i][j] = min; 
        }
    }


/* ---- erosion in y direction ---- */

for (i=1; i<=nx; i++)
    {
    /* copy column in vector with refelcted boundary layer */
    for (j=1; j<=m; j++)
        f[j] = u[i][m+1-j];
    for (j=1; j<=ny; j++)
        f[m+j] = u[i][j];
    for (j=1; j<=m; j++)
        f[m+ny+j] = u[i][ny-j+1];

    /* perform dilation for each pixel */
    for (j=1; j<=ny; j++) 
        {
        min = f[j];
        for (k=j+1; k<=j+2*m; k++)
            if (f[k] < min) min = f[k];
        u[i][j] = min;
        }
    }


/* ---- disallocate storage for f ---- */

disalloc_vector (f, n+2*m+1);
return;

} /* erosion */

/*--------------------------------------------------------------------------*/

void closing 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Closing with a square of size (2m + 1) * (2m + 1) as structuring element.
*/

{
    dilation(nx, ny, m, u);
    erosion(nx, ny, m, u);
    
/*
 INSERT CODE
*/
return;
} 

/*--------------------------------------------------------------------------*/

void opening 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Opening with a square of size (2m + 1) * (2m + 1) as structuring element.
*/

{
    erosion(nx, ny, m, u);
    dilation(nx, ny, m, u);

/*
 INSERT CODE
*/
return;
} 

/*--------------------------------------------------------------------------*/

void white_top_hat 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 White top hat with a square of size (2m + 1) * (2m + 1) as structuring 
 element.
*/

{
long   i, j;       /* loop variables */
float  **uo;       /* opening of input image */

alloc_matrix (&uo, nx+2, ny+2);


for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        uo[i][j] = u[i][j];
    }
}

opening(nx, ny, m, uo);

for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        u[i][j] = u[i][j] - uo[i][j];
    }
}

/*
 INSERT CODE
*/

disalloc_matrix (uo, nx+2, ny+2);
return;

} /* white_top_hat */ 

/*--------------------------------------------------------------------------*/

void black_top_hat 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Black top hat with a square of size (2m + 1) * (2m + 1) as structuring 
 element.
*/

{
long   i, j;       /* loop variables */
float  **uc;       /* closing of input image */

alloc_matrix (&uc, nx+2, ny+2);

for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        uc[i][j] = u[i][j];
    }
}


closing(nx, ny, m, uc);

for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        u[i][j] = uc[i][j] - u[i][j];
    }
}
/*
 INSERT CODE
*/

disalloc_matrix (uc, nx+2, ny+2);
return;

} /* black_top_hat */ 

/*--------------------------------------------------------------------------*/

void selfdual_top_hat 

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of structuring element: (2m+1)*(2m+1) */
      float   **u)        /* input: original image; output: processed */

/*
 Selfdual top hat with a square of size (2m + 1) * (2m + 1) as structuring 
 element.
*/

{
long   i, j;       /* loop variables */
float  **uc;       /* copy of input image */
alloc_matrix (&uc, nx+2, ny+2);

for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        uc[i][j] = u[i][j];
    }
}


white_top_hat(nx, ny, m, u);
black_top_hat(nx, ny, m, uc);



for(int i = 1; i <= nx; i++)
{
    for(int j = 1; j <= ny; j++)
    {
        u[i][j] = u[i][j] + uc[i][j];
    }
}

/*
 INSERT CODE
*/

disalloc_matrix (uc, nx+2, ny+2);
return;

} /* selfdual_top_hat */ 

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
long   goal;                 /* filter type */
long   m;                    /* size of structuring element: (2m+1)*(2m+1) */
float  factor;               /* enhancement factor */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("MORPHOLOGY WITH SQUARE-SHAPED STRUCTURING ELEMENT\n\n");
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

printf ("input image (pgm):                          ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("\n");
printf ("available filters:                        \n");
printf ("  (0) dilation                            \n");
printf ("  (1) erosion                             \n");
printf ("  (2) closing                             \n");
printf ("  (3) opening                             \n");
printf ("  (4) white top hat                       \n");
printf ("  (5) black top hat                       \n");
printf ("  (6) selfdual top hat                    \n");
printf ("your choice:                                ");
read_long (&goal);

printf ("\n");
printf ("size m of (2m+1) * (2m+1) mask (integer):   ");
read_long (&m);

printf ("output image (pgm):                         ");
read_string (out);
printf ("\n");


/* ---- process image ---- */

if (goal == 0) 
   dilation (nx, ny, m, u);
if (goal == 1) 
   erosion (nx, ny, m, u);
if (goal == 2)  
   closing (nx, ny, m, u);
if (goal == 3)  
   opening (nx, ny, m, u);
if (goal == 4)  
   white_top_hat (nx, ny, m, u);
if (goal == 5)  
   black_top_hat (nx, ny, m, u);
if (goal == 6)  
   selfdual_top_hat (nx, ny, m, u);


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
if (goal == 0) comment_line (comments, "# dilation\n");
if (goal == 1) comment_line (comments, "# erosion\n");
if (goal == 2) comment_line (comments, "# closing\n");
if (goal == 3) comment_line (comments, "# opening\n");
if (goal == 4) comment_line (comments, "# white top hat\n");
if (goal == 5) comment_line (comments, "# black top hat\n");
if (goal == 6) comment_line (comments, "# selfdual top hat\n");
comment_line (comments, "# struct. element: square\n");
comment_line (comments, "# size:            %8ld * %8ld\n", 2*m+1, 2*m+1);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);

return(0);
}
