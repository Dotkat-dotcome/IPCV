#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                              QUANTISATION                                */
/*                                                                          */
/*  (Copyright by Markus Mainberger, 10/2009 and Joachim Weickert, 8/2014)  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 quantisation
*/


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

void quantisation 

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      long     q,         /* number of bits used to represent a value
                             in the output image */
      float    **u)       /* input: original image; output: quantised */
      
/*
 quantisation
*/

{
long    i, j;             /* loop variables */

/* quantise the input image */
float d = powf(2.0f,8-q);
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
    {
    /*
    SUPPLEMENT MISSING CODE HERE
    */
    u[i][j] = (floor((double)u[i][j]/(double)d) + 0.5)*d;
    }

return;

} /* quantisation */

/*--------------------------------------------------------------------------*/


void quantisation_with_noise

     (long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      long     q,         /* number of bits used to represent a value
                             in the output image */
      float    **u)       /* input: original image; output: quantised */
      
/*
 quantisation with uniformly distributed noise
*/

{
long    i, j;             /* loop variables */
double  n = 0.0f;	  /* variable for uniformly distributed noise */

/* quantise the input image */
float d = powf(2.0f,8-q);
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
    {
    /*
    SUPPLEMENT MISSING CODE HERE
    */
    //uniform_noise = ((float)rand()/RAND_MAX) - 0.5;
    u[i][j] = (floor((double)u[i][j]/(double)d + (((float)rand()/RAND_MAX) - 0.5)) + 0.5)*d;      
    }
    
return;

} /* quantisation_with_noise */

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
long   q;                    /* number of bits used to represent a value
                                in the output image */                                
long   flag;                 /* option flag */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("QUANTISATION\n\n");
printf ("***************************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert                         \n");
printf ("    and 2009 by Markus Mainberger                              \n");
printf ("    Dept. of Mathematics and Computer Science                  \n");
printf ("    Saarland University, Saarbruecken, Germany                 \n\n");
printf ("    All rights reserved. Unauthorized usage,                   \n");
printf ("    copying, hiring, and selling prohibited.                   \n\n");
printf ("    Send bug reports to                                        \n");
printf ("    weickert@mia.uni-saarland.de                               \n\n");
printf ("***************************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                                     ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &u);


/* ---- read parameters ---- */

printf ("\nnumber of bits to represent quantised value (integer): ");
read_long (&q);


printf ("\n\nyou have the following options:     \n"); 
printf ("\n    (1) quantisation without noise    ");
printf ("\n    (2) quantisation with uniformly   ");
printf ("\n        distributed noise             ");
printf ("\nchoose the quantisation option (1-2): ");

read_long (&flag);
printf("\n\n");

printf ("output image (pgm):                                    ");
read_string (out);
printf ("\n");


/* ---- quantise image ---- */
switch(flag)
{
  case 1 :
    /* perform quantisation without noise */
    quantisation (nx, ny, q, u);
    break;
  case 2 :
    /* perform quantisation with uniformly distributed noise */
    quantisation_with_noise (nx, ny, q, u);
    break;
  default :
    printf ("option (%ld) not available! \n\n\n",flag);
    disalloc_matrix (u, nx+2, ny+2);
    return(0);
}

/* ---- analyse filtered image ---- */

analyse (u, nx, ny, &min, &max, &mean, &std);
printf ("quantised image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# quantisation\n");
comment_line (comments, "# q: %8ld\n", q);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);

return(0);
}
