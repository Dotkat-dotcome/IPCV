#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                          RGB TO YCBCR CONVERSION                         */
/*                                                                          */
/*    (Copyright by Andres Bruhn, 11/2007 and Joachim Weickert, 8/2014)     */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 RGB to YCBCR conversion
*/


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

void read_pgm_or_ppm_and_allocate_memory

     (const char  *file_name,    /* name of pgm file */ 
      long        *nc,           /* number of colour channels */
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      float       ****u)         /* image, output */   

/* 
  reads a colour image that has been encoded in ppm format P5;
  allocates memory for the image u; 
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
FILE   *inimage;    /* input file */
char   row[80];     /* for reading data */
long   i, j, m;     /* loop variables */

/* open file */
inimage = fopen (file_name, "rb");
if (NULL == inimage) 
   {
   printf ("could not open file '%s' for reading, aborting.\n", file_name);
   exit (1);
   }

/* read header */
fgets (row, 80, inimage);                        /* image type: P5 or P6 */
if ((row[0]=='P') && (row[1]=='5'))
   *nc = 1;                                      /* P5: grey scale image */
else if ((row[0]=='P') && (row[1]=='6'))
   *nc = 3;                                      /* P6: colour image */
else
   {
   printf ("unknown image format");
   exit(0);
   }
fgets (row, 80, inimage);        
while (row[0]=='#')                /* skip comments */
      fgets (row, 80, inimage);
sscanf (row, "%ld %ld", nx, ny);   /* read image size */
fgets (row, 80, inimage);          /* read maximum grey value */

/* allocate memory */
alloc_cubix (u, (*nc), (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++) 
 for (i=1; i<=(*nx); i++) 
  for (m=0; m<=(*nc)-1; m++) 
      (*u)[m][i][j] = (float) getc(inimage);

/* close file */
fclose(inimage);

return;

} /* read_pgm_or_ppm_and_allocate_memory */

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

void write_pgm_or_ppm

     (float ***u,         /* colour image, unchanged */ 
      long  nc,           /* number of channels */
      long  nx,           /* size in x direction */
      long  ny,           /* size in y direction */
      char  *file_name,   /* name of ppm file */
      char  *comments)    /* comment string (set 0 for no comments) */

/* 
  writes an image into a pgm P5 (greyscale) or ppm P6 (colour) file;
*/

{
FILE           *outimage;  /* output file */
long           i, j, m;    /* loop variables */
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
if (nc == 1)
   fprintf (outimage, "P5\n");                  /* greyscale format */
else if (nc == 3)
   fprintf (outimage, "P6\n");                  /* colour format */
else
   {
   printf ("unsupported number of channels\n");
   exit (0);
   }
if (comments != 0)
   fprintf (outimage, comments);             /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
  for (m=0; m<=nc-1; m++)
     {
     aux = u[m][i][j] + 0.499999;    /* for correct rounding */
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

void RGB_to_YCbCr

     (float  ***u_RGB,     /* RGB image */
      float  ***u_YCbCr,   /* YCbCr image */  
      long   nx,           /* pixel number in x-direction */
      long   ny)           /* pixel number in y-direction */      

/*
  converts RGB to YCbCr colour space
*/

{
long    i, j;        /* loop variables */

/* computes the YCbCr values */ 
for (i=1;i<=nx;i++)
  for (j=1;j<=ny;j++)
      {
        /*
	 * SUPPLEMENT CODE HERE
	 */
      u_YCbCr[0][i][j] = 0.0 + 0.2990   *  u_RGB[0][i][j]
                   +   0.5870   * u_RGB[1][i][j]
                   +   0.1140 * u_RGB[2][i][j];

      u_YCbCr[1][i][j] = 127.5 + (-0.1687   *  u_RGB[0][i][j])
                   +   (-0.3313   * u_RGB[1][i][j])
                   +   0.5 * u_RGB[2][i][j];

      u_YCbCr[2][i][j] = 127.5 + 0.5   *  u_RGB[0][i][j]
                   +   (-0.4187   * u_RGB[1][i][j])
                   +   (-0.0813 * u_RGB[2][i][j]);
      }

return;

} /* RGB_to_YCbCr */


/*--------------------------------------------------------------------------*/

void YCbCr_to_RGB

     (float  ***u_YCbCr,   /* YCbCr image */  
      float  ***u_RGB,     /* RGB image */      
      long   nx,           /* pixel number in x-direction */
      long   ny)           /* pixel number in y-direction */      

/*
  converts YCbCr to RGB colour space
*/

{
long    i, j;        /* loop variables */

/* computes the RGB values */ 
for (i=1;i<=nx;i++)
  for (j=1;j<=ny;j++)
      {
      u_RGB[0][i][j] = 1.0   *  u_YCbCr[0][i][j]
                   +   0.0   * (u_YCbCr[1][i][j] - 127.5)
                   +   1.402 * (u_YCbCr[2][i][j] - 127.5);

      u_RGB[1][i][j] = 1.0   *  u_YCbCr[0][i][j]
                   -   0.344 * (u_YCbCr[1][i][j] - 127.5)
                   -   0.714 * (u_YCbCr[2][i][j] - 127.5);

      u_RGB[2][i][j] = 1.0   *  u_YCbCr[0][i][j]
                   +   1.773 * (u_YCbCr[1][i][j] - 127.5)
                   +   0.0   * (u_YCbCr[2][i][j] - 127.5);
      }

return;

} /* YCbCr_to_RGB */


/*--------------------------------------------------------------------------*/

void subsample_channel

     (float  **c,          /* image channel, changed */           
      long   nx,           /* pixel number in x-direction */
      long   ny,           /* pixel number in y-direction */    
      long   S)            /* subsample factor */


/*
  reduce resolution by averaging blocks of SxS neighbouring pixels
*/

{
long    i, j;        /* loop variables */
long    k, l;        /* loop variables */
float   sum;         /* summation variable */

/* replace SxS block by block average */ 
for (i=1;i<=nx;i+=S)
 for (j=1;j<=ny;j+=S)
     {
     /* initialise sum */
     sum = 0.0;

     /* compute block average */
     for (k=0;k<S;k++)
      for (l=0;l<S;l++)
          {
          sum += c[i+k][j+l];
          }
     sum=sum/(S*S);   

     /* set all block entries to average */
     for (k=0;k<S;k++)
      for (l=0;l<S;l++)
          {
          c[i+k][j+l] = sum;
          }
     }

return;

} /* subsample_channel */

/*--------------------------------------------------------------------------*/

void analyse_RGB

     (float   ***u,        /* image, unchanged */
      long    nc,          /* number of channels in the image */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *std)        /* standard deviation, output */

/*
 computes minimum, maximum, mean, and standard deviation of an RGB image u
*/

{
long    i, j, k;    /* loop variables */
double  help1;      /* auxiliary variable */
float   help2;      /* auxiliary variable */

*min  = u[0][1][1];
*max  = u[0][1][1];
help1 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
  for (k=0; k<nc; k++)
      {
      if (u[k][i][j] < *min) *min = u[k][i][j];
      if (u[k][i][j] > *max) *max = u[k][i][j];
      help1 = help1 + (double)u[k][i][j];
      }
*mean = (float)help1 / (nx * ny * nc);

*std = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
  for (k=0; k<nc; k++)
      {
      help2  = u[k][i][j] - *mean;
      *std = *std + help2 * help2;
      }
*std = sqrt(*std / (nx * ny * nc));

return;

} /* analyse_RGB */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  ***u_RGB;             /* RGB image */
float  ***u_YCbCr;           /* YCbCr image */
long   nx, ny;               /* image size in x, y direction */ 
long   nc;                   /* number of channels in the image */
long   S;                    /* subsampling factor */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("RGB TO YCBCR CONVERSION\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert            \n");
printf ("    and 2007 by Andres Bruhn                      \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (ppm format P6) ---- */

printf ("input image (ppm):                ");
read_string (in);
read_pgm_or_ppm_and_allocate_memory (in, &nc, &nx, &ny, &u_RGB);


/* ---- read parameters ---- */

printf ("subsampling factor: (integer):    ");
read_long (&S);

printf ("output image (ppm):               ");
read_string (out);
printf ("\n");


/* ---- analyse input image ---- */

analyse_RGB (u_RGB, nc, nx, ny, &min, &max, &mean, &std);
printf ("input image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- check if image can be downsampled by a factor of S ----*/

if ((nx%S!=0)||(ny%S!=0))
   {
   printf("\n\n Image size does not allow downsampling by factor %d! \n\n",S);
   return(0);
   }


/* ---- allocate storage for YCbCr image ---- */

alloc_cubix (&u_YCbCr, nc, nx+2, ny+2);


/* ---- process image ---- */

RGB_to_YCbCr(u_RGB, u_YCbCr, nx, ny);
subsample_channel(u_YCbCr[1], nx, ny, S);
subsample_channel(u_YCbCr[2], nx, ny, S);
YCbCr_to_RGB(u_YCbCr, u_RGB, nx, ny);


/* ---- analyse filtered image ---- */

analyse_RGB (u_RGB, nc, nx, ny, &min, &max, &mean, &std);
printf ("processed image:\n");
printf ("minimum:       %8.2f \n", min);
printf ("maximum:       %8.2f \n", max);
printf ("mean:          %8.2f \n", mean);
printf ("standard dev.: %8.2f \n\n", std);


/* ---- write output image (ppm format P6) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# RGB to YCbCr conversion\n");
comment_line (comments, "# S: %8ld\n", S);

/* write image */
write_pgm_or_ppm (u_RGB, nc, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_cubix (u_RGB,   nc, nx+2, ny+2);
disalloc_cubix (u_YCbCr, nc, nx+2, ny+2);

return(0);
}
