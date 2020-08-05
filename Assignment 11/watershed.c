#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                     TOBOGGAN WATERSHED SEGMENTATION                      */
/*                                                                          */
/*       (Copyright by Markus Mainberger and Joachim Weickert, 8/2014)      */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Toboggan Watershed Algorithm
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
      

/* ------------------------ convolution in x direction -------------------- */

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


/* ------------------------ convolution in y direction -------------------- */

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
 
     (float  **f,         /* original image */
      float  **u,         /* smoothed version */
      long   nx,          /* size in x-direction */
      long   ny,          /* size in y-direction */
      float  sigma)       /* std. dev. of Gaussian */

/* smoothes an image f (f has to have a border of 1 pixel) */

{
long  i, j;   /* loop variables */
 
/* copy image */
for (i=0; i<nx+2; ++i)
 for (j=0; j<ny+2; ++j)
     u[i][j] = f[i][j];
   
if (sigma > 0.0)
   gauss_conv (sigma, nx, ny, 1.0, 1.0, 5.0, 1, u);

return;
}

/* ------------------ Stack container using void pointers ------------------ */

typedef struct StackElement StackElement;

struct StackElement
{
    void         *data;
    StackElement *next;
};

typedef struct
{
    StackElement *start;
    StackElement *end;
    unsigned int size;
} Stack;

/*--------------------------------------------------------------------------*/

StackElement* new_StackElement

     (void *data)         /* data of the new stack element */
     
/* creates new stack element */
     
{
StackElement *e = (StackElement*) malloc (sizeof (StackElement));
e->next = NULL;
e->data = data;

return e;
}

/*--------------------------------------------------------------------------*/

void free_StackElement

     (StackElement *q)    /* stack element to free */
     
/* frees stuck element */
     
{
free(q->data);
free(q);

return;
}

/*--------------------------------------------------------------------------*/

Stack* new_Stack

     ()
     
/* creates new stack */
     
{
Stack *s = (Stack*) malloc (sizeof (Stack));
s->start = NULL;
s->end   = NULL;
s->size  = 0;

return s;
}

/*--------------------------------------------------------------------------*/

void free_Stack 
  
     (Stack *s)           /* stack to free */
     
/* frees stack */
     
{
StackElement *e = s->start;
while (e != NULL)
      {
      s->start = e->next;
      free_StackElement (e);
      e = s->start;
      }
free(s);

return;
}

/*--------------------------------------------------------------------------*/

void push_Stack

     (Stack *s,           /* input stack */
      void  *data)        /* data to push to stack */
      
/* pushes to stack */
      
{
if (s->end == NULL)
   {
   s->end = new_StackElement (data);
   s->start = s->end;
   }
else
   {
   StackElement *e = new_StackElement (data);
   e->next = s->start;
   s->start = e;
   }
s->size++;

return;
}

/*--------------------------------------------------------------------------*/

void pop_Stack

     (Stack* s)           /* input stack */
     
/* pops from stack */
     
{
if (s->start != NULL)
   {
   StackElement *e = s->start;
   if (s->start == s->end)
      s->end = NULL;
   s->start = e->next;
   free_StackElement (e);
   s->size--;
   }
   
return;
}

/*--------------------------------------------------------------------------*/

void* front_Stack

     (Stack* s)           /* input stack */
     
/* returns front of the stack*/
     
{
if (s->start == NULL)
   {
   printf("\nError: front_Stack(): Stack is empty! size: %i\n\n",s->size);
   exit(0);
   }

return (s->start->data);
}

/*--------------------------------------------------------------------------*/

void* back_Stack

     (Stack* s)           /* input stack */

/* returns back of the stack */

{
if (s->end == NULL)
   {
   printf("\nError: back_Stack(): Stack is empty!\n\n");
   exit(0);
   }

return (s->end->data);
}

/*--------------------------------------------------------------------------*/

long empty_Stack

     (Stack* s)           /* input stack */

/* empties stack */

{
return (s->size == 0);
}

/*--------------------------------------------------------------------------*/

typedef struct
{
    long  x;
    long  y;
} intPair;

/*--------------------------------------------------------------------------*/

inline intPair* new_intPair

    (long  x,             /* first element of int pair */
     long  y)             /* second element of int pair */

/* creates intPair */

{
intPair *d = (intPair*) malloc (sizeof (intPair));
d->x = x;
d->y = y;

return d;
}

/*--------------------------------------------------------------------------*/

void set_bounds

     (float  **f,         /* image */
      long   nx,          /* size in x direction */
      long   ny,          /* size in y direction */
      float  a)           /* set boundaries to a */

/* set boundaries of f to value a */

{
long  i, j;           /* loop variables  */

/* upper and lower boundary */
for (i=0; i<nx+2; i++)
    {
    f[i][0]    = a;
    f[i][nx+1] = a;
    }

/* left and right boundary */
for (j=0; j<ny+2; j++)
    {
    f[0][j]    = a;
    f[nx+1][j] = a;
    }

return;
}

/*--------------------------------------------------------------------------*/

void watershed 

     (long     nx,                  /* image dimension in x direction */
      long     ny,                  /* image dimension in y direction */
      float    sigma,               /* std. dev. of Gaussian */
      long     show_segment_lines,  /* draw white segment boundaries? */
      float    **f,                 /* input: original; output: processed */
      float    **u)                 /* output image */

/*
 toboggan watershed algorithm
*/

{
float  **gradmagn;           /* averaged gradient magnitude */
float  **visited;            /* log visited points */
long   i, j;                 /* loop variables */

/* memory allocations */
alloc_matrix(&gradmagn, nx+2, ny+2);
alloc_matrix(&visited,  nx+2, ny+2);

/* smooth image */
if (sigma > 0.0f)
   gauss_conv(sigma, nx, ny, 1, 1, 3, 1, f);

/* mirror boundaries */
dummies(f, nx, ny);

/* compute gradient magnitude using the sobel operator */
for (i=1; i<nx+1; ++i)
 for (j=1; j<ny+1; ++j)
     {
     float dx = (-        f[i-1][j-1] +        f[i+1][j-1]
                 - 2.0f * f[i-1][j]   + 2.0f * f[i+1][j]
                 -        f[i-1][j+1] +        f[i+1][j+1])
                 /8.0f;

     float dy = (-        f[i-1][j-1] -  2.0f * f[i][j-1]
                 -        f[i+1][j-1] +         f[i-1][j+1]
                 + 2.0f * f[i][j+1]   +         f[i+1][j+1])
                 /8.0f;

     gradmagn[i][j] = sqrtf(dx * dx + dy * dy);
     }

/* set boundaries of gradient magnitude to infinity */
set_bounds (gradmagn, nx, ny, FLT_MAX);

/* create an empty stack */
Stack *s = new_Stack ();

/* initialise visited-array (visited pixels will be marked with 1.0f) */
for (i=0; i<nx+2; ++i)
 for (j=0; j<ny+2; ++j)
     visited[i][j] = 0.0f;

 
/* ---- perform watershed segmentation ---- */

long  x, y;     /* coordinates giving the current position of a "water drop" */
float color;    /* colour found at minimum */
long  next_x;   /* x-coordinate of the minimum in the 8-neighbourhood */
long  next_y;   /* y-coordinate of the minimum in the 8-neighbourhood */
float min;      /* auxillary variable giving the minimum in the
                   8-neighbourhood of pixel (x,y) */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* if pixel was not visited yet */
     if (visited[i][j] != 1.0f)
        {
        /* put pixel on stack and mark as visited */
        x = i, y = j;
        push_Stack (s, new_intPair(x,y));
        visited[i][j] = 1.0f;

        /* while minimum not reached follow the water drop on its path
           down to minimum */
        long minimum_reached = 0;
        
         while (!minimum_reached)
            {

            if (   gradmagn[x][y]>gradmagn[x+1][y  ]
                || gradmagn[x][y]>gradmagn[x-1][y  ]
                || gradmagn[x][y]>gradmagn[x  ][y+1]
                || gradmagn[x][y]>gradmagn[x  ][y-1]
                || gradmagn[x][y]>gradmagn[x+1][y+1]
                || gradmagn[x][y]>gradmagn[x+1][y-1]
                || gradmagn[x][y]>gradmagn[x-1][y+1]
                || gradmagn[x][y]>gradmagn[x-1][y-1])
                {
                /* weight for diagonal distance */
                const float sqrt2 = 1.0f/sqrtf(2.0f);

                /* initialise with one neighbour */
                min = (gradmagn[x+1][y+1]-gradmagn[x][y])*sqrt2;
                next_x = x+1;
                next_y = y+1;

                /* check for steepest descent
                   and update next_x, next_y and min */
                
               /*
                 SUPPLEMENT CODE
               */
               if ((gradmagn[x+1][y]-gradmagn[x][y])<min)
               {
                   min = (gradmagn[x+1][y]-gradmagn[x][y]);
                   next_x = x+1;
                   next_y = y;
               }
               if ((gradmagn[x+1][y-1]-gradmagn[x][y])*sqrt2<min)
               {
                   min = (gradmagn[x+1][y-1]-gradmagn[x][y])*sqrt2;
                   next_x = x+1;
                   next_y = y-1;
               }
               if ((gradmagn[x][y+1]-gradmagn[x][y])<min)
               {
                   min = (gradmagn[x][y+1]-gradmagn[x][y]);
                   next_x = x;
                   next_y = y+1;
               }
               if ((gradmagn[x][y-1]-gradmagn[x][y])<min)
               {
                   min = (gradmagn[x][y-1]-gradmagn[x][y]);
                   next_x = x;
                   next_y = y-1;
               }
               if ((gradmagn[x-1][y+1]-gradmagn[x][y])*sqrt2<min)
               {
                   min = (gradmagn[x-1][y+1]-gradmagn[x][y])*sqrt2;
                   next_x = x-1;
                   next_y = y+1;
               }
               if ((gradmagn[x-1][y]-gradmagn[x][y])<min)
               {
                   min = (gradmagn[x-1][y]-gradmagn[x][y]);
                   next_x = x-1;
                   next_y = y;
               }
               if ((gradmagn[x-1][y-1]-gradmagn[x][y])*sqrt2<min)
               {
                   min = (gradmagn[x-1][y-1]-gradmagn[x][y])*sqrt2;
                   next_x = x-1;
                   next_y = y-1;
               }

                /* set new current coordinates of waterdrop */
                x = next_x;
                y = next_y;
                }

                if (visited[x][y] != 1.0f)
                  {     
                  /* if pixel was not visited yet, mark pixel as visited 
                     and put it on stack */
                  push_Stack(s,new_intPair(x,y));
                  visited[x][y] = 1.0f;
                  }
                else  
                  /* if pixel was already visited stop since we know then
                     the colour of the minimum */
                  minimum_reached = 1;
                }

        /* get colour of minimum */
        color = f[x][y];

        /* mark all pixels on stack using the colour of the minimum */
        while (!(empty_Stack(s)))
              {
              x = ((intPair*)front_Stack(s))->x;
              y = ((intPair*)front_Stack(s))->y;
              pop_Stack(s);
              f[x][y] = color;
              }
        }
     }


/* ---- prepare output image ---- */

for (i=1; i<nx+1; i++)
 for (j=1; j<ny+1; j++)
     {
     if (show_segment_lines && 
            (((i < nx) && ((int)f[i][j] != (int)f[i+1][j]))
         ||  ((j < ny) && ((int)f[i][j] != (int)f[i][j+1]))))
        u[i][j] =  255.0f;
     else
        u[i][j] = f[i][j];
     }
    
/* free memory */
free_Stack(s);
disalloc_matrix (gradmagn, nx+2, ny+2);
disalloc_matrix (visited,  nx+2, ny+2);

return;

} /* watershed */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **f;                  /* image */
float  **u;                  /* image */
long   nx, ny;               /* image size in x, y direction */ 
float  sigma;                /* Gaussian standard deviation */
long   show_segment_lines;   /* draw white segment boundaries? */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("TOBOGGAN WATERSHED SEGMENTATION\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by                             \n");
printf ("    Joachim Weickert and Markus Mainberger        \n");
printf ("    Dept. of Mathematics and Computer Science     \n");
printf ("    Saarland University, Saarbruecken, Germany    \n\n");
printf ("    All rights reserved. Unauthorized usage,      \n");
printf ("    copying, hiring, and selling prohibited.      \n\n");
printf ("    Send bug reports to                           \n");
printf ("    weickert@mia.uni-saarland.de                  \n\n");
printf ("**************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                   ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &f);
alloc_matrix (&u, nx+2, ny+2);

/* ---- read parameters ---- */

printf ("sigma (presmoothing) (float):        ");
read_float (&sigma);

printf ("show segment lines? (0: no, 1: yes): ");
read_long (&show_segment_lines);

printf ("output image (pgm):                  ");
read_string (out);
printf ("\n");


/* ---- segment image ---- */

watershed (nx, ny, sigma, show_segment_lines, f, u);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# toboggan watershed algorithm\n");
comment_line (comments, "# sigma: %8.4f\n", sigma);
comment_line (comments, "# show lines: %8ld\n", show_segment_lines);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);
disalloc_matrix (f, nx+2, ny+2);
return(0);
}
