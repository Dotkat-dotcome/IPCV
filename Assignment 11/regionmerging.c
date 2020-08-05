#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                             REGION MERGING                               */
/*                                                                          */
/*       (Copyright by Markus Mainberger and Joachim Weickert, 8/2014)      */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 Region Merging Algorithm
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

void region_merging 

     (long     nx,                  /* image dimension in x direction */
      long     ny,                  /* image dimension in y direction */
      float    T,                   /* merging threshold */
      long     show_segment_lines,  /* draw white segment boundaries? */
      float    **f,                 /* input: original; output: processed */
      float    **u)                 /* output image */

/* perform region merging */

{
long  x, y;             /* coordinates giving the current position */
long  i, j;             /* loop variables */
float **visited;        /* log visited points */
Stack *segment;         /* contains pixels that belomg to current segment */
float seg_mean;         /* mean of currently considered segment */
long  sum_pixel;        /* number of pixels of currently considered segment */
Stack *to_be_visited;   /* contains all pixels that have still to be visited
                           for current segment */

/* allocate memory */
alloc_matrix(&visited, nx+2, ny+2);

/* initialise visited-array (visited pixels will be marked with 1.0f) */
for (i=0; i<nx+2; ++i)
 for (j=0; j<ny+2; ++j)
     visited[i][j] = 0.0f;

/* set boundary as visited */
set_bounds(visited, nx, ny, 1.0f);

/* create empty stacks */
segment = new_Stack ();
to_be_visited = new_Stack ();

/* start region merging */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* if pixel was not visited yet start looking for a new segment */
     if (visited[i][j] != 1.0f)
        {
        seg_mean = 0.0f;
        sum_pixel = 0;

        /* add pixel to the set of pixels which have still to be visited
           for the currently considered segment */
        push_Stack(to_be_visited,new_intPair(i,j));

        while (!empty_Stack(to_be_visited))
              { 
              /* get a pixel which has still to be visited */
              x = ((intPair*) front_Stack (to_be_visited))->x;
              y = ((intPair*) front_Stack (to_be_visited))->y;

              /* remove it from "to_be_visited" stack */
              pop_Stack (to_be_visited);
              if (visited[x][y] != 1.0f)
                 {
                 /* mark pixel as visited 
                    and add it to the current segment */
                 visited[x][y] = 1.0f;
                 push_Stack (segment, new_intPair (x,y));
                 sum_pixel++;
                 seg_mean += f[x][y];

                 /* if neighbour has not been visited yet and belongs to the 
                    same segment add to "to_be_visited" stack */
                 /*
                   SUPPLEMENT CODE HERE
                 */
                 if ((visited[x+1][y] != 1.0f) && (abs(f[x+1][y]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x+1,y));
                 }
                 if ((visited[x-1][y] != 1.0f) && (abs(f[x-1][y]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x-1,y));
                 }
                 if ((visited[x][y+1] != 1.0f) && (abs(f[x][y+1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x,y+1));
                 }
                 if ((visited[x][y-1] != 1.0f) && (abs(f[x][y-1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x,y-1));
                 }
                 if ((visited[x+1][y+1] != 1.0f) && (abs(f[x+1][y+1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x+1,y+1));
                 }
                 if ((visited[x+1][y-1] != 1.0f) && (abs(f[x+1][y-1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x+1,y-1));
                 }
                 if ((visited[x-1][y+1] != 1.0f) && (abs(f[x-1][y+1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x-1,y+1));
                 }
                 if ((visited[x-1][y-1] != 1.0f) && (abs(f[x-1][y-1]-f[x][y])<T))
                 {
                     push_Stack(to_be_visited,new_intPair(x-1,y-1));
                 }
                 }
              }

        /* compute mean value of segment */
        seg_mean /= (float) sum_pixel;

        /* set value of segment pixels to mean value of segment */
        while (!empty_Stack(segment))
              {
               x = ((intPair*) front_Stack (segment))->x;
               y = ((intPair*) front_Stack (segment))->y;
               pop_Stack (segment);
               u[x][y] = seg_mean;
              }
        }
     }


/* ---- prepare output image ---- */

for (i=1; i<nx+1; i++)
 for (j=1; j<ny+1; j++)
     {
     if (show_segment_lines && 
           (((i < nx) && ((int)u[i][j] != (int)u[i+1][j]))
         || ((j < ny) && ((int)u[i][j] != (int)u[i][j+1]))))
        u[i][j] =  255.0f;
     }
    
/* free memory */
free_Stack (segment);
free_Stack (to_be_visited);
disalloc_matrix (visited, nx+2, ny+2);

return;

} /* region_merging */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **f;                  /* image */
float  **u;                  /* image */
long   nx, ny;               /* image size in x, y direction */ 
float  T;                    /* merging threshold */
long   show_segment_lines;   /* draw white segment boundaries? */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  std;                  /* standard deviation */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("REGION MERGING\n\n");
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

printf ("threshold (float):                   ");
read_float (&T);

printf ("show segment lines? (0: no, 1: yes): ");
read_long (&show_segment_lines);

printf ("output image (pgm):                  ");
read_string (out);
printf ("\n");


/* ---- segment image ---- */

region_merging (nx, ny, T, show_segment_lines, f, u);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# region merging algorithm\n");
comment_line (comments, "# T: %8.4f\n", T);
comment_line (comments, "# show lines: %8ld\n", show_segment_lines);

/* write image */
write_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (u, nx+2, ny+2);
disalloc_matrix (f, nx+2, ny+2);
return(0);
}
