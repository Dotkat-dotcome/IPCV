#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#define NONE -1
#define OCCLUSION -2
#define RECONSTRUCTED -3
#define TRUE 1
#define FALSE 0
#define CENTRALPOINT 1
#define PATCHFILLING 2
#define CENTRAL 1
#define FORWARD 2
#define BACKWARD 3
#define NODERIVATIVE 0

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                         TEXTURE RECONSTRUCTION                           */
/*                                                                          */
/* (Copyright by Markus Schwinn, Pascal Peter and Joachim Weickert, 8/2014) */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 texture reconstruction
*/


/* -------------------- global variables -- structures -------------------- */

struct longTupel
{
    long x;
    long y;
};

struct doubleTupel
{
    double x;
    double y;
};

typedef struct longTupel Pixel;
typedef struct doubleTupel Vector;

long   modus = PATCHFILLING;  /* mode: CENTRALPOINT or  PATCHFILLING */
long   greyValues = 255;      /* color space range */
float  hx = 1.0f;             /* spatial grid size in x direction */
float  hy = 1.0f;             /* spatial grid size in y direction */
double **alphachannel;        /*extra channel for defining the defect */
long   occludedPixel;         /*extra channel for defining the occlusion */
double **priorities;          /* priorities */
double **confidenceValues;    /* confidence values */
double **copyAlphachannel;    /* copy of alphachannel */
long   leftX, rightX;         /* original boundaries in X direction */
long   upperY, lowerY;        /* original boundaries in Y direction */

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

void disalloc_matrix

     (double **matrix,   /* matrix */
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

     (double **u,          /* image, unchanged */ 
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

void dummies
 
     (double **u,        /* image matrix */
      long   m,          /* patch size */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/* creates dummy boundaries with Dirichlet boundary condition */

{
long i, j;  /* loop variables */

for (i=0; i<m; i++) 
 for (j=0; j<(2*m)+ny; j++) 
     {
     u[i][j] = 0.0f;
     u[m + nx + i][j] = 0.0f;
     }

for (j=0; j<m; j++) 
 for (i=0; i<(2*m)+nx; i++) 
     {
     u[i][j] = 0.0f;
     u[i][m + ny + j] = 0.0f;
     }
  
return;
}  

/*--------------------------------------------------------------------------*/

long is_border

     (double  **u,        /* input image */
      Pixel   p)          /* pixel to check */

/* checks if p is the border of u */
      
{
long isB = FALSE;

if (u[p.x][p.y] == OCCLUSION) 
   {
   if (  (u[p.x + 1][p.y]     != OCCLUSION)
      || (u[p.x - 1][p.y]     != OCCLUSION)
      || (u[p.x][p.y + 1]     != OCCLUSION)
      || (u[p.x][p.y - 1]     != OCCLUSION)
      || (u[p.x - 1][p.y + 1] != OCCLUSION)
      || (u[p.x + 1][p.y + 1] != OCCLUSION)
      || (u[p.x - 1][p.y - 1] != OCCLUSION)
      || (u[p.x - 1][p.y + 1] != OCCLUSION))
         isB = TRUE;
   }

return isB;
}

/*--------------------------------------------------------------------------*/

double** create_patch

     (double  **u,        /* processed image */
      long    posX,       /* central x coordinate */
      long    posY,       /* central y coordinate */
      long    m)          /* patch size */
      
 /* creates patch */
 
{
long    i, j, k, l;     /* loop variables */
double  ** patch;       /* patch */
long    mHalf = m / 2;  /* time saving variable */

alloc_matrix (&patch, m, m);

k = -1;
for (i=posX-mHalf; i<=posX+mHalf; i++)
    {
    k++;
    l = -1;
    for (j=posY-mHalf; j<=posY+mHalf; j++) 
        {
        l++;
        if (alphachannel[i][j] != OCCLUSION)
           patch[k][l] = u[i][j];
        else 
           patch[k][l] = alphachannel[i][j];
        }
    }

return patch;
}

/*--------------------------------------------------------------------------*/

double** create_alpha_patch

     (double  **u,        /* processed image */
      long    posX,       /* central x coordinate */
      long    posY,       /* central y coordinate */
      long    m)          /* patch size */

/*
  equivalent to createPatch, but with information if non-occluded
  pixel are reconstructed or original
*/

{
long    i, j, k, l;     /* loop variables */
double  ** patch;       /* patch */
long    mHalf = m / 2;  /* time saving variable */
    
alloc_matrix (&patch, m, m);

k = -1;
for (i=posX-mHalf; i<=posX+mHalf; i++) 
    {
    k++;
    l = -1;
    for (j=posY-mHalf; j<=posY+mHalf; j++) 
        {
        l++;
        if (alphachannel[i][j] == NONE)
           patch[k][l] = u[i][j];
        else
           patch[k][l] = alphachannel[i][j];
        }
    }

return patch;
}

/*--------------------------------------------------------------------------*/

char patch_equivalence

     (double  **targetPatch,  /* target patch */
      double  **sourcePatch,  /* source patch */
      long    m)              /* length of the patch */

/*
  additional testing, if source-patch fits to target-patch
  (change this method if implentation becomes rotationally invariant)
*/

{
long  i, j;   /* loop variables */

for (i=0; i<m-1; i++) 
 for (j=0; j<m-1; j++)
     if (sourcePatch[i][j] == OCCLUSION)
        if (targetPatch[i][j] != OCCLUSION)
           return FALSE;

return TRUE;
}

/*--------------------------------------------------------------------------*/

void clear_patch

     (double  **patch,    /* input matrix */
      long    m)          /* size of matrix in each dimension */
      
/* disallocates the memory for patch */
      
{
disalloc_matrix (patch, m, m);
return;
}

/*--------------------------------------------------------------------------*/

Pixel high_priory_pixel

     () 
     
/* finds high priory pixel */     
     
{
Pixel   p = {-1, -1};      /* high priory pixel */
double  max = 0.0f;        /* max value of priorities */
double  tmpMax = 0.0f;     /* helping variable */
long    k, l;              /* loop variables */

for (k=leftX; k<rightX; k++) 
    {
    for (l=upperY; l<lowerY; l++) 
        {
        tmpMax = priorities[k][l];
        Pixel dummypixel = {k, l};
        if ((max <= tmpMax)
             && (alphachannel[dummypixel.x][dummypixel.y] == OCCLUSION)
             && (is_border (copyAlphachannel, dummypixel))) 
	   {
           max = tmpMax;
           p.x = k;
           p.y = l;
           }
        }
    }

return p;
}

/*--------------------------------------------------------------------------*/

double inner_prodcut

     (Vector v1,          /* first vector */
      Vector v2)          /* second vector */

/* applies inner product */

{
return (v1.x * v2.x) + (v1.y * v2.y);
}

/*--------------------------------------------------------------------------*/

double L2Norm

     (Vector v)           /* input vector */

/* finds L2-norm of vector v */

{
return sqrtf (inner_prodcut (v, v));
}

/*--------------------------------------------------------------------------*/

double confidence_term

     (double  ** confCopy, /* input matrix */
      Pixel   p,           /* input pixel */
      long    m)           /* patch size */
      
/* finds the confidence term */

{
double  **patch;           /* patch of confCopy and p */
double  confidence = 0.0f; /* confidence term */
long    i, j;              /* loop variables */

patch = create_patch (confCopy, p.x, p.y, m);

for (i=0; i<m; i++) 
 for (j=0; j<m; j++)
     confidence += (patch[i][j] >= 0.0f) ? patch[i][j] : 0.0f;

confidence /= (double)(m * m);

confidenceValues[p.x][p.y] = confidence;

if (confidence < 0.0f) 
   {
   printf("confidence is lower zero\n");
   exit(0);
   }

clear_patch (patch, m);

return confidence;
}

/*--------------------------------------------------------------------------*/

Vector image_gradient

     (Pixel   p,          /* input pixel */
      double  ** patch,   /* input patch */
      long    m)          /* patch size */

/*
  computing image gradient with forward derivative
 */
{
long i, j;                                /* loop variables */
Vector gradient = {0, 0};                 /* relative values inside patch */
long countValuesX = 0, countValuesY = 0;  /* relative values inside patch */

for (i=0; i<m-1; i++) 
    {
    for (j=0; j<m-1; j++) 
        {
        /* derivative in x direction */
        if ((patch[i][j] != OCCLUSION) && (patch[i + 1][j] != OCCLUSION)) 
           {
           gradient.x += (patch[i + 1][j] - patch[i][j]) / hx;
           countValuesX++;
           }

           /* Derivative in y Direction */
           if ((patch[i][j] != OCCLUSION) && (patch[i][j + 1] != OCCLUSION)) 
              {
              gradient.y += (patch[i][j + 1] - patch[i][j]) / hy;
              countValuesY++;
              }
        }
    }

if (countValuesX > 0)
   gradient.x /= countValuesX;

if (countValuesY > 0)
   gradient.y /= countValuesY;

return gradient;
}

/*--------------------------------------------------------------------------*/

Vector normal_vector

     (double   **u,       /* input matrix */
      Pixel    p)         /* pixel to create patch */

/*
  Returns the normalized normal vector with respect to border front,
*/

{
Vector n;       /* normal vector */
long   i, j;    /* loop variables */
double **edgePattern = create_patch(u, p.x, p.y, 3);

/* compute the derivatives for channel m using the Sobel operator */
for (i=0; i<3; i++)
 for (j=0; j<3; j++) 
     edgePattern[i][j] = (edgePattern[i][j] != OCCLUSION) ? 0.0f : 1.0f;

n.x = (   - edgePattern[0][0] +     edgePattern[2][0]
      - 2 * edgePattern[0][1] + 2 * edgePattern[2][1]
          - edgePattern[0][2] +     edgePattern[2][2]) / (8.0 * hx);

n.y = (   - edgePattern[0][0] - 2 * edgePattern[1][0]
          - edgePattern[2][0] +     edgePattern[0][2]
      + 2 * edgePattern[1][2] +     edgePattern[2][2]) / (8.0 * hy);

double norm = L2Norm (n);

if (norm != 0.0f) 
   {
   n.x = n.x / norm;
   n.y = n.y / norm;
   }

clear_patch (edgePattern, 3);

return n;
}

/*--------------------------------------------------------------------------*/

double data_term

     (double   **u,       /* input image */
      Pixel    p,         /* input pixel */
      long     m)         /* patch size */
      
/* returns data term */      
      
{
double d = 0.0f;    /* data term */

double **patch = create_patch (u, p.x, p.y, m);
Vector imageGrad = {0.0, 0.0};
Vector normal = {0.0, 0.0};
Vector orthogonalImageGrad = {0, 0};

imageGrad = image_gradient (p, patch, m);
normal = normal_vector (u, p);
orthogonalImageGrad.x =   imageGrad.y;
orthogonalImageGrad.y = - imageGrad.x;

/* compute data term */
d = fabs (orthogonalImageGrad.x * normal.x + orthogonalImageGrad.y * normal.y);
d /= (double)greyValues;
clear_patch(patch, m);

return d;
}

/*--------------------------------------------------------------------------*/

long priority

     (double   **u,       /* input image */
      long     nx,        /* image size in x direction */
      long     ny,        /* image size in y direction */
      Pixel    p,         /* input pixel */
      long     m)         /* patch size */
      
/* computes priority */      
      
{
long   i, j;         /* loop variables */
long   pixels = 0;   /* number of pixels */
long   nonZero = 0;  /* number of pixels with zero prio */
double prio;         /* priority */
double **copyConf;   /* copy of confidence */

if (!copyAlphachannel) 
   alloc_matrix (&copyAlphachannel, nx, ny);

alloc_matrix (&copyConf, nx, ny);

/*copy confidences */
for (i=0; i<nx; i++) 
 for (j=0; j<ny; j++) 
     {
     copyConf[i][j] = confidenceValues[i][j];
     copyAlphachannel[i][j] = alphachannel[i][j];
     }

for (i=leftX; i<rightX; i++) 
 for (j=upperY; j<lowerY; j++)
     {
     p.x = i;
     p.y = j;
     if (is_border (copyAlphachannel, p)) 
        {
        double conf = confidence_term (copyConf, p, m);
        double d = data_term (u, p, m);
        
        prio = conf * d;
        priorities[p.x][p.y] = prio;
        pixels++;
        /* if (prio > 0.0001f)
	      {
              nonZero++;
              }
        */
        }
     }

disalloc_matrix (copyConf, nx, ny);

return (nonZero != 0) ? nonZero : pixels;
}

/*--------------------------------------------------------------------------*/

double ASAD

     (double   **patch,       /* target-patch */
      double   **compPatch,   /* source-patch */
      double   **alphaPatch,  /* alpha information */
      long     m)             /* patch size */

/*
  computes adapted SAD for m * m patch
*/

{
long   i, j;                 /* loop variables */
double sim = 0.0f;           /* similarity */ 
long   pixels = m * m;       /* number of pixels in the patch */
double theta = 0.25 * sqrtf( 2.0f * m * m);  /* theta */
long   mHalf = ceil(m / 2);  /* half of the patch size */
double dist, distX, distY;   /* distance in x, y directions */
double weightFactor;         /* weighting factor */

/* weighting distance in the patch. central will be higher weighted
   theta is equal to one forth of patch diagonal */
for (i=0; i<m; i++) 
 for (j=0; j<m; j++) 
     {
     distX = abs (mHalf - i);
     distY = abs (mHalf - j);
     dist  = sqrtf ((distX * distX) + (distY * distY));
     weightFactor = expf (-(dist * dist) / (2.0f * theta * theta));
     if ((patch[i][j] != OCCLUSION) && (compPatch[i][j] != OCCLUSION)) 
        {
        /* trusting own reconstruction half than original pixels */
        if (alphaPatch[i][j] == RECONSTRUCTED) 
           sim = sim + (weightFactor * 2 * abs(patch[i][j] - compPatch[i][j]));
        else 
           sim = sim + (weightFactor * abs(patch[i][j] - compPatch[i][j]));
        } 
     else 
        pixels--;
     }

/* normalize sim to compare patch with different size of non-occluded pixels */
sim = sim / (double)(pixels);

return sim;
}

/*--------------------------------------------------------------------------*/

long sliceCounter = 0;    /* number of slices */
long frontPixel = 0;      /* number of front pixels */

Pixel find_next_patch

     (double  **u,        /* image */
      long    nx,         /* image width */
      long    ny,         /* image height */
      long    m)          /* patch size */

/*
  returns the position of defect central pixel which patch has 
  the highest information
*/

{
Pixel p;

p.x = -1;
p.y = -1;
if (modus == PATCHFILLING) 
   {
   /* compute priority for next decision */
   priority (u, nx + m - 1, ny + m - 1, p, m);

   /* plot intermediate steps */
   sliceCounter++;
   return high_priory_pixel();
   }

/* modus == CENTRALPIXEL */

/* recompute priority only for a complete boundary and process each 
   boundary pixel with fixed prio for a slice */
if (frontPixel == 0) 
   {
   /* returns the total amount of boundary pixels */
   frontPixel = priority (u, nx + m - 1, ny + m - 1, p, m);

   /* plot intermediate steps */
   sliceCounter++;
   }

frontPixel--;
p = high_priory_pixel();

return p;

} /* find_next_patch */

/*--------------------------------------------------------------------------*/

int replace_central_pixel 
   
     (double  **u,        /* grey value image */
      Pixel   p,          /* pixel to replace */
      long    matchPosX,  /* source position in x */
      long    matchPosY,  /* source position in y */
      long    m)          /* patch size*/
                
/* replaces central pixel */

{
long  counter = 0;   /* number of reconstructed pixels */
  
u[p.x][p.y] = u[matchPosX][matchPosY];
counter++;
alphachannel[p.x][p.y] = RECONSTRUCTED;

return counter;
}

/*--------------------------------------------------------------------------*/

long replace_occluded_patch_parts

     (double  **u,        /* grey value image */
      Pixel   p,          /* central pixel of target-patch */
      long    matchPosX,  /* central position of source-patch in x */
      long    matchPosY,  /* central position of source-patch in y */
      long    m)          /* patch size*/

/* replaces occluded patch parts */

{
typedef struct 
    {
    Pixel globalPos;
    Pixel localPos;
    } Cache;
    
long   i, j;               /* loop variables */
long   xOffset, yOffset;   /* offsets in x and y direction */
long   pixelCounter = 0;   /* number of pixels replaced */
long   counter = 0;        /* counter */
Cache  pixels[m * m];      /* cache of pixels */
double **copyC;            /* copy pixels */
double **sourcePatch = create_patch (u, matchPosX, matchPosY, m);
double **targetPatch = create_patch (u, p.x, p.y, m);

/* collect pixels for replacement */
xOffset = - (m / 2) - 1;
for (i=0; i<m; i++) 
    {
    xOffset++;
    yOffset = - (m / 2) - 1;
    for (j=0; j<m; j++) 
        {
        yOffset++;

        /* skip all pixels which are not occluded */
        if (targetPatch[i][j] != OCCLUSION)
           continue;

        /* skip if no source pixel for target available */
        if (sourcePatch[i][j] == OCCLUSION) 
           continue;

        /* save local position within the target patch */
        pixels[pixelCounter].localPos.x = i;
        pixels[pixelCounter].localPos.y = j;

        /* save global position wrt to the complete image */
        pixels[pixelCounter].globalPos.x = p.x + xOffset;
        pixels[pixelCounter].globalPos.y = p.y + yOffset;
        pixelCounter++;
        }
    }

/* copy pixels to target patch and update confidences and alpha channel
   confidence has to be set only for boundary pixels to ensure that 
   confidence is descending by distance to original data */
alloc_matrix(&copyC, (rightX + (m / 2)), (lowerY + (m / 2)));

long filled = FALSE;

while (!filled) 
      {
      /* copy confidences */
      for (i=leftX; i<rightX; i++) 
       for (j=upperY; j<lowerY; j++) 
           copyC[i][j] = confidenceValues[i][j];

      filled = TRUE;
      for (i=0; i<pixelCounter; i++) 
          {
          if (is_border (alphachannel, pixels[i].globalPos)) 
             {
             filled = FALSE;
             long globX = pixels[i].globalPos.x;
             long globY = pixels[i].globalPos.y;
             long locX = pixels[i].localPos.x;
             long locY = pixels[i].localPos.y;
          
             u[globX][globY] = sourcePatch[locX][locY];
             alphachannel[globX][globY] = RECONSTRUCTED;
             confidenceValues[globX][globY] =
             confidence_term(copyC, pixels[i].globalPos, m);
             counter++;
             }
          }
      }

disalloc_matrix (copyC, (rightX + (m / 2)), (lowerY + (m / 2)));
clear_patch (sourcePatch, m);
clear_patch (targetPatch, m);

return counter;

} /* replace_occluded_patch_parts */

/*--------------------------------------------------------------------------*/

void init_prio

     (long nx,            /* image size in x direction */
      long ny)            /* image size in y direction */
      
/* 
  initializes confidenceValues with 1 if there is original image information
  or 0 if occluded.
*/
      
{
long  k, l;   /* loop variables */

alloc_matrix (&priorities, nx, ny);
alloc_matrix (&confidenceValues, nx, ny);
for (l=0; l<ny; l++)
 for (k=0; k<nx; k++) 
     {
     priorities[k][l] = 0;
     confidenceValues[k][l] = 0;
     }

for (k=0; k<nx; k++) 
 for (l=0; l<ny; l++) 
     {
     priorities[k][l] = 0.0f;
     if ((k < leftX) || (k >= rightX) || (l < upperY) || (l >= lowerY))
        confidenceValues[k][l] = 0.0f;
     else 
        if (alphachannel[k][l] == NONE) 
           confidenceValues[k][l] = 1.0f;
        else 
           confidenceValues[k][l] = 0.0f;
     }
     
return;
}

/*--------------------------------------------------------------------------*/

void reconstruction

     (long    nx,         /* image dimension in x direction */
      long    ny,         /* image dimension in y direction */
      long    m,          /* size of pattern-mask: m * m */
      double  **u)        /* input: occluded image; output: processed */
    
/*
  Inpainting missing texture parts with a m * m mask
*/
    
{          
long    i, j;                 /* loop variables */
double  **patch;              /* patch which has to be restored */
double  **compPatch;          /* patch which has to be compared */
double  **alphaPatch;         /* patch with information from alpha channel */
double  asad;                 /* similarity value for compared patch */
double  asadMin;              /* maximal similarity value */
long    matchPosX, matchPosY; /* position for central pixel to be replaced */
Pixel   p;                    /* pixel */
long    counter = 0;          /* counts occluded pixels*/
long    printCounter = 0;     /* variables to steer the output frequency */
long    printDelta = 100;     /* variables to steer the output frequency */

/* replace all occluded pixel until each pixel is occluded */
while (counter < occludedPixel) 
      {
      if (!((modus == CENTRALPOINT) || (modus == PATCHFILLING)))
         exit( printf ("Illegal Mode"));

      /* get central position of next patch to dissoclude */
      p = find_next_patch (u, nx, ny, m);

      /* create next patch for disoccluding */
      patch = create_patch(u, p.x, p.y, m);
      asadMin = INFINITY;

      /* search nearest match */
      for (i=leftX+(m/2); i<rightX-(m/2); i++) 
       for (j=upperY+(m/2); j<lowerY-(m/2); j++) 
           {
           /* skip patch with central pixels which are part of the occlusion */
           if (alphachannel[i][j] == OCCLUSION)
              continue;
           compPatch = create_patch (u, i, j, m);
           alphaPatch = create_alpha_patch (u, i, j, m);

           /* skip patchs which are occluded
              not covers the relative information */
           if (!patch_equivalence (patch, compPatch, m))
              {
              clear_patch (compPatch, m);
              clear_patch (alphaPatch, m);
              continue;
              }

              /* compute asad */
              asad = ASAD(patch, compPatch, alphaPatch, m);
              if (asad < asadMin) 
                 {
                 asadMin   = asad;
                 matchPosX = i;
                 matchPosY = j;
                 }

              clear_patch (compPatch, m);
              clear_patch (alphaPatch, m);
            }
        
      /* replace central pixel or occluded patch parts by best fitting */
      counter += (modus == CENTRALPOINT) ?
                  replace_central_pixel (u, p, matchPosX, matchPosY, m) :
                  replace_occluded_patch_parts (u, p, matchPosX, matchPosY, m);

      printCounter++;
      if (printCounter >= printDelta) 
         {
         printf("pixel restored - %i/%ld \n", counter, occludedPixel);
         printCounter = 0;
         }
      clear_patch (patch, m);
      }
printf ("\n");      
      
return;

} /* reconstruction */

/*--------------------------------------------------------------------------*/

int main ()

{
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
double **u;                  /* image */
double **tmp;                /* help matrix */
double **mask;               /* mask */
long   nx, ny;               /* image size in x, y direction */ 
long   m;                    /* size for pattern-mask */
long   imax, jmax;           /* help indices */
long   i, j, k, l;           /* loop variables */
char   comments[1600];       /* string for comments */

printf ("\n");
printf ("TEXTURE RECONSTRUCTION\n\n");
printf ("**************************************************\n\n");
printf ("    Copyright 2014 by Joachim Weickert,           \n");
printf ("    Markus Schwinn and Pascal Peter               \n");
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
read_pgm_and_allocate_memory (in, &nx, &ny, &tmp);

printf ("input mask (pgm):                 ");
read_string (in);
read_pgm_and_allocate_memory (in, &nx, &ny, &mask);


/* ---- read parameters ---- */

printf ("patch size m (integer):           ");
read_long (&m);

printf ("output image (pgm):               ");
read_string (out);
printf ("\n");


/* ---- set boundaries of original image ---- */

leftX  = m;
rightX = m + nx;
upperY = m;
lowerY = m + ny;


/* ---- memory allocations ---- */

alloc_matrix (&u, nx+2*m, ny+2*m);
alloc_matrix (&alphachannel, nx+2*m, ny+2*m);


/* ---- initialize image u ---- */

imax=nx+m-1;
jmax=ny+m-1;
for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     u[i][j] = tmp[k][l];
  
/* use Dirichlet boundary contitions */
dummies (u, m, nx, ny);


/* ---- process mask ---- */

occludedPixel = 0;
for (j=0; j<ny+2*m; j++)
 for (i=0; i<nx+2*m; i++)
     alphachannel[i][j] = NONE;
for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     if (mask[k][l] < 0.5) 
        {
        alphachannel[i][j] = OCCLUSION;
        occludedPixel++;
        }
     else 
        alphachannel[i][j] = NONE;


/* ---- process image ---- */

init_prio (nx+2*m, ny+2*m);
reconstruction (nx, ny, 2*m+1, u);


/* ---- write solution back ---- */

for (j=m, l=1; j<=jmax; j++, l++)
 for (i=m, k=1; i<=imax; i++, k++)
     tmp[k][l] = u[i][j];


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# texture reconstruction\n");
comment_line (comments, "# m: %8ld\n", m);

/* write image */
write_pgm (tmp, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

disalloc_matrix (tmp,  nx+2,   ny+2);
disalloc_matrix (u,    nx+2*m, ny+2*m);
disalloc_matrix (mask, nx+2,   ny+2);
disalloc_matrix (alphachannel, nx+2*m, ny+2*m);

return(0);
}
