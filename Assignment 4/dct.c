#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                        DISCRETE COSINE TRANSFORM                         */
/*                                                                          */
/*    (Copyright by Andres Bruhn, 11/2007 and Joachim Weickert, 8/2014)     */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/* 
 Discrete Cosine Transform
*/

/*--------------------------------------------------------------------------*/

void alloc_vector

    (float **vector, /* vector */
     long n1)        /* size */

/* allocates memory for a vector of size n1 */

{
    *vector = (float *)malloc(n1 * sizeof(float));
    if (*vector == NULL)
    {
        printf("alloc_vector: not enough memory available\n");
        exit(1);
    }
    return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

    (float ***matrix, /* matrix */
     long n1,         /* size in direction 1 */
     long n2)         /* size in direction 2 */

/* allocates memory for matrix of size n1 * n2 */

{
    long i;

    *matrix = (float **)malloc(n1 * sizeof(float *));
    if (*matrix == NULL)
    {
        printf("alloc_matrix: not enough memory available\n");
        exit(1);
    }
    for (i = 0; i < n1; i++)
    {
        (*matrix)[i] = (float *)malloc(n2 * sizeof(float));
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

    (float *vector, /* vector */
     long n1)       /* size */

/* disallocates memory for a vector of size n1 */

{
    free(vector);
    return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

    (float **matrix, /* matrix */
     long n1,        /* size in direction 1 */
     long n2)        /* size in direction 2 */

/* disallocates memory for matrix of size n1 * n2 */

{
    long i;

    for (i = 0; i < n1; i++)
        free(matrix[i]);

    free(matrix);

    return;
}

/*--------------------------------------------------------------------------*/

void read_string

    (char *v) /* string to be read */

/*
 reads a long value v
*/

{
    fgets(v, 80, stdin);
    if (v[strlen(v) - 1] == '\n')
        v[strlen(v) - 1] = 0;
    return;
}

/*--------------------------------------------------------------------------*/

void read_long

    (long *v) /* value to be read */

/*
 reads a long value v
*/

{
    char row[80]; /* string for reading data */

    fgets(row, 80, stdin);
    if (row[strlen(row) - 1] == '\n')
        row[strlen(row) - 1] = 0;
    sscanf(row, "%ld", &*v);
    return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

    (const char *file_name, /* name of pgm file */
     long *nx,              /* image size in x direction, output */
     long *ny,              /* image size in y direction, output */
     float ***u)            /* image, output */

/* 
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u; 
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
    FILE *inimage; /* input file */
    char row[80];  /* for reading data */
    long i, j;     /* loop variables */

    /* open file */
    inimage = fopen(file_name, "rb");
    if (NULL == inimage)
    {
        printf("could not open file '%s' for reading, aborting.\n", file_name);
        exit(1);
    }

    /* read header */
    fgets(row, 80, inimage); /* skip format definition */
    fgets(row, 80, inimage);
    while (row[0] == '#') /* skip comments */
        fgets(row, 80, inimage);
    sscanf(row, "%ld %ld", nx, ny); /* read image size */
    fgets(row, 80, inimage);        /* read maximum grey value */

    /* allocate memory */
    alloc_matrix(u, (*nx) + 2, (*ny) + 2);

    /* read image data row by row */
    for (j = 1; j <= (*ny); j++)
        for (i = 1; i <= (*nx); i++)
            (*u)[i][j] = (float)getc(inimage);

    /* close file */
    fclose(inimage);

    return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void comment_line

    (char *comment,    /* comment string (output) */
     char *lineformat, /* format string for comment line */
     ...)              /* optional arguments */

/* 
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",float_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
    char line[80];
    va_list arguments;

    /* get list of optional function arguments */
    va_start(arguments, lineformat);

    /* convert format string and arguments to plain text line string */
    vsprintf(line, lineformat, arguments);

    /* add line to total commentary string */
    strncat(comment, line, 80);

    /* add line break if input string does not end with one */
    if (line[strlen(line) - 1] != '\n')
        sprintf(comment, "%s\n", comment);

    /* close argument list */
    va_end(arguments);

    return;

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_pgm

    (float **u,       /* image, unchanged */
     long nx,         /* image size in x direction */
     long ny,         /* image size in y direction */
     char *file_name, /* name of pgm file */
     char *comments)  /* comment string (set 0 for no comments) */

/* 
  writes a greyscale image into a pgm P5 file;
*/

{
    FILE *outimage;     /* output file */
    long i, j;          /* loop variables */
    float aux;          /* auxiliary variable */
    unsigned char byte; /* for data conversion */

    /* open file */
    outimage = fopen(file_name, "wb");
    if (NULL == outimage)
    {
        printf("Could not open file '%s' for writing, aborting\n", file_name);
        exit(1);
    }

    /* write header */
    fprintf(outimage, "P5\n"); /* format */
    if (comments != 0)
        fprintf(outimage, comments);        /* comments */
    fprintf(outimage, "%ld %ld\n", nx, ny); /* image size */
    fprintf(outimage, "255\n");             /* maximal value */

    /* write image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
        {
            aux = u[i][j] + 0.499999; /* for correct rounding */
            if (aux < 0.0)
                byte = (unsigned char)(0.0);
            else if (aux > 255.0)
                byte = (unsigned char)(255.0);
            else
                byte = (unsigned char)(aux);
            fwrite(&byte, sizeof(unsigned char), 1, outimage);
        }

    /* close file */
    fclose(outimage);

    return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

void analyse

    (float **u,   /* image, unchanged */
     long nx,     /* pixel number in x direction */
     long ny,     /* pixel number in x direction */
     float *min,  /* minimum, output */
     float *max,  /* maximum, output */
     float *mean, /* mean, output */
     float *std)  /* standard deviation, output */

/*
 computes minimum, maximum, mean, and standard deviation of an image u
*/

{
    long i, j;    /* loop variables */
    double help1; /* auxiliary variable */
    float help2;  /* auxiliary variable */

    *min = u[1][1];
    *max = u[1][1];
    help1 = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            if (u[i][j] < *min)
                *min = u[i][j];
            if (u[i][j] > *max)
                *max = u[i][j];
            help1 = help1 + (double)u[i][j];
        }
    *mean = (float)help1 / (nx * ny);

    *std = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
        {
            help2 = u[i][j] - *mean;
            *std = *std + help2 * help2;
        }
    *std = sqrt(*std / (nx * ny));

    return;

} /* analyse */

/*--------------------------------------------------------------------------*/

void jpeg_multiply_block

    (float **c_block) /* coefficients of the DCT */

/*
  weights 8x8 coefficient block with JPEG weighting matrix
*/

{
    long i, j; /* loop variables */

    c_block[0][0] *= 10;
    c_block[0][1] *= 15;
    c_block[0][2] *= 25;
    c_block[0][3] *= 37;
    c_block[0][4] *= 51;
    c_block[0][5] *= 66;
    c_block[0][6] *= 82;
    c_block[0][7] *= 100;

    c_block[1][0] *= 15;
    c_block[1][1] *= 19;
    c_block[1][2] *= 28;
    c_block[1][3] *= 39;
    c_block[1][4] *= 52;
    c_block[1][5] *= 67;
    c_block[1][6] *= 83;
    c_block[1][7] *= 101;

    c_block[2][0] *= 25;
    c_block[2][1] *= 28;
    c_block[2][2] *= 35;
    c_block[2][3] *= 45;
    c_block[2][4] *= 58;
    c_block[2][5] *= 72;
    c_block[2][6] *= 88;
    c_block[2][7] *= 105;

    c_block[3][0] *= 37;
    c_block[3][1] *= 39;
    c_block[3][2] *= 45;
    c_block[3][3] *= 54;
    c_block[3][4] *= 66;
    c_block[3][5] *= 79;
    c_block[3][6] *= 94;
    c_block[3][7] *= 111;

    c_block[4][0] *= 51;
    c_block[4][1] *= 52;
    c_block[4][2] *= 58;
    c_block[4][3] *= 66;
    c_block[4][4] *= 76;
    c_block[4][5] *= 89;
    c_block[4][6] *= 103;
    c_block[4][7] *= 119;

    c_block[5][0] *= 66;
    c_block[5][1] *= 67;
    c_block[5][2] *= 72;
    c_block[5][3] *= 79;
    c_block[5][4] *= 89;
    c_block[5][5] *= 101;
    c_block[5][6] *= 114;
    c_block[5][7] *= 130;

    c_block[6][0] *= 82;
    c_block[6][1] *= 83;
    c_block[6][2] *= 88;
    c_block[6][3] *= 94;
    c_block[6][4] *= 103;
    c_block[6][5] *= 114;
    c_block[6][6] *= 127;
    c_block[6][7] *= 142;

    c_block[7][0] *= 100;
    c_block[7][1] *= 101;
    c_block[7][2] *= 105;
    c_block[7][3] *= 111;
    c_block[7][4] *= 119;
    c_block[7][5] *= 130;
    c_block[7][6] *= 142;
    c_block[7][7] *= 156;

    return;

} /* jpeg_multiply_block */

/*--------------------------------------------------------------------------*/

void equal_multiply_block

    (float **c_block, /* coefficients of the DCT */
     long factor)     /* factor to multiply */

/*
  multiplies an 8x8 coefficient block with a given factor 
*/

{
    long i, j; /* loop variables */

    for (i = 0; i <= 7; i++)
        for (j = 0; j <= 7; j++)
            c_block[i][j] *= factor;

    return;
}

/*--------------------------------------------------------------------------*/

void jpeg_divide_block

    (float **c_block) /* coefficients of the DCT */

/*
  weights 8x8 coefficient block with inverse JPEG weighting matrix
*/

{
    long i, j; /* loop variables */

    c_block[0][0] /= 10;
    c_block[0][1] /= 15;
    c_block[0][2] /= 25;
    c_block[0][3] /= 37;
    c_block[0][4] /= 51;
    c_block[0][5] /= 66;
    c_block[0][6] /= 82;
    c_block[0][7] /= 100;

    c_block[1][0] /= 15;
    c_block[1][1] /= 19;
    c_block[1][2] /= 28;
    c_block[1][3] /= 39;
    c_block[1][4] /= 52;
    c_block[1][5] /= 67;
    c_block[1][6] /= 83;
    c_block[1][7] /= 101;

    c_block[2][0] /= 25;
    c_block[2][1] /= 28;
    c_block[2][2] /= 35;
    c_block[2][3] /= 45;
    c_block[2][4] /= 58;
    c_block[2][5] /= 72;
    c_block[2][6] /= 88;
    c_block[2][7] /= 105;

    c_block[3][0] /= 37;
    c_block[3][1] /= 39;
    c_block[3][2] /= 45;
    c_block[3][3] /= 54;
    c_block[3][4] /= 66;
    c_block[3][5] /= 79;
    c_block[3][6] /= 94;
    c_block[3][7] /= 111;

    c_block[4][0] /= 51;
    c_block[4][1] /= 52;
    c_block[4][2] /= 58;
    c_block[4][3] /= 66;
    c_block[4][4] /= 76;
    c_block[4][5] /= 89;
    c_block[4][6] /= 103;
    c_block[4][7] /= 119;

    c_block[5][0] /= 66;
    c_block[5][1] /= 67;
    c_block[5][2] /= 72;
    c_block[5][3] /= 79;
    c_block[5][4] /= 89;
    c_block[5][5] /= 101;
    c_block[5][6] /= 114;
    c_block[5][7] /= 130;

    c_block[6][0] /= 82;
    c_block[6][1] /= 83;
    c_block[6][2] /= 88;
    c_block[6][3] /= 94;
    c_block[6][4] /= 103;
    c_block[6][5] /= 114;
    c_block[6][6] /= 127;
    c_block[6][7] /= 142;

    c_block[7][0] /= 100;
    c_block[7][1] /= 101;
    c_block[7][2] /= 105;
    c_block[7][3] /= 111;
    c_block[7][4] /= 119;
    c_block[7][5] /= 130;
    c_block[7][6] /= 142;
    c_block[7][7] /= 156;

    return;

} /* jpeg_divide_block */

/*--------------------------------------------------------------------------*/

void round_block_coeff

    (float **c_block) /* coefficients of the DCT */

/*
  round entries of a 8x8 coefficient block
*/

{
    long i, j; /* loop variables */

    for (i = 0; i <= 7; i++)
        for (j = 0; j <= 7; j++)
            c_block[i][j] = rintf(c_block[i][j]);

    return;
}

/*--------------------------------------------------------------------------*/

void equal_divide_block

    (float **c_block, /* coefficients of the DCT */
     long factor)     /* factor to divide */

/*
  divides an 8x8 coefficient block by a given factor 
*/

{
    long i, j; /* loop variables */

    for (i = 0; i <= 7; i++)
        for (j = 0; j <= 7; j++)
            c_block[i][j] /= factor;

    return;
}

/*--------------------------------------------------------------------------*/

void DCT_2d

    (float **u, /* image, unchanged */
     float **c, /* coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */

/*
  calulates DCT of input image
*/

{
    long i, j, m, p; /* loop variables */
    float nx_1;      /* time saver */
    float ny_1;      /* time saver */
    float pi;        /* variable pi */
    float **tmp;     /* temporary image */
    float *cx, *cy;  /* arrays for coefficients */

    /* ---- derive pi ---- */

    pi = 2.0 * asinf(1.0);

    /* ---- set time savers ---- */

    nx_1 = pi / (2.0 * nx);
    ny_1 = pi / (2.0 * ny);

    /* ---- allocate memory ---- */

    alloc_matrix(&tmp, nx, ny);
    alloc_vector(&cx, nx);
    alloc_vector(&cy, ny);

    /* ---- set up coefficients ---- */

    cx[0] = sqrt(1.0 / nx);
    for (p = 1; p < nx; p++)
        cx[p] = sqrt(2.0 / nx);

    cy[0] = sqrt(1.0 / ny);
    for (p = 1; p < ny; p++)
        cy[p] = sqrt(2.0 / ny);

    /* ---- DCT in y-direction ---- */

    for (i = 0; i < nx; i++)
        for (p = 0; p < ny; p++)
        {
            float sum = 0;
            for (m = 0; m < ny; m++)
            {
                sum += u[i][m] * cy[p] * cos(ny_1 * (2 * m + 1) * p);
            }
            tmp[i][p] = sum;
        }

    /* ---- DCT in x-direction ---- */

    for (p = 0; p < nx; p++)
        for (j = 0; j < ny; j++)
        {
            float sum = 0;
            for (m = 0; m < ny; m++)
            {
                sum += tmp[m][j] * cx[p] * cos(nx_1 * (2 * m + 1) * p);
            }
            c[p][j] = sum;
        }

    /* ---- free memory ---- */

    disalloc_matrix(tmp, nx, ny);
    disalloc_vector(cx, nx);
    disalloc_vector(cy, ny);

    return;

} /* DCT_2d */

/*--------------------------------------------------------------------------*/

void IDCT_2d

    (float **u, /* image, unchanged */
     float **c, /* coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */

/*
  calulates inverse DCT of input image
*/

{
    long i, j, m, p; /* loop variables */
    float nx_1;      /* time saver */
    float ny_1;      /* time saver */
    float pi;        /* variable pi */
    float **tmp;     /* temporary image */
    float *cx, *cy;  /* arrays for coefficients */

    /* ---- derive pi ---- */

    pi = 2.0 * asinf(1.0);

    /* ---- set time savers ---- */

    nx_1 = pi / (2.0 * nx);
    ny_1 = pi / (2.0 * ny);

    /* ---- allocate memory ---- */

    alloc_matrix(&tmp, nx, ny);
    alloc_vector(&cx, nx);
    alloc_vector(&cy, ny);

    /* ---- set up coefficients ---- */

    cx[0] = sqrt(1.0 / nx);
    for (m = 1; m < nx; m++)
        cx[m] = sqrt(2.0 / nx);

    cy[0] = sqrt(1.0 / ny);
    for (m = 1; m < ny; m++)
        cy[m] = sqrt(2.0 / ny);

    /* ---- DCT in y-direction ---- */

    for (i = 0; i < nx; i++)
        for (m = 0; m < ny; m++)
        {
            float sum = 0;
            for (p = 0; p < ny; p++)
            {
                sum += c[i][p] * cy[p] * cos(ny_1 * (2 * m + 1) * p);
            }
            tmp[i][m] = sum;
        }

    /* ---- DCT in x-direction ---- */

    for (m = 0; m < nx; m++)
        for (j = 0; j < ny; j++)
        {
            float sum = 0;
            for (p = 0; p < ny; p++)
            {
                sum += tmp[p][j] * cx[p] * cos(nx_1 * (2 * m + 1) * p);
            }
            u[m][j] = sum;
        }

    /* ---- free memory ---- */

    disalloc_matrix(tmp, nx, ny);
    disalloc_vector(cx, nx);
    disalloc_vector(cy, ny);

    return;

} /* IDCT_2d */

/*--------------------------------------------------------------------------*/

void remove_freq_2d

    (float **c, /* in and out: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */

/*
 remove frequencies
*/

{
    long i, j; /* loop variables */

    /* sets frequencies to zero */
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            if ((i >= nx / sqrt(10)) || (j >= ny / sqrt(10)))
                c[i][j] = 0;

    return;
}

/*--------------------------------------------------------------------------*/

void blockwise_DCT_2d

    (float **u, /* in: image */
     float **c, /* out: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */

/*
 calculates DCT in image blocks of size 8x8
*/

{
    long i, j, k, l; /* loop variables */
    float **u_block; /* 8x8 block */
    float **c_block; /* 8x8 block */

    /* allocate memory for 8x8 blocks */
    alloc_matrix(&u_block, 8, 8);
    alloc_matrix(&c_block, 8, 8);

    for (i = 1; i <= nx; i += 8)
        for (j = 1; j <= ny; j += 8)
        {
            /* copy 8x8 block */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    u_block[k][l] = u[i + k][j + l];

            /* DCT of 8x8 block */
            DCT_2d(u_block, c_block, 8, 8);

            /* copy back coefficients */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c[i + k][j + l] = c_block[k][l];
        }

    /* free memory for 8x8 block */
    disalloc_matrix(u_block, 8, 8);
    disalloc_matrix(c_block, 8, 8);

    return;

} /* blockwise_DCT_2d */

/*--------------------------------------------------------------------------*/

void blockwise_IDCT_2d

    (float **u, /* out: image */
     float **c, /* in: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */
/*
 calculates inverse DCT in image blocks of size 8x8
*/

{
    long i, j, k, l; /* loop variables */
    float **u_block; /* 8x8 block */
    float **c_block; /* 8x8 block */

    /* allocate memory for 8x8 blocks */
    alloc_matrix(&u_block, 8, 8);
    alloc_matrix(&c_block, 8, 8);

    for (i = 1; i <= nx; i += 8)
        for (j = 1; j <= ny; j += 8)
        {
            /* copy 8x8 block */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c_block[k][l] = c[i + k][j + l];

            /* inverse DCT of 8x8 block */
            IDCT_2d(u_block, c_block, 8, 8);

            /* copy back coefficients */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    u[i + k][j + l] = u_block[k][l];
        }

    /* free memory for 8x8 block */
    disalloc_matrix(u_block, 8, 8);
    disalloc_matrix(c_block, 8, 8);

    return;

} /* blockwise_IDCT_2d */

/*--------------------------------------------------------------------------*/

void blockwise_remove_freq_2d

    (float **c, /* in and out: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */
/*
 removes frequencies within 8x8 block
*/

{
    long i, j, k, l; /* loop variables */
    float **c_block; /* 8x8 block */

    /* allocate memory for 8x8 block */
    alloc_matrix(&c_block, 8, 8);

    /* scale coefficients blockwise */
    for (i = 1; i <= nx; i += 8)
        for (j = 1; j <= ny; j += 8)
        {
            /* copy 8x8 block */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c_block[k][l] = c[i + k][j + l];

            /* set frequencies to zero */
            remove_freq_2d(c_block, 8, 8);

            /* count zeros */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    if (c_block[k][l] == 0)

                        /* copy back coefficients */
                        for (k = 0; k <= 7; k++)
                            for (l = 0; l <= 7; l++)
                                c[i + k][j + l] = c_block[k][l];
        }

    /* free memory for 8x8 block */
    disalloc_matrix(c_block, 8, 8);

    return;

} /* blockwise_remove_freq_2d */

/*--------------------------------------------------------------------------*/

void blockwise_quantisation_jpeg_2d

    (float **c, /* in and out: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */
/*
 quantises coefficients of 8x8 block
*/

{
    long i, j, k, l; /* loop variables */
    float **c_block; /* 8x8 block */

    /* allocate memory for 8x8 block */
    alloc_matrix(&c_block, 8, 8);

    /* scale coefficients blockwise */
    for (i = 1; i <= nx; i += 8)
        for (j = 1; j <= ny; j += 8)
        {
            /* copy 8x8 block */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c_block[k][l] = c[i + k][j + l];

            /* scale coefficients of 8x8 block */
            jpeg_divide_block(c_block);

            /* round coefficients */
            round_block_coeff(c_block);

            /* rescale coefficients of 8x8 block */
            jpeg_multiply_block(c_block);

            /* copy back coefficients */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c[i + k][j + l] = c_block[k][l];
        }

    /* free memory for 8x8 block */
    disalloc_matrix(c_block, 8, 8);

    return;

} /* blockwise_quantisation_jpeg_2d */

/*--------------------------------------------------------------------------*/

void blockwise_quantisation_equal_2d

    (float **c, /* in and out: coefficients of the DCT */
     long nx,   /* pixel number in x-direction */
     long ny)   /* pixel number in y-direction */
/*
 quantises coefficients of 8x8 block
*/

{
    long i, j, k, l; /* loop variables */
    long count;      /* counter */
    float **c_block; /* 8x8 block */

    /* allocate memory for 8x8 block */
    alloc_matrix(&c_block, 8, 8);

    /* scale coefficients blockwise */
    for (i = 1; i <= nx; i += 8)
        for (j = 1; j <= ny; j += 8)
        {
            /* copy 8x8 block */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c_block[k][l] = c[i + k][j + l];

            /* scale coefficients of 8x8 block */
            equal_divide_block(c_block, 40);

            /* round coefficients */
            round_block_coeff(c_block);

            /* rescale coefficients of 8x8 block */
            equal_multiply_block(c_block, 40);

            /* copy back coefficients */
            for (k = 0; k <= 7; k++)
                for (l = 0; l <= 7; l++)
                    c[i + k][j + l] = c_block[k][l];
        }

    /* free memory for 8x8 block */
    disalloc_matrix(c_block, 8, 8);

    return;

} /* blockwise_quantisation_equal_2d */

/*--------------------------------------------------------------------------*/

int main()

{
    char in[80];         /* for reading data */
    char out1[80];       /* for reading data */
    char out2[80];       /* for reading data */
    float **u;           /* image */
    float **c;           /* DCT coefficients */
    long nx, ny;         /* image size in x, y direction */
    long i, j;           /* loop variables */
    long flag;           /* processing flag */
    long counter;        /* count variable */
    float sigma;         /* standard deviation */
    float max, min;      /* largest, smallest grey value */
    float mean;          /* average grey value */
    float std;           /* standard deviation */
    char comments[1600]; /* string for comments */

    printf("\n");
    printf("DISCRETE COSINE TRANSFORM\n\n");
    printf("**************************************************\n\n");
    printf("    Copyright 2014 by Joachim Weickert            \n");
    printf("    and 2007 by Andres Bruhn                      \n");
    printf("    Dept. of Mathematics and Computer Science     \n");
    printf("    Saarland University, Saarbruecken, Germany    \n\n");
    printf("    All rights reserved. Unauthorized usage,      \n");
    printf("    copying, hiring, and selling prohibited.      \n\n");
    printf("    Send bug reports to                           \n");
    printf("    weickert@mia.uni-saarland.de                  \n\n");
    printf("**************************************************\n\n");

    /* ---- read input image (pgm format P5) ---- */

    printf("input image (pgm):                ");
    read_string(in);
    read_pgm_and_allocate_memory(in, &nx, &ny, &u);

    /* check if image can be devided in blocks of size 8x8 */
    if ((nx % 8 != 0) || (ny % 8 != 0))
    {
        printf("\n\n Image size does not allow decomposition! \n\n");
        return (0);
    }

    /* ---- read parameters ---- */

    printf("spectrum image (pgm):             ");
    read_string(out1);

    printf("output image (pgm):               ");
    read_string(out2);

    printf("\n\nyou have the following options:      \n");
    printf("\n    (1) DCT/IDCT of the whole image   ");
    printf("\n    (2) DCT/IDCT of 8x8 blocks        ");
    printf("\n    (3) DCT/IDCT of the whole image   ");
    printf("\n        with removal of frequencies   ");
    printf("\n    (4) DCT/IDCT of 8x8 blocks        ");
    printf("\n        with removal of frequencies   ");
    printf("\n    (5) DCT/IDCT of 8x8 blocks        ");
    printf("\n        with equal quantisation       ");
    printf("\n    (6) DCT/IDCT of 8x8 blocks        ");
    printf("\n        with JPEG quantisation        \n");
    printf("\nchose the processing mode (1-6):      ");

    read_long(&flag);
    printf("\n\n");

    /* ---- allocate storage for c ---- */

    alloc_matrix(&c, nx + 2, ny + 2);

    /* ---- analyse input image ---- */

    analyse(u, nx, ny, &min, &max, &mean, &std);
    printf("input image:\n");
    printf("minimum:       %8.2f \n", min);
    printf("maximum:       %8.2f \n", max);
    printf("mean:          %8.2f \n", mean);
    printf("standard dev.: %8.2f \n\n", std);

    /* ---- process image ---- */

    switch (flag)
    {
    case 1:
        /* perform DCT and IDCT for the whole image */
        DCT_2d(u, c, nx, ny);
        IDCT_2d(u, c, nx, ny);
        break;
    case 2:
        /* perform DCT and IDCT in 8x8 blocks */
        blockwise_DCT_2d(u, c, nx, ny);
        blockwise_IDCT_2d(u, c, nx, ny);
        break;
    case 3:
        /* perform DCT and IDCT for the whole image */
        /* remove frequencies */
        DCT_2d(u, c, nx, ny);
        remove_freq_2d(c, nx, ny);
        IDCT_2d(u, c, nx, ny);
        break;
    case 4:
        /* perform DCT and IDCT in 8x8 blocks */
        /* remove frequencies */
        blockwise_DCT_2d(u, c, nx, ny);
        blockwise_remove_freq_2d(c, nx, ny);
        blockwise_IDCT_2d(u, c, nx, ny);
        break;
    case 5:
        /* perform DCT and IDCT in 8x8 blocks */
        /* and use equal quantisation */
        blockwise_DCT_2d(u, c, nx, ny);
        blockwise_quantisation_equal_2d(c, nx, ny);
        blockwise_IDCT_2d(u, c, nx, ny);
        break;
    case 6:
        /* perform DCT and IDCT in 8x8 blocks */
        /* and use JPEG quantisation */
        blockwise_DCT_2d(u, c, nx, ny);
        blockwise_quantisation_jpeg_2d(c, nx, ny);
        blockwise_IDCT_2d(u, c, nx, ny);
        break;
    default:
        printf("option (%ld) not available! \n\n\n", flag);
        return (0);
    }

    /* ---- count zero entries in the DCT coefficients ---- */

    counter = 0;
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            if (c[i][j] == 0)
                counter++;

    printf("percentage of zero coefficients %f\n\n", (float)counter / (nx * ny));

    /* ---- compute logarithmised spectrum of c ---- */

    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            c[i][j] = logf(1.0f + fabsf(c[i][j]));

    /* ---- normalize spectrum of c ---- */

    analyse(c, nx, ny, &min, &max, &mean, &std);
    if (max != 0)
    {
        for (j = 1; j <= ny; j++)
            for (i = 1; i <= nx; i++)
                c[i][j] = c[i][j] * 255.0 / max;
    }

    /* ---- analyse filtered image ---- */

    analyse(u, nx, ny, &min, &max, &mean, &std);
    printf("filtered image:\n");
    printf("minimum:       %8.2f \n", min);
    printf("maximum:       %8.2f \n", max);
    printf("mean:          %8.2f \n", mean);
    printf("standard dev.: %8.2f \n\n", std);

    /* ---- write output image 1 (spectrum) (pgm format P5) ---- */

    /* generate comment string */
    comments[0] = '\0';
    comment_line(comments, "# Discrete Cosine Transform (spectrum)\n");
    comment_line(comments, "# menu option: %8ld\n", flag);

    /* write image */
    write_pgm(c, nx, ny, out1, comments);
    printf("output image %s successfully written\n\n", out1);

    /* ---- write output image 2 (image) (pgm format P5) ---- */

    /* generate comment string */
    comments[0] = '\0';
    comment_line(comments, "# Discrete Cosine Transform (image)\n");
    comment_line(comments, "# menu option: %8ld\n", flag);

    /* write image */
    write_pgm(u, nx, ny, out2, comments);
    printf("output image %s successfully written\n\n", out2);

    /* ---- free memory  ---- */

    disalloc_matrix(u, nx + 2, ny + 2);
    disalloc_matrix(c, nx + 2, ny + 2);

    return (0);
}
