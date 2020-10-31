#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k;
    double* tmpr = (double*)malloc(sizeof(double) * n);
    for (i = 0; i < n; i++)
    {
        int maxidx = i;
        double max = fabs(A[i * n + i]);
        for (j = i + 1; j < n; j++)
        {
            double tmp = fabs(A[j * n + i]);
            if (tmp - max > 1e-6)
            {
                maxidx = j;
                max = tmp;
            }
        }

        //too small pivot is also unacceptable
        if (fabs(max - 0.0) < 1e-3) 
            return -1;

        if (maxidx != i)
        {
            ipiv[maxidx] = ipiv[maxidx] ^ ipiv[i];
            ipiv[i] = ipiv[maxidx] ^ ipiv[i];
            ipiv[maxidx] = ipiv[maxidx] ^ ipiv[i];

            swap(A, tmpr, n, i, maxidx);
        }

        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            double A_j = A[j * n + i];
            for (k = i + 1; k < n; k++)
            {
                A[j * n + k] -= A_j * A[i * n + k];
            }
        }
    }
    free(tmpr);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double *newB = (double*) malloc(n * sizeof(double));
    int i, j;
    if (UPLO == 'L')
    {
        for (i = 0; i < n; i++)
        {
            newB[i] = B[ipiv[i]];
        }

        for (i = 0; i < n; i++)
        {
            double sub = newB[i];
            for (j = 0; j < i; j++)
            {
                sub -= B[j] * A[i * n + j];
            }
            B[i] = sub;
        }
    }
    else
    {//avoid cache miss
        for (i = n-1; i >=0; i--)
        {
            double sum = 0;
            for (j = i+1; j < n; j++)
            {
                sum += B[j] * A[i * n + j];
            }
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    free(newB);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    return 0;
}

void swap(double* A, double* tmpr, int n, int r1, int r2)
{
    memcpy(tmpr, A + r1 * n, n * sizeof(double));
    memcpy(A + r1 * n, A + r2 * n, n * sizeof(double));
    memcpy(A + r2 * n, tmpr, n * sizeof(double));
}