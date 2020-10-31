#include "../include/for_you_to_do.h"
#include <math.h>
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
void swap(double* A, double* tmpr, int n, int r1, int r2)
{
    memcpy(tmpr, A + r1 * n, n * sizeof(double));
    memcpy(A + r1 * n, A + r2 * n, n * sizeof(double));
    memcpy(A + r2 * n, tmpr, n * sizeof(double));
}

int mydgetrf(double* A, int* ipiv, int n)
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
void mydtrsv(char UPLO, double* A, double* B, int n, int* ipiv)
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
void mydgemm(const double* A, const double* B, double* C, const int m, const int p, const int n, const int b)
{//mxp pxn = mxn
    int i = 0;
    for (i = 0; i < m; i += b)
    {
        int j = 0;
        for (j = 0; j < n; j += b)
        {
            int k = 0;
            for (k = 0; k < p; k += b)
            {
                int i1 = 0;
                for (i1 = i; i1 < (i + b > m? m : (i + b)); i1 += 3)
                {
                    int j1 = 0;
                    for (j1 = j; j1 < (j + b > n? n : (j + b)); j1 += 3)
                    {
                        register double C_0_0 = C[i1 * n + j1];
                        register double C_1_0 = C[(i1 + 1) * n + j1];
                        register double C_2_0 = C[(i1 + 2) * n + j1];

                        register double C_0_1 = C[i1 * n + (j1 + 1)];
                        register double C_1_1 = C[(i1 + 1) * n + (j1 + 1)];
                        register double C_2_1 = C[(i1 + 2) * n + (j1 + 1)];

                        register double C_0_2 = C[i1 * n + (j1 + 2)];
                        register double C_1_2 = C[(i1 + 1) * n + (j1 + 2)];
                        register double C_2_2 = C[(i1 + 2) * n + (j1 + 2)];

                        int k1 = 0;
                        for (k1 = k; k1 < (k + b > p? p : (k + b)); k1++)
                        {
                            register double A_0_M = A[i1 * p + k1];
                            register double A_1_M = A[(i1 + 1) * p + k1];
                            register double A_2_M = A[(i1 + 2) * p + k1];

                            register double B_M =  B[k1 * n + j1];
                            C_0_0 += A_0_M * B_M;
                            C_1_0 += A_1_M * B_M;
                            C_2_0 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 1)];
                            C_0_1 += A_0_M * B_M;
                            C_1_1 += A_1_M * B_M;
                            C_2_1 += A_2_M * B_M;

                            B_M = B[k1 * n + (j1 + 2)];
                            C_0_2 += A_0_M * B_M;
                            C_1_2 += A_1_M * B_M;
                            C_2_2 += A_2_M * B_M;

                        }
                        C[i1 * n + j1] = C_0_0;
                        C[(i1 + 1) * n + j1] = C_1_0;
                        C[(i1 + 2) * n + j1] = C_2_0;

                        C[i1 * n + (j1 + 1)] = C_0_1;
                        C[(i1 + 1) * n + (j1 + 1)] = C_1_1;
                        C[(i1 + 2) * n + (j1 + 1)] = C_2_1;

                        C[i1 * n + (j1 + 2)] = C_0_2;
                        C[(i1 + 1) * n + (j1 + 2)] = C_1_2;
                        C[(i1 + 2) * n + (j1 + 2)] = C_2_2;
                    
                    }
                }
            }
        }
    }
}


void mydgemm_sub(double *A, double *B, double *C, int m, int p, int n, int rowsize, int b)
{
int i = 0;
    for (i = 0; i < m; i += b)
    {
        int j = 0;
        for (j = 0; j < n; j += b)
        {
            int k = 0;
            for (k = 0; k < p; k += b)
            {
                int i1 = 0;
                for (i1 = i; i1 < (i + b > m? m : (i + b)); i1 += 3)
                {
                    int j1 = 0;
                    for (j1 = j; j1 < (j + b > n? n : (j + b)); j1 += 3)
                    {
                        register double C_0_0 = C[i1 * rowsize + j1];
                        register double C_1_0 = C[(i1 + 1) * rowsize + j1];
                        register double C_2_0 = C[(i1 + 2) * rowsize + j1];

                        register double C_0_1 = C[i1 * rowsize + (j1 + 1)];
                        register double C_1_1 = C[(i1 + 1) * rowsize + (j1 + 1)];
                        register double C_2_1 = C[(i1 + 2) * rowsize + (j1 + 1)];

                        register double C_0_2 = C[i1 * rowsize + (j1 + 2)];
                        register double C_1_2 = C[(i1 + 1) * rowsize + (j1 + 2)];
                        register double C_2_2 = C[(i1 + 2) * rowsize + (j1 + 2)];

                        int k1 = 0;
                        for (k1 = k; k1 < (k + b > p? p : (k + b)); k1++)
                        {
                            register double A_0_M = A[i1 * rowsize + k1];
                            register double A_1_M = A[(i1 + 1) * rowsize + k1];
                            register double A_2_M = A[(i1 + 2) * rowsize + k1];

                            register double B_M =  B[k1 * rowsize + j1];
                            C_0_0 -= A_0_M * B_M;
                            C_1_0 -= A_1_M * B_M;
                            C_2_0 -= A_2_M * B_M;

                            B_M = B[k1 * rowsize + (j1 + 1)];
                            C_0_1 -= A_0_M * B_M;
                            C_1_1 -= A_1_M * B_M;
                            C_2_1 -= A_2_M * B_M;

                            B_M = B[k1 * rowsize + (j1 + 2)];
                            C_0_2 -= A_0_M * B_M;
                            C_1_2 -= A_1_M * B_M;
                            C_2_2 -= A_2_M * B_M;
                        }

                        C[i1 * rowsize + j1] = C_0_0;
                        C[(i1 + 1) * rowsize + j1] = C_1_0;
                        C[(i1 + 2) * rowsize + j1] = C_2_0;

                        C[i1 * rowsize + (j1 + 1)] = C_0_1;
                        C[(i1 + 1) * rowsize + (j1 + 1)] = C_1_1;
                        C[(i1 + 2) * rowsize + (j1 + 1)] = C_2_1;

                        C[i1 * rowsize + (j1 + 2)] = C_0_2;
                        C[(i1 + 1) * rowsize + (j1 + 2)] = C_1_2;
                        C[(i1 + 2) * rowsize + (j1 + 2)] = C_2_2;                
                    }
                }
            }
        }
    }
    return;
}


inline void mydgemm_sub_best(double *ptr, int m, int rowsize, int b)
{
    register int i, j, i1, j1, k1;
    double *A = ptr + b * rowsize;
    double *B = ptr + b;
    double *C = ptr + b * rowsize + b;
    for (i = 0; i < m; i += b)
    {
        register int i1bound = i + b > m? m : (i + b);
        for (j = 0; j < m; j += b)
        {
            register int j1bound = j + b > m? m : (j + b);
            for (i1 = i; i1 < i1bound; i1 += 3)
            {
                for (j1 = j; j1 < j1bound; j1 += 3)
                {
                    register int i00 = i1 * rowsize + j1;
                    register int i10 = i00 + rowsize;
                    register int i20 = i10 + rowsize;

                    register double C_0_0 = C[i00];
                    register double C_1_0 = C[i10];
                    register double C_2_0 = C[i20];

                    register double C_0_1 = C[i00 + 1];
                    register double C_1_1 = C[i10 + 1];
                    register double C_2_1 = C[i20 + 1];

                    register double C_0_2 = C[i00 + 2];
                    register double C_1_2 = C[i10 + 2];
                    register double C_2_2 = C[i20 + 2];

                    for (k1 = 0; k1 < b; k1++)
                    {
                        register int j00 = i1 * rowsize + k1;
                        register int k00 = k1 * rowsize + j1;

                        register double A_0_M = A[j00];
                        register double A_1_M = A[j00 + rowsize];
                        register double A_2_M = A[j00 + 2 * rowsize];

                        register double B_M =  B[k00];
                        C_0_0 -= A_0_M * B_M;
                        C_1_0 -= A_1_M * B_M;
                        C_2_0 -= A_2_M * B_M;

                        B_M = B[k00 + 1];
                        C_0_1 -= A_0_M * B_M;
                        C_1_1 -= A_1_M * B_M;
                        C_2_1 -= A_2_M * B_M;

                        B_M = B[k00 + 2];
                        C_0_2 -= A_0_M * B_M;
                        C_1_2 -= A_1_M * B_M;
                        C_2_2 -= A_2_M * B_M;
                    }

                    C[i00] = C_0_0;
                    C[i10] = C_1_0;
                    C[i20] = C_2_0;

                    C[i00 + 1] = C_0_1;
                    C[i10 + 1] = C_1_1;
                    C[i20 + 1] = C_2_1;

                    C[i00 + 2] = C_0_2;
                    C[i10 + 2] = C_1_2;
                    C[i20 + 2] = C_2_2;                
                }
            }
        }
    }
    return;
}

/**
 *  
 * this function transposes a square matrix
 * 
 * syntax 
 *  
  *  input : 
 *      A     n by n , a square matrix
 *      n            , length of the whole vector / size of the whole matrix
 * 
 **/

void transpose(double* A, int m, int n)
{
    int i, j;
    double *tmp = (double*)malloc(sizeof(double) * n * m);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            tmp[i * n + j] = A[j * n + i];
        }
    }
    memcpy(A, tmp, sizeof(double) * n * m);
    free(tmp);
}

/**
 *  
 * this function computes LU factorization
 * for a non-square matrix
 * 
 * syntax 
 *  
  *  input : 
 *      A     n by n , a non-square part of the whole matrix
 *      pos          , start position of the whole matrix, both rows and columns.
 *      ipiv  1 by n , vector
 *      n            , length of the whole vector / size of the whole matrix
 *      bm           , number of rows
 *      nm           , number of columns
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_non_squrare_naive(double* A, int pos, int* ipiv, int n, int bm, int bn, int b)
{
    /* add your code here */
    int i, j, k, i1, j1;
    int bn2 = bm - bn;
    double* tmpr = (double*)malloc(sizeof(double) * n);
    double* LLT = (double*)malloc(sizeof(double) * bn * bn);
    double* AUR = (double*)malloc(sizeof(double) * bn * bn2);
    double* AURD = (double*)malloc(sizeof(double) * bn * bn2);
    double* ALLD = (double*)malloc(sizeof(double) * bn2 * bn);
    double* LL = (double*)malloc(sizeof(double) * bn * bn);
    int* ipivl = (int*)malloc(sizeof(int) * bn);

    for (i = 0; i < bn; i++)
    {
        int maxidx = i;
        double max = fabs(A[i * n + i]);
        for (j = i + 1; j < bm; j++)
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
            int newMaxidx = pos + maxidx;
            int newI      = pos + i;
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newI] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];

            swap(A-pos, tmpr, n, i, maxidx);
        }

        for (j = i + 1; j < bm; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            double A_j = A[j * n + i];
            for (k = i + 1; k < bn; k++)
            {
                A[j * n + k] -= A_j * A[i * n + k];
            }
        }
    }

    if (bn2 > 0)
    {
        memset(LLT, 0, bn * bn * sizeof(double));
        memset(LL, 0, bn * bn * sizeof(double));
        memset(ipivl, 0, bn * sizeof(int));
        memset(AUR, 0, bn * bn2 * sizeof(double));
        memset(AURD, 0, bn * bn2 * sizeof(double));
        memset(ALLD, 0, bn * bn2 * sizeof(double));

        for (i = 0; i < bn; i++)
        {
            LLT[i * bn + i] = 1;
            ipivl[i] = i;
            
            LL[i * bn + i] = 1;
            for (j = 0; j < i; j++)
            {
                LL[i * bn + j] = A[i * n + j];
            }
        }
            
        for (i = 0; i < bn; i++)
        {
            memcpy(AUR + i * bn2, A + i * n + bn, bn2 * sizeof(double));
        }

        //get LL inverse and store in LLT
        for (i = 0; i < bn; i++)
        {
            mydtrsv('L', LL, LLT + i * bn, bn, ipivl);
        }

        transpose(LLT, bn, bn);

        //A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n)
        mydgemm(LLT, AUR, AURD, bn, bn, bn2, b);

        for (i = 0; i < bn; i++)
        {
            memcpy(A + i * n + bn, AURD + i * bn2, bn2 * sizeof(double));
        }
        //A(end+1:n , end+1:n )-= A(end+1:n , ib:end) * A(ib:end , end+1:n)    
        mydgemm_sub(A + bn * n, A + bn, A + bn * n + bn, bn2, bn, bn2, n, b);
    }

    free(LLT);
    free(LL);
    free(ipivl);
    free(AUR);
    free(AURD);
    free(tmpr);
    return 0;
}

int mydgetrf_non_squrare(double* A, int pos, int* ipiv, int n, int bm, int bn, int b)
{
    /* add your code here */
    register int i, j, k, i1, j1;
    register int bn2 = bm - bn;
    double* tmpr = (double*)malloc(sizeof(double) * n);

    for (i = 0; i < bn; i++)
    {
        register int maxidx = i;
        register double max = fabs(A[i * n + i]);
        for (j = i + 1; j < bm; j++)
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
            int newMaxidx = pos + maxidx;
            int newI      = pos + i;
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newI] = ipiv[newMaxidx] ^ ipiv[newI];
            ipiv[newMaxidx] = ipiv[newMaxidx] ^ ipiv[newI];

            swap(A-pos, tmpr, n, i, maxidx);
        }

        for (j = i + 1; j < bm; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
        }
        register int biasj = bm - (bm - i - 1) % 3;
        register int biask = bn - (bn - i - 1) % 3;
        for (j = i + 1; j < biasj; j += 3)
        {
            register int i00 = j * n + i;
            register int i10 = i00 + n;
            register int i20 = i10 + n;
            register double L1 = A[i00];
            register double L2 = A[i10];
            register double L3 = A[i20];

            for (k = i + 1; k < biask; k += 3)
            {
                register int j00 = i * n + k;
                register int k00 = j * n + k;
                register int k10 = k00 + n;
                register int k20 = k10 + n;
                register double R1 = A[j00];
                register double R2 = A[j00 + 1];
                register double R3 = A[j00 + 2];

                A[k00] -= L1 * R1;
                A[k10] -= L2 * R1;
                A[k20] -= L3 * R1;

                A[k00 + 1] -= L1 * R2;
                A[k10 + 1] -= L2 * R2;
                A[k20 + 1] -= L3 * R2;

                A[k00 + 2] -= L1 * R3;
                A[k10 + 2] -= L2 * R3;
                A[k20 + 2] -= L3 * R3;
            }
            for (k = biask; k < bn; k++)
            {
                register int k00 = j * n + k;
                register double R1 = A[i * n + k];
                A[k00] -= L1 * R1;
                A[k00 + n] -= L2 * R1;
                A[k00 + n + n] -= L3 * R1;
            }   
        }
        for (j = biasj; j < bm; j++)
        {
            register int jn = j * n;
            register int in = i * n;
            register double A_j = A[jn + i];
            for (k = i + 1; k < bn; k++)
            {
                A[jn + k] -= A_j * A[in + k];
            }
        }
    }

    if (bn2 > 0)
    {
        int blocksize = bn;
        for (j = bn; j < bm; j += blocksize)
        {
            for (i = 0; i < bn; i++)
            {
                register int jbound = j + blocksize > bm ? bm : j + blocksize;
                for (j1 = j; j1 < jbound; j1++)
                {
                    register int in = i * n;
                    register double A_i_j = A[in + j1];
                    for (k = 0; k < i; k++)
                    {
                        A_i_j -= A[in + k] * A[k * n + j1];
                    }
                    A[in + j1] = A_i_j;
                }
            }
        }
        mydgemm_sub_best(A, bn2, n, b);
    }
    free(tmpr);
    return 0;
}
/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
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
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int i, j, k;

    double* Aptr = A;
    for (i = 0; i < n - b; i += b)
    {
        mydgetrf_non_squrare(Aptr, i, ipiv, n, n - i, b, b);
        Aptr += b * n + b;
    }
    int blocksize = n % b > 0 ? n % b : b;
    int bias = n - blocksize;
    mydgetrf_non_squrare(Aptr, bias, ipiv, n, blocksize, blocksize, blocksize);
    return 0;
}

int mydgetrf_block_naive(double *A, int *ipiv, int n, int b) 
{
    int i, j, k;

    double* Aptr = A;
    for (i = 0; i < n - b; i += b)
    {
        mydgetrf_non_squrare_naive(Aptr, i, ipiv, n, n - i, b, b);
        Aptr += b * n + b;
    }
    int blocksize = n % b > 0 ? n % b : b;
    int bias = n - blocksize;
    mydgetrf_non_squrare_naive(Aptr, bias, ipiv, n, blocksize, blocksize, blocksize);
    return 0;
}

