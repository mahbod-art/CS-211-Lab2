#include "../include/for_you_to_do.h"
#include <math.h>

int get_block_size()
{
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
    int i, j, k, max_of_index, var;
    for (i = 0; i < n - 1; i++)
    {
        max_of_index = i;
        double max = fabs (A[i * n + i]);
        int m;
        for (m = i + 1; m < n; m++)
        {
            if (fabs(A[i + m * n]) > max)
            {
                max_of_index = m;
                max = fabs(A[m * n + i]);
            }
        }

        if (max == 0) 
            return -1;
        else if (max_of_index != i)
        {
            //swap rows
            var = ipiv[i];
            ipiv[i] = ipiv[max_of_index];
            ipiv[max_of_index] = var;
            
            int l;
            for (l = 0; l < n; l++)
            {
                double tempv;
                tempv = A[i * n + l];
                A[i * n + l] = A[max_of_index * n + l];
                A[max_of_index * n + l] = tempv;
            }
        }

        //factorization part
        int z;
        for (z = i + 1; z < n; z++)
        {
            A[z * n + i] = A[z * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++)
                A[z * n + k] -= A[z * n + i] * A[i * n + k];
        }
    }
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
    int i, j;
    double *new_var = (double*) malloc(n * sizeof(double));
    if (UPLO == 'L')
    {
        for (i = 0; i < n; i++)
            new_var[i] = B[ipiv[i]];

        for (i = 0; i < n; i++)
        {
            double sub = new_var[i];
            for (j = 0; j < i; j++)
                sub -= B[j] * A[i * n + j];
            B[i] = sub;
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            double sum = 0;
            for (j = i+1; j < n; j++)
                sum = sum + B[j] * A[i * n + j];
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    free(new_var);
    return;
}

/**
 *
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 *
 **/


void mydgemm(double *A, double *B, double *C, int n, int matx, int maty, int b)
{
	int i , j , k , i_BLOCK, j_BLOCK , k_BLOCK;
    for( k = 0;k < maty ; k += b)
    {
        for( i = 0;i < matx;i += b)
        {
            for( j = 0;j < matx ;j += b)
            {
                for( k_BLOCK = k;k_BLOCK < k + b && k_BLOCK < maty;k_BLOCK += 2)
                {
                    for( i_BLOCK = i;i_BLOCK <i + b && i_BLOCK < matx;i_BLOCK += 2)
                    {
                        register int A1 = i_BLOCK *n + k_BLOCK;
                        register int A2 = A1 + n;
                        register double a0 = A[A1],a1 = A[A1 + 1],a2 = A[A2],a3 = A[A2 + 1];
                        for(j_BLOCK = j;j_BLOCK < j + b && j_BLOCK < matx;j_BLOCK += 2)
                        {
                            register int B1 = k_BLOCK * n + j_BLOCK,C1 = i_BLOCK * n + j_BLOCK;
                            register int B2 = B1 + n,C2 = C1 + n;
                            register double b0 = B[B1],b1 = B[B1 + 1],b2 = B[B2],b3 = B[B2 + 1];
                            register double c0 = C[C1],c1 = C[C1 + 1],c2 = C[C2],c3 = C[C2 + 1];
                            C[C1] -= a0 * b0 + a1 * b2;
                            C[C1+1] -= a0 * b1 + a1 * b3;
                            C[C2] -= a2 * b0 + a3 * b2;
                            C[C2+1] -= a2 * b1 + a3 * b3;

                        }
                    }
                }
            }
        }
    }
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
		int i, j, k, e, max_index, i_BLOCK, var, t;
		double max;
		for(i_BLOCK = 0; i_BLOCK < n - 1; i_BLOCK += b)
		{
			e = ((n-1) > (i_BLOCK + b -1)) ? (i_BLOCK + b - 1) : n-1;
			for(i = i_BLOCK; i <= e; i++)
			{
				max_index = i;
				max = fabs(A[i*n+i]);
				for(k = i+1; k < n; k++)
				{
					if(fabs(A[k * n + i]) > max)
					{
						max_index = k;
						max = fabs(A[k * n + i]);
					}
				}
				if(max == 0) return -1;
				else if (max_index !=i)
				{

					var = ipiv[i];
					ipiv[i] = ipiv[max_index];
					ipiv[max_index] = var;
					for(j = 0; j < n; j++)
					{
						double tempv;
						tempv = A[i * n + j];
						A[i * n + j] = A[max_index * n + j];
						A[max_index * n + j] = tempv;
					}
				}

				for(j = i + 1; j < n; j++)
				{
					A[j * n + i] = (double)A[j * n + i] / A[i * n + i];
					for(t = i + 1;t <= e; t++)
					{
						A[j*n+t] = A[j*n+t] - A[j*n+i] * A[i*n+t];
					}
				}
			}

			for(i = i_BLOCK; i <= e; i++)
			{
				for(k = e +1; k < n; k++)
				{
					double sum = 0;
					for(j = i_BLOCK; j < i; j++)
					{
						sum += A[i * n + j] * A [j * n + k];
					}
					A[i * n + k] -= sum;
				}
			}
			mydgemm(&A[(e+1) * n + i_BLOCK], &A[i_BLOCK * n + e +1], &A[(e+1) * n + (e + 1)], n , (n - e - 1) , (e-i_BLOCK+1), 32);
		}
    return 0;
}


