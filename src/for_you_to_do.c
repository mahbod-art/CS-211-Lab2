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


void mydgemm(double *A, double *B, double *C, int n, int z, int m, int b)
{
    int i, j, k;
    int i_BLOCK, j_BLOCK, k_BLOCK;
    for (k = 0; k < m; k += b)
    {
        for (i = 0; i < z; i += b)
        {
            for (j = 0; j < z; j += b)
            {
                for (k_BLOCK = k; k_BLOCK < k + b && k_BLOCK < m; k_BLOCK += 2)
                {
                    for (i_BLOCK = i; i_BLOCK < i + b && i_BLOCK < z; i_BLOCK += 2)
                    {
                        register int RA1 = i_BLOCK * n + k_BLOCK;
                        register int RA2 = RA1 + n;
                        register double a_1 = A[RA1];
                        register double a_2 = A[RA1 + 1];
                        register double a_3 = A[RA2];
                        register double a_4 = A[RA2 + 1];

                        for (j_BLOCK = j; j_BLOCK < j + b && j_BLOCK < z; j_BLOCK += 2)
                        {
                            register int RB1 = k_BLOCK * n + j_BLOCK;
                            register int RC1 = i_BLOCK * n + j_BLOCK;
                            register int RB2 = RB1 + n; 
                            register int RC2 = RC1 + n;
                            register double b_0 = B[RB1];
                            register double b_1 = B[RB1 + 1];
                            register double b_2 = B[RB2]; 
                            register double b_3 = B[RB2 + 1];
                            register double c_0 = C[RC1];
                            register double c_1 = C[RC1 + 1];
                            register double c_2 = C[RC2];
                            register double c_3 = C[RC2 + 1];

                            C[RC1] = C[RC1] - a_1 * b_0 + a_2 * b_2;
                            C[RC1 + 1] = C[RC1 + 1] - a_1 * b_1 + a_2 * b_3;
                            C[RC2] = C[RC2] - a_3 * b_0 + a_4 * b_2;
                            C[RC2 + 1] = C[RC2 + 1] - a_3 * b_1 + a_4 * b_3;
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
		int i, j, k, p, e, m_index, iblock, vars;
		double max;
        
		for (iblock = 0; iblock < n - 1; iblock += b)
		{
			e = ((n-1) > (iblock + b -1)) ? (iblock + b - 1) : n-1;
			for (i = iblock; i <= e; i++)
			{
				m_index = i;
				max = fabs(A[i*n+i]);
				for (k = i + 1; k < n; k++)
				{
					if (fabs(A[k * n + i]) > max)
					{
						m_index = k;
						max = fabs(A[k * n + i]);
					}
				}
				if (max == 0) return -1;
				else if (m_index !=i)
				{
					vars = ipiv[i];
					ipiv[i] = ipiv[m_index];
					ipiv[m_index] = vars;
					for (j = 0; j < n; j++)
					{
						double tempv;
						tempv = A[i * n + j];
						A[i * n + j] = A[m_index * n + j];
						A[m_index * n + j] = tempv;
					}
				}

				for (j = i + 1; j < n; j++)
				{
					A[j * n + i] = (double)A[j * n + i] / A[i * n + i];
					for (p = i + 1; p <= e; p++)
						A[j * n + p] = A[j * n + p] - A[j * n + i] * A[i * n + p];
				}
			}

			for (i = iblock; i <= e; i++)
			{
				for (k = e +1; k < n; k++)
				{
					double sum = 0;
					for(j = iblock; j < i; j++)
						sum = sum + A[i * n + j] * A [j * n + k];
					A[i * n + k] -= sum;
				}
			}
			mydgemm(&A[(e+1) * n + iblock], &A[iblock * n + e +1], &A[(e+1) * n + (e + 1)], n , (n - e - 1) , (e-iblock+1), 32);
		}
    return 0;
}


