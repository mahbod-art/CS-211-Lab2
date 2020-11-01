#include "../include/for_you_to_do.h"
#include <math.h>

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

#include "../include/for_you_to_do.h"
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
//we have it in the lab2 it has just been changed to C
int i, maxind, t,p , temps,j,k;
double max;
for(i=0; i<n-1; i++)
{
maxind = i;
max = fabs(A[i*n+i]);
for(t=i+1; t<n; t++)
{
if(fabs(A[i+t*n])>max)
{
maxind=t;
max=fabs(A[t*n+i]);
}
}
if(max==0)
{
return -1;
}
else if(maxind != i)
{
temps =ipiv[i];
ipiv[i] = ipiv[maxind];
ipiv[maxind] = temps;
for(p=0;p<n;p++)
{
double tempv;
tempv = A[i*n+p];
A[i*n+p] = A[maxind*n+p];
A[maxind*n+p] = tempv;
}
}for(p=i+1;p<n;p++)
{
A[p*n+i] = A[p*n+i]/A[i*n+i];
for(k=i+1;k<n;k++){
A[p*n+k] = A[p*n+k] - A[p*n+i] * A[i*n+k];
}
}
}
return 0;}



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
{int i,j;
double y[n];
//since we call this function with 'L' and 'U'we need an 'if' which we do not have in the code which is given
if(UPLO =='L')
{
y[0] = B[ipiv[0]];
for (i=1; i<n; i++)
{
double sum=0;
for (j = 0; j < i; j++)
sum += y[j] * A[i*n+j];
y[i] = B[ipiv[i]]-sum;
}
}
  else if(UPLO =='U')
{
y[n-1] =B[n-1] / A[(n - 1) * n + (n - 1)];
for (i = n-2; i >= 0; i--)
{
double sum =0;
for (j = i + 1; j < n; j++)
sum += y[j] * A[i * n + j];
y[i] = (B[i] - sum) / A[i * n + i];}
}else{
printf("something is wrong");
}
for(i = 0; i < n; i++){
B[i] = y[i];
}
return;
}
/**
 *
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 *
 **/
void mydgemm(double *A, double *B, double *C, int n, int i1, int j1, int b)
{
// my code from lab1 with some changes
if(n%2==0){
       register int reg1,reg2;
        int i , j , k , iblock, jblock , kblock;
        for( k = 0;k < j1 ; k += b)
                for( i = 0;i < i1;i += b)
                        for( j = 0;j < i1;j += b)
                        {
                                for( kblock = k;kblock < k + b && kblock <j1 ;kblock += 2)
                                        for( iblock = i;iblock <i + b && iblock < i1;iblock += 2)
                                              {
                                                register int reg1 = iblock *n + kblock;
                                                register int reg2 = reg1 + n;
                                                register double a00 = A[reg1],a01 = A[reg1 + 1],a10 = A[reg2],a11 = A[reg2 + 1];
                                                for(jblock = j;jblock < j + b && jblock < i1;jblock += 2)
                                                {
                                                        register int pb = kblock * n + jblock;
                                                        register int pb1 = n+pb ;
                                                        register int pc = iblock * n + jblock;                                                                                                  register int pc1 = n+ pc ;
                                                        register double b00 = B[pb],b10 = B[pb1];
                                                        register double c00 = C[pc],c01 = C[pc + 1],c10 = C[pc1],c11 = C[pc1 + 1];
                                                        c00 -= a00 * b00 + a01 * b10;
                                                      c10 -= a10 * b00 + a11 * b10;                                                                                                       b10= B[pb + 1];
                                                           b00=B[pb1 + 1];                                                                                                              c01 -= a00 * b10 + a01 * b00;
                                                        c11-= a10 * b10 + a11 * b00;
                                                         C[pc]=c00;C[pc + 1]=c01; C[pc1]=c10;C[pc1 + 1]=c11;

                                                }
                                        }


}}
// handling if n%2!=0 as the above code dose not cover it corectly.
else if(n%2==1){
   int i,j,k,iblock,kblock,jblock;
    for( k=0;k<j1;k+=b)
        for( j=0;j<i1;j+=b)
            for( i=0;i<i1;i+=b)
            {
 for( kblock=k;kblock<k+b && kblock<j1;kblock++)
                    for( jblock=j;jblock<j+b && jblock<i1;jblock++)
                    {
                        register double r=B[kblock*n+jblock];
                        for( iblock=i;iblock<i+b && iblock<i1;iblock++)
                            C[iblock*n+jblock]-=r*A[iblock*n+kblock];
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
int mydgetrf_block(double *A, int *ipiv, int n, int b) {
// the site "http://people.eecs.berkeley.edu/~demmel/cs267/lecture12/lecture12.html" was used .it has the explanations like the course slides and the question given and was found when trying to learn more.
int i,j,k, maxind, ib,t ,end;
double max, addition=0;
double p, change ;
for(ib = 0; ib < n - 1; ib =ib + b)
{
if(ib+b-1>n-1){
end=n-1;
}else if(ib+b-1<n-1){
end=ib+b-1;
}
for(i = ib; i <= end; i++)
{
maxind = i;
max = fabs(A[i*n+i]);
for(k = i+1; k < n; k++)
{
if(fabs(A[k * n + i]) > max)
{
maxind = k;
max = fabs(A[k * n + i]);
}
}
if(max == 0){
return -1;}
else if (maxind !=i)
{
p = ipiv[i];
ipiv[i] = ipiv[maxind];
ipiv[maxind] = p;

for(j = 0; j < n; j++)
{
change = A[i * n + j];
A[i*n+j] = A[maxind*n+j];
A[maxind*n+j] = change;
}}
for(j = i + 1; j < n; j++)
{
A[j * n + i] = A[j * n + i] / A[i * n + i];
for(t = i + 1;t <= end; t++)
{
A[j*n+t] = A[j*n+t] - A[j*n+i] * A[i*n+t];
}
}
}
for(i = ib; i <= end; i++)
{
for(k = end +1; k < n; k++)
{
 addition=0;
for(j = ib; j < i; j++)
{
addition =addition + A[i * n + j] * A [j * n + k];
}
A[i * n + k] = A[i * n + k]-addition;
}
}int num=end+1;
 mydgemm(&A[((num)*n)+ ib], &A[num +n*ib],&A[((num)*n) +(num)], n ,(n-end-1),num-ib, 64);
}
return 0;
}



