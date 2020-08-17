#include <iostream>
#include <cmath>
using namespace std;
/*Solving 1D steady state Adv-Diffusion
  equation, we have got a set of 4 linear
  equations
 */
const int n=4;
const double eps = 1.e-15;
const int variable = 2;
/*Preliminary pivoting strategy
  Pivoting function
 */
 double pivot(double a[][n],double b[],int i)
 {
     int j=i;
     double t=0;
     
     for(int k=i; k<n; k+=1)
     {
         double aki = fabs(a[k][i]);
         if(aki>t)
         {
             t=aki;
             j=k;
         }
     }
     if(j>i)
     {
         double dummy;
         for(int L=0; L<n; L+=1)
         {
             dummy = a[i][L];
             a[i][L]= a[j][L];
             a[j][L]= dummy;
         }
         double temp = b[j];
         b[i]=b[j];
         b[j]=temp;
     }
     return a[i][i];
 }        
 
/* Forward elimination */
void triang(double a[][n],double b[])
{
    for(int i=0; i<n-1; i+=1)
    {
        double diag = pivot(a,b,i);
        if(fabs(diag)<eps)
        {
            cout<<"zero det"<<endl;
            return;
        }
        for(int j=i+1; j<n; j+=1)
        {
            double mult = a[j][i]/diag;
            for(int k = i+1; k<n; k+=1)
            {
                a[j][k]-=mult*a[i][k];
            }
            b[j]-=mult*b[i];
        }
    }
}
/*
DOT PRODUCT OF TWO VECTORS
*/
double dotProd(double u[], double v[], int k1, int k2)
{
  double sum = 0;
  for(int i = k1; i <= k2; i += 1)
  {
     sum += u[i] * v[i];
  }
  return sum;
}
/*
BACK SUBSTITUTION STEP
*/
void backSubst(double a[][n], double b[], double x[])
{
  for(int i =  n-1; i >= 0; i -= 1)
  {
    x[i] = (b[i] - dotProd(a[i], x, i + 1,  n-1))/ a[i][i];
  }
}
/*
REFINED GAUSSIAN ELIMINATION PROCEDURE
*/
void gauss(double a[][n ], double b[], double x[])
{
   triang(a, b);
   backSubst(a, b, x);
}                               
// EXAMPLE MAIN PROGRAM
int main()
{
  double a[n][n];
  double b[n];
  double x[n];
  //hard code matrix for clarity
  a[0][0] = 1.0;
  a[0][1] = 0.0;
  a[0][2] = 0.0;
  a[0][3] = 0.0;
  
  a[1][0] = 0.0;
  a[1][1] = 2.0;
  a[1][2] = -1.0;
  a[1][3] = 0.0;
  
  a[2][0] = 0.0;
  a[2][1] = -1.0;
  a[2][2] = 2.0;
  a[2][3] = -1.0;
  
  a[3][0] = 0.0;
  a[3][1] = 0.0;
  a[3][2] = -1.0;
  a[3][3] = 1.0;
  
  
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
  b[3] = 80.0;
  
  gauss(a, b, x);
  for( int i = 0; i < n; i += 1)
  {
     cout << x[i]/variable<< "   " << endl;
  
  }
  std::cin.get();
   return 0;
}
