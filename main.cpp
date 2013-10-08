#include <iostream>
#include <math.h>
#include <armadillo>
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

int main()
{
    // initialization of matrices

    int n = 400;
    clock_t start, finish;
    mat A = mat(n,n);
    double rho_max =10;
    double rho_min =0;
    double h = (rho_max - rho_min)/n;                   // step length
    double w = 5;                                       // gauge of potential strength

    for(int i=0;i<n;i++){
       for(int j=0;j<n;j++){
           if(i==j)                { A(i,j) = (2/(h*h)) + pow((i+1)*h,2)*pow(w,2) + 1/((i+1)*h); }
                             else if(j==i+1||j==i-1) { A(i,j) = -1/pow(h,2); }
                             else                    { A(i,j) = 0; }}}

// using armadillo to get eigenvalues and eigenvectors
    vec eigenval;
    mat eigenvec;
    eig_sym(eigenval,eigenvec,A);
    cout << eigenval(0) << endl;

    int k,l;
    double e = pow(10,-8); //  tolerance for the maximum value of offdiagonal elements
    int iterations = 0;
    double max_iter = pow(n,3);
    double max_off = 1;
    double s,c,t,T;
    double A_kk, A_ll, A_ik, A_il;

    start = clock();

while(fabs(max_off) > e && (double)iterations < max_iter){

        // finding the max offdiagonal element

        max_off=0.0;
        for(int i=0;i < n;i++){
            for(int j=i+1;j<n;j++){
                if(fabs(max_off) < fabs(A(i,j)))
                        { max_off = fabs(A(i,j));
                          l = i;
                          k = j; }
            }
        }


        // calculation of T,t,c and s.


        if(A(k,l)!= 0.0)
        {
            T = (A(l,l)-A(k,k))/(2*A(k,l));
            if(T >= 0.0)  {t = (-T + sqrt(1+pow(T,2)));}
            else          {t = (-T - sqrt(1+pow(T,2)));}
            c = 1/sqrt(1+pow(t,2));
            s = t*c;
        }
        else  { c = 1.0 ; s = 0.0; }

        // updating the elements of the new matrix after the similarity transform

        A_kk = A(k,k);
        A_ll = A(l,l);
        A(k,k) = pow(c,2)*A_kk - 2*c*s*A(k,l) + pow(s,2)*A_ll;
        A(l,l) = pow(s,2)*A_kk + 2*c*s*A(k,l) + pow(c,2)*A_ll;
        A(k,l) = 0;
        A(l,k) = 0;


        for(int i=0;i<n;i++){
            if(i!=k && i!=l){
                A_ik = A(i,k);
                A_il = A(i,l);
                A(i,k) = c*A_ik - s*A_il;
                A(k,i) = A(i,k);
                A(i,l) = c*A_il + s*A_ik;
                A(l,i) = A(i,l);
            }
        }

        iterations++;
}

// finding all eigenvalues with value less than 12
for(int i=0;i < n;i++)
{
    if(A(i,i) < 12)
    cout << A(i,i)<< endl;
}
}
