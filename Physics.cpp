#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;


// extern double h, w,r0,sig,epsilon;

double sigma(VectorXd X)
{
    return sqrt(2);
}

double sigma_C(VectorXd X)
{
    return sigma(X);
}



VectorXd gradVF(VectorXd X)
{
    VectorXd  GradV(X.rows());

    double x_carre=pow(X(0),2);
    GradV(0) = 6*(X(0))*(x_carre-1)*(x_carre+3)/pow((x_carre+1),2);

    return GradV;
}
//
// VectorXd gradVC(VectorXd X)
// {
//
//     VectorXd GradV(X.rows());
//
//
//     double x=X(0);
//     double x_carre=pow(x,2);
//
//     if ((-2.26 <x) && (x<0) )
//         GradV(0)= -8*pow((-x-1),3) - 4*(-x-1);
//     else
//     {
//         if ((0<x) && (x < 2.26 ))
//             GradV(0)=8*pow((x-1),3) + 4*(x-1);
//         else
//             GradV(0)=6*X(0)*(x_carre-1)*(x_carre+3)/pow((x_carre+1),2);
//     }
//
//     return GradV;
// }


VectorXd gradVC(VectorXd X)
{

    VectorXd GradV(X.rows());

    double x=X(0);
    if (abs(x)>0.3)
        GradV(0)=gradVF(X)[0];

    else
        GradV(0)=0;

    return GradV;
}



double Psi(vector<VectorXd> X, double dt)
{

    double mean_time_in_wheel=0;
    int number_of_wheel=0;
    int signe=X[0](1);
    int n=1;
    int vec_size=X.size();
    int temp=0;

    while( n<vec_size && signe*(X[n](1)) >=0 )
    { // on attends de quitter le premier puits
        n=n+1;
    }

    signe=-signe;
    temp=1;
    // cout << "temps avant début " << (n-1)*dt << endl;

    for (int k=n+1; k<vec_size; k++)
    {
        // si on est toujours dans le même puits
        if (signe*(X[k](0))>=0)
        {
            temp+=1;
        }
        else // on a change de puits
        {
            // if (abs(L2_norm(X[k].row(1)-X[k].row(0))-(r0+w))< 0.1*w)
            // {
            //     // on a a peine change de  (moins de 10%), on neglige
            //     temp+=1;
            // }
            // else
            {
                // on a vraiment change de puits
                mean_time_in_wheel+=temp;
                temp=1;
                number_of_wheel+=1;
                signe=-signe;
            }
        }
    }

    // cout << "=) " << number_of_wheel<< endl;

    if (number_of_wheel==0)
    {
        cout << "Pas de puits" << endl;
        return (vec_size-1)*dt;
    }

    if (n==vec_size)
    {
            // cout << "Pas de puits " << (vec_size-1)*dt << endl;
            return (vec_size-1)*dt;
    }

    return mean_time_in_wheel*dt/number_of_wheel;



    return n*dt;
}



double Bond_distance(VectorXd X)
{
    return X(0);
}
