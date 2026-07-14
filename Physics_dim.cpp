#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include <vector>

#include "UsefuleFunctions.h"


using namespace std;
using namespace Eigen;


extern double h, w,r0,sig,epsilon, L;


double sigma(VectorXd X)
{
     return sqrt(2./0.25);
}

double sigma_C(VectorXd X)
{
    return sigma(X);
}



double GradV_S_of_r(double r)
{
    return  -4*h * (1-pow((r-r0-w)/w,2)) * (r-r0-w)/pow(w,2);
}

double GradV_S_of_r_coarse(double r)
{
    if (r<r0 + w*(1-1./sqrt(3)))
        return -8*h/(3*pow(w,2))*(r-r0)/(1-sqrt(3));

    else if (r<r0 + w*(1+1./sqrt(3)))
        return -8*h/(3*pow(w,2))*(r-r0-w);

    else
        return -8*h/(3*pow(w,2))*(r-r0-2*w)/(1-sqrt(3));
}



double GradV_LJ_of_r(double r)
{
    return 24*epsilon/sig*(-2*pow(sig/r,13) + pow(sig/r,7));
}

double GradV_WCA_of_r(double r)
{
    if (r<=r0)
        return GradV_LJ_of_r(r);
    else
        return 0;
}


VectorXd Periodic_vector(VectorXd Z)
{

    VectorXd v(2);
    double a;

    a=fmod(Z(0),L);
    if (abs(a)>abs(a-sign(a)*L))
        v[0]=(a-sign(a)*L);
    else
        v[0]=a;

    a=fmod(Z(1),L);
    if (abs(a)>abs(a-sign(a)*L))
        v[1]=(a-sign(a)*L);
    else
        v[1]=a;

    return v;
}



VectorXd gradVF(VectorXd X)
{
    int dim=X.rows();
    int Nb_atoms=dim/2;
    VectorXd GradV(dim);

    GradV.setZero(dim);

    VectorXd V_periodic=X.segment(0,2)-X.segment(2,2);
    double r = L2_norm(V_periodic);
    double Sr = GradV_S_of_r(r);

    GradV.segment(0,2) = Sr/r*V_periodic;
    GradV.segment(2,2)= - GradV.segment(0,2);


    double S_LJ;

    vector<vector<VectorXd>>
        M(Nb_atoms, vector<VectorXd>(Nb_atoms, VectorXd(2)));

    for (int s=0; s<Nb_atoms; s++)
    {
        M[s][s](0)=0; M[s][s](1)=0;
        for (int t=max(2,s+1) ; t<Nb_atoms; t++) // max(2,s+1) pour eviter d'avoir l'interaction entre 0 et 1
        {
            V_periodic=Periodic_vector(X.segment(2*s,2)-X.segment(2*t,2));
            r = L2_norm(V_periodic);
            S_LJ = GradV_LJ_of_r(r);


            M[s][t] = S_LJ/r*V_periodic;
            M[t][s] = - M[s][t];
        }

        for (int t=0 ; t<Nb_atoms; t++)
            GradV.segment(2*s,2) += M[s][t];

    }

    return GradV;
}






VectorXd gradVC(VectorXd X)
{
    int dim=X.rows();
    int Nb_atoms=dim/2;
    VectorXd GradV(dim);

    GradV.setZero(dim);

    VectorXd V_periodic=X.segment(0,2)-X.segment(2,2);
    double r = L2_norm(V_periodic);
    double Sr = GradV_S_of_r_coarse(r);

    GradV.segment(0,2) = Sr/r*V_periodic;
    GradV.segment(2,2)= - GradV.segment(0,2);


    double S_LJ;

    vector<vector<VectorXd>>
        M(Nb_atoms, vector<VectorXd>(Nb_atoms, VectorXd(2)));

    for (int s=0; s<Nb_atoms; s++)
    {
        M[s][s](0)=0; M[s][s](1)=0;
        for (int t=max(2,s+1) ; t<Nb_atoms; t++) // max(2,s+1) pour eviter d'avoir l'interaction entre 0 et 1
        {
            V_periodic=Periodic_vector(X.segment(2*s,2)-X.segment(2*t,2));
            r = L2_norm(V_periodic);
            S_LJ = GradV_WCA_of_r(r);

            M[s][t] = S_LJ/r*V_periodic;
            M[t][s] = - M[s][t];
        }

        for (int t=0 ; t<Nb_atoms; t++)
            GradV.segment(2*s,2) += M[s][t];

    }

    return GradV;
}












double Psi(vector<VectorXd> X, double dt)
{

    int signe=sign(L2_norm(X[0].segment(0,2)-X[0].segment(2,2))-(r0+w));
    // double mean_time_in_wheel=0;
    // int number_of_wheel=0;
    int n=1;
    int vec_size=X.size();
    // int temp=0;
    while(n<vec_size && signe*(L2_norm(X[n].segment(0,2)-X[n].segment(2,2))-(r0+w) )>=0 )
    { // on attends de quitter le premier puits
        n=n+1;
    }

    // signe=-signe;
    // temp=1;
    // cout << "temps avant début " << (n-1)*dt << endl;
    //
    // for (int k=n+1; k<vec_size; k++)
    // {
    //     // si on est toujours dans le même puits
    //     if (signe*(L2_norm(X[k].row(1)-X[k].row(0))-(r0+w))>=0)
    //     {
    //         temp+=1;
    //     }
    //     else // on a change de puits
    //     {
    //         // if (abs(L2_norm(X[k].row(1)-X[k].row(0))-(r0+w))< 0.1*w)
    //         // {
    //         //     // on a a peine change de  (moins de 10%), on neglige
    //         //     temp+=1;
    //         // }
    //         // else
    //         {
    //             // on a vraiment change de puits
    //             mean_time_in_wheel+=temp;
    //             temp=1;
    //             number_of_wheel+=1;
    //             signe=-signe;
    //         }
    //     }
    // }

    // cout << "=) " << number_of_wheel<< endl;

    // if (number_of_wheel==0)
    // {
    //     cout << "Pas de puits" << endl;
    //     return (vec_size-1)*dt;
    // }

    if (n==vec_size)
    {
            // cout << "Pas de puits " << (vec_size-1)*dt << endl;
            return (vec_size-1)*dt;
    }


    return n*dt;
}


double Bond_distance(VectorXd X)
{
    return L2_norm(X.segment(0,2) - X.segment(2,2)) ;
}
