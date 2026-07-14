
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include "Eigen/Dense"
#include "Physics.h"
#include "UsefuleFunctions.h"


using namespace std;
using namespace Eigen;


vector<VectorXd> F_fine(vector<VectorXd> G, VectorXd X0, double dt, string EDP)
{
    int N = G.size();
    int d = X0.size();
    vector<VectorXd> X_fine(N+1, VectorXd(d));
    X_fine[0]=X0;


    if (EDP=="Overdamped")
    {
        for (unsigned int n=0; n<N; n++)
            X_fine[n+1] = X_fine[n] - dt*gradVF(X_fine[n])  + sigma(X_fine[n])*sqrt(dt)*G[n];
    }

    // 
    // else if (EDP=="Langevin")
    // {
    //
    //     double gamma=1;
    //     double beta_moins_1=1;
    //     double m_moins_1=1;
    //
    //     VectorXd q(d), p(d);
    //
    //     p=X0.segment(0,d);
    //     q=X0.segment(d,d);
    //
    //     VectorXd Gradient_force = gradVF(q);
    //     for (unsigned int n=0; n<G.size(); n++)
    //     {
    //
    //         p = p - dt*gradVF(q) - dt*gamma*p + sqrt(2*gamma*beta_moins_1*dt)*G[n+1];
    //
    //         q = q + dt*m_moins_1*p_1_demi;
    //
    //         X_fine[n+1].segment(0,d) = p;
    //         X_fine[n+1].segment(d,d) = q;
    //     }
    //
    //
    // }

    return X_fine ;
}




vector<VectorXd> F_moins_1(vector<VectorXd> X, double dt)
{
    vector<VectorXd> G(X.size()-1);

    for (unsigned int n=0; n<X.size()-1; n++)
    {
        G[n] = 1./(sqrt(dt)*sigma(X[n])) * (X[n+1] - X[n] + dt*gradVF(X[n]));
    }

    return G;
}


vector<VectorXd> Gaussian_vector(int N, int dimension)
{

    random_device rd;
    mt19937 generator(rd());
    // default_random_engine generator;

    normal_distribution<double> distribution(0.0, 1.0);

    vector<VectorXd> G(N, VectorXd(dimension));

    for (int n=0; n<N; n++)
    {
        for (int s=0; s< dimension; s ++)
        {
            G[n](s) = distribution(generator);
        }
    }

    return G;
}


double weight(vector<VectorXd> X, double NormeG_carre, double dt)
{
    return exp(0.5*( NormeG_carre - pow( L2_norm( F_moins_1(X,dt) ), 2) ));
}



double weight(vector<VectorXd> X, vector<VectorXd> G, double dt, int N_star)
{

    vector<VectorXd> Gcal= F_moins_1(X,dt);

    double delta=0;
    for (unsigned int n=0; n< N_star; n++)
    {
        delta+= pow(L2_norm(Gcal[n]),2) -  pow(L2_norm(G[n]),2);
    }

    return exp(-0.5*delta);
}
