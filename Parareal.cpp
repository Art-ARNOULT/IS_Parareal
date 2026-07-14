#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

#include "UsefuleFunctions.h"
#include "Functions.h"
#include "Physics.h"


extern double TailleFinale;

vector<VectorXd> Parareal(
        vector<VectorXd> G,
        VectorXd X0,
        double dt,
        int K
    )
{

    int N = G.size();

    vector<vector<VectorXd>> X(K+1, vector<VectorXd>(N+1, VectorXd(X0.rows())));
    vector<VectorXd> Gold(N+1, VectorXd(X0.rows()));

  //fisrt coarse Solver
  X[0][0] = X0;
  for (int n=0; n<N; n++)
  {
     Gold[n]= - dt*gradVC(X[0][n]) ;
     X[0][n+1] = X[0][n] + Gold[n] + sigma(X[0][n])*sqrt(dt)*G[n];

     // X[0][n+1] = X[0][n] - dt*gradVC(X[0][n]) + sqrt(2*dt)*G[n];
 }


  // Parareal Solver
  for (int k=0; k<K; k++)
  {
    X[k+1][0] = X0;

    for (int n=0; n<N; n++)
    {

        X[k+1][n+1] = X[k+1][n] - Gold[n] + sigma(X[k+1][n])*sqrt(dt)*G[n];
        Gold[n] = - dt *gradVC(X[k+1][n]);
        X[k+1][n+1] += Gold[n] - dt *(gradVF(X[k][n]));

        // X[k+1][n+1] = X[k+1][n] - dt *(gradVC(X[k+1][n]) + gradVF(X[k][n]) - gradVC(X[k][n]) )
        //                                           //+ sqrt(dt)*(sigma_C(X[k+1][n]) +  sigma(X[k][n])  - sigma_C(X[k][n]))*G[n];

    }

  }


  return X[K];
}
