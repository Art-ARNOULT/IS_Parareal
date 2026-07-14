#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

#include "UsefuleFunctions.h"
#include "Functions.h"
#include "Physics.h"


vector<VectorXd> Parareal_multiTS(
        vector<VectorXd> G,
        VectorXd X0,
        double dt,
        int K
    )
{
    // N : nb de pas de temps total
    int N = G.size();
    int d = X0.size();

    vector<vector<VectorXd>> X(K+1, vector<VectorXd>(N+1, VectorXd(d)));

    int S=10; // nb de pas de temps par fenêtre
    int R=N/S; // nb de fenêtres
    // cout << "\t Nb de fenetres " << R << endl;

    vector<VectorXd> Gold_k(R, VectorXd(d));
    vector<VectorXd> Gold_k1(R, VectorXd(d));


  //fisrt coarse Solver
  X[0][0] = X0;
  for (int r=0; r<R; r++) // boucle sur les fenetres
  {
      Gold_k[r] = - dt*gradVC(X[0][r*S]) ;
      for (int s=0; s<S; s++) // boucle sur les pas de temps dans la fenetre
      {
          X[0][r*S+s+1] = X[0][r*S+s] +Gold_k[r] + sqrt(dt)*sigma(X[0][r*S])*G[r*S+s];
      }
 }


  // Parareal Solver
  for (int k=0; k<K; k++)
  {
      // cout << "k = " << k+1 << endl;
    X[k+1][0] = X0;
    for (int r=0; r<R; r++)
    {

        Gold_k1[r] = - dt*gradVC(X[k+1][r*S]);

        for (int s=0; s<S; s++)
        {

            X[k+1][r*S+s+1] = X[k+1][r*S+s] + Gold_k1[r] - dt*gradVF(X[k][r*S+s]) - Gold_k[r] + sqrt(dt)*sigma(X[k+1][r*S])*G[r*S+s];

        }

        Gold_k[r] = Gold_k1[r];

    }

  }

  return X[K];
}
