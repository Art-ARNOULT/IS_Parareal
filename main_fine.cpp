#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include "Eigen/Dense"


#include "data.h"
#include "Functions.h"
#include "Physics.h"
#include "Parareal.h"
#include "Parareal_multiTS.h"
#include "UsefuleFunctions.h"

using namespace std;
double h=20;
double w=0.7;

double epsilon=1;
double sig=1;
double r0=pow(2,1./6)*sig;

double L=5;





int main(int argc, char** argv)
{
    // int Monte_Carlo_frequence_save=1;
    // cout << "Monte_Carlo_frequence_save = " <<  Monte_Carlo_frequence_save << endl;

    double T;
    int N,K,I;
    VectorXd X0;
    string file_name_parameter;
    double Aplus, Amoins, Bplus, Bmoins, Sigma;
    string MODEL;

    Parameter(N,K,T,X0,I,file_name_parameter, atoi(argv[1]), MODEL, Amoins, Aplus, Bmoins, Bplus, Sigma);
    int d=X0.size();


    double dt=T/N;
    cout << "dt = " << dt << endl;

    cout << X0 << endl;


    vector<VectorXd> G(N, VectorXd(d));
    vector<VectorXd> X_fine(N, VectorXd(d));

    ofstream Fichier_Fine;
    Fichier_Fine.open("RES/Fine_Psi_"+file_name_parameter);

    vector<vector<VectorXd>> Xnk(K+1, vector<VectorXd>(N+1, VectorXd(d)));


    for (int i=1 ; i< I+1 ; i++)
    {
        // if (i%10==0)
            cout << i << endl;

        G=Gaussian_vector(N,d);
        X_fine = F_fine(G, X0, dt);
        cout << "Fin : " << X_fine[N] << endl;
        // // Fichier_Fine << Psi(X_fine, dt) << endl;
        // // cout << Psi(X_fine, dt) << endl;
        // cout << "FIN" << X_fine[N] << endl;
        //
        Xnk=Parareal(X0,dt,N,K,G);
        cout<< "Para : " <<Xnk[K][N]  << endl;

        Xnk=Parareal_multiTS(X0,dt,N,K,G);
        cout<< "Para MS : " <<Xnk[K][N]  << endl;

        cout << " - - - - " << endl;
        // cout << X_fine[N] << endl;
        // cout << 100*abs(Bond_distance(X_fine[N]) - Bond_distance(Xnk[K][N]))/Bond_distance(X_fine[N]) << endl;


    }
    // cout << "COUCOU" << X_fine[N] << endl;
    Fichier_Fine.close();

    cout << "|| FINI ||" << endl;

}
