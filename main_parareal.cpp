#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include "Eigen/Dense"


#include "data.h"
#include "Functions.h"
#include "Physics.h"
#include "Parareal.h"
#include "UsefuleFunctions.h"


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

    Parameter(N,K,T,X0,I,file_name_parameter,atoi(argv[1]));

    int d=X0.size();
    double dt=T/N;

    cout << "dt = " << dt << endl;


    vector<VectorXd> G(N, VectorXd(d));
    // vector<VectorXd> X_fine(N, VectorXd(d));
    vector<vector<VectorXd>> Xnk(K+1, vector<VectorXd>(N+1, VectorXd(d)));


    double NormeG_carre;
    double w,P;

    // vector<double> mean_time(K+1);
    // vector<double> weightk(K+1);
    // vector<double> weighted_mean_time(K+1);
    // mean_time.assign(K+1,0);
    // weightk.assign(K+1,0);
    // weighted_mean_time.assign(K+1,0);

    ofstream Fichier_CV_PSI, Fichier_CV_WEIGHT, Fichier_CV_PSI_WEIGHT;
    Fichier_CV_PSI.open("RES/MC_Parareal_PSI_"+file_name_parameter);
    Fichier_CV_WEIGHT.open("RES/MC_Parareal_Weight_"+file_name_parameter);
    Fichier_CV_PSI_WEIGHT.open("RES/MC_Parareal_Psi_Weight_"+file_name_parameter);



    for (int i=1 ; i< I+1 ; i++)
    {
        if (i%100==0)
            cout << i << endl;

        G=Gaussian_vector(N,d);
        NormeG_carre=pow(L2_norm(G),2);

        // double G2=0;
        // for (int n =0; n< N; n++)
        //     G2+=pow(G[n](0),2)  + pow(G[n](1),2);


        Xnk=Parareal(X0,dt,N,K,G);

        for (int k=0; k<K+1; k++)
        {
            w = weight(Xnk[k], NormeG_carre, dt);
            // cout << w << endl;

            // vector<VectorXd> Gcal= F_moins_1(Xnk[k],dt) ;
            //
            // double Gcal2=0;
            // for (int n =0; n< N; n++)
            //     Gcal2+=pow(Gcal[n](0),2) + pow(Gcal[n](1),2);
            //
            // w = exp(0.5*(Gcal2 - G2));
            // cout << w << endl;

            P = Psi(Xnk[k], dt);

            Fichier_CV_PSI << P << "\t";
            Fichier_CV_WEIGHT << w << "\t";
            Fichier_CV_PSI_WEIGHT << w*P << "\t";
        }

        Fichier_CV_PSI << endl;
        Fichier_CV_WEIGHT << endl;
        Fichier_CV_PSI_WEIGHT << endl;


    }
    Fichier_CV_PSI.close();
    Fichier_CV_WEIGHT.close();
    Fichier_CV_PSI_WEIGHT.close();

    cout << "|| FINI ||" << endl;

}
