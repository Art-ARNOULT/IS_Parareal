#include <sstream>
#include "UsefuleFunctions.h"
#include <iostream>
#include <fstream>

#include "Eigen/Dense"


using namespace std;
using namespace Eigen;

extern double r0, L;

void Parameter(int &N,
    int &K,
    double &T,
    VectorXd &X0,
    int &I,
    string &file_name_parameter,
    int Numero_fichier,
    string &MODEL,
    string &EDP,
    double &Amoins,
    double &Aplus,
    double &Bmoins,
    double &Bplus,
    double &Sigma)
    {

        MODEL="1D";
        EDP="Overdamped";

        N=200;
        K=10;
        T=1;
        I=1250;

        int Nb_atoms=1;
        X0.resize(Nb_atoms);

        if (EDP=="Overdamped")
            X0.resize(Nb_atoms);
        else if (EDP=="Langevin")
            X0.resize(2*Nb_atoms);



        Amoins=-1.2;
        Aplus=-0.8;

        Bmoins=0.8;
        Bplus=1.2;

        Sigma=-0.75;




        cout << endl;
        cout << "\t \t -----------------------" << endl;
        cout << "\t \t ||\t" << MODEL << endl;
        cout << "\t \t ||\t T = " << T << endl;
        cout << "\t \t ||\t N = " << N << endl;
        cout << "\t \t ||\t K = " << K << endl;
        cout << "\t \t ||\t I = "<< I << endl;
        cout << "\t \t -----------------------" << endl;



        cout << "A = [" <<  Amoins << "," << Aplus<< "] " << endl;
        cout << "B = [" <<  Bmoins << "," << Bplus<< "] " << endl;

        cout << "Sigma = " << Sigma << endl;



        file_name_parameter="T="+ str(T)
        +"_N="+str(N)
        // + "_X0="+str(X0)
        + "_K="+str(K)
        + "_I="+str(I)
        +"_Nb_atoms="+str(Nb_atoms)
        +"_" + str(Numero_fichier)
        +".txt";

    }
