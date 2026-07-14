#include <sstream>
#include "UsefuleFunctions.h"
#include <iostream>
#include <fstream>

#include "Eigen/Dense"


using namespace std;
using namespace Eigen;

extern double r0, w, L;


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

        MODEL="Dimere_L="+str(L);
        EDP = "Overdamped";
        // cout << MODEL<< endl;


        // N=200;
        // K=2;
        // T=0.1;
        // I=30;

        // Nb_atoms=4;
        // X0.resize(Nb_atoms,2);
        //
        // X0(0,0)=0; X0(0,1)=0;
        // X0(1,0)=1; X0(1,1)=1;
        // X0(2,0)=0; X0(2,1)=1.1;
        // X0(3,0)=1.1; X0(3,1)=0;



        N = 100;
        K = 0;
        T = 0.005;
        I = 3;

        int Nb_atoms=16;
        X0.resize(2*Nb_atoms);

        // Etat compact
        Amoins = r0 - 0.25*r0;
        Aplus  = r0 + 0.25*r0;

        // Etat etendu
        Bmoins = 2*w + r0 - 0.25*r0;
        Bplus  = 2*w + r0 + 0.25*r0;

        Sigma  = r0 + 0.3*r0;



        double size=L/4;
        cout << "size = " << size << endl;
        X0(0)=size; X0(1)=size;
        X0(2)=size; X0(3)=2*size;
        int offset=0;
        for (int i = 0; i < Nb_atoms-2; i++)
        {
              if(i==5)
                  offset=2;
              // if(i==6)
              //     offset=2;
              double x=(i+offset)/4;
              double y=(i+offset)%4;
              X0(2*(i+2))=size*x; X0(2*(i+2)+1)=size*y;
              // cout << i << "\t" << x <<"\t" << y << endl;

        }

        // cout << X0 << endl;

        //
        // ofstream Fichier_Initial_Positions;
        // Fichier_Initial_Positions.open("POSITION/Position_initiale_L="+str(L));
        // for (int s=0; s<Nb_atoms; s++)
        // {
        //     Fichier_Initial_Positions << X0(s, 0) << "\t" <<  X0(s, 1) << endl;
        // }
        // Fichier_Initial_Positions.close();




        // ifstream Fichier_Initial_Positions;
        // Fichier_Initial_Positions.open("POSITION/Position_initiale_L="+str(L));
        // for (unsigned int s=0; s<Nb_atoms; s++)
        // {
        //     Fichier_Initial_Positions >> X0(2*s) >> X0(2*s+1);
        // }
        // Fichier_Initial_Positions.close();

        cout << "\t \t -----------------------" << endl;
        cout << "\t \t ||\t T = " << T << endl;
        cout << "\t \t ||\t N = " << N << endl;
        cout << "\t \t ||\t K = " << K << endl;
        cout << "\t \t ||\t I = "<< I << endl;
        cout << "\t \t -----------------------" << endl;




        file_name_parameter="T="+ str(T)
        +"_N="+str(N)
        // + "_X0="+str(X0)
        + "_K="+str(K)
        + "_I="+str(I)
        +"_Nb_atoms="+str(Nb_atoms)
        +"_" + str(Numero_fichier)
        +".txt";



        cout << "A = [" <<  Amoins << "," << Aplus<< "] " << endl;
        cout << "B = [" <<  Bmoins << "," << Bplus<< "] " << endl;

        cout << "Sigma = " << Sigma << endl;





    }
