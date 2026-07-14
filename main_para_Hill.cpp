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
    string MODEL, EDP;

    Parameter(N,K,T,X0,I,file_name_parameter, atoi(argv[1]), MODEL, EDP, Amoins, Aplus, Bmoins, Bplus, Sigma);
    int d=X0.size();

    int Number_of_Initial_Conditions=1000;
    uniform_int_distribution<mt19937::result_type> dist_I(1,Number_of_Initial_Conditions);

    double dt=T/N;
    cout << "dt = " << dt << endl;

    random_device rd;
    random_device dev;
    mt19937 generator(rd());
    mt19937 rng(dev());

    // default_random_engine generator;
    // normal_distribution<double> distribution(0.0, 1.0);


    int Sigma_atteint=0;

    double tau_1(0);


          //////////////
         // Etape 2  //
        //////////////


        vector<VectorXd> G(N, VectorXd(d));
        vector<VectorXd> Xnk(vector<VectorXd>(N+1, VectorXd(d)));

        double Y_1_atteint;
        int Pb=0;
        double NormeG_carre;

        ofstream Fichier_Hill;
        Fichier_Hill.open("RES/HILL_para_K="+str(K)+"_"+str(atoi(argv[1]))+".txt");
        int IA=-1;

        int j=0;
        while (Pb<I)
        {
            // if(Pb==I) // si on fait I passages dans B, on s'arrete
            //     break;

            X0=Init_X0(MODEL,dist_I(rng),d, EDP);

            Sigma_atteint=0;
            Y_1_atteint=0;
            tau_1=0;

            double w_temp=1;

            while (Y_1_atteint==0)
            {
                G=Gaussian_vector(N,d);
                Xnk=Parareal_multiTS(G,X0,dt,K);

                for (int n=1; n<=N; n++)
                {
                    tau_1+=dt;
                    double X=Bond_distance(Xnk[n]);

                    if (X<Bplus && X> Bmoins) // si on est dans B
                    {
                        w_temp *= weight(Xnk, G, dt, n);
                        Pb+=1;
                        Y_1_atteint=1;
                        if (Pb%1==0)
                            cout << Pb << endl;

                        Fichier_Hill << tau_1  << "\t" << j-IA << "\t" << w_temp << endl;

                        IA=j;
                        break;
                    }

                    else if (Sigma_atteint==0)
                    {
                        if (X>Sigma) // si on a dépasse sigma
                            Sigma_atteint=1;
                    }

                    else if (Sigma_atteint==1)
                    {
                        if (X>Amoins && X<Aplus) // si on est de retour dans A
                        {
                            Y_1_atteint=1;
                            break;
                        }
                    }
                } // end time steps


                // si on arrive jusque là c'est qu'on a pas eu de Y0
                // on va refaire une simulation de temps T en partant du dernier pas de temps
                if (Y_1_atteint==0)
                {
                    NormeG_carre=pow(L2_norm(G),2);
                    w_temp *= weight(Xnk, NormeG_carre, dt);
                    X0=Xnk[N];
                }
            }

            j++;
        } // end loop on I

        Fichier_Hill.close();

    cout << "|| FINI ||" << endl;

}
