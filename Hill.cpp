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

using namespace std;
void Hill(int N, int d, string MODEL, string EDP, string type, double dt, int I, double Amoins, double Aplus, double Bmoins, double Bplus, double Sigma, int numero_fichier, int K)
{
        vector<VectorXd> G(N, VectorXd(d));
        vector<VectorXd> X_n(N, VectorXd(d));
        VectorXd X0;

        int Sigma_atteint=0;
        int CI_sortie_A=0;
        double Y_1_atteint;
        double tau_1;
        int nb_non_atteint=0;
        int Pb=0;

        int Number_of_Initial_Conditions=1000;
        uniform_int_distribution<mt19937::result_type> dist_I(1,Number_of_Initial_Conditions);
        random_device dev; mt19937 rng(dev());


        ofstream Fichier_Hill;
        Fichier_Hill.open("RES/HILL_"+type+ "_"+str(numero_fichier)+".txt");
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
                if (type=="fine")
                    X_n = F_fine(G, X0, dt, EDP);
                if (type=="para")
                    X_n = Parareal(G, X0, dt, K);

                for (int n=1; n<=N; n++)
                {
                    tau_1+=dt;
                    double X=Bond_distance(X_n[n]);

                    if (X<Bplus && X> Bmoins) // si on est dans B
                    {

                        Pb+=1;
                        Y_1_atteint=1;
                        if (Pb%10==0)
                            cout << Pb << endl;

                        if (type=="fine")
                            Fichier_Hill << tau_1  << "\t" << j-IA << endl;
                        if (type=="para")
                        {
                            w_temp *= weight(X_n, G, dt, n);
                            Fichier_Hill << tau_1  << "\t" << j-IA << "\t" << w_temp << endl;
                        }
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
                    if (type=="para")
                    {
                        w_temp *= weight(X_n, G, dt, N);
                    }
                    X0=X_n[N];
                }
            }

            j++;
        } // end loop on I

        Fichier_Hill.close();

    }
