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
#include "Hill.h"

using namespace std;
double h=20;
double w=0.7;

double epsilon=1;
double sig=1;
double r0=pow(2,1./6)*sig;


double L=5;


int main(int argc, char** argv)
{
    // INTIALISATION DES PARAMETRES
    double T;
    int N,K,I;
    VectorXd X0;
    string file_name_parameter;
    double Aplus, Amoins, Bplus, Bmoins, Sigma;
    string MODEL, EDP;

    Parameter(N,K,T,X0,I,file_name_parameter, atoi(argv[1]), MODEL, EDP, Amoins, Aplus, Bmoins, Bplus, Sigma);

    int d=X0.size();

    cout << "d = " << d << endl;
    double dt=T/N;
    int Number_of_Initial_Conditions=1000;
    uniform_int_distribution<mt19937::result_type> dist_I(1,Number_of_Initial_Conditions);


    // cout << "dt = " << dt << endl;

    random_device dev; mt19937 rng(dev());

    int i=0;
    vector<VectorXd> Gn(1,VectorXd(d));
    VectorXd X_fine_n(d);

    int Sigma_atteint=0;
    int CI_sortie_A=0;

    double tau_A(0), tau_B(0), tau_0(0), tau_1(0);
    int I_B(0);

    // //   /////////////
    // //   Etape 0  //
    // // /////////////
    // // Trajectoire directe
    //
    // ofstream Fichier_Fine;
    // Fichier_Fine.open("RES/Direct"+str(atoi(argv[1]))+".txt");
    //
    // X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    // tau_1=0;
    // while (i<I)
    // {
    //     Gn = Gaussian_vector(1,d);
    //     X_fine_n = F_fine(Gn, X_fine_n, dt, EDP)[1];
    //
    //     tau_1 += dt;
    //
    //     double X = Bond_distance(X_fine_n);
    //
    //     if (X<Bplus && X> Bmoins) // si on est dans B
    //     {
    //
    //         Fichier_Fine << tau_1 << endl;
    //
    //         tau_1=0;
    //         i+=1;
    //
    //         cout << i << endl;
    //
    //         X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    //     }
    // }
    // Fichier_Fine.close();
    //
    //   ///////////////
    //  //  Etape 0' //
    // ///////////////
    // // Trajectoire directe
    // ofstream Fichier_Fine;
    // Fichier_Fine.open("RES/Direct_"+str(atoi(argv[1]))+"_prime_a_jeter.txt");

    // X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    // cout <<X_fine_n << endl;
    //
    //
    //
    // tau_1=0;
    // int fromA=1;
    // while (i<I)
    // {
    //     Gn = Gaussian_vector(1,d);
    //     X_fine_n = F_fine(Gn, X_fine_n, dt, EDP)[1];
    //     tau_1+= dt;
    //
    //     double X = Bond_distance(X_fine_n);
    //     cout <<X << endl;
    //
    //
    //     // if (fromA==1)
    //     // {
    //     //     if (X<Bplus && X> Bmoins) // si on est dans B
    //     //     {
    //     //         Fichier_Fine << tau_1 << endl;
    //     //         tau_A+=tau_1;
    //     //
    //     //         tau_1=0; i+=1; fromA=0;
    //     //     }
    //     // }
    //     // else
    //     {
    //         if (X<Aplus && X> Amoins) // si on est dans A
    //         {
    //             Fichier_Fine << tau_1 << endl;
    //             tau_A+=tau_1;
    //
    //             tau_1=0; i+=1; fromA=1;
    //
    //             if (i%1==0)
    //                 cout << "Direct prime" << tau_A/i << endl;
    //             Save_position(MODEL, X_fine_n, i/2);
    //         }
    //     }
    // }
    // // cout << "T_AB = " << tau_A/I << endl;
    //
    // Fichier_Fine.close();


    //
    // // /////////////
    // //  Etape 1  //
    // // ///////////
    // // tau_A
    //
    // i=0;
    // X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    // cout<< X_fine_n << endl;
    //
    // while (i<I)
    // {
    //     Gn= Gaussian_vector(1,d);
    //     X_fine_n = F_fine(Gn, X_fine_n, dt, EDP)[1];
    //     tau_1 += dt;
    //
    //     double X=Bond_distance(X_fine_n);
    //     // cout << "\t " << X << endl;
    //
    //     if (X<Bplus && X> Bmoins) // si on est dans B
    //     {
    //         X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    //         tau_1=0; Sigma_atteint=0;
    //     }
    //
    //     else if (Sigma_atteint==0)
    //     {
    //         if (X>Sigma) {Sigma_atteint=1;}
    //     }
    //
    //     else if (Sigma_atteint==1)
    //     {
    //         if (X>Amoins && X<Aplus)  // si on est dans A
    //         {
    //             tau_A += tau_1;
    //             tau_1 = 0;
    //
    //             if (i%1==0)
    //                 cout << i << "\t" << tau_A/(i+1) << endl;
    //             // X_fine_n=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    //             i+=1; Sigma_atteint=0;
    //
    //             Save_position(MODEL, X_fine_n, i);
    //
    //         }
    //     }
    // }


    //   //////////////
    //  // Etape 2  //
    // //////////////
    //
    // vector<VectorXd> G(N, VectorXd(d));
    // vector<VectorXd> X_fine(N, VectorXd(d));
    //
    // double Y_1_atteint;
    // int nb_non_atteint=0;
    // int Pb=0;
    //
    // ofstream Fichier_Hill;
    // Fichier_Hill.open("RES/HILL_"+str(atoi(argv[1]))+".txt");
    // int IA=-1;
    //
    // for (int j=0; j<10000*I; j++)
    // {
    //     if(Pb==I) // si on fait I passages dans B, on s'arrete
    //         break;
    //
    //     X0=Init_X0(MODEL,dist_I(rng),d, EDP);
    //
    //     Sigma_atteint=0;
    //     Y_1_atteint=0;
    //     tau_1=0;
    //
    //     while (Y_1_atteint==0)
    //     {
    //         G=Gaussian_vector(N,d);
    //         X_fine = F_fine(G, X0, dt, EDP);
    //
    //         for (int n=1; n<=N; n++)
    //         {
    //             tau_1+=dt;
    //             double X=Bond_distance(X_fine[n]);
    //
    //             if (X<Bplus && X> Bmoins) // si on est dans B
    //             {
    //
    //                 Pb+=1;
    //                 Y_1_atteint=1;
    //                 if (Pb%10==0)
    //                     cout << Pb << endl;
    //
    //
    //                 Fichier_Hill << tau_1  << "\t" << j-IA << endl;
    //                 IA=j;
    //                 break;
    //             }
    //
    //             else if (Sigma_atteint==0)
    //             {
    //                 if (X>Sigma) // si on a dépasse sigma
    //                     Sigma_atteint=1;
    //             }
    //
    //             else if (Sigma_atteint==1)
    //             {
    //                 if (X>Amoins && X<Aplus) // si on est de retour dans A
    //                 {
    //                     Y_1_atteint=1;
    //                     break;
    //                 }
    //             }
    //         } // end time steps
    //
    //
    //         // si on arrive jusque là c'est qu'on a pas eu de Y0
    //         // on va refaire une simulation de temps T en partant du dernier pas de temps
    //         if (Y_1_atteint==0)
    //         {
    //             // nb_non_atteint+=1;
    //             X0=X_fine[N];
    //         }
    //     }
    // } // end loop on I
    //
    // Fichier_Hill.close();



    string type="para";
    Hill(N, d, MODEL, EDP, type, dt, I, Amoins, Aplus, Bmoins, Bplus, Sigma, atoi(argv[1]), K);


    cout << "|| FINI ||" << endl;

}
