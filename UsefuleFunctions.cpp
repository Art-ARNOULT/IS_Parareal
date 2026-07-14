#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
// #include <sstream>
#include <iomanip> // pour ficxer le nombre de chiffres significatifs

#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace std;
using namespace Eigen;


void print(vector<double> vec)
{
    cout << fixed << showpos << setprecision(5);

    for (int k=0; k<vec.size(); k++)
      cout << "\t" << vec[k] << endl;
}
void print(VectorXd vec)
{
    cout << fixed << showpos << setprecision(0);

    for (int k=0; k<vec.size(); k++)
      cout << "\t" << vec(k) << endl;
}

void print(vector<VectorXd> mat)
{

    cout << fixed << showpos << setprecision(5);

    for (int n=0; n<mat.size(); n++)
    {
        for (int s =0; s < mat[0].size(); s++)
        {
            cout << mat[n](s) << "\t";
        }
        cout <<  "\n";
    }
}


double mean(vector<double> vec)
{
    double m=0;//
    for (unsigned int k=0; k<vec.size(); k++)
      m+=vec[k];
     m=m/vec.size();

     return m;
}




double L2_norm(vector<double> vec)
{
    double m=0;
    for (unsigned int k=0; k<vec.size(); k++)
        m+=vec[k]*vec[k];

    return sqrt(m);
}




double L2_norm(VectorXd vec)
{
    double m=0;
    for (unsigned int k=0; k<vec.size(); k++)
        m+=vec(k)*vec(k);

    return sqrt(m);
}


// double L2_norm(MatrixXd M)
// {
//     double m=0;
//     for (unsigned int k=0; k<M.rows(); k++)
//     {
//         for (unsigned int l=0; l<M.cols(); l++)
//         {
//             // cout << "!" << endl;
//             m+=pow(M(k,l),2);
//
//         }
//     }
//
//     return sqrt(m);
// }
//
// double L2_norm(vector<MatrixXd> M)
// {
//     double m=0;
//     for (unsigned int n=0; n<M.size(); n++)
//     {
//             m+=pow(L2_norm(M[n]),2);
//     }
//
//     return sqrt(m);
// }
//


double L2_norm(vector<VectorXd> M)
{
    double m=0;
    for (unsigned int n=0; n<M.size(); n++)
    {
        for (unsigned int p=0; p<M[0].rows(); p++)
        {
            m+=pow(M[n](p),2);
        }
    }

    return sqrt(m);
}




string str(double value) {
    ostringstream oss;
    oss << value;
    return oss.str();
}

double sign(double x)
{
    if (x >= 0) return 1;
    else return -1;
}



VectorXd Init_X0(string MODEL, int rand, int d, string EDP)
{
    VectorXd X0(d);
    if (EDP=="Overdamped")
    {
        ifstream Fichier_Initial_Positions;
        Fichier_Initial_Positions.open("../IS_parareal/POSITION_HILL_"+MODEL+"/"+MODEL+"_"+str(rand)+".txt");
        for (int s=0; s<d; s++)
            Fichier_Initial_Positions >> X0(s);
        Fichier_Initial_Positions.close();
    }
    else { cout << "PB CONDITION INITIALE" << endl;}

    // X0(0)=-1.1;

    return X0;
}




void Save_position(string MODEL, VectorXd X0, int indice)
{
    ofstream Fichier_Initial_Positions;
    Fichier_Initial_Positions.open("POSITION_HILL_"+MODEL+"/"+MODEL+"_"+str(indice)+".txt");
    Fichier_Initial_Positions << X0 << endl;
    Fichier_Initial_Positions.close();

}
