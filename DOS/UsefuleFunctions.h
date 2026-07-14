#include "Eigen/Dense"
#include <vector>

using namespace std;
using namespace Eigen;


using Eigen::VectorXd;
using Eigen::MatrixXd;

void print(vector<double>);
void print(VectorXd vec);
void print(vector<VectorXd> mat);
void print(MatrixXd mat);

double mean(vector<double>);

double L2_norm(vector<double>);
double L2_norm(VectorXd);
// double L2_norm(MatrixXd);
double L2_norm(vector<MatrixXd>);
double L2_norm(vector<VectorXd>);


string str(double);

double sign(double);

VectorXd Init_X0(string MODEL ,int rand, int d, string EDP);
void Save_position(string MODEL, VectorXd X0, int indice);
