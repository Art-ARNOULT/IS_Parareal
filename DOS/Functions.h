#include "Eigen/Dense"
#include <random>
using namespace std;
using namespace Eigen;



vector<VectorXd> F_fine(vector<VectorXd> G, VectorXd X0, double dt, string EDP);
vector<VectorXd> F_moins_1(vector<VectorXd>, double dt);


vector<VectorXd> Gaussian_vector(int, int);


double weight(vector<VectorXd> X, double NormeG_carre, double dt);
double weight(vector<VectorXd> X, vector<VectorXd> G, double dt, int N_star);
