#include "Eigen/Dense"
using namespace std;
using namespace Eigen;

VectorXd gradVF(VectorXd);
VectorXd gradVC(VectorXd);

double Psi(vector<VectorXd>, double);

double sigma(VectorXd);
double sigma_C(VectorXd);


double GradV_S_of_r(double r);
double GradV_S_of_r_coarse(double r);

double GradV_LJ_of_r_(double r);
double GradV_WCA_of_r(double r);

double Periodic_distance(VectorXd v);
VectorXd Periodic_vector(VectorXd v);


double Bond_distance(VectorXd X);
