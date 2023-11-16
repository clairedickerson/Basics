#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <vector>
#include <cmath>
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using std::cout; 
using std::endl;

//Get distance matrix
double distance(double x1, double y1, 
            double z1, double x2, 
            double y2, double z2)
{
    double d = sqrt(pow(x2 - x1, 2) + 
                pow(y2 - y1, 2) + 
                pow(z2 - z1, 2) * 1.0);
    std::cout << std::fixed;
    return d;
}

//Get bond angles
void bonds(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector <std::vector<double> > R, int atoms)
{
    std::vector <std::vector<double> > ex(atoms, std::vector<double>(atoms));
    std::vector <std::vector<double> > ey(atoms, std::vector<double>(atoms));
    std::vector <std::vector<double> > ez(atoms, std::vector<double>(atoms));
    for (int i=0; i < atoms; i++) {
	for (int j=0; j < i; j++) {
 		for (int k=0; k < j; k++) {
        		double norm;
        		norm = sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2) + pow(z[j] - z[i], 2) * 1.0);
      			ex[j][i] = (x[j] - x[i])/norm; 
        		ey[j][i] = (y[j] - y[i])/norm;
        		ez[j][i] = (z[j] - z[i])/norm;
			double norm2;
  			norm2 = sqrt(pow(x[j] - x[k], 2) + pow(y[j] - y[k], 2) + pow(z[j] - z[k], 2) * 1.0);
                        ex[j][k] = (x[j] - x[k])/norm2;
                        ey[j][k] = (y[j] - y[k])/norm2;
                        ez[j][k] = (z[j] - z[k])/norm2;
        	}
	 }
    }
    //Only print interesting bond angles (< 4 bohr dist.) 
    double ang; 
    for(int i=0; i < atoms; i++) {
      for(int j=0; j < i; j++) {
        for(int k=0; k < j; k++) {
        	if(R[i][j] < 4.0 && R[j][k] < 4.0){
            		ang = acos((ex[j][i] * ex[j][k])+(ey[j][i] * ey[j][k])+(ez[j][i] * ez[j][k])); 
            		cout << i <<"-"<< j <<"-"<< k << ": " << ang*(180.0/acos(-1.0)) << endl;}
        }
      }
    } 
    
   return;
}

//Calculate rotational constants from moments of inertia
void rotational(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<int> atmtype, int atoms)
{
    std::vector<double> m(20); //mass of elements up to Ca
    m[1]=1.00784;
    m[6]=12.01;
    // if m=0, add error message: Add atom mass! 
    //Find COM and translate coords
    double M = 0.0;
    for(int i=0; i < atoms; i++) M += m[atmtype[i]];
    double xcm=0;
    double ycm=0;
    double zcm=0;
    double mi;
    for(int i=0; i < atoms; i++) {
    	mi = m[atmtype[i]];
    	xcm += mi * x[i];
    	ycm += mi * y[i];
    	zcm += mi * z[i];
    }

    xcm /= M;
    ycm /= M;
    zcm /= M; 

    //Translate coords to CM
    for(int i=0; i < atoms; i++) {
        x[i] = x[i]-xcm;
        y[i] = y[i]-ycm;
        z[i] = z[i]-zcm;
    }

    double ixx;
    for(int i=0; i < atoms; i++){        	
        ixx += m[atmtype[i]] * (y[i]*y[i] + z[i]*z[i]);
    }
    double iyy;
    for(int i=0; i < atoms; i++){
        iyy += m[atmtype[i]] * (x[i]*x[i] + z[i]*z[i]);
    }
    double izz;
    for(int i=0; i < atoms; i++){
        izz += m[atmtype[i]] * (x[i]*x[i] + y[i]*y[i]);
    }
    double ixy;
    for(int i=0; i < atoms; i++){
        ixy += (m[atmtype[i]] * x[i] * y[i]);
    }
    double ixz;
    for(int i=0; i < atoms; i++){
        ixz += (m[atmtype[i]] * x[i] * z[i]);
    }
    double iyz;
    for(int i=0; i < atoms; i++){
        iyz += (m[atmtype[i]] * y[i] * z[i]);
    }
    Matrix I(3,3);
    I(0,0) = ixx;
    I(1,0) = I(0,1) = -1*ixy;
    I(1,2) = I(2,1) = -1*iyz;
    I(2,0) = I(0,2) = -1*ixz;
    I(1,1) = iyy;
    I(2,2) = izz;
    cout << "\nMoment of Inertia Tensor:\n" << I << endl;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();
    cout << "\nPrincipal moments of inertia (amu * bohr^2):" << endl;
    cout <<  evals(0,0) << ", " << evals(1,0) << ", " << evals(2,0) << endl;
    //Rotational constants
    double h = 6.626070040e-34;
    double c = 299792458;
    double pi = acos(-1.0);
    double hz;
    hz = 6.6260755E-34/(8.0 * pi * pi);
    hz /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
    hz *= 1e-6;
    cout << "\nRotational Constants (MHz):" << endl << hz/evals(0,0) << ", " << hz/evals(1,0) << ", " << hz/evals(2,0) << endl;
    return;
}


int main()
//Read in variables from .xyz file
{
  std::vector<double> xvec;
  std::vector<double> yvec;
  std::vector<double> zvec;
  std::vector<int> typvec;
  std::ifstream input("benzene.xyz");
  int atoms;
  input >> atoms;
  cout << "Number of atoms: " << atoms << endl;
  cout << "Input Coords (bohr):" << endl;
  double x, y, z;
  int atomtype;
  int count=0;
  while ( input >> atomtype >> x >> y >> z)
  {
  cout << count << " " << atomtype << " " << x << ", " << y << ", " << z << "\n" ;
  //Check if need to allocate memory for array?
  xvec.push_back(x);
  yvec.push_back(y);
  zvec.push_back(z);
  typvec.push_back(atomtype);
  count++;
  }
  input.close();
  // delete/clear up variables (to-do)  
  //Check array formed correctly
  //for (int i = 0; i < xvec.size(); i++)
  //  {cout << xvec[i] << endl;} 
  //Print Distance Matrix
  std::vector <std::vector<double> > R(atoms, std::vector<double>(atoms));
  double rr;
  cout << "\nDistance Matrix (bohr):" << endl;
  for (int i = 0; i < xvec.size(); i++){
     for(int j=0; j < i; j++) {
       rr = distance(xvec[i],yvec[i],zvec[i],xvec[j],yvec[j],zvec[j]);
       R[i][j] = rr;
       cout << " R(" << i << "," << j << "): " << R[i][j] << endl;}}
  cout << "\nBond Angles (degrees):" << endl;
  bonds(xvec,yvec,zvec,R,atoms);
  rotational(xvec,yvec,zvec,typvec,atoms);
  return 0;
}


