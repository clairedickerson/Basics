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

Matrix parser1e(std::string file, int col){
  //Read in 1e integral files of structure row, col, value
  //col = which column is variable to store, determines size of matrix, make dynamic with fscan later
  std::ifstream input;
  input.open(file.c_str());
  double val;
  std::vector<double> inp;
  while (true) {
    input >> val;
    if( input.eof() ) break;
    inp.push_back(val);
  }
  input.close();
  int n = inp.size();
  int m = inp[n-col]; //size of matrix 
  Matrix Val(m,m);
  Val.setZero();
  int el=0;
  for(int i=0; i < m; i++){
        for(int j=0; j <= i; j++){
                Val(i,j) = Val(j,i) = inp[el+2];
                el=el+3;
        }
  }
 
  return Val;
}

std::vector<double> parser2e(std::string file){
  int ioff[3000] = {};
  ioff[0]=0;
  for(int z=1; z < 3000; z++){
    ioff[z] = ioff[z-1] + z;
    }
  #define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)
  std::vector<double> TEI(3000); //make this size dynamic later
  std::ifstream input;
  input.open(file.c_str());
  int i, j, k, l, ij, kl, ijkl;
  double val;
  while (true) {
    input >> i >> j >> k >> l >> val;
    //store via compound index, Yoshimine sort
    ij = INDEX(i,j);
    kl = INDEX(k,l);
    ijkl = INDEX(ij,kl);
    if( input.eof() ) break;
    TEI[ijkl] = val; 
    }
  input.close();
  return TEI;
}

Matrix orthogonalize(Matrix S){
  Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();
  Matrix lambda=(evals.col(0)).cwiseSqrt().asDiagonal(); // Construct the diagonal matrix from the square root of eigenvalues
  Matrix Ssym;
  Ssym = evecs * lambda.inverse() * evecs.transpose(); // Raise S^(-1/2) using eigenvectors and eigenvalues
  return Ssym;
}

std::tuple<Matrix, Matrix> initialD(Matrix symS, Matrix Hcore){
  Matrix F0;
  F0 = symS.transpose() * Hcore * symS;
  Eigen::SelfAdjointEigenSolver<Matrix> solver(F0);
  Matrix evecs = solver.eigenvectors();
  Matrix evals = solver.eigenvalues();
  Matrix C;
  int m=evals.size(); 
  Matrix D(m,m);
  D.setZero();
  //can exchange columns of occ MOs (first 5) or virt (last 2) and won't change E. also, signs of eigvecs doesn't matter as long as it's consistent
  C = symS * evecs;  
  for (int i=0; i < m; i++){
    for (int j = 0; j < m; ++j) {
       //only 5 occ orbitals
       for (int k = 0; k < 5; ++k) {
               D(i,j) += C(i,k) * C(j,k);
   }}}
 
  return std::make_tuple(D,F0);
}

std::tuple<double, double> energy(Matrix D,Matrix H, Matrix F, double nuc, std::vector<double> TEI){
  double el_e = (D.cwiseProduct(H+F).sum());
  double final_E = el_e + nuc;
  return std::make_tuple(el_e,final_E);
}

Matrix fock(Matrix D,Matrix H, Matrix F, double nuc, std::vector<double> TEI){
  F.setZero();
  int ioff[3000] = {};
  ioff[0]=0;
  for(int z=1; z < 3000; z++){
    ioff[z] = ioff[z-1] + z;
    }
  #define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)
  int nao = D.cols();
  int ij,kl,ijkl,ik,jl,ikjl;
  for(int i=0; i < nao; i++){
    for(int j=0; j < nao; j++) {
      F(i,j) = H(i,j);
      for(int k=0; k < nao; k++){
        for(int l=0; l < nao; l++) {
          ij = INDEX(i+1,j+1);
          kl = INDEX(k+1,l+1);
          ijkl = INDEX(ij,kl);
          ik = INDEX(i+1,k+1);
          jl = INDEX(j+1,l+1);
          ikjl = INDEX(ik,jl);  
          F(i,j) = F(i,j) + D(k,l)*(2.0*TEI[ijkl] - TEI[ikjl]);
     }
  }}}

  return F;
}

int main()
{
  //Read in integrals
  double enuc;
  std::ifstream input;
  input.open("enuc.dat");
  input >> enuc;
  input.clear();
  input.close(); 

  Matrix S,T,V;
  int size;
  size=3;
  S=parser1e("s.dat",size);
  T=parser1e("t.dat",size);
  V=parser1e("v.dat",size);
  Matrix Hcore(size,size);
  Hcore=T+V;
  //cout << endl << "Initial Core Hamiltonian:\n" << Hcore << endl;

  //Read in 2e integral file
  std::vector<double> TEI;
  TEI=parser2e("eri.dat");

  //Orthogonalize Matrix
  Matrix symS;
  symS=orthogonalize(S);
  //Initialize density, fock
  Matrix Din, Fin;
  std::tie(Din,Fin) = initialD(symS,Hcore);
  double e0, E0;
  //SCF Energy
  std::tie(e0,E0) = energy(Din,Hcore,Hcore,enuc,TEI);
  cout << "\nBeginning Hartree-Fock calculation..." << endl << "-----------------------------------------------------------------------------" << endl;
  cout << std::setw(5) << "Iter" << std::setw(18) << "El. E" << std::setw(18) << "Total E" << std::setw(18) << "Delta(E)" << std::setw(18) << "RMS(D)" << endl;
  cout << std::setw(5) << std::setprecision(12) << std::fixed << "0" << std::setw(18) << e0 << std::setw(18) << E0 << std::setw(18) << "--" << std::setw(18) << "--" << endl;

  Matrix F10;
  F10=fock(Din,Hcore,Fin,enuc,TEI);
  //SCF Cycle
  double rms_tol, delE_tol;
  rms_tol=8e-12;
  delE_tol=2e-12;
  for(int i=1; i < 40; i++){
  Matrix D,F;
  std::tie(D,F) = initialD(symS,F10);
  double e, E;
  std::tie(e,E) = energy(D,Hcore,F10,enuc,TEI);
  //Check convergence
  double del_E;
  del_E = E-E0;
  double rms;
  Matrix RMS;
  RMS=Din;
  RMS.setZero();
  for (int u=0; u < D.cols(); u++){
    for (int v=0; v < D.rows(); ++v) {
        RMS(u,v) = (D(u,v)-Din(u,v))*(D(u,v)-Din(u,v));
  }}
  rms = sqrt(RMS.sum());
  //Print results
  cout << std::setw(5) << i << std::setprecision(12) << std::setw(18) << std::fixed << e << std::setw(18) << E << std::setw(18) << del_E << std::setw(18) << rms << endl;
  if (rms < rms_tol and del_E < delE_tol){
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << std::setw(21) << "Results Converged!" << std::setprecision(2) << std::scientific << std::setw(18) << "Delta(E) Conv: " << delE_tol << std::setw(18) << "RMS(D) Conv: " << rms_tol << endl << endl;
    break;
  }

  F10.setZero();
  F10=fock(D,Hcore,F,enuc,TEI);
  E0=E;
  Din=D;

  }
 

  return 0;
}


