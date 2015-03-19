#ifndef TWODWIND_H
#define TWODWIND_H

#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

class pauli //defines a pauli matrix based on pauli_type_ = {1,2,3}
{
    private:
    Eigen::Matrix2cd emptyp;

    public:
    pauli( int pauli_type_ );
    Eigen::Matrix2cd output() { return emptyp;}
};

class sys //defines the hamiltonian and computes the eigenvectors, eigenvalues and winding vector
{
    private:
    Eigen::Matrix4cd hamiltonian;
    Eigen::Vector4d eigenvalue;
    Eigen::Vector4cd eigenvector1;
    Eigen::Vector4cd eigenvector2;
    double gapdiff;

    public:
    sys(double pX_, double pY_, double mu, double delta);
    double getgap() { return gapdiff; }
    Eigen::Vector4d eigenvaloutput() {return eigenvalue;}
    Eigen::Vector4cd vecoutput1() { return eigenvector1;}
    Eigen::Vector4cd vecoutput2() { return eigenvector2;}
};


class chern
{
    private:
    std::vector<std::vector<sys>> momSpace;
    std::vector<double> pVec;
    double ChernN;
    double gapmin;

    public:
    chern( double mu_, double delta_, int pint_ );
    void linspace(double pmin, double pmax, int pint);
    double getChern() { return ChernN; }
    double getGapC() { return gapmin; }
};

class phasespace
{
    private:
    std::vector<std::vector<chern>> space;
    std::vector<double> muVec;
    std::vector<double> deltaVec;

    public:
    phasespace(double mumin_, double mumax_, int muint_, double deltamin_, double deltamax_, int deltaint_, int pint_);
    void mulinspace(double mumin_, double mumax_, int muint_);
    void deltalinspace(double deltamin_, double deltamax_, int deltaint_);
    double getmu( int munum_ ) {return muVec[munum_];}
    double getdelta( int deltanum_ ) {return deltaVec[deltanum_];}
    double getWinding(int munum_, int deltanum_);
};

#endif
