#include "twodwind.h"
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <eigen3/Eigen/Eigenvalues>

sys::sys(double pX_, double pY_, double mu_, double delta_)
{
    /*
    std::complex<double> eps( 2 * (cos( pX_ ) + cos( pY_)) - mu_, 0.0 );
    std::complex<double> gap1( -delta_ * sin( pX_ ), delta_ * sin( pY_ ) );
    std::complex<double> gap2( delta_ * sin( pX_ ), delta_ * sin( pY_ ) );


    hamiltonian(0,0) = eps;
    hamiltonian(0,1) = 0.0;
    hamiltonian(0,2) = gap1;
    hamiltonian(0,3) = 0.0;
    hamiltonian(1,1) = eps;
    hamiltonian(1,2) = 0.0;
    hamiltonian(1,3) = gap2;
    hamiltonian(2,2) = -eps;
    hamiltonian(2,3) = 0.0;
    hamiltonian(3,3) = -eps;

    hamiltonian(1,0) = std::conj(hamiltonian(0,1));
    hamiltonian(2,0) = std::conj(hamiltonian(0,2));
    hamiltonian(3,0) = std::conj(hamiltonian(0,3));
    hamiltonian(2,1) = std::conj(hamiltonian(1,2));
    hamiltonian(3,1) = std::conj(hamiltonian(1,3));
    hamiltonian(3,2) = std::conj(hamiltonian(2,3));
    */

    std::complex<double> ipx(0.0,pX_);
    std::complex<double> e2ipx = exp(ipx);
    std::complex<double> e2mipx = exp(-ipx);
    std::complex<double> one(1.0,0.0);
    std::complex<double> I(0.0,1.0);

    std::complex<double> eps00( mu_ + delta_ + 2.0 * cos( pY_ ), 0.0 );
    std::complex<double> eps11( mu_ - delta_ + 2.0 * cos( pY_ ), 0.0 );
    std::complex<double> eps01 = I * (1.0 + e2mipx);
    std::complex<double> del02( 0.0, 4 * sin( pY_) );
    std::complex<double> del03 = 2.0 * (1.0 - e2mipx);
    std::complex<double> del12 = -2.0 * (1.0 - e2ipx);

    hamiltonian(0,0) = eps00;
    hamiltonian(0,1) = eps01;
    hamiltonian(0,2) = del02;
    hamiltonian(0,3) = del03;
    hamiltonian(1,1) = eps11;
    hamiltonian(1,2) = del12;
    hamiltonian(1,3) = del02;
    hamiltonian(2,2) = -eps00;
    hamiltonian(2,3) = eps01;
    hamiltonian(3,3) = -eps11;

    hamiltonian(1,0) = std::conj(hamiltonian(0,1));
    hamiltonian(2,0) = std::conj(hamiltonian(0,2));
    hamiltonian(3,0) = std::conj(hamiltonian(0,3));
    hamiltonian(2,1) = std::conj(hamiltonian(1,2));
    hamiltonian(3,1) = std::conj(hamiltonian(1,3));
    hamiltonian(3,2) = std::conj(hamiltonian(2,3));

    hamiltonian = 0.5 * hamiltonian;

    Eigen::Matrix4cd HermT;

    HermT = hamiltonian - hamiltonian.adjoint();

    //std::cout << HermT;

    Eigen::ComplexEigenSolver<Eigen::Matrix4cd> ces(hamiltonian);

    Eigen::Vector4cd eigenvaltemp;
    eigenvaltemp = ces.eigenvalues();
    Eigen::Vector4d sorteigenvaltemp;

    Eigen::Vector4cd temp;

    sorteigenvaltemp[0] = std::real(eigenvaltemp[0]);
    sorteigenvaltemp[1] = std::real(eigenvaltemp[1]);
    sorteigenvaltemp[2] = std::real(eigenvaltemp[2]);
    sorteigenvaltemp[3] = std::real(eigenvaltemp[3]);
    //std::cout << sorteigenvaltemp;

    Eigen::Vector4d refer(0,1,2,3);
    int flag = 1;
    double tempeva, temprefer;
    for (int i = 0; (i < 4) && flag; i++)
    {
        flag = 0;
        for (int j = 0; j < 3; j++)
        {
            if (sorteigenvaltemp[j+1] < sorteigenvaltemp[j])
            {
                tempeva = sorteigenvaltemp[j];
                sorteigenvaltemp[j] = sorteigenvaltemp[j+1];
                sorteigenvaltemp[j+1] = tempeva;
                temprefer = refer[j];
                refer[j] = refer[j+1];
                refer[j+1] = temprefer;

                flag = 1;
            }
        }
    }

    //std::cout << refer << std::endl;

    for (int k = 0; k < 4; k++)
    {
        if ((refer[k] = 0))
        {
            eigenvector1 = ces.eigenvectors().col(k);
        }
    }

    for (int k = 0; k < 4; k++)
    {
        if ((refer[k] = 1))
        {
            eigenvector2 = ces.eigenvectors().col(k);
        }
    }

    eigenvalue = sorteigenvaltemp;
    //std::cout << eigenvalue << std::endl;
    gapdiff = std::real(eigenvalue[2] - eigenvalue[1]);

    //std::cout << gapdiff << " ";
}

