#include "twodwind.h"
#include <complex>
#include <eigen3/Eigen/Dense>


pauli::pauli( int pauli_type_ )
{
    std::complex<double> zero(0.0,0.0);
    std::complex<double> one(1.0,0.0);
    std::complex<double> I(0.0,1.0);

    if ( pauli_type_ == 1 )
    {
        emptyp(0,0) = zero;
        emptyp(0,1) = one;
        emptyp(1,0) = one;
        emptyp(1,1) = zero;
    }
    else if ( pauli_type_ == 2 )
    {
        emptyp(0,0) = zero;
        emptyp(0,1) = -I;
        emptyp(1,0) = I;
        emptyp(1,1) = zero;
    }
    else if ( pauli_type_ == 3 )
    {
        emptyp(0,0) = one;
        emptyp(0,1) = zero;
        emptyp(1,0) = zero;
        emptyp(1,1) = -one;
    }
    else
    {
        emptyp(0,0) = one;
        emptyp(0,1) = zero;
        emptyp(1,0) = zero;
        emptyp(1,1) = one;
    }
}


