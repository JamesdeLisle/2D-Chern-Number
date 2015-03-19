#include "twodwind.h"
#define PI 3.141592

void chern::linspace(double pmin, double pmax, int pint)
{
    double step = (pmax-pmin) / (pint);
    while(pmin <= pmax)
    {
        pVec.push_back(pmin);
        pmin += step;
    }
}

chern::chern( double mu_, double delta_, int pint_ )
{
    int i,j,k,l;
    linspace(-PI,PI,pint_);
    momSpace.resize( pint_ );
    for ( auto& row : momSpace ) {
        row.reserve( pint_ );
    }

    for(i = 0; i < pint_; i++)
    {
        for( j = 0; j < pint_; j++ )
        {
            momSpace[i][j] = sys(pVec[i],pVec[j], mu_, delta_);
        }
    }
    ChernN = 0.0;
    Eigen::Matrix2cd C1, C2, C3, C4;
    double cherntemp;
    for(i = 0; i < pint_-1; i++)
    {
        for( j = 0; j < pint_-1; j++ )
        {
            C1(0,0) = momSpace[i][j].vecoutput1().adjoint().dot(momSpace[i+1][j].vecoutput1());
            C1(0,1) = momSpace[i][j].vecoutput1().adjoint().dot(momSpace[i+1][j].vecoutput2());
            C1(1,0) = momSpace[i][j].vecoutput2().adjoint().dot(momSpace[i+1][j].vecoutput1());
            C1(1,1) = momSpace[i][j].vecoutput2().adjoint().dot(momSpace[i+1][j].vecoutput2());
            C2(0,0) = momSpace[i+1][j].vecoutput1().adjoint().dot(momSpace[i+1][j+1].vecoutput1());
            C2(0,1) = momSpace[i+1][j].vecoutput1().adjoint().dot(momSpace[i+1][j+1].vecoutput2());
            C2(1,0) = momSpace[i+1][j].vecoutput2().adjoint().dot(momSpace[i+1][j+1].vecoutput1());
            C2(1,1) = momSpace[i+1][j].vecoutput2().adjoint().dot(momSpace[i+1][j+1].vecoutput2());
            C3(0,0) = momSpace[i+1][j+1].vecoutput1().adjoint().dot(momSpace[i][j+1].vecoutput1());
            C3(0,1) = momSpace[i+1][j+1].vecoutput1().adjoint().dot(momSpace[i][j+1].vecoutput2());
            C3(1,0) = momSpace[i+1][j+1].vecoutput2().adjoint().dot(momSpace[i][j+1].vecoutput1());
            C4(1,1) = momSpace[i+1][j+1].vecoutput2().adjoint().dot(momSpace[i][j+1].vecoutput2());
            C4(0,0) = momSpace[i][j+1].vecoutput1().adjoint().dot(momSpace[i][j].vecoutput1());
            C4(0,1) = momSpace[i][j+1].vecoutput1().adjoint().dot(momSpace[i][j].vecoutput2());
            C4(1,0) = momSpace[i][j+1].vecoutput2().adjoint().dot(momSpace[i][j].vecoutput1());
            C4(1,1) = momSpace[i][j+1].vecoutput2().adjoint().dot(momSpace[i][j].vecoutput2());

            ChernN += std::arg(C1.determinant()*C2.determinant()*C3.determinant()*C4.determinant())/(2*PI);
        }
    }

    double gaptemp = 100.0;
    for (i = 0; i < pint_-1; i++)
    {
        for( j = 0; j<pint_ -1; j++ )
        {

            if ( momSpace[i][j].getgap() < gaptemp )
            {
                gaptemp = momSpace[i][j].getgap();
                //std::cout << momSpace[i][j].getgap() << " ";
            }
        }
    }
    //std::cout << gaptemp << std::endl;
    gapmin = gaptemp;
    //std::cout << gapmin << std::endl << std::endl;
}


