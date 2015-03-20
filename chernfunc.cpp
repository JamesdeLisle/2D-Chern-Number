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
    for(i = 0; i < pint_-1; i++)
    {
        for( j = 0; j < pint_-1; j++ )
        {
            for(k = 0; k < 2; k++)
            {
                for(l = 0; l < 2; l++)
                {
                    C1(k,l) = momSpace[i][j].eigenvectorsoutput().col(k).dot(momSpace[i+1][j].eigenvectorsoutput().col(l));
                    C2(k,l) = momSpace[i+1][j].eigenvectorsoutput().col(k).dot(momSpace[i+1][j+1].eigenvectorsoutput().col(l));
                    C3(k,l) = momSpace[i+1][j+1].eigenvectorsoutput().col(k).dot(momSpace[i][j+1].eigenvectorsoutput().col(l));
                    C4(k,l) = momSpace[i][j+1].eigenvectorsoutput().col(k).dot(momSpace[i][j].eigenvectorsoutput().col(l));
                }
            }

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


