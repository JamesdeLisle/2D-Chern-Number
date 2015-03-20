
#include "twodwind.h"
#include <iostream>
#include <fstream>

void phasespace::mulinspace(double mumin_, double mumax_, int muint_)
{
    double step = (mumax_-mumin_) / (muint_-1);
    while(mumin_ <= mumax_)
    {
        muVec.push_back(mumin_);
        mumin_ += step;
    }
}

void phasespace::deltalinspace(double deltamin_, double deltamax_, int deltaint_)
{
    double step = (deltamax_-deltamin_) / (deltaint_-1);
    while(deltamin_ <= deltamax_)
    {
        deltaVec.push_back(deltamin_);
        deltamin_ += step;
    }
}

phasespace::phasespace(double mumin_, double mumax_, int muint_, double deltamin_, double deltamax_, int deltaint_, int pint_)
{

    mulinspace( mumin_, mumax_, muint_);
    deltalinspace( deltamin_, deltamax_, deltaint_ );

    space.resize(muint_);
    for ( auto& row : space ) {
        row.reserve(deltaint_);
    }

    std::ofstream myfile;
    myfile.open ("wind.dat");

    int i,j;
    double mu_, delta_;

    for( i = 0; i < muint_; i++)
    {
        for( j = 0; j < deltaint_; j++)
        {
            mu_ = muVec[i];
            delta_ = deltaVec[j];
            //std::cout << mu_ << delta_;
            space[i][j] = chern(mu_, delta_, pint_);
            myfile <<  getdelta(j) << " " << getmu(i) << " " << space[i][j].getChern() << " " <<  space[i][j].getGapC() << std::endl;

            //std::cout << space[i][j].getGapC() << std::endl;
        }
        std::cout << "\r" << muint_-i << "  " << std::flush;
    }
    myfile.close();

}

double phasespace::getWinding( int munum_, int deltanum_ )
{
    double out;
    out = space[munum_][deltanum_].getChern();
    //std::cout << out;

    return out;
}
