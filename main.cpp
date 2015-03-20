#include "twodwind.h"
#include <iostream>

int main()
{
    double mumin_   = 0.0;     //minimum chemical potential
    double mumax_   = 6.0;      //maximum chemical potential
    int muint_      = 21;       //descritisation of chemical potential
    double deltamin_= 0.0;     //minimum pairing coupling
    double deltamax_= 4.0;      //maximum pairing coupling
    int deltaint_   = 21;       //descretisation of pairing coupling
    int pint_       = 41;      //descretisation of momentum space


    //wind test(1.0,1.0,150);
    //std::cout << test.getWind();

    //chern test(-1.0,3.0,200);
    //std::cout << test.getChern();

    //sys test(2,2,1.0,1.0);
    //std::cout << test.getgap();
    phasespace PS(mumin_, mumax_, muint_, deltamin_, deltamax_, deltaint_, pint_);
    return 0;
}
