// this example shows how to compute first order, second order of tensor function

#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>

#include "MCSFD/MCSFDCore.h"

LOBO_MAKE_TYPEDEFS(double, t);

int main()
{
    double h = 1e-30;
    lobo::LoboComplexMatrix3<LoboComplext, double> F;
    lobo::LoboComplexMatrix3<LoboComplext, double> R;

    F.setIdentity();
    //compute dR_dF_11
    F.data[0].image_ = h;

    R = F.transpose()*F;
    std::cout<< R*(1/h)<<std::endl;

    //compute dE_dF_11

    //3x3 inverse
    R = R.inverse();
    LoboComplext E = R.trace()+ lobo::inner_product(R,R);

    std::cout<< E.image_/h <<std::endl;

    return 0;
}
