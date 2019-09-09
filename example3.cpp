// this example shows how to compute first order, second order and third order of function
// y = exp(x)/(x^4+x^2+1.0)

#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>

#include "MCSFD/MCSFDCore.h"

LOBO_MAKE_TYPEDEFS(double, t);

int main()
{
    double h = 1e-30;
    LoboComplext result;
    LoboComplext x;
    x.real_ = 4.0;
    x.image_ = h;

    result = lobo::exp(x)/(x*x*x*x+x*x+1.0);

    std::cout<<"first order"<<std::endl;
    std::cout<<result.image_/h<<std::endl;

    LoboComplexDualt result_d;
    LoboComplexDualt x_d;
    x_d.real_.real_ = 4.0;
    x_d.real_.image_ = h;
    x_d.image_.real_ = h;

    result_d = lobo::exp(x_d)/(x_d*x_d*x_d*x_d+x_d*x_d+1.0);
    std::cout<<"second order"<<std::endl;
    std::cout<<result_d.image_.image_/h/h<<std::endl;

    LoboComplexTrit result_t;
    LoboComplexTrit x_t;
    x_t.real_.real_.real_ = 4.0;
    x_t.image_.real_.real_ = h;
    x_t.real_.real_.image_ = h;
    x_t.real_.image_.real_ = h;

    result_t = lobo::exp(x_t)/(x_t*x_t*x_t*x_t+x_t*x_t+1.0);
    std::cout<<"third order"<<std::endl;
    std::cout<<result_t.image_.image_.image_/h/h/h<<std::endl;

    return 0;
}
