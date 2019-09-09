// this example shows how to handle composite binary operators

#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>

#include "MCSFD/MCSFDCore.h"

LOBO_MAKE_TYPEDEFS(double, t);

void function_value_acc(LoboComplext &result)
{
    LoboComplext f1;
    LoboComplext f2;
    std::vector<LoboComplext> f_list(5);
    LoboComplext x;
    LoboComplext tmp;

    double h = 1e-30;
    x.real_ = 1.5;
    x.image_ = h;
    int numtest = 1000000;

    for (int j = 0; j < f_list.size(); j++)
    {
        f_list[j] = x;
        for (int k = 0; k < j; k++)
        {
            //lobo::multi<LoboComplexDualt,double>(f_list[j],x);
            lobo::add_a(f_list[j], x);
        }
    }

    for (int i = 0; i < numtest; i++)
    {
        f_list[0].real_ += 1e-6 * i;
        lobo::multi_composite(f_list, result.image_);
    }
}

void function_value(LoboComplext &result)
{
    LoboComplext f1;
    LoboComplext f2;
    std::vector<LoboComplext> f_list(5);
    LoboComplext x;
    LoboComplext tmp;

    double h = 1e-30;
    x.real_ = 1.5;
    x.image_ = h;
    int numtest = 1000000;
    for (int j = 0; j < f_list.size(); j++)
    {
        f_list[j] = x;
        for (int k = 0; k < j; k++)
        {
            //lobo::multi<LoboComplexDualt,double>(f_list[j],x);
            lobo::add_a(f_list[j], x);
        }
    }

    for (int i = 0; i < numtest; i++)
    {
        f_list[0].real_ += 1e-6 * i;

        double real = 1.0;
        result.image_ = 0;
        for (int k = 0; k < f_list.size(); k++)
        {
            real = 1.0;
            for (int j = 0; j < f_list.size(); j++)
            {
                if (j == k)
                    real *= f_list[j].image_;
                else
                {
                    real *= f_list[j].real_;
                }
            }
            result.image_ += real;
        }
    }
}

int main()
{
    //lobo::LoboComplexMatrix3<LoboComplext, double> mat;
    //lobo::LoboComplexMatrix3<LoboComplext, double> result;
    LoboComplext result;
    clock_t t1 = clock();
    function_value(result);
    clock_t t2 = clock() - t1;
    std::cout << "MSCFD " << (double)t2 / CLOCKS_PER_SEC << "s" << std::endl;
    std::cout << result << std::endl;

    clock_t t3 = clock();
    function_value_acc(result);
    clock_t t4 = clock() - t3;
    std::cout << "acc " << (double)t4 / CLOCKS_PER_SEC << "s" << std::endl;

    std::cout << result << std::endl;

    return 0;
}
