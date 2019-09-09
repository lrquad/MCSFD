#pragma once
#include "MCSFD/LoboComplex.h"
#include "MCSFD/LoboComplexMatrix.h"

namespace lobo
{

inline void multi_composite(const std::vector<LoboComplexS> &a,double &result)
{
    double real = a[0].real_;
    result = 0;
    for(int i=1;i<a.size();i++)
    {
        real*=a[i].real_;
    }
    for(int i=0;i<a.size();i++)
    {
        result+=real/a[i].real_*a[i].image_;
    }
}

inline void multi_all(const LoboComplexD &a, const LoboComplexD &b, LoboComplexD &c)
{
    c.real_.real_ += a.real_.real_ * b.real_.real_;
    c.real_.image_ += a.real_.image_ * b.real_.real_ + a.real_.real_ * b.real_.image_;
    c.image_.real_ += a.image_.real_ * b.real_.real_ + a.real_.real_ * b.image_.real_;
    c.image_.image_ += a.image_.image_ * b.real_.real_ + a.real_.image_ * b.image_.real_ + a.image_.real_ * b.real_.image_ + a.real_.real_ * b.image_.image_;
}

inline void multi_image(const LoboComplexD &a, const LoboComplexD &b, LoboComplexD &c)
{
    c.image_.image_ += a.image_.image_ * b.real_.real_ + a.real_.image_ * b.image_.real_ + a.image_.real_ * b.real_.image_ + a.real_.real_ * b.image_.image_;
}

inline void multi_image(const LoboComplexS &a, const LoboComplexS &b, LoboComplexS &c)
{
    c.image_ += a.image_ * b.real_ + a.real_ * b.image_;
}


inline void multi_real_image(const LoboComplexD &a, const LoboComplexD &b, LoboComplexD &c)
{
    c.image_.real_ += a.image_.real_ * b.real_.real_ + a.real_.real_ * b.image_.real_;
    c.image_.image_ += a.image_.image_ * b.real_.real_ + a.real_.image_ * b.image_.real_ + a.image_.real_ * b.real_.image_ + a.real_.real_ * b.image_.image_;
}

inline void multiMTd_all(const LoboComplexMatrix3dd &lhs, const LoboComplexMatrix3dd &rhs, LoboComplexMatrix3dd &result)
{
    //LoboComplexMatrix<COMPLEX_TYPE, TYPE> result(row, col);

    memset(result.data.data(), 0, sizeof(LoboComplexDualt) * 3 * 3);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                lobo::multi_all(lhs.data[i * 3 + k], rhs.data[j * 3 + k], result.data[j * 3 + i]);
                //result.data[j*row + i] += lhs.data[k*lhs.row_ + i] * rhs.data[j*rhs.row_ + k];
                //lobo::add_a(result.data[j * 3 + i], tmp, flag);
            }
        }
    }
}

inline void multiMTd_image(const LoboComplexMatrix3dd &lhs, const LoboComplexMatrix3dd &rhs, LoboComplexMatrix3dd &result)
{
    //LoboComplexMatrix<COMPLEX_TYPE, TYPE> result(row, col);
    for (int i = 0; i < 9; i++)
    {
        result.data[i].image_.image_ = 0.0;
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                lobo::multi_image(lhs.data[i * 3 + k], rhs.data[j * 3 + k], result.data[j * 3 + i]);
            }
        }
    }
}

inline void multiMTd_real_image(const LoboComplexMatrix3dd &lhs, const LoboComplexMatrix3dd &rhs, LoboComplexMatrix3dd &result)
{
    //LoboComplexMatrix<COMPLEX_TYPE, TYPE> result(row, col);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result.data[j * 3 + i].image_.real_ = 0.0;
            result.data[j * 3 + i].image_.image_ = 0.0;
        }
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                lobo::multi_real_image(lhs.data[i * 3 + k], rhs.data[j * 3 + k], result.data[j * 3 + i]);
            }
        }
    }
}

} // namespace lobo