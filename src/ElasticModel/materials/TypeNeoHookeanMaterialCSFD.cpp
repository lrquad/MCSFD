#include "TypeNeoHookeanMaterialCSFD.h"
#include <complex>
#include <omp.h>
#include "Functions/LoboMacros.h"
#include "MCSFD/MatrixOp.h"

LOBO_TEMPLATE_INSTANT(TypeNeoHookeanMaterialCSFD)

template <class TYPE>
TYPE TypeNeoHookeanMaterialCSFD<TYPE>::ComputeEnergy(int elementIndex, TYPE *invariants)
{
    TYPE IC = invariants[0];
    TYPE IIIC = invariants[2];
    TYPE J = std::sqrt(IIIC);
    TYPE logJ = std::log(J);

    TYPE energy = 0.5 * muLame[elementIndex] * (IC - 3.0) - muLame[elementIndex] * logJ + 0.5 * lambdaLame[elementIndex] * logJ * logJ;

    this->AddCompressionResistanceEnergy(elementIndex, invariants, &energy);
    return energy;
}

template <class TYPE>
void TypeNeoHookeanMaterialCSFD<TYPE>::ComputeEnergyGradient(int elementIndex, TYPE *invariants, TYPE *gradient)
{
    this->element_dE_dI[elementIndex * 3 + 0].real_ = invariants[0];
    this->element_dE_dI[elementIndex * 3 + 0].image_ = this->h_CSFD;
    this->element_dE_dI[elementIndex * 3 + 2].real_ = invariants[2];
    this->element_dE_dI[elementIndex * 3 + 2].image_ = this->h_CSFD;

    LoboComplext J = lobo::sqrt(this->element_dE_dI[elementIndex * 3 + 2]);
    LoboComplext logJ = lobo::log(J);
    LoboComplext tmp = -muLame[elementIndex] * logJ + 0.5 * lambdaLame[elementIndex] * logJ * logJ;

    gradient[0] = 0.5 * muLame[elementIndex];
    gradient[1] = 0.0;
    gradient[2] = tmp.image_ / this->h_CSFD;

    this->AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

template <class TYPE>
void TypeNeoHookeanMaterialCSFD<TYPE>::ComputeEnergyHessian(int elementIndex, TYPE *invariants, TYPE *hessian)
{
    this->element_dE_dIdI[elementIndex*3+2].real_.real_ = invariants[2];
    this->element_dE_dIdI[elementIndex*3+2].image_.real_ = this->h_CSFD;
    this->element_dE_dIdI[elementIndex*3+2].real_.image_ = this->h_CSFD;
    this->element_dE_dIdI[elementIndex*3+2].image_.image_ = 0;

    LoboComplexDualt J = lobo::sqrt(this->element_dE_dIdI[elementIndex * 3 + 2]);
    LoboComplexDualt logJ = lobo::log(J);
    LoboComplexDualt tmp = -muLame[elementIndex] * logJ + 0.5 * lambdaLame[elementIndex] * logJ * logJ;

    // 11
    hessian[0] = 0.0;
    // 12
    hessian[1] = 0.0;
    // 13
    hessian[2] = 0.0;
    // 22
    hessian[3] = 0.0;
    // 23
    hessian[4] = 0.0;
    // 33

    //get hessian value from dual complex
    hessian[5] = tmp.image_.image_/this->h_CSFD/this->h_CSFD;

    this->AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}