#include "TypeStVKMaterialCSFD.h"
#include <complex>
#include <omp.h>
#include "Functions/LoboMacros.h"
#include "MCSFD/MatrixOp.h"

LOBO_TEMPLATE_INSTANT(TypeStVKMaterialCSFD)
template <class TYPE>
TYPE TypeStVKMaterialCSFD<TYPE>::ComputeEnergy(int elementIndex, TYPE *invariants)
{
    TYPE IC = invariants[0];
    TYPE IIC = invariants[1];
    //double IIIC = invariants[2]; // not needed for StVK

    TYPE energy = 0.125 * lambdaLame[elementIndex] * (IC - 3.0) * (IC - 3.0) + 0.25 * muLame[elementIndex] * (IIC - 2.0 * IC + 3.0);

    this->AddCompressionResistanceEnergy(elementIndex, invariants, &energy);

    return energy;
}
template <class TYPE>
void TypeStVKMaterialCSFD<TYPE>::ComputeEnergyGradient(int elementIndex, TYPE *invariants, TYPE *gradient)
{
    this->element_dE_dI[elementIndex*3+0].real_ = invariants[0];
    this->element_dE_dI[elementIndex*3+0].image_ = this->h_CSFD;
    this->element_dE_dI[elementIndex*3+1].real_ = invariants[1];
    this->element_dE_dI[elementIndex*3+1].image_ = this->h_CSFD;

    //I1-3.0
    LoboComplext tmp = (this->element_dE_dI[elementIndex*3+0]-3.0);
    LoboComplext tmp2;
    //(I1-3.0)*(I1-3.0)
    lobo::multi_image(tmp,tmp,tmp2);
    tmp2.image_*=0.125*lambdaLame[elementIndex];
    tmp2.image_-=0.25* muLame[elementIndex] * 2.0*this->element_dE_dI[elementIndex*3+0].image_;
    gradient[0] = tmp2.image_/this->h_CSFD;

    gradient[1] = 0.25 * muLame[elementIndex] * this->element_dE_dI[elementIndex*3+1].image_/this->h_CSFD;
    gradient[2] = 0.0; // no I3

    this->AddCompressionResistanceGradient(elementIndex, invariants, gradient);
}

template <class TYPE>
void TypeStVKMaterialCSFD<TYPE>::ComputeEnergyHessian(int elementIndex, TYPE *invariants, TYPE *hessian)
{
    this->element_dE_dIdI[elementIndex*3+0].real_.real_ = invariants[0];
    this->element_dE_dIdI[elementIndex*3+0].image_.real_ = this->h_CSFD;
    this->element_dE_dIdI[elementIndex*3+0].real_.image_ = this->h_CSFD;
    this->element_dE_dIdI[elementIndex*3+0].image_.image_ = 0;
    LoboComplexDualt tmp = (this->element_dE_dIdI[elementIndex*3+0]-3.0);
    tmp = tmp*tmp*0.125 * lambdaLame[elementIndex];

    //get hessian value from dual complex
    hessian[0] = tmp.image_.image_/this->h_CSFD/this->h_CSFD;
	// 12
	hessian[1] = 0.0;
	// 13
	hessian[2] = 0.0;
	// 22
	hessian[3] = 0.0;
	// 23
	hessian[4] = 0.0;
	// 33
	hessian[5] = 0.0;

	this->AddCompressionResistanceHessian(elementIndex, invariants, hessian);
}