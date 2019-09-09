#pragma once
#include "TypeNeoHookeanMaterial.h"

template <class TYPE>

class TypeNeoHookeanMaterialCSFD : public TypeNeoHookeanMaterial<TYPE>
{
    typedef Eigen::Matrix<TYPE, 3, 3> Matrix3;
    typedef Eigen::Matrix<TYPE, -1, -1> MatrixX;

    LOBO_MAKE_TYPEDEFS(TYPE, t);

public:
    using TypeTetElementMaterial<TYPE>::lambdaLame;
    using TypeTetElementMaterial<TYPE>::muLame;
    using TypeTetElementMaterial<TYPE>::compressionResistance;
    using TypeTetElementMaterial<TYPE>::EdivNuFactor;
    using TypeTetElementMaterial<TYPE>::tetmesh;
    using TypeTetElementMaterial<TYPE>::element_dF_du;
    using TypeTetElementMaterial<TYPE>::element_dF_dudu;
    using TypeTetElementMaterial<TYPE>::element_dF_dudv;

    TypeNeoHookeanMaterialCSFD(Lobo::LoboTetMesh *tetmesh, int enableCompressionResistance = 0, TYPE compressionResistance = 0.0) : TypeNeoHookeanMaterial<TYPE>(tetmesh, enableCompressionResistance, compressionResistance)
    {
    }

    virtual TYPE ComputeEnergy(int elementIndex, TYPE *invariants);
    virtual void ComputeEnergyGradient(int elementIndex, TYPE *invariants, TYPE *gradient);
    virtual void ComputeEnergyHessian(int elementIndex, TYPE *invariants, TYPE *hessian);
};
