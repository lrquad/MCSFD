#pragma once
#include "LoboVolumtricMesh/LoboTetMesh.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Utils/pugixml/pugixml.hpp"
#include "Utils/SparseMatrix/SparseMatrixTopology.h"

template <class TYPE>
class TypeIsotropicMaterial;



namespace Lobo
{

enum ComputationFlag {
    Computeflags_energy = 1 << 0,
    Computeflags_fisrt = 1 << 1,  // generted mesh by tetgen
    Computeflags_second = 1 << 2,
    Computeflags_reset = 1 << 3,
    };


class HyperelasticModel 
{
private:
    /* data */
public:
    HyperelasticModel( LoboTetMesh *tetmesh_);
    ~HyperelasticModel();

    virtual void precompute();
    virtual void computeStiffnessMatrixTopology(Eigen::SparseMatrix<double> *K);

    virtual void runXMLscript(pugi::xml_node &xml_node);

    virtual void computeEnergySparse(Eigen::VectorXd *free_variables, double *energy, Eigen::VectorXd *jacobi, Eigen::SparseMatrix<double> *hessian, int computationflags);
    virtual void computeEnergyDense(Eigen::VectorXd *free_variables, double *energy, Eigen::VectorXd *jacobi, Eigen::MatrixXd *hessian, int computationflags){};

    std::string materialtype;
    bool isinvertible;
    double inversion_Threshold;
    bool useMCSFD;

    virtual void computeAccelerationIndices(SparseMatrixTopologyTYPE<double> *sparsetopology);



    virtual void precomputedFdU();
    virtual void initMultiThreadBuffer();

    virtual void computeElementDeformationshapeMatrix(int eleid, Eigen::Matrix3d &Ds);
    virtual void diffDeformationGradient(int elementid, std::vector<Eigen::Matrix3d> &dF);

    virtual void getTetForceMatrixCSFD(int eleid, Eigen::VectorXd *internalforce, Eigen::MatrixXd *stiffness);
    virtual void getTetForceCSFD(int eleid, Eigen::VectorXd *internalforce);

    virtual void getTetForceMatrix(int eleid, Eigen::VectorXd *internalforce, Eigen::MatrixXd *stiffness, int computationflags);


	virtual void assignTetElementForceAndMatrix(int eleid,double* energy, Eigen::VectorXd *internalforce, Eigen::SparseMatrix<double>* sparseMatrix, int flags, double weights = 1.0);


    virtual void ComputeDiagonalPFromStretches(int eleid, double* lambda, double* PDiag, double &energy);
    virtual void compute_dPdF(int eleid, Eigen::Matrix3d dPdF[9], double* fhat, Eigen::Matrix3d &U, Eigen::Matrix3d &V);
	virtual void compute_dfdF(int eleid, Eigen::Matrix3d (&dP_dF)[9], Eigen::MatrixXd &dfdF, std::vector<Eigen::Matrix3d> &dF);


    inline double gammaValue(int i, int j, double sigma[3], double invariants[3], double gradient[3], double hessian[6]);
	int tensor9x9Index(int i, int j, int m, int n);

    LoboTetMesh *tetmesh;
    TypeIsotropicMaterial<double> *elastic_material;

    Eigen::VectorXd currentdisplacement; //buffer

    //precompute
    std::vector<std::vector<Eigen::Matrix3d>> dF; // 12* numelements
    std::vector<Eigen::MatrixXd> dF_du;           // 9X12
    std::vector<Eigen::MatrixXd> dF_du_9X9;       // 9X12

    std::vector<Eigen::VectorXd *> internalforce_list;
    std::vector<Eigen::MatrixXd *> stiffness_list;
    std::vector<double> energy_list;

    //
    int colMajorMatrixToTeran[9];
	int teranToColMajorMatrix[9];

    int num_DOFs;
    
    int ** row_;
	int ** column_;
    int *diagonal_;
};
} // namespace Lobo