// this example shows how to compute elastic model energy, energy gradient and energy hessian

#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>

#include "MCSFD/MCSFDCore.h"
#include "LoboVolumtricMesh/LoboTetMesh.h"
#include "ElasticModel/HyperelasticModel.h"

#include "Utils/pugixml/pugixml.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Utils/SparseMatrix/SparseMatrixTopology.h"
#include <omp.h>

LOBO_MAKE_TYPEDEFS(double, t);

using namespace Lobo;

void loadmaterial(Lobo::HyperelasticModel *elasticmodel, const char *xmlfile)
{
    pugi::xml_document xml_doc;
    pugi::xml_parse_result result = xml_doc.load_file(xmlfile);
    pugi::xml_node model_node = xml_doc.child("Scene").child("HyperelasticModel");
    elasticmodel->runXMLscript(model_node);
    elasticmodel->precompute();
    elasticmodel->useMCSFD = false;
    elasticmodel->isinvertible = false;
}

int main()
{
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(
        1); // Use 4 threads for all consecutive parallel regions

    Lobo::LoboTetMesh *tetmesh = new Lobo::LoboTetMesh();
    tetmesh->loadTetMeshAscii("2tet");
    std::cout << tetmesh->tet_indices.size() << std::endl;
    Lobo::HyperelasticModel *stvkmodel = new Lobo::HyperelasticModel(tetmesh);
    loadmaterial(stvkmodel, "FEM_stvk.xml");
    Lobo::HyperelasticModel *stvkmodel_csfd = new Lobo::HyperelasticModel(tetmesh);
    loadmaterial(stvkmodel_csfd, "FEM_stvkCSFD.xml");

    Lobo::HyperelasticModel *neohookeanmodel = new Lobo::HyperelasticModel(tetmesh);
    loadmaterial(neohookeanmodel, "FEM_neohookean.xml");
    Lobo::HyperelasticModel *neohookeanmodel_csfd = new Lobo::HyperelasticModel(tetmesh);
    loadmaterial(neohookeanmodel_csfd, "FEM_neohookeanCSFD.xml");

    int num_dofs = tetmesh->getNumVertex() * 3;
    Eigen::VectorXd u(num_dofs);
    Eigen::VectorXd internalforce(num_dofs);

    int flags_all = 0;
    flags_all |= Computeflags_energy | Computeflags_fisrt | Computeflags_second;

    srand(time(NULL));
    u.setRandom();

    stvkmodel->currentdisplacement = u;
    stvkmodel_csfd->currentdisplacement = u;
    neohookeanmodel->currentdisplacement = u;
    neohookeanmodel_csfd->currentdisplacement = u;

    int numtest = 10000;
    double energy;
    clock_t t1 = clock();

    //compute energy, energy gradient, energy hessian with analytical function
    for (int i = 0; i < numtest; i++)
    {
        stvkmodel->getTetForceMatrix(0, stvkmodel->internalforce_list[0], stvkmodel->stiffness_list[0], flags_all);
    }

    clock_t t2 = clock();
    std::cout << energy << std::endl;
    //std::cout << *stvkmodel->stiffness_list[0] << std::endl;
    std::cout << "no MSCFD stvk " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;

    t1 = clock();

    //compute energy, energy gradient, energy hessian with MSCFD
    for (int i = 0; i < numtest; i++)
    {
        stvkmodel_csfd->getTetForceMatrix(0, stvkmodel_csfd->internalforce_list[0], stvkmodel_csfd->stiffness_list[0], flags_all);
    }

    t2 = clock();

    std::cout << energy << std::endl;
    //std::cout << *stvkmodel_csfd->stiffness_list[0] << std::endl;
    std::cout << "MSCFD stvk " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;

    t1 = clock();

    //compute energy, energy gradient, energy hessian with MSCFD
    for (int i = 0; i < numtest; i++)
    {
        stvkmodel_csfd->getTetForceMatrixCSFD(0, stvkmodel_csfd->internalforce_list[0], stvkmodel_csfd->stiffness_list[0]);
    }

    t2 = clock();

    std::cout << energy << std::endl;
    //std::cout << *stvkmodel_csfd->stiffness_list[0] << std::endl;
    std::cout << "no image MSCFD stvk " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;


    t1 = clock();

    //compute energy, energy gradient, energy hessian with MSCFD
    for (int i = 0; i < numtest; i++)
    {
        neohookeanmodel->getTetForceMatrix(0, neohookeanmodel->internalforce_list[0], neohookeanmodel->stiffness_list[0], flags_all);
    }

    t2 = clock();

    std::cout << energy << std::endl;
    //std::cout << *neohookeanmodel->stiffness_list[0] << std::endl;
    std::cout << "no MSCFD neohookean " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;


    t1 = clock();

    //compute energy, energy gradient, energy hessian with MSCFD
    for (int i = 0; i < numtest; i++)
    {
        neohookeanmodel_csfd->getTetForceMatrix(0, neohookeanmodel_csfd->internalforce_list[0], neohookeanmodel_csfd->stiffness_list[0], flags_all);
    }

    t2 = clock();

    std::cout << energy << std::endl;
    //std::cout << *neohookeanmodel_csfd->stiffness_list[0] << std::endl;
    std::cout << "MSCFD neohookean " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << std::endl;


    //delete tetmesh;
    //delete hyperelasticmodel_csfd;
    //delete hyperelasticmodel;

    return 0;
}
