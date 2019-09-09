#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <string>
#include "ObjLoader/tiny_obj_loader.h"

enum TetMeshStatusFlags_
{
    TetMeshStatusFlags_None = 0,
    TetMeshStatusFlags_tetgened = 1 << 0, // generted mesh by tetgen
    TetMeshStatusFlags_loadtet = 1 << 1,
    TetMeshStatusFlags_initialGL = 1 << 2,
    TetMeshStatusFlags_datasizeUpdated = 1 << 3,
    TetMeshStatusFlags_precomputed = 1 << 4
};

namespace Lobo
{

struct Material
{
    /* data */
    double density;
    double youngsmodulus;
    double possionratio;

    double getLambda()
    {
        return (possionratio * youngsmodulus) / ((1 + possionratio) * (1 - 2 * possionratio));
    }

    double getMu()
    {
        return youngsmodulus / (2 * (1 + possionratio));
    }

    double getE(){return youngsmodulus;}
    double getNu(){return possionratio;}
};

struct TetElementData
{
    Eigen::Matrix3d Dm;
    Eigen::Matrix3d Dm_inverse;
    Eigen::Matrix4d shape_function_inv;
    Eigen::Vector3d center_p;

    double volume;
    /* data */
};

class LoboTetMesh
{
public:
    LoboTetMesh();
    ~LoboTetMesh();

    virtual void reinitialTetMesh();

    virtual void updateTetVertices(Eigen::VectorXd *u);

    virtual void setInputPolygon(Eigen::VectorXd *vertices,
                                 Eigen::VectorXi *Faces);

    virtual void generateTet(const char *tetgen_command = NULL);

    virtual void exportTetMesh();
    virtual void loadTetMesh();

    virtual void loadTetMeshBinary(const char *filename);
    virtual void loadTetMeshAscii(const char *filebase);

    virtual void exportTetMeshBinary(const char *filename);
    virtual void exportTetMeshAscii(const char *filebase);

    virtual void exportConstrainedVertices(const char *filename);

    //for simulation
    virtual void setAllMaterial(double density, double youngsmodulus, double possionratio);

    virtual void computeDiagMassMatrix(Eigen::SparseMatrix<double>* mass);

    virtual void precomputeElementData();


    //adapter for trimesh
    virtual void generateBarycentricCoordinate();
    virtual void computeBarycentricWeights(int eleid,Eigen::Vector3d&position,Eigen::Vector4d& weights);
    virtual int getContainedElement(Eigen::Vector3d &position);
    virtual int getCloesetElement(Eigen::Vector3d &position);

    virtual bool containsVertex(int eleid, Eigen::Vector3d &pos);

    // generate tetgen
    std::string filebase;
    std::string tetgen_command;
    bool usebinary;
    bool initializedGL;
    int status_flags;

    // Input polygon
    Eigen::MatrixXd tri_vertices;
    Eigen::MatrixXi tri_faces;

    // Tetrahedralized interior
    int numElementVertices;
    Eigen::VectorXd tet_vertice;       //#numVertices*3 X 1
    Eigen::VectorXd ori_tet_vertice;       //#numVertices*3 X 1

    Eigen::VectorXf tet_vertice_attri; //#numVertices*3 X1
    Eigen::VectorXi tet_indices;       //#numTet*4 X 1
    Eigen::VectorXi tet_faces;         //#numface*3 // for rendering X 1

    Eigen::MatrixXd tet_vertice_col;
    Eigen::MatrixXi tet_faces_col;
    Eigen::MatrixXi tet_indices_col;


    //for test
    std::vector<unsigned int> tet_faces_glint;
    std::vector<unsigned int> vertices_flags;


    unsigned int VAO; // vertex array object
    unsigned int VBO; // vertex buffer object
    unsigned int EBO; // element buffer object

    //interface
    int getNumElements() { return tet_indices.size() / 4; }
    int getNumVertex(){return tet_vertice.size()/3;}
    Material *getElementMaterial(int elementid) { return &materials[materialid[elementid]]; }
    TetElementData* getTetElement(int elementid){return &elements_data[elementid];};

    //Eigen interface
    void getNodeRestPosition(int nodeid,Eigen::Vector3d &p);
    Eigen::Vector3d getNodeRestPosition(int nodeid);


    //for triangle mesh binding
    Eigen::VectorXi tri_ele_idl; //tri vertex to tet ele
    Eigen::VectorXi tri_vertices_idl; //tri vertex to tet node
    Eigen::VectorXd tri_ele_weights;

protected:

    virtual void updateTetAttri(Eigen::VectorXd &inputattri, int offset, int attrisize, int totalsize);
    virtual void setTetAttriColor(double r, double g, double b, int offset = 8, int totalsize = 11);
    virtual void setTetVetAttriColor(int vid, double r, double g, double b, int offset = 8, int totalsize = 11);


    //tetmesh precompute function
    virtual void correctElementNodeOrder(int elementid);
    virtual void computeElementVolume(int elementid);
    virtual void computeElementShapeFunctionDerivate(int elementid);

    tinyobj::material_t default_material;

    std::vector<Material> materials;
    std::vector<int> materialid;
    std::vector<TetElementData> elements_data;

    

    double mesh_total_volume;

};
} // namespace Lobo