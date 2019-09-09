#include "LoboTetMesh.h"
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/project.h>
#include <igl/barycentric_coordinates.h>
#include <igl/in_element.h>
#include <igl/AABB.h>

#include "Functions/EigenMatrixIO.h"
#include "Functions/computeTriangle.h"
#include "Utils/glmEigenConverter.h"
#include "Utils/glmMyFunctions.h"


Lobo::LoboTetMesh::LoboTetMesh() {
    initializedGL = false;
    tetgen_command = "pq1.414";
    status_flags = 0;

    mesh_total_volume = 0.0;
    numElementVertices = 4;
}

Lobo::LoboTetMesh::~LoboTetMesh() {}

void Lobo::LoboTetMesh::reinitialTetMesh() {
    vertices_flags.resize(tet_vertice.size() / 3);
    std::fill(vertices_flags.begin(), vertices_flags.end(), 0);

    tet_vertice_col =
        Lobo::eigen_vec_2_mat(tet_vertice, tet_vertice.size() / 3, 3);
    tet_faces_col = Lobo::eigen_vec_2_mat(tet_faces, tet_faces.size() / 3, 3);
    tet_indices_col = Lobo::eigen_vec_2_mat(tet_indices, tet_indices.size() / 4, 4);


    tet_vertice_attri.resize(tet_vertice.size() / 3 * 11);
    tet_vertice_attri.setZero();
    setTetAttriColor(0.8, 0.8, 0.8);
    tet_faces_glint.resize(tet_faces.size());
    for (int i = 0; i < tet_faces.size(); i++) {
        tet_faces_glint[i] = tet_faces[i];
    }
    updateTetAttri(tet_vertice, 0, 3, 11);

    // default material
    setAllMaterial(1.0, 1000.0, 0.4);

    ori_tet_vertice = tet_vertice;

    status_flags &= ~TetMeshStatusFlags_precomputed;
    precomputeElementData();
}



void Lobo::LoboTetMesh::updateTetVertices(Eigen::VectorXd *u) {
    tet_vertice = *u + ori_tet_vertice;
    tet_vertice_col =
        Lobo::eigen_vec_2_mat(tet_vertice, tet_vertice.size() / 3, 3);
    updateTetAttri(tet_vertice, 0, 3, 11);

    // update tri mesh
}

void Lobo::LoboTetMesh::generateTet(const char *tetgen_command) {
    std::string command_ = "pq1.414Y";
    if (tetgen_command != NULL) {
        command_ = tetgen_command;
    }
    Eigen::MatrixXd TV;
    Eigen::MatrixXi TT;

    std::cout << "tetgen check" << std::endl;
    std::cout << "tri_vertices" << tri_vertices.rows() << " "
              << tri_vertices.cols() << std::endl;
    std::cout << "tri_faces" << tri_faces.rows() << " " << tri_faces.cols()
              << std::endl;

    int result = igl::copyleft::tetgen::tetrahedralize(
        tri_vertices, tri_faces, command_.c_str(), TV, TT, tet_faces_col);

    // test
    // Lobo::exportSimpleObj("test.obj",TV,TF);

    // copy data

    if (result == 0) {
        tet_vertice.resize(TV.rows() * TV.cols());
        for (int i = 0; i < TV.rows(); i++) {
            for (int j = 0; j < TV.cols(); j++) {
                tet_vertice.data()[i * TV.cols() + j] =
                    TV.data()[j * TV.rows() + i];
            }
        }
        tet_indices.resize(TT.rows() * TT.cols());
        for (int i = 0; i < TT.rows(); i++) {
            for (int j = 0; j < TT.cols(); j++) {
                // one tet has 4 vertices
                tet_indices.data()[i * TT.cols() + j] =
                    TT.data()[j * TT.rows() + i];
            }
        }

        tet_faces.resize(tet_faces_col.rows() * tet_faces_col.cols());
        for (int i = 0; i < tet_faces_col.rows(); i++) {
            for (int j = 0; j < tet_faces_col.cols(); j++) {
                tet_faces.data()[i * tet_faces_col.cols() + j] =
                    tet_faces_col.data()[j * tet_faces_col.rows() + i];
            }
        }

        status_flags |= TetMeshStatusFlags_datasizeUpdated;
        status_flags |= TetMeshStatusFlags_tetgened;
    } else {
        std::cout << "tetgen failed" << command_ << filebase << std::endl;
    }

    reinitialTetMesh();
}



void Lobo::LoboTetMesh::setInputPolygon(Eigen::VectorXd *vertices,
                                        Eigen::VectorXi *faces) {
    tri_vertices.resize(vertices->rows(), vertices->cols());
    memcpy(tri_vertices.data(), vertices->data(),
           sizeof(double) * vertices->rows() * vertices->cols());
    tri_faces.resize(faces->rows(), faces->cols());
    memcpy(tri_faces.data(), faces->data(),
           sizeof(int) * faces->rows() * faces->cols());
}

void Lobo::LoboTetMesh::exportTetMesh() {
    if (usebinary) {
        exportTetMeshBinary(filebase.c_str());
    } else {
        exportTetMeshAscii(filebase.c_str());
    }
}
void Lobo::LoboTetMesh::loadTetMesh() {
    if (usebinary) {
        loadTetMeshBinary(filebase.c_str());
    } else {
        loadTetMeshAscii(filebase.c_str());
    }
}

void Lobo::LoboTetMesh::loadTetMeshBinary(const char *filebase_) {
    std::ostringstream stringStream;
    stringStream << filebase_ << ".tet";
    std::string filename = stringStream.str();

    std::cout << "loadTetMeshBinary " << filename << std::endl;

    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in.good()) {
        std::cout << filename << "file not open" << std::endl;
        return;
    }
    EigenMatrixIO::read_binary(in, tet_vertice);
    EigenMatrixIO::read_binary(in, tet_indices);
    EigenMatrixIO::read_binary(in, tet_faces);
    in.close();
    status_flags |= TetMeshStatusFlags_loadtet;
    status_flags |= TetMeshStatusFlags_datasizeUpdated;
    reinitialTetMesh();
}

void Lobo::LoboTetMesh::exportTetMeshBinary(const char *filebase_) {
    if (!(status_flags &
          (TetMeshStatusFlags_tetgened | TetMeshStatusFlags_loadtet))) {
        return;
    }
    std::ostringstream stringStream;
    stringStream << filebase_ << ".tet";
    std::string filename = stringStream.str();

    std::cout << "exportTetMeshBinary " << filename << std::endl;

    std::ofstream out(filename,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    EigenMatrixIO::write_binary(out, tet_vertice);
    EigenMatrixIO::write_binary(out, tet_indices);
    EigenMatrixIO::write_binary(out, tet_faces);
    out.close();
}

void Lobo::LoboTetMesh::loadTetMeshAscii(const char *filebase_) {
    std::cout << "loadTetMeshAscii " << filebase_ << std::endl;

    std::ostringstream stringStream;
    stringStream << filebase_ << ".ele";
    std::string elementfile = stringStream.str();
    stringStream.str("");
    stringStream.clear();
    stringStream << filebase_ << ".node";
    std::string nodefile = stringStream.str();
    stringStream.str("");
    stringStream.clear();
    stringStream << filebase_ << ".face";
    std::string facefile = stringStream.str();

    int tmp;
    int numele, numvet, numface;
    std::ifstream inputstream(elementfile);
    inputstream >> numele >> tmp >> tmp;
    tet_indices.resize(numele * 4);
    for (int i = 0; i < numele; i++) {
        inputstream >> tmp;
        for (int j = 0; j < 4; j++) {
            inputstream >> tet_indices.data()[i * 4 + j];
        }
    }
    inputstream.close();

    inputstream.open(nodefile);
    inputstream >> numvet >> tmp >> tmp >> tmp;
    tet_vertice.resize(numvet * 3);
    for (int i = 0; i < numvet; i++) {
        inputstream >> tmp;
        for (int j = 0; j < 3; j++) {
            inputstream >> tet_vertice.data()[i * 3 + j];
        }
    }
    inputstream.close();

    inputstream.open(facefile);
    inputstream >> numface;
    tet_faces.resize(numface * 3);
    for (int i = 0; i < numface; i++) {
        inputstream >> tmp;
        for (int j = 0; j < 3; j++) {
            inputstream >> tet_faces.data()[i * 3 + j];
        }
    }
    inputstream.close();
    status_flags |= TetMeshStatusFlags_datasizeUpdated;
    status_flags |= TetMeshStatusFlags_loadtet;
    reinitialTetMesh();
}
void Lobo::LoboTetMesh::exportTetMeshAscii(const char *filebase_) {
    if (!(status_flags &
          (TetMeshStatusFlags_tetgened | TetMeshStatusFlags_loadtet))) {
        return;
    }
    std::cout << "exportTetMeshAscii " << filebase_ << std::endl;
    std::ostringstream stringStream;
    stringStream << filebase_ << ".ele";
    std::string elementfile = stringStream.str();
    stringStream.str("");
    stringStream.clear();
    stringStream << filebase_ << ".node";
    std::string nodefile = stringStream.str();
    stringStream.str("");
    stringStream.clear();
    stringStream << filebase_ << ".face";
    std::string facefile = stringStream.str();

    std::ofstream outstream(elementfile);
    outstream << tet_indices.size() / 4 << " 4 0 " << std::endl;
    for (int i = 0; i < tet_indices.rows() / 4; i++) {
        outstream << i << " ";
        for (int j = 0; j < 4; j++) {
            outstream << tet_indices.data()[i * 4 + j] << " ";
        }
        outstream << std::endl;
    }
    outstream.close();

    outstream.open(nodefile);
    outstream << tet_vertice.size() / 3 << " 3 0 0 " << std::endl;
    for (int i = 0; i < tet_vertice.rows() / 3; i++) {
        outstream << i << " ";
        for (int j = 0; j < 3; j++) {
            outstream << tet_vertice.data()[i * 3 + j] << " ";
        }
        outstream << std::endl;
    }
    outstream.close();
    outstream.open(facefile);
    outstream << tet_faces.size() / 3 << std::endl;
    for (int i = 0; i < tet_faces.rows() / 3; i++) {
        outstream << i << " ";
        for (int j = 0; j < 3; j++) {
            outstream << tet_faces.data()[i * 3 + j] << " ";
        }
        outstream << std::endl;
    }
    outstream.close();
}

void Lobo::LoboTetMesh::exportConstrainedVertices(const char *filename) {
    std::vector<unsigned int> constrained_DoFs;
    for (int i = 0; i < vertices_flags.size(); i++) {
        if (vertices_flags[i] == 1) {
            constrained_DoFs.push_back(i * 3);
            constrained_DoFs.push_back(i * 3 + 1);
            constrained_DoFs.push_back(i * 3 + 2);
        }
    }
    std::ofstream output(filename);
    output << constrained_DoFs.size() << std::endl;
    for (int i = 0; i < constrained_DoFs.size(); i++) {
        output << constrained_DoFs[i] << std::endl;
    }
    output.close();
}

void Lobo::LoboTetMesh::updateTetAttri(Eigen::VectorXd &inputattri, int offset,
                                       int attrisize, int totalsize) {
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < inputattri.size() / attrisize; i++) {
            for (int j = 0; j < attrisize; j++)
                tet_vertice_attri.data()[i * totalsize + offset + j] =
                    inputattri.data()[i * attrisize + j];
        }
    }
}

void Lobo::LoboTetMesh::setTetAttriColor(double r, double g, double b,
                                         int offset, int totalsize) {
    for (int i = 0; i < tet_vertice_attri.size() / totalsize; i++) {
        tet_vertice_attri.data()[i * totalsize + offset + 0] = r;
        tet_vertice_attri.data()[i * totalsize + offset + 1] = g;
        tet_vertice_attri.data()[i * totalsize + offset + 2] = b;
    }
}

void Lobo::LoboTetMesh::setTetVetAttriColor(int vid, double r, double g,
                                            double b, int offset,
                                            int totalsize) {
    tet_vertice_attri.data()[vid * totalsize + offset + 0] = r;
    tet_vertice_attri.data()[vid * totalsize + offset + 1] = g;
    tet_vertice_attri.data()[vid * totalsize + offset + 2] = b;
}

void Lobo::LoboTetMesh::setAllMaterial(double density, double youngsmodulus,
                                       double possionratio) {
    materials.clear();
    materialid.resize(getNumElements());
    std::fill(materialid.begin(), materialid.end(), 0);

    Material m;
    m.density = density;
    m.youngsmodulus = youngsmodulus;
    m.possionratio = possionratio;
    materials.push_back(m);
}

void Lobo::LoboTetMesh::computeDiagMassMatrix(
    Eigen::SparseMatrix<double> *mass) {
    int r = getNumVertex() * 3;
    mass->resize(r, r);
    typedef Eigen::Triplet<double> EIGEN_TRI_T;
    std::vector<EIGEN_TRI_T> coefficients;
    int numElements = getNumElements();
    for (int i = 0; i < numElements; i++) {
        double vol = elements_data[i].volume;
        double density = materials[materialid[i]].density;
        for (int j = 0; j < 4; j++) {
            double mass_value = 0;
            mass_value = density * vol / 4.0;
            int row = tet_indices[i * 4 + j] * 3;
            for (int d = 0; d < 3; d++) {
                coefficients.push_back(
                    EIGEN_TRI_T(row + d, row + d, mass_value));
            }
        }
    }
    mass->setFromTriplets(coefficients.begin(), coefficients.end());
}

void Lobo::LoboTetMesh::precomputeElementData() {
    if (status_flags & TetMeshStatusFlags_precomputed) {
        // already precomputed
        return;
    }
    int numElements = getNumElements();
    elements_data.resize(numElements);

    mesh_total_volume = 0.0;
    for (int i = 0; i < numElements; i++) {
        correctElementNodeOrder(i);
        computeElementVolume(i);
        computeElementShapeFunctionDerivate(i);  // precompute dF_du
        mesh_total_volume += elements_data[i].volume;

        // update element center
        elements_data[i].center_p.setZero();
        for (int j = 0; j < 4; j++) {
            int nodeid = tet_indices.data()[i * 4 + j];
            elements_data[i].center_p.data()[0] +=
                tet_vertice.data()[nodeid * 3 + 0];
            elements_data[i].center_p.data()[1] +=
                tet_vertice.data()[nodeid * 3 + 1];
            elements_data[i].center_p.data()[2] +=
                tet_vertice.data()[nodeid * 3 + 2];
        }
        elements_data[i].center_p /= 4;
    }

    generateBarycentricCoordinate();

    status_flags |= TetMeshStatusFlags_precomputed;
}

void Lobo::LoboTetMesh::getNodeRestPosition(int nodeid, Eigen::Vector3d &p) {
    for (int j = 0; j < 3; j++) {
        p.data()[j] = ori_tet_vertice.data()[nodeid * 3 + j];
    }
}

Eigen::Vector3d Lobo::LoboTetMesh::getNodeRestPosition(int nodeid) {
    Eigen::Vector3d tmp;
    getNodeRestPosition(nodeid, tmp);
    return tmp;
}

void Lobo::LoboTetMesh::correctElementNodeOrder(int elementid) {
    int ni[4];
    Eigen::Vector3d node_p[4];
    for (int i = 0; i < 4; i++) {
        ni[i] = tet_indices[elementid * 4 + i];
        for (int j = 0; j < 3; j++) {
            getNodeRestPosition(ni[i], node_p[i]);
        }
    }
    Eigen::Vector3d direction_v;
    Lobo::computeTriangleNorm(node_p[0], node_p[1], node_p[2], direction_v);
    Eigen::Vector3d n3n0 = node_p[3] - node_p[0];
    if (n3n0.dot(direction_v) < 0) {
        // element->node_indices[1] = n2;
        // element->node_indices[2] = n1;
        tet_indices[elementid * 4 + 1] = ni[2];
        tet_indices[elementid * 4 + 2] = ni[1];
        std::cout << "bad order element" << std::endl;
    }
}

void Lobo::LoboTetMesh::computeElementVolume(int elementid) {
    Eigen::Vector3d a =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 0]);
    Eigen::Vector3d b =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 1]);
    Eigen::Vector3d c =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 2]);
    Eigen::Vector3d d =
        this->getNodeRestPosition(tet_indices[elementid * 4 + 3]);

    elements_data[elementid].volume = Lobo::computeTetVolumeABS(a, b, c, d);
}

void Lobo::LoboTetMesh::computeElementShapeFunctionDerivate(int elementid) {
    int ni[4];
    Eigen::Vector3d node_p[4];
    for (int i = 0; i < 4; i++) {
        ni[i] = tet_indices[elementid * 4 + i];
        for (int j = 0; j < 3; j++) {
            getNodeRestPosition(ni[i], node_p[i]);
        }
    }

    TetElementData *te = &elements_data[elementid];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            te->Dm.data()[j * 3 + i] =
                node_p[j].data()[i] - node_p[3].data()[i];
        }
    }
    te->Dm_inverse = te->Dm.inverse();

    Eigen::Matrix4d referenceShape;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            referenceShape.data()[i * 4 + j] = node_p[i].data()[j];
        }
        referenceShape.data()[i * 4 + 3] = 1;
    }

    // need test
    Eigen::Matrix4d inverseShapefunction = referenceShape.inverse();
    te->shape_function_inv = inverseShapefunction;
}

void Lobo::LoboTetMesh::generateBarycentricCoordinate() {
    
}

void Lobo::LoboTetMesh::computeBarycentricWeights(int eleid,
                                                  Eigen::Vector3d &pos,
                                                  Eigen::Vector4d &weights) {
    Eigen::Vector4d point;
    point.data()[0] = pos.data()[0];
    point.data()[1] = pos.data()[1];
    point.data()[2] = pos.data()[2];
    point.data()[3] = 1;

    if (eleid == -1) {
        weights.setZero();
        return;
    }

    weights = elements_data[eleid].shape_function_inv * point;
}

int Lobo::LoboTetMesh::getContainedElement(Eigen::Vector3d &position) {
    int numElements = this->getNumElements();
    for (int element = 0; element < numElements; element++) {
        if (containsVertex(element, position)) return element;
    }
    return -1;
}

int Lobo::LoboTetMesh::getCloesetElement(Eigen::Vector3d &position) {
    double distance = DBL_MAX;
    int eleid = -1;
    int numElements = this->getNumElements();
    for (int element = 0; element < numElements; element++) {
       double dis = (position-elements_data[element].center_p).norm();
       if(dis<distance)
       {
           eleid = element;
           distance = dis;
       }
    }
    return eleid;
}

bool Lobo::LoboTetMesh::containsVertex(int eleid, Eigen::Vector3d &pos) {
    Eigen::Vector4d weight;
    this->computeBarycentricWeights(eleid, pos, weight);

    return ((weight.data()[0] >= -1e-2) && (weight.data()[1] >= -1e-2) &&
            (weight.data()[2] >= -1e-2) && (weight.data()[3] >= -1e-2));
}
