#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "Eigen/Dense"
#include "cube.hpp"

class HBSurface {
public:
    HBSurface(int npx, int npy, int res);
    glm::vec3* get_vertices();
    glm::vec3* get_normals();
    unsigned int get_vertices_size();
    unsigned int get_n_vertices();
    unsigned int get_cp_vertices_size();
    glm::vec3* get_cp_vertices();
    unsigned int get_n_cp_vertices();

private:
    glm::vec3 eval_point(int x, int y, float u, float v);
    void init_test();
    unsigned int fxy(int x, int y);

    int npx;
    int npy;
    int res;
    glm::vec3* Pbuffer;
    glm::vec3* Nbuffer;
    int npatches;
    int ncpx;
    int ncpy;
    int k = 4;
    int l = 4;
    Eigen::MatrixXf* cpsx;
    Eigen::MatrixXf* cpsy;
    Eigen::MatrixXf* cpsz;
    Eigen::Matrix4f B;
    glm::vec3* vert;
    glm::vec3* norm;
    glm::vec3* cp_verts;
    unsigned int nvert = 0;
    bool refined = true;
};
