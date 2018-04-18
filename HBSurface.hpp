#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "Eigen/Dense"

class HBSurface {
public:
    HBSurface(int npx, int npy, int res);
    glm::vec3* get_vertices();
    unsigned int get_vertices_size();
    unsigned int get_n_vertices();

private:
    glm::vec3 eval_point(int x, int y, float u, float v);
    void init_test();
    unsigned int fxy(int x, int y);

    int npx;
    int npy;
    int res;
    glm::vec3* Pbuffer;
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
    unsigned int nvert = 0;
    bool refined = true;
};
