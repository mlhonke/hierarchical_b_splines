#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "cs488-framework/CS488Window.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

#include <thread>
#include <mutex>
#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Dense"
#include "cube.hpp"

class A1;

typedef Eigen::Matrix<float, 5, 5> Matrix5f;
typedef Eigen::Matrix<float, 7, 7> Matrix7f;

class HBSurface {
public:
    HBSurface(A1* GLapp, ShaderProgram* m_shader, ShaderProgram* b_shader, int npx, int npy, int res);
    void render_points(glm::mat4 W, glm::mat4 proj, glm::mat4 view, bool do_picking, int level);
    void render_surface(glm::mat4 W, glm::mat4 proj, glm::mat4 view, bool do_picking);
    glm::vec3* get_vertices();
    glm::vec3* get_normals();
    unsigned int get_vertices_size();
    unsigned int get_n_vertices();
    unsigned int get_cp_vertices_size();
    glm::vec3* get_cp_vertices();
    unsigned int get_n_cp_vertices();
    unsigned int get_n_cps();
    void select_cp(int idx, bool second, int level);
    glm::vec3 get_cp_col(int idx);
    bool is_cp_selected();
    void move_selected_cp(glm::vec3 delta);
    void split_selected_cp();
    unsigned int get_selected_cp_idx();
    void update_cps();
    glm::vec3 get_selected_cp_coords();
    void update_references();
    void set_mode_edit_cp(bool val);
    void get_colour(int level, std::vector<float *> &col_list);

    // Public variables
    int ncpx;
    int ncpy;
    int npx;
    int npy;
    Eigen::MatrixXf* Ocpsx;
    Eigen::MatrixXf* Ocpsy;
    Eigen::MatrixXf* Ocpsz;
    Eigen::MatrixXi* cpmask;
    Eigen::MatrixXf* Rcpsx;
    Eigen::MatrixXf* Rcpsy;
    Eigen::MatrixXf* Rcpsz;
    Eigen::MatrixXf* cpsx;
    Eigen::MatrixXf* cpsy;
    Eigen::MatrixXf* cpsz;
    GLFWwindow * m_window;
    static int idx_start;
    int parent_sel_cp_i_2, parent_sel_cp_j_2, parent_sel_cp_j, parent_sel_cp_i;

private:
    void init_render();
    A1* GLapp;
    void split_patch(int i, int j, Matrix5f& X, Matrix5f& Y, Matrix5f& Z);
    glm::vec3 eval_point(int x, int y, float u, float v);
    glm::vec3 eval_normal(int x, int y, float u, float v);
    void init_test();
    unsigned int fxy(int x, int y);
    void split();
    bool is_cp_mine(int idx);
    void get_split_reference(Eigen::MatrixXf& RX, Eigen::MatrixXf& RY, Eigen::MatrixXf& RZ, int cp_i, int cp_j, int cp_i_2, int cp_j_2);
    glm::vec3 get_edit_point(int i, int j);
    glm::vec3* get_edit_vertices();
    unsigned int get_n_eps();
    unsigned int get_n_ep_vertices();
    unsigned int get_ep_vertices_size();
    glm::vec3 get_ep_col(int idx);
    void gen_vertex_patch(int& elem, int i, int j);

    float* colour;
    std::mutex mtx;
    bool edit_cp_mode = false;
    int my_idx_start = idx_start;
    bool has_parent = false;
    bool has_children = false;
    std::vector<HBSurface*> child_list;
    Eigen::MatrixXi* render_patch;
    HBSurface* parent;
    bool selected = false;
    int res;
    //glm::vec3* Pbuffer;
    //glm::vec3* Nbuffer;
    std::thread* threads;
    int npatches;
    int sel_cp_i;
    int sel_cp_j;
    int sel_cp_i_2;
    int sel_cp_j_2;
    int sel_idx = -1;
    int sel_idx_2 = -1;
    int k = 4;
    int l = 4;
    Eigen::Matrix4f B;
    Eigen::Matrix4f alphaA1;
    Eigen::Matrix4f alphaA2;
    glm::vec3* vert;
    glm::vec3* norm;
    glm::vec3* cp_verts;
    glm::vec3* ep_verts;
    unsigned int nvert = 0;
    unsigned int n_cp_verts = 0;
    unsigned int ncps = 0;
    bool refined = true;
    bool refinedcps = true;
    bool refinedeps = true;

    //OpenGL parameters
    ShaderProgram* m_shader;
    ShaderProgram* b_shader;
    GLuint vaos[2];
    GLuint vbos[3];
    GLuint m_surface_vao; // Vertex Array Object
    GLuint m_surface_vbo; // Vertex Buffer Object
    GLuint m_surface_normals_vbo;
    GLuint m_cp_vao;
    GLuint m_cp_vbo;
    GLint m_normalAttribLocation;
    GLint kd_location;
    GLint ks_location;
    GLint sh_location;
    GLint norm_location;
    GLint Pos;

    // Set up the uniforms
    GLint Pers;
    GLint Model;

    //Set up basic shader uniforms (non-phong)
    GLint P_uni;
    GLint V_uni;
    GLint M_uni;
    GLint col_uni;
    GLint posAttrib2;
};
