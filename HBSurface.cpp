#include "HBSurface.hpp"

HBSurface::HBSurface(int npx, int npy, int res):npx(npx), npy(npy), res(res){
    Pbuffer = new glm::vec3[(res+1)*(res+1)];
    Nbuffer = new glm::vec3[(res+1)*(res+1)];
    npatches = npx*npy;
    int m = npx + k - 2;
    int n = npy + l - 2;
    ncpx = m + 1;
    ncpy = n + 1;
    ncps = ncpx*ncpy;
    cp_verts = new glm::vec3[36*ncpx*ncpy];
    cpsx = new Eigen::MatrixXf(ncpx, ncpy);
    cpsy = new Eigen::MatrixXf(ncpx, ncpy);
    cpsz = new Eigen::MatrixXf(ncpx, ncpy);
    //std::cout << *cps << std::endl;

    //Initialize all the control points
    for(int i = 0; i<ncpx; i++){
        for(int j = 0; j<ncpy; j++){
            (*cpsx)(i,j) = i - (ncpx-1.0f)/2.0f;
            (*cpsy)(i,j) = 0.0f;
            (*cpsz)(i,j) = j - (ncpy-1.0f)/2.0f;
        }
    }

    (*cpsy)(3,3) = 3.0f;

    //Initialize the BSpline operator matrix
    B <<    1, 4, 1, 0,
            -3, 0, 3, 0,
            3, -6, 3, 0,
            -1, 3, -3, 1;
    B = B*(1.0f/6.0f);

    init_test();
}

void HBSurface::init_test(){
    eval_point(0,0,1,1);
}

glm::vec3 HBSurface::eval_point(int x, int y, float u, float v){
    Eigen::Matrix4f Vx = (*cpsx).block<4,4>(x, y);
    Eigen::Matrix4f Vy = (*cpsy).block<4,4>(x, y);
    Eigen::Matrix4f Vz = (*cpsz).block<4,4>(x, y);

    Eigen::Vector4f U(1, u, u*u, u*u*u);
    Eigen::Vector4f V(1, v, v*v, v*v*v);

    //std::cout << U.transpose() << std::endl;
    float Px = U.transpose()*B*Vx*B.transpose()*V;
    float Py = U.transpose()*B*Vy*B.transpose()*V;
    float Pz = U.transpose()*B*Vz*B.transpose()*V;

    glm::vec3 P(Px, Py, Pz);

    return P;
}

glm::vec3* HBSurface::get_vertices(){
    if (refined){
        refined = false;
        if (nvert != 0){
            delete [] vert;
        }

        vert = new glm::vec3[get_n_vertices()];
        norm = new glm::vec3[get_n_vertices()];
    }

    //Iterate through the patches
    float h = 0.001f;
    for (int i = 0; i < npx; i++){
        for (int j = 0; j < npy; j++){
            //Iterate in the patch
            for (int xp = 0; xp <= res; xp++){
                float u = (float) xp / (float) res;
                for (int yp = 0; yp <= res; yp++){
                    float v = (float) yp / (float) res;
                    glm::vec3 C = eval_point(i, j, u, v);
                    Pbuffer[fxy(xp, yp)] =  C;

                    // Find the normal
                    glm::vec3 A0, A1;
                    glm::vec3 B0, B1;
                    float sign = -1.0f;
                    if (xp == res && i != npx){
                        A0 = eval_point(i+1,j, h, v);
                        A1 = eval_point(i, j, u-h, v);
                        //sign *= -1.0f;
                    } else {
                        A0 = eval_point(i,j, u+h, v);
                        A1 = eval_point(i, j, u-h, v);
                    }
                    if (yp == res && j != npy){
                        B0 = eval_point(i,j+1, u, h);
                        B1 = eval_point(i, j, u, v-h);
                        //sign *= -1.0f;
                    } else {
                        B0 = eval_point(i,j, u, v+h);
                        B1 = eval_point(i, j, u, v-h);
                    }
                    Nbuffer[fxy(xp, yp)] = sign*glm::normalize(glm::cross(A0-A1, B0-B1));
                    //std::cout << Pbuffer[fxy(xp, yp)].x << std::endl;
                }
            }

            for (int xp = 0; xp < res; xp++){
                glm::vec3 P0 = Pbuffer[fxy(xp, 0)];
                glm::vec3 P1 = Pbuffer[fxy(xp+1, 0)];
                glm::vec3 N0 = Nbuffer[fxy(xp, 0)];
                glm::vec3 N1 = Nbuffer[fxy(xp+1, 0)];
                for (int yp = 0; yp < res; yp++){
                    glm::vec3 P2 = Pbuffer[fxy(xp+1, yp+1)];
                    glm::vec3 P3 = Pbuffer[fxy(xp, yp+1)];
                    glm::vec3 N2 = Nbuffer[fxy(xp+1, yp+1)];
                    glm::vec3 N3 = Nbuffer[fxy(xp, yp+1)];

                    //Find the current square element of the surface (ugly)
                    int elem = (i*npx+j)*res*res*2*3+xp*res*2*3+yp*2*3;
                    //std::cout << P0.x << " " << P0.y << " " << P0.z << std::endl;
                    vert[elem + 0] = P2;
                    vert[elem + 1] = P1;
                    vert[elem + 2] = P0;
                    vert[elem + 3] = P0;
                    vert[elem + 4] = P3;
                    vert[elem + 5] = P2;

                    norm[elem + 0] = N2;
                    norm[elem + 1] = N1;
                    norm[elem + 2] = N0;
                    norm[elem + 3] = N0;
                    norm[elem + 4] = N3;
                    norm[elem + 5] = N2;

                    //Shift over a square
                    P1 = P2;
                    P0 = P3;
                    N1 = N2;
                    N0 = N3;
                }
            }
        }
    }

    return vert;
}

glm::vec3* HBSurface::get_cp_vertices(){
    for (int i = 0; i < ncpx; i++){
        for (int j = 0; j < ncpy; j++){
            float x = (*cpsx)(i,j);
            float y = (*cpsy)(i,j);
            float z = (*cpsz)(i,j);
            //float x,y,z;
            //std::cout << x << " " << y << " " << z << std::endl;
            Cube CP(x, y, z);
            for (int k = 0; k < 36; k++){
                cp_verts[36*(i*ncpx + j) + k] = CP.get_vertices()[k];
            }
        }
    }

    return cp_verts;
}

glm::vec3 HBSurface::get_cp_col(int idx){
    if (idx == sel_idx){
        //std::cout << "what I think is selected " << idx << std::endl;
        return glm::vec3(1.0f, 0.0f, 0.0f);
    } else {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

void HBSurface::select_cp(int idx){
    if (idx < ncps){
        sel_idx = idx;
        sel_cp_i = idx/ncpx;
        sel_cp_j = idx%ncpx;
    }
}

unsigned int HBSurface::get_n_cps(){
    return ncps;
}

unsigned int HBSurface::get_cp_vertices_size(){
    return 3*get_n_cp_vertices()*sizeof(float);
}

unsigned int HBSurface::get_n_cp_vertices(){
    return 36*ncpx*ncpy;
}

glm::vec3* HBSurface::get_normals(){
    return norm;
}

unsigned int HBSurface::fxy(int x, int y){
    return x*(res+1) + y;
}

unsigned int HBSurface::get_vertices_size(){
    return 3 * get_n_vertices() * sizeof(float);
}

unsigned int HBSurface::get_n_vertices(){
    return npatches * res * res * 2 * 3;
}
