#include "HBSurface.hpp"
#include "A1.hpp"

int HBSurface::idx_start = 0;

HBSurface::HBSurface(A1* GLapp, ShaderProgram* m_shader, ShaderProgram* b_shader, int npx, int npy, int res)
    : npx(npx), npy(npy), res(res), GLapp(GLapp), m_shader(m_shader), b_shader(b_shader){
    children = new HBSurface**[npx];
    for (int i = 0; i<npx; i++){
        children[i] = new HBSurface*[npy];
        for (int j=0; j<npy; j++){
            children[i][j] = this;
        }
    }

    npatches = npx*npy;
    int m = npx + k - 2;
    int n = npy + l - 2;
    ncpx = m + 1;
    ncpy = n + 1;
    ncps = ncpx*ncpy;
    cpsx = new Eigen::MatrixXf(ncpx, ncpy);
    cpsy = new Eigen::MatrixXf(ncpx, ncpy);
    cpsz = new Eigen::MatrixXf(ncpx, ncpy);
    Rcpsx = new Eigen::MatrixXf(ncpx, ncpy);
    Rcpsy = new Eigen::MatrixXf(ncpx, ncpy);
    Rcpsz = new Eigen::MatrixXf(ncpx, ncpy);
    (*Rcpsx) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    (*Rcpsy) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    (*Rcpsz) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    Ocpsx = new Eigen::MatrixXf(ncpx, ncpy);
    Ocpsy = new Eigen::MatrixXf(ncpx, ncpy);
    Ocpsz = new Eigen::MatrixXf(ncpx, ncpy);
    (*Ocpsx) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    (*Ocpsy) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    (*Ocpsz) = Eigen::MatrixXf::Zero(ncpx, ncpy);
    cpmask = new Eigen::MatrixXi(ncpx, ncpy);
    *cpmask = Eigen::MatrixXi::Constant(ncpx, ncpy, 1);

    threads = new std::thread[npx*npy];

    //std::cout << my_idx_start << std::endl;
    //std::cout << *cps << std::endl;

    //Initialize all the control points only if first surface.
    if (idx_start == 0){
        for(int i = 0; i<ncpx; i++){
            for(int j = 0; j<ncpy; j++){
                (*Ocpsx)(i,j) = i;
                (*Ocpsy)(i,j) = 0.0f;
                (*Ocpsz)(i,j) = j;
            }
        }
    }

    idx_start += ncps;

    update_cps();

    //(*cpsy)(3,3) = 3.0f;

    //Initialize the BSpline operator matrix
    B <<    1, 4, 1, 0,
            -3, 0, 3, 0,
            3, -6, 3, 0,
            -1, 3, -3, 1;
    B = B*(1.0f/6.0f);

    alphaA1 <<   1.0f/2.0f, 1.0f/2.0f, 0, 0,
            1.0f/8.0f, 3.0f/4.0f, 1.0f/8.0f, 0,
            0, 1.0f/2.0f, 1.0f/2.0f, 0,
            0, 1.0/8.0f, 3.0f/4.0f, 1.0f/8.0f;

    alphaA2 <<   1.0f/8.0f, 3.0f/4.0f, 1.0f/8.0f, 0,
            0, 1.0f/2.0f, 1.0f/2.0f, 0,
            0, 1.0/8.0f, 3.0f/4.0f, 1.0f/8.0f,
            0, 0, 1.0f/2.0f, 1.0f/2.0f;
    init_test();
    init_render();
}

void HBSurface::update_cps(){
    *cpsx = *Ocpsx + *Rcpsx;
    *cpsy = *Ocpsy + *Rcpsy;
    *cpsz = *Ocpsz + *Rcpsz;
}

void HBSurface::init_render(){
    GLapp->gen_buffers(vaos, vbos);
    m_surface_vao = vaos[0];
	m_cp_vao = vaos[1];

	m_surface_vbo = vbos[0];
	m_surface_normals_vbo = vbos[1];
	m_cp_vbo = vbos[2];

    m_normalAttribLocation = m_shader->getAttribLocation("normal");
    kd_location = m_shader->getUniformLocation("material.kd");
    ks_location = m_shader->getUniformLocation("material.ks");
    sh_location = m_shader->getUniformLocation("material.shininess");
    norm_location = m_shader->getUniformLocation("NormalMatrix");
    Pos = m_shader->getAttribLocation( "position" );

    // Set up the uniforms
    Pers = m_shader->getUniformLocation( "Perspective" );
    Model = m_shader->getUniformLocation( "ModelView" );

    //Set up basic shader uniforms (non-phong)
    P_uni = b_shader->getUniformLocation( "P" );
    V_uni = b_shader->getUniformLocation( "V" );
    M_uni = b_shader->getUniformLocation( "M" );
    col_uni = b_shader->getUniformLocation( "colour" );
    posAttrib2 = b_shader->getAttribLocation( "position");

    CHECK_GL_ERRORS;
}

void HBSurface::render_points(glm::mat4 W, glm::mat4 proj, glm::mat4 view, bool do_picking, int level){
    if (level <= 0){
        int n_points;
        //CP vertices
        glBindVertexArray(m_cp_vao);
        glBindBuffer(GL_ARRAY_BUFFER, m_cp_vbo);
        if (edit_cp_mode && level == 0){
            glBufferData(GL_ARRAY_BUFFER, get_ep_vertices_size(), get_edit_vertices(), GL_DYNAMIC_DRAW);
            n_points = get_n_eps();
        } else {
            glBufferData(GL_ARRAY_BUFFER, get_cp_vertices_size(), get_cp_vertices(), GL_DYNAMIC_DRAW);
            n_points = get_n_cps();
        }
        glEnableVertexAttribArray( posAttrib2 );
        glVertexAttribPointer( posAttrib2, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

        glBindVertexArray( 0 );
        glBindBuffer( GL_ARRAY_BUFFER, 0 );
        glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

        CHECK_GL_ERRORS;

        b_shader->enable();
            glEnable( GL_DEPTH_TEST );
            glEnable( GL_CULL_FACE );
            glUniformMatrix4fv( P_uni, 1, GL_FALSE, value_ptr( proj ) );
            glUniformMatrix4fv( V_uni, 1, GL_FALSE, value_ptr( view ) );
            glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W ) );
            glBindVertexArray( m_cp_vao );

            if (do_picking){
                for (int idx = my_idx_start; idx < (n_points + my_idx_start); idx++){
                    float r = float(idx&0xff) / 255.0f;
                    float g = float((idx>>8)&0xff) / 255.0f;
                    float b = float((idx>>16)&0xff) / 255.0f;
                    glUniform3f( col_uni, r, g, b);
                    glDrawArrays(GL_TRIANGLES, (idx-my_idx_start)*36, 36);
                }
            } else {
                for (int idx = my_idx_start; idx < (n_points + my_idx_start); idx++){
                    glm::vec3 cp_col;
                    if (edit_cp_mode && level == 0){
                        cp_col = get_ep_col(idx);
                    } else {
                        cp_col = get_cp_col(idx);
                    }
                    glUniform3f( col_uni, cp_col.x, cp_col.y, cp_col.z );
                    glDrawArrays( GL_TRIANGLES, (idx-my_idx_start)*36, 36);
                }
            }
            CHECK_GL_ERRORS;
        b_shader->disable();
    }

    if (level != 0){
        level--;
        for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
            (*child)->render_points(W, proj, view, do_picking, level);
        }
    }
}

void HBSurface::render_surface(glm::mat4 W, glm::mat4 proj, glm::mat4 view, bool do_picking){
    //Render surface vertices
    glBindVertexArray(m_surface_vao);

    //Position data
    glBindBuffer(GL_ARRAY_BUFFER, m_surface_vbo);
    glBufferData(GL_ARRAY_BUFFER, get_vertices_size(), get_vertices(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray( Pos );
    glVertexAttribPointer( Pos, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

    //Normal data
    glBindBuffer(GL_ARRAY_BUFFER, m_surface_normals_vbo);
    glBufferData(GL_ARRAY_BUFFER, get_vertices_size(), get_normals(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(m_normalAttribLocation);
    glVertexAttribPointer(m_normalAttribLocation, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    glBindVertexArray( 0 );
    glBindBuffer( GL_ARRAY_BUFFER, 0 );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

    CHECK_GL_ERRORS;

    m_shader->enable();
        glEnable( GL_DEPTH_TEST );
        glDisable( GL_CULL_FACE );

        glUniformMatrix4fv( Pers, 1, GL_FALSE, value_ptr( proj ) );
        mat4 M = view * W;
        glUniformMatrix4fv( Model, 1, GL_FALSE, value_ptr( M ) );
        mat3 normalMatrix = glm::transpose(glm::inverse(mat3(M)));
        glUniformMatrix3fv(norm_location, 1, GL_FALSE, value_ptr(normalMatrix));
        CHECK_GL_ERRORS;

        // Set Material values:
        glm::vec3 kd(0.5f, 0.7f, 0.5f);
        glUniform3fv(kd_location, 1, value_ptr(kd));
        glm::vec3 ks(0.5f, 0.5f, 0.5f);
        glUniform3fv(ks_location, 1, value_ptr(ks));
        float sh = 100.0f;
        glUniform1f(sh_location, sh);
        CHECK_GL_ERRORS;

        // Draw the surface.
        glBindVertexArray( m_surface_vao );
        glDrawArrays( GL_TRIANGLES, 0, get_n_vertices());
    m_shader->disable();

    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        (*child)->render_surface(W, proj, view, do_picking);
    }
}

void HBSurface::init_test(){
    eval_point(0,0,1,1);
    eval_normal(1,1,0,0);
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

glm::vec3 HBSurface::eval_normal(int x, int y, float u, float v){
    Eigen::Matrix4f Vx = (*cpsx).block<4,4>(x, y);
    Eigen::Matrix4f Vy = (*cpsy).block<4,4>(x, y);
    Eigen::Matrix4f Vz = (*cpsz).block<4,4>(x, y);

    // Take derivative along U.
    Eigen::Vector4f U(0.0f, 1.0f, 2.0f*u, 3.0f*u*u);
    Eigen::Vector4f V(1.0f, v, v*v, v*v*v);

    float Px = U.transpose()*B*Vx*B.transpose()*V;
    float Py = U.transpose()*B*Vy*B.transpose()*V;
    float Pz = U.transpose()*B*Vz*B.transpose()*V;

    glm::vec3 Du(Px, Py, Pz);

    // Take derivative along V.
    V = Eigen::Vector4f(0.0f, 1.0f, 2.0f*v, 3.0f*v*v);
    U = Eigen::Vector4f(1.0f, u, u*u, u*u*u);

    //std::cout << U << std::endl;

    Px = U.transpose()*B*Vx*B.transpose()*V;
    Py = U.transpose()*B*Vy*B.transpose()*V;
    Pz = U.transpose()*B*Vz*B.transpose()*V;

    glm::vec3 Dv(Px, Py, Pz);

    // Calculate the normal from derivatives.
    glm::vec3 P = eval_point(x, y, u, v);
    glm::vec3 A = Du;
    glm::vec3 B = Dv;
    glm::vec3 N = -glm::normalize(glm::cross(A, B));

    //std::cout << "P " << glm::to_string(P) << std::endl;
    //std::cout << "N " << glm::to_string(N) << std::endl << std::endl;

    return N;
}

glm::vec3* HBSurface::get_vertices(){
    if (refined){
        refined = false;
        if (nvert != 0){
            delete [] vert;
            delete [] norm;
        }

        nvert = 1;
        vert = new glm::vec3[get_n_vertices()];
        norm = new glm::vec3[get_n_vertices()];
    }

    //Iterate through the patches
    int elem = -6;
    for (int i = 0; i < npx; i++){
    for (int j = 0; j < npy; j++){
        //std::cout << (children[i][j] == this) << std::endl;
        threads[i*npy + j] = std::thread(&HBSurface::gen_vertex_patch, this, std::ref(elem), i, j);
        //gen_vertex_patch(elem, i, j);
    }
    }

    for (int i = 0; i < npx*npy; i++){
        threads[i].join();
    }

    return vert;
}

void HBSurface::gen_vertex_patch(int& elem, int i, int j){
    glm::vec3 Pbuffer[(res+1)*(res+1)];
    glm::vec3 Nbuffer[(res+1)*(res+1)];
    float h = 0.001f;
    if (children[i][j] == this){
        //Iterate in the patch
        for (int xp = 0; xp <= res; xp++){
            float u = (float) xp / (float) res;
            for (int yp = 0; yp <= res; yp++){
                float v = (float) yp / (float) res;
                glm::vec3 C = eval_point(i, j, u, v);
                Pbuffer[fxy(xp, yp)] =  C;

                //Find the normal
                glm::vec3 A0, A1;
                glm::vec3 B0, B1;
                float sign = -1.0f;
                if (xp == res && i < (npx-1)){
                    A0 = eval_point(i+1,j, h, v);
                    A1 = eval_point(i, j, u-h, v);
                    //sign *= -1.0f;
                } else if (xp == res){
                    A0 = eval_point(i, j, u, v);
                    A1 = eval_point(i, j, u-h, v);
                    //sign *= -1.0f;
                } else {
                    A0 = eval_point(i,j, u+h, v);
                    A1 = eval_point(i, j, u-h, v);
                }
                if (yp == res && j < (npy-1)){
                    B0 = eval_point(i, j+1, u, h);
                    B1 = eval_point(i, j, u, v-h);
                    //sign *= -1.0f;
                } else if (yp == res){
                    B0 = eval_point(i, j, u, v);
                    B1 = eval_point(i, j, u, v-h);
                    //sign *= -1.0f;
                } else {
                    B0 = eval_point(i,j, u, v+h);
                    B1 = eval_point(i, j, u, v-h);
                }
                //Nbuffer[fxy(xp, yp)] = sign*glm::normalize(glm::cross(A0-A1, B0-B1));
                //glm::vec3 result = sign*glm::normalize(glm::cross(A0-A1, B0-B1));
                //std::cout << "NA " << glm::to_string(result) << std::endl;
                Nbuffer[fxy(xp, yp)] = eval_normal(i, j, u, v);
                //std::cout << Pbuffer[fxy(xp, yp)].x << std::endl;
            }
        }
        mtx.lock();
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

                elem += 6;

                //Find the current square element of the surface (ugly)
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
        mtx.unlock();
    }
}

void HBSurface::set_mode_edit_cp(bool val){
    if (val == true){
        edit_cp_mode = true;
    } else {
        edit_cp_mode = false;
    }

    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        (*child)->set_mode_edit_cp(val);
    }
}

glm::vec3 HBSurface::get_edit_point(int i, int j){
    float P[3];
    Eigen::MatrixXf* CP [3] = {cpsx, cpsy, cpsz};

    for (int k = 0; k < 3; k++){
        P[k] =  (1.0f/6.0f)*((1.0f/6.0f)*(*CP[k])(i-1, j-1) + (4.0f/6.0f)*(*CP[k])(i, j-1) + (1.0f/6.0f)*(*CP[k])(i+1, j-1)) +
                (4.0f/6.0f)*((1.0f/6.0f)*(*CP[k])(i-1, j) + (4.0f/6.0f)*(*CP[k])(i, j) + (1.0f/6.0f)*(*CP[k])(i+1, j)) +
                (1.0f/6.0f)*((1.0f/6.0f)*(*CP[k])(i-1, j+1) + (4.0f/6.0f)*(*CP[k])(i, j+1) + (1.0f/6.0f)*(*CP[k])(i+1, j+1));
    }

    return glm::vec3(P[0], P[1], P[2]);
}

glm::vec3* HBSurface::get_edit_vertices(){
    if (refinedeps){
        refinedeps = false;
        ep_verts = new glm::vec3[get_n_ep_vertices()];
    }

    for (int i = 1; i < (ncpx-1); i++){
        for (int j = 1; j < (ncpy-1); j++){
            glm::vec3 P = get_edit_point(i,j);
            Cube EP(P.x, P.y, P.z);
            for (int k = 0; k < 36; k++){
                ep_verts[36*((i-1)*(ncpx-2) + (j-1)) + k] = EP.get_vertices()[k];
            }
        }
    }

    return ep_verts;
}

glm::vec3* HBSurface::get_cp_vertices(){
    if (refinedcps){
        if (n_cp_verts != 0){
            delete [] cp_verts;
        }
        refinedcps = false;

        n_cp_verts = 1;
        cp_verts = new glm::vec3[get_n_cp_vertices()];
    }

    for (int i = 0; i < ncpx; i++){
        for (int j = 0; j < ncpy; j++){
            float x = (*cpsx)(i,j);
            float y = (*cpsy)(i,j);
            float z = (*cpsz)(i,j);
            //float x,y,z;
            //std::cout << x << " " << y << " " << z << std::endl;
            Cube CP(x, y, z);
            for (int k = 0; k < 36; k++){
                cp_verts[36*(i*ncpy + j) + k] = CP.get_vertices()[k];
            }
        }
    }

    return cp_verts;
}

bool HBSurface::is_cp_selected(){
    return selected;
}

glm::vec3 HBSurface::get_ep_col(int idx){
    int x = (idx - my_idx_start)/(ncpx-2) + 1;
    int y = (idx - my_idx_start)%(ncpx-2) + 1;

    if (idx == sel_idx && selected == true){
        return glm::vec3(1.0f, 0.0f, 0.0f);
    } else if ((*cpmask)(x,y) != 1){
        return glm::vec3(0.5f, 0.5f, 0.5f);
    } else {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

glm::vec3 HBSurface::get_cp_col(int idx){
    int x = (idx - my_idx_start)/ncpy;
    int y = (idx - my_idx_start)%ncpy;

    if (idx == sel_idx && selected == true){
        //std::cout << "what I think is selected " << idx << std::endl;
        return glm::vec3(1.0f, 0.0f, 0.0f);
    } else if (idx == sel_idx_2 && selected == true){
        return glm::vec3(0.0f, 0.0f, 1.0f);
    } else if ((*cpmask)(x,y) != 1){
        return glm::vec3(0.5f, 0.5f, 0.5f);
    } else {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

glm::vec3 HBSurface::get_selected_cp_coords(){
    return glm::vec3((*cpsx)(sel_cp_i, sel_cp_j), (*cpsy)(sel_cp_i, sel_cp_j), (*cpsz)(sel_cp_i, sel_cp_j));
}

void HBSurface::select_cp(int idx, bool second, int level){
    if (level <= 0){
        if (is_cp_mine(idx)){
            std::cout << "running select idx " << idx << std::endl;
            selected = true;
            if (second){
                sel_idx_2 = idx;
                sel_cp_i_2 = (idx - my_idx_start)/ncpy;
                sel_cp_j_2 = (idx - my_idx_start)%ncpy;
            } else {
                if(edit_cp_mode){
                    sel_idx = idx;
                    sel_cp_i = (idx - my_idx_start)/(ncpy-2);
                    sel_cp_j = (idx - my_idx_start)%(ncpy-2);
                } else {
                    sel_idx = idx;
                    sel_cp_i = (idx - my_idx_start)/ncpy;
                    sel_cp_j = (idx - my_idx_start)%ncpy;
                }
            }
        } else {
            selected = false;
            sel_idx = -1;
            sel_idx_2 = -1;
        }
    } else {
        selected = false;
        sel_idx = -1;
        sel_idx_2 = -1;
    }

    if (level != 0){
    level -= 1;
        for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
            (*child)->select_cp(idx, second, level);
        }
    }
}

bool HBSurface::is_cp_mine(int idx){
    return (idx >= my_idx_start && idx < (my_idx_start + get_n_cps()));
}

void HBSurface::move_selected_cp(glm::vec3 delta){
    //Require delta to be CP's coordinates
    if(selected){
        if (edit_cp_mode){
            int i = sel_cp_i+1;
            int j = sel_cp_j+1;
            if ( glm::length(delta) < 1000 && ((*cpmask)(i, j) == 1) ){
                glm::vec3 P = get_edit_point(i, j);
                glm::vec3 V( (*cpsx)(i, j), (*cpsy)(i, j), (*cpsz)(i, j) );
                glm::vec3 D = (36.0f/16.0f)*(delta - P) + V;
                (*Ocpsx)(i, j) = D.x - (*Rcpsx)(i, j);
                (*Ocpsy)(i, j) = D.y - (*Rcpsy)(i, j);
                (*Ocpsz)(i, j) = D.z - (*Rcpsz)(i, j);
            }
        } else {
            if ( glm::length(delta) < 1000  && ((*cpmask)(sel_cp_i, sel_cp_j) == 1) ){
                (*Ocpsx)(sel_cp_i, sel_cp_j) = delta.x - (*Rcpsx)(sel_cp_i, sel_cp_j);
                (*Ocpsy)(sel_cp_i, sel_cp_j) = delta.y - (*Rcpsy)(sel_cp_i, sel_cp_j);
                (*Ocpsz)(sel_cp_i, sel_cp_j) = delta.z - (*Rcpsz)(sel_cp_i, sel_cp_j);
                update_cps();
            }
        }

        update_references();
    } else {
        for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
            (*child)->move_selected_cp(delta);
        }
    }
}

void HBSurface::update_references(){
    int dim_x, dim_y;
    //std::cout << "Updating References" << std::endl;

    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        dim_x = 7 + std::abs((*child)->parent_sel_cp_i - (*child)->parent_sel_cp_i_2)*2;
        dim_y = 7 + std::abs((*child)->parent_sel_cp_j - (*child)->parent_sel_cp_j_2)*2;
        //std::cout << "Updating a child" << std::endl;

        Eigen::MatrixXf RX(dim_x, dim_y);
        Eigen::MatrixXf RY(dim_x, dim_y);
        Eigen::MatrixXf RZ(dim_x, dim_y);

        get_split_reference(RX, RY, RZ, (*child)->parent_sel_cp_i, (*child)->parent_sel_cp_j, (*child)->parent_sel_cp_i_2, (*child)->parent_sel_cp_j_2);

        *((*child)->Rcpsx) = RX;
        *((*child)->Rcpsy) = RY;
        *((*child)->Rcpsz) = RZ;
        (*child)->update_cps();

        //std::cout << "Done updating a child" << std::endl;
    }

    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        (*child)->update_references();
    }
}

void HBSurface::split_patch(int i, int j, Matrix5f& X, Matrix5f& Y, Matrix5f& Z){
    update_cps();
    Eigen::Matrix4f Vx = (*cpsx).block<4,4>(i-2, j-2);
    Eigen::Matrix4f Vy = (*cpsy).block<4,4>(i-2, j-2);
    Eigen::Matrix4f Vz = (*cpsz).block<4,4>(i-2, j-2);

    Eigen::Matrix4f W1x = alphaA1 * Vx * alphaA1.transpose();
    Eigen::Matrix4f W1y = alphaA1 * Vy * alphaA1.transpose();
    Eigen::Matrix4f W1z = alphaA1 * Vz * alphaA1.transpose();

    Eigen::Matrix4f W2x = alphaA2 * Vx * alphaA1.transpose();
    Eigen::Matrix4f W2y = alphaA2 * Vy * alphaA1.transpose();
    Eigen::Matrix4f W2z = alphaA2 * Vz * alphaA1.transpose();

    Eigen::Matrix4f W3x = alphaA1 * Vx * alphaA2.transpose();
    Eigen::Matrix4f W3y = alphaA1 * Vy * alphaA2.transpose();
    Eigen::Matrix4f W3z = alphaA1 * Vz * alphaA2.transpose();

    Eigen::Matrix4f W4x = alphaA2 * Vx * alphaA2.transpose();
    Eigen::Matrix4f W4y = alphaA2 * Vy * alphaA2.transpose();
    Eigen::Matrix4f W4z = alphaA2 * Vz * alphaA2.transpose();

    // std::cout << W1z << std::endl << std::endl;
    // std::cout << W2z << std::endl << std::endl;
    // std::cout << W3z << std::endl << std::endl;
    // std::cout << W4z << std::endl << std::endl;

    X.block<4,4>(0,0) = W1x;
    Y.block<4,4>(0,0) = W1y;
    Z.block<4,4>(0,0) = W1z;

    X.block<4,4>(1,0) = W2x;
    Y.block<4,4>(1,0) = W2y;
    Z.block<4,4>(1,0) = W2z;

    X.block<4,4>(0,1) = W3x;
    Y.block<4,4>(0,1) = W3y;
    Z.block<4,4>(0,1) = W3z;

    X.block<4,4>(1,1) = W4x;
    Y.block<4,4>(1,1) = W4y;
    Z.block<4,4>(1,1) = W4z;

    //std::cout << X << std::endl;
    //std::cout << Z << std::endl;
}

void HBSurface::split_selected_cp(){
    if(selected){
        split();
    } else {
        for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
            (*child)->split_selected_cp();
        }
    }
}

void HBSurface::get_split_reference(Eigen::MatrixXf& RX, Eigen::MatrixXf& RY, Eigen::MatrixXf& RZ, int cp_i, int cp_j, int cp_i_2, int cp_j_2){
    Matrix5f X, Y, Z; // individual refined patches

    int i_start, i_end, j_start, j_end;

    if ( cp_i < cp_i_2 ){
        i_start = cp_i;
        i_end = cp_i_2;
    } else {
        i_start = cp_i_2;
        i_end = cp_i;
    }

    if ( cp_j < cp_j_2 ){
        j_start = cp_j;
        j_end = cp_j_2;
    } else {
        j_start = cp_j_2;
        j_end = cp_j;
    }

    for (int i = i_start; i <= (i_end+1); i++){
        for (int j = j_start; j <= (j_end+1); j++){
            split_patch(i, j, X, Y, Z);
            (RX).block<5,5>((i-i_start)*2, (j-j_start)*2) = X;
            (RY).block<5,5>((i-i_start)*2, (j-j_start)*2) = Y;
            (RZ).block<5,5>((i-i_start)*2, (j-j_start)*2) = Z;
        }
    }
}

void HBSurface::split(){
    if (sel_idx_2 < 0){
        sel_cp_i_2 = sel_cp_i;
        sel_cp_j_2 = sel_cp_j;
    }

    int i_start, i_end, j_start, j_end;

    if ( sel_cp_i < sel_cp_i_2 ){
        i_start = sel_cp_i;
        i_end = sel_cp_i_2;
    } else {
        i_start = sel_cp_i_2;
        i_end = sel_cp_i;
    }

    if ( sel_cp_j < sel_cp_j_2 ){
        j_start = sel_cp_j;
        j_end = sel_cp_j_2;
    } else {
        j_start = sel_cp_j_2;
        j_end = sel_cp_j;
    }


    // Go through all children and find conflicting patches on this level.
    std::vector<HBSurface*> adjacent_children;
    bool adjacent_exists = false;
    int max_ci, max_cj, min_ci, min_cj;
    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        bool adjacent = false;
        int ci_start = (*child)->parent_sel_cp_i;
        int cj_start = (*child)->parent_sel_cp_j;
        int ci_end = (*child)->parent_sel_cp_i_2;
        int cj_end = (*child)->parent_sel_cp_j_2;

        if (child == child_list.begin()){
            max_ci = i_end;
            max_cj = j_end;
            min_ci = i_start;
            min_cj = j_start;
        }

        for (int i = i_start; i <= i_end; i++){
            for (int j = j_start; j <= j_end; j++){
                for (int k = ci_start; k <= ci_end; k++){
                    for (int l = cj_start; l <= cj_end; l++){
                        if ( std::abs(k - i) <= 1 && std::abs(l - j) <= 1 ){
                            adjacent = true;
                            adjacent_exists = true;
                            std::cout << "Found adjacent surface to combine with." << std::endl;
                        }
                    }
                }
            }
        }

        if(adjacent){
            adjacent_children.push_back(*child);

            if (ci_start < min_ci){
                min_ci = ci_start;
            }
            if (cj_start < min_cj){
                min_cj = cj_start;
            }
            if (ci_end > max_ci){
                max_ci = ci_end;
            }
            if (cj_end > max_cj){
                max_cj = cj_end;
            }
        }
    }

    if( adjacent_exists){
        std::cout << "Creating new combined surface " << min_ci << "," << min_cj << " to " << max_ci << "," << max_cj << std::endl;
        sel_cp_i = min_ci;
        sel_cp_j = min_cj;
        sel_cp_i_2 = max_ci;
        sel_cp_j_2 = max_cj;
    }

    int dim_x, dim_y;
    dim_x = 7 + std::abs(sel_cp_i - sel_cp_i_2)*2;
    dim_y = 7 + std::abs(sel_cp_j - sel_cp_j_2)*2;

    std::cout << "Adding new surface." << std::endl;
    std::cout << "Dim " << dim_x << " by " << dim_y << " control points." << std::endl;

    Eigen::MatrixXf RX(dim_x, dim_y);
    Eigen::MatrixXf RY(dim_x, dim_y);
    Eigen::MatrixXf RZ(dim_x, dim_y);

    get_split_reference(RX, RY, RZ, sel_cp_i, sel_cp_j, sel_cp_i_2, sel_cp_j_2);

    int pdim_x = std::abs(sel_cp_i - sel_cp_i_2)+2;
    int pdim_y = std::abs(sel_cp_j - sel_cp_j_2)+2;
    std::cout << "Dim " << pdim_x << " by " << pdim_y << " patches." << std::endl;
    npatches -= pdim_x * pdim_y; //Needs to be updated when deleting a sub surface.

    // Set up the child.
    has_children = true;
    HBSurface* new_surface = new HBSurface(GLapp, m_shader, b_shader, 2*pdim_x, 2*pdim_y, res);
    *(new_surface->Rcpsx) = RX;
    *(new_surface->Rcpsy) = RY;
    *(new_surface->Rcpsz) = RZ;
    new_surface->parent_sel_cp_i = sel_cp_i;
    new_surface->parent_sel_cp_j = sel_cp_j;
    new_surface->parent_sel_cp_i_2 = sel_cp_i_2;
    new_surface->parent_sel_cp_j_2 = sel_cp_j_2;
    new_surface->update_cps();
    *new_surface->cpmask = Eigen::MatrixXi::Constant(dim_x, dim_y, 0);
    //std::cout << *new_surface->cpmask << std::endl;
    (*new_surface->cpmask).block(3, 3, dim_x - 6, dim_y - 6) = Eigen::MatrixXi::Constant(dim_x - 6, dim_y - 6, 1);
    //std::cout << *new_surface->cpmask << std::endl;
    new_surface->has_parent = true;
    new_surface->parent = this;
    std::cout << "New surface has " << new_surface->get_n_cps() << " number of cps." << std::endl;
    //std::cout << *new_surface->cpsx << std::endl;
    child_list.push_back(new_surface);


    //This for checking patches
    for (int i = -2; i <= -1 + (i_end-i_start); i++){
        for (int j = -2; j <= -1 + (j_end - j_start); j++){
            children[i_start+i][j_start+j] = new_surface;
        }
    }

    for (std::vector<HBSurface*>::iterator child = adjacent_children.begin(); child != adjacent_children.end(); child++){
        std::cout << "Adding back " << (((*child)->npx)/2) * (((*child)->npy)/2) << " patches." << std::endl;
        npatches += (((*child)->npx)/2) * (((*child)->npy)/2);
        int i = 2*((*child)->parent_sel_cp_i - new_surface->parent_sel_cp_i);
        int j = 2*((*child)->parent_sel_cp_j - new_surface->parent_sel_cp_j);

        (*new_surface->Ocpsx).block(i, j, (*child)->ncpx, (*child)->ncpy) = *(*child)->Ocpsx;
        (*new_surface->Ocpsy).block(i, j, (*child)->ncpx, (*child)->ncpy) = *(*child)->Ocpsy;
        (*new_surface->Ocpsz).block(i, j, (*child)->ncpx, (*child)->ncpy) = *(*child)->Ocpsz;

        for (std::vector<HBSurface*>::iterator childp = child_list.begin(); childp != child_list.end(); childp++){
            if( *childp == *child ){
                child_list.erase(childp);
            }
        }
    }
    if (adjacent_exists){
        npatches++;
    }

    std::cout << "Completed Split" << std::endl;
    // std::cout << RX << std::endl << std::endl;
    // std::cout << RZ << std::endl << std::endl;
}

unsigned int HBSurface::get_selected_cp_idx(){
    if (sel_idx >= 0){
        return sel_idx;
    } else {
        int child_idx;
        for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
            child_idx = (*child)->get_selected_cp_idx();
            if (child_idx >= 0)
                return child_idx;
        }
    }
}


unsigned int HBSurface::get_n_eps(){
    return (ncpy-2)*(ncpx-2);
}

unsigned int HBSurface::get_n_ep_vertices(){
    return 36*get_n_eps();
}

unsigned int HBSurface::get_ep_vertices_size(){
    return 3*get_n_ep_vertices()*sizeof(float);
}

unsigned int HBSurface::get_n_cps(){
    return ncpy*ncpx;
}

unsigned int HBSurface::get_cp_vertices_size(){
    return 3*get_n_cp_vertices()*sizeof(float);
}

unsigned int HBSurface::get_n_cp_vertices(){
    return 36*get_n_cps();
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
