#include "HBSurface.hpp"
#include "A1.hpp"

HBSurface::HBSurface(A1* GLapp, int npx, int npy, int res):npx(npx), npy(npy), res(res), GLapp(GLapp){
    children = new HBSurface**[npx];
    for (int i = 0; i<npx; i++){
        children[i] = new HBSurface*[npy];
        for (int j=0; j<npy; j++){
            children[i][j] = NULL;
        }
    }

    Pbuffer = new glm::vec3[(res+1)*(res+1)];
    Nbuffer = new glm::vec3[(res+1)*(res+1)];
    npatches = npx*npy;
    int m = npx + k - 2;
    int n = npy + l - 2;
    ncpx = m + 1;
    ncpy = n + 1;
    ncps = ncpx*ncpy;
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

    GLapp->gen_buffers(vaos, vbos);
    m_surface_vao = vaos[0];
	m_cp_vao = vaos[1];

	m_surface_vbo = vbos[0];
	m_surface_normals_vbo = vbos[1];
	m_cp_vbo = vbos[2];
    //init_test();
}

void HBSurface::render(ShaderProgram& m_shader, ShaderProgram& b_shader, glm::mat4 W, glm::mat4 proj, glm::mat4 view, bool do_picking){
    GLint m_normalAttribLocation = m_shader.getAttribLocation("normal");
	GLint kd_location = m_shader.getUniformLocation("material.kd");
	GLint ks_location = m_shader.getUniformLocation("material.ks");
	GLint sh_location = m_shader.getUniformLocation("material.shininess");
	GLint norm_location = m_shader.getUniformLocation("NormalMatrix");
    GLint Pos = m_shader.getAttribLocation( "position" );

    // Set up the uniforms
    GLint Pers = m_shader.getUniformLocation( "Perspective" );
    GLint Model = m_shader.getUniformLocation( "ModelView" );

    //Set up basic shader uniforms (non-phong)
    GLint P_uni = b_shader.getUniformLocation( "P" );
    GLint V_uni = b_shader.getUniformLocation( "V" );
    GLint M_uni = b_shader.getUniformLocation( "M" );
    GLint col_uni = b_shader.getUniformLocation( "colour" );
    GLint posAttrib2 = b_shader.getAttribLocation( "position");

    CHECK_GL_ERRORS;

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

	m_shader.enable();
		glEnable( GL_DEPTH_TEST );
		glEnable( GL_CULL_FACE );
		glfwSwapInterval(1);

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
	m_shader.disable();

    //CP vertices
	glBindVertexArray(m_cp_vao);
	glBindBuffer(GL_ARRAY_BUFFER, m_cp_vbo);
	glBufferData(GL_ARRAY_BUFFER, get_cp_vertices_size(), get_cp_vertices(), GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray( posAttrib2 );
	glVertexAttribPointer( posAttrib2, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

	glBindVertexArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

	CHECK_GL_ERRORS;

	b_shader.enable();
		glEnable( GL_DEPTH_TEST );
		glEnable( GL_CULL_FACE );
		glUniformMatrix4fv( P_uni, 1, GL_FALSE, value_ptr( proj ) );
		glUniformMatrix4fv( V_uni, 1, GL_FALSE, value_ptr( view ) );
		glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W ) );
		glBindVertexArray( m_cp_vao );
		if (do_picking){
			for (int idx = 0; idx < get_n_cps(); idx++){
				float r = float(idx&0xff) / 255.0f;
				float g = float((idx>>8)&0xff) / 255.0f;
				float b = float((idx>>16)&0xff) / 255.0f;
				glUniform3f( col_uni, r, g, b);
				glDrawArrays(GL_TRIANGLES, idx*36, 36);
			}
		} else {
			for (int idx = 0; idx < get_n_cps(); idx++){
				glm::vec3 cp_col = get_cp_col(idx);
				glUniform3f( col_uni, cp_col.x, cp_col.y, cp_col.z );
				glDrawArrays( GL_TRIANGLES, idx*36, 36);
			}
		}
		CHECK_GL_ERRORS;
	b_shader.disable();

    for (std::vector<HBSurface*>::iterator child = child_list.begin(); child != child_list.end(); child++){
        (*child)->render(m_shader, b_shader, W, proj, view, do_picking);
    }
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
            delete [] norm;
        }

        nvert = 1;
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
                cp_verts[36*(i*ncpx + j) + k] = CP.get_vertices()[k];
            }
        }
    }

    return cp_verts;
}

bool HBSurface::is_cp_selected(){
    return selected;
}

glm::vec3 HBSurface::get_cp_col(int idx){
    if (idx == sel_idx && selected == true){
        //std::cout << "what I think is selected " << idx << std::endl;
        return glm::vec3(1.0f, 0.0f, 0.0f);
    } else {
        return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

glm::vec3 HBSurface::get_selected_cp_coords(){
    return glm::vec3((*cpsx)(sel_cp_i, sel_cp_j), (*cpsy)(sel_cp_i, sel_cp_j), (*cpsz)(sel_cp_i, sel_cp_j));
}

void HBSurface::select_cp(int idx){
    if (idx < ncps){
        selected = true;
        sel_idx = idx;
        sel_cp_i = idx/ncpx;
        sel_cp_j = idx%ncpx;
    }
}

void HBSurface::move_selected_cp(glm::vec3 delta){
    //Require delta to be CP's coordinates
    (*cpsx)(sel_cp_i, sel_cp_j) = delta.x;
    (*cpsy)(sel_cp_i, sel_cp_j) = delta.y;
    (*cpsz)(sel_cp_i, sel_cp_j) = delta.z;
}

void HBSurface::split_patch(int i, int j, Matrix5f& X, Matrix5f& Y, Matrix5f& Z){
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
    std::cout << "Adding new surface." << std::endl;
    Matrix5f X, Y, Z; // individual refined patches
    Matrix7f RX, RY, RZ; // refined points about central point
    split_patch(sel_cp_i, sel_cp_j, X, Y, Z);
    RX.block<5,5>(0,0) = X;
    RY.block<5,5>(0,0) = Y;
    RZ.block<5,5>(0,0) = Z;
    split_patch(sel_cp_i+1, sel_cp_j, X, Y, Z );
    RX.block<5,5>(2,0) = X;
    RY.block<5,5>(2,0) = Y;
    RZ.block<5,5>(2,0) = Z;
    split_patch(sel_cp_i, sel_cp_j+1, X, Y, Z);
    RX.block<5,5>(0,2) = X;
    RY.block<5,5>(0,2) = Y;
    RZ.block<5,5>(0,2) = Z;
    split_patch(sel_cp_i+1, sel_cp_j+1, X, Y, Z);
    RX.block<5,5>(2,2) = X;
    RY.block<5,5>(2,2) = Y;
    RZ.block<5,5>(2,2) = Z;

    has_children = true;
    HBSurface* new_surface = new HBSurface(GLapp, 4,4, res);
    (*new_surface->cpsx) = RX;
    (*new_surface->cpsy) = RY;
    (*new_surface->cpsz) = RZ;
    child_list.push_back(new_surface);

    //This for checking patches
    for (int i = -2; i <= -1; i++){
        for (int j = -2; j <= -1; j++){
            children[sel_cp_i+i][sel_cp_j+j] = new_surface;
        }
    }

    // std::cout << RX << std::endl << std::endl;
    // std::cout << RZ << std::endl << std::endl;
}

unsigned int HBSurface::get_selected_cp_idx(){
    return sel_idx;
}

unsigned int HBSurface::get_n_cps(){
    return ncps;
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
