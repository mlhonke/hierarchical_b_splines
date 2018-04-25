#pragma once

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "cs488-framework/CS488Window.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"

#include "Eigen/Dense"
#include "trackball.hpp"
#include <fstream>

class HBSurface;

class A1 : public CS488Window {
public:
	A1();
	virtual ~A1();
	void gen_buffers(GLuint* vaos, GLuint* vbos);

protected:
	virtual void init() override;
	virtual void appLogic() override;
	virtual void guiLogic() override;
	virtual void draw() override;
	virtual void cleanup() override;

	virtual bool cursorEnterWindowEvent(int entered) override;
	virtual bool mouseMoveEvent(double xPos, double yPos) override;
	virtual bool mouseButtonInputEvent(int button, int actions, int mods) override;
	virtual bool mouseScrollEvent(double xOffSet, double yOffSet) override;
	virtual bool windowResizeEvent(int width, int height) override;
	virtual bool keyInputEvent(int key, int action, int mods) override;

private:
	void restart_app();
	void updateLighting();
	void rotate_surface(float xPos, float yPos);
	glm::mat4 update_W();
	int pick_object();
	glm::vec3 GetOGLPos(float x, float y, float depth);
	float getDepth(float x, float y);
	float depth_val;
	bool depth_set = false;
	bool drag_cp = false;
	bool edit_cp = false;
	glm::mat4 W;
	glm::mat4 W_rot;
	glm::mat4 W_rot_old;
	glm::mat4 W_trans;
	glm::mat4 W_trans_init;
	glm::mat4 W_trans_old;
	glm::mat4 W_trans_center;
	glm::mat4 W_scale;

	//Saving/loading variables
	template<typename Derived>
	void write_matrix(std::ofstream output_file, const Eigen::MatrixBase<Derived>& M);
	void load_matrix();
	void save();
	std::fstream save_file;


	//Fields for object picking
	bool do_picking = false;
	bool lctrl = false;
	bool hide_surface = false;

	// Fields related to surface properties
	int npx = 4;
	int npy = 4;
	int level = -1;

	// Fields related to the shader and uniforms.
	ShaderProgram m_shader;
	ShaderProgram b_shader;

	GLint Pers; // Uniform location for Projection matrix.
	GLint Model; // Uniform location for Model matrix.

	GLint P_uni;
	GLint V_uni;
	GLint M_uni;
	GLint col_uni;   // Uniform location for cube colour.
	GLint posAttrib2;

	// Fields related to grid geometry.
	HBSurface* surface;

	// Fields related to movement.
	bool dragging = false;
	bool dragging_right = false;
	float old_x = 0;
	float old_y = 0;
	float rot_rads = 0.0f;
	float old_rot_rads = 0.0f;
	float zoom = 1.0f;
	float drag_x = 0;
	float drag_y = 0;

	// Matrices controlling the camera and projection.
	glm::mat4 proj;
	glm::mat4 view;

	GLint light_position;
	GLint light_colour;
	GLint light_ambient;
	GLint Pos;
	GLint m_normalAttribLocation;
	GLint ks_location;
	GLint kd_location;
	GLint sh_location;
	GLint norm_location;
};
