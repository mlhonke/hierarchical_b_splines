#pragma once

#include <glm/glm.hpp>

#include "cs488-framework/CS488Window.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"

#include "HBSurface.hpp"

class A1 : public CS488Window {
public:
	A1();
	virtual ~A1();

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
	void initBuffers();
	void initSurface();
	void updateSurface(glm::mat4 W);
	void updateCPs(glm::mat4 W);
	void updateLighting();
	glm::mat4 get_W();

	//Fields for object picking
	bool do_picking = false;

	// Fields related to surface properties
	int npx = 4;
	int npy = 4;

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
	GLuint m_surface_vao; // Vertex Array Object
	GLuint m_surface_vbo; // Vertex Buffer Object
	GLuint m_surface_normals_vbo;
	GLuint m_cp_vao;
	GLuint m_cp_vbo;

	// Fields related to movement.
	bool dragging = false;
	float old_x = 0;
	float rot_rads = 0.0f;
	float old_rot_rads = 0.0f;
	float zoom = 1.0f;

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
