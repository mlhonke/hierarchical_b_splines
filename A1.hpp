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
	void updateSurface();

	// Fields related to surface properties
	int npx = 4;
	int npy = 4;

	// Fields related to the shader and uniforms.
	ShaderProgram m_shader;
	GLint P_uni; // Uniform location for Projection matrix.
	GLint V_uni; // Uniform location for View matrix.
	GLint M_uni; // Uniform location for Model matrix.
	GLint col_uni;   // Uniform location for cube colour.

	// Fields related to grid geometry.
	HBSurface* surface;
	GLuint m_surface_vao; // Vertex Array Object
	GLuint m_surface_vbo; // Vertex Buffer Object
	GLuint m_surface_normals_vbo;

	// Fields related to movement.
	bool dragging;
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
};
