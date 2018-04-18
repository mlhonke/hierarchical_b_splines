#include "A1.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

#include <iostream>
#include <algorithm>
#include <math.h>

#include <imgui/imgui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;
using namespace std;
int DIM = 2;

//----------------------------------------------------------------------------------------
// Constructor
A1::A1()
{
	initSurface();
}

void A1::initSurface(){
	surface = new HBSurface(npx, npy, 4);
}

//----------------------------------------------------------------------------------------
// Destructor
A1::~A1()
{
	delete surface;
}

//----------------------------------------------------------------------------------------
/*
 * Called once, at program start.
 */
void A1::init()
{
	// Set the background colour.
	glClearColor( 0.3, 0.5, 0.7, 1.0 );

	// Build the shader
	m_shader.generateProgramObject();
	m_shader.attachVertexShader( getAssetFilePath( "VertexShader.vs" ).c_str() );
	m_shader.attachFragmentShader( getAssetFilePath( "FragmentShader.fs" ).c_str() );
	m_shader.link();

	// Set up the uniforms
	P_uni = m_shader.getUniformLocation( "Perspective" );
	M_uni = m_shader.getUniformLocation( "ModelView" );

	//Set up lighting
	light_position = m_shader.getUniformLocation("light.position");
	light_colour = m_shader.getUniformLocation("light.rgbIntensity");
	light_ambient = m_shader.getUniformLocation("ambientIntensity");

	//Initialize the vertex buffers
	initBuffers();

	//Initialize the B-Spline surface
	updateSurface();

	// Set up initial view and projection matrices (need to do this here,
	// since it depends on the GLFW window being set up correctly).
	view = glm::lookAt(
		glm::vec3( -2.0f, 4.0f, -2.0f ),
		glm::vec3( 0.0f, 0.0f, 0.0f ),
		glm::vec3( 0.0f, 1.0f, 0.0f ) );

	proj = glm::perspective(
		glm::radians( 45.0f ),
		float( m_framebufferWidth ) / float( m_framebufferHeight ),
		1.0f, 1000.0f );
}

/*
 * initBuffers: Initializes vertex arrays and buffers for the grid,
 * cursor and cubes.
 */
void A1::initBuffers()
{
	GLuint vaos[1];
	GLuint vbos[2];

	// Create the vertex array to record buffer assignments.
	glGenVertexArrays( 1, vaos );

	// Create the buffer.
	glGenBuffers( 2, vbos );

	//Assign vertex arrays and buffers to named variables for easier tracking.
	m_surface_vao = vaos[0];

	m_surface_vbo = vbos[0];
	m_surface_normals_vbo = vbos[1];
}

/*
initSurface: Updates the surface vertices in the buffer.
*/
void A1::updateSurface(){
	glBindVertexArray(m_surface_vao);

	//Position data
	glBindBuffer(GL_ARRAY_BUFFER, m_surface_vbo);
	glBufferData(GL_ARRAY_BUFFER, surface->get_vertices_size(), surface->get_vertices(), GL_DYNAMIC_DRAW);
	GLint posAttrib2 = m_shader.getAttribLocation( "position" );
	glEnableVertexAttribArray( posAttrib2 );
	glVertexAttribPointer( posAttrib2, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

	//Normal data
	glBindBuffer(GL_ARRAY_BUFFER, m_surface_normals_vbo);
	glBufferData(GL_ARRAY_BUFFER, surface->get_vertices_size(), surface->get_normals(), GL_DYNAMIC_DRAW);
	GLint m_normalAttribLocation = m_shader.getAttribLocation("normal");
	glEnableVertexAttribArray(m_normalAttribLocation);
	glVertexAttribPointer(m_normalAttribLocation, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

	glBindVertexArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, before guiLogic().
 */
void A1::appLogic()
{

}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after appLogic(), but before the draw() method.
 */
void A1::guiLogic()
{
	// We already know there's only going to be one window, so for
	// simplicity we'll store button states in static local variables.
	// If there was ever a possibility of having multiple instances of
	// A1 running simultaneously, this would break; you'd want to make
	// this into instance fields of A1.
	static bool showTestWindow(false);
	static bool showDebugWindow(true);

	ImGuiWindowFlags windowFlags(ImGuiWindowFlags_AlwaysAutoResize);
	float opacity(0.5f);

	ImGui::Begin("Debug Window", &showDebugWindow, ImVec2(100,100), opacity, windowFlags);
		if( ImGui::Button( "Quit Application" ) ) {
			glfwSetWindowShouldClose(m_window, GL_TRUE);
		}
		if( ImGui::Button( "Restart" ) ) {
			restart_app();
		}
/*
		// For convenience, you can uncomment this to show ImGui's massive
		// demonstration window right in your application.  Very handy for
		// browsing around to get the widget you want.  Then look in
		// shared/imgui/imgui_demo.cpp to see how it's done.
		if( ImGui::Button( "Test Window" ) ) {
			showTestWindow = !showTestWindow;
		}
*/
		ImGui::Text( "Framerate: %.1f FPS", ImGui::GetIO().Framerate );

	ImGui::End();

	if( showTestWindow ) {
		ImGui::ShowTestWindow( &showTestWindow );
	}
}

/*
 * restart_app: Clear all cubes from screen, reset all values to default.
 */
void A1::restart_app(){
	zoom = 1.0f;
	rot_rads = 0.0f;
	old_rot_rads = 0.0f;
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after guiLogic().
 */
void A1::draw()
{
	updateSurface();
	// Create a global transformation for the model (centre it).
	mat4 W;
	W = glm::rotate(W, rot_rads, vec3(0, 1, 0));
	W = glm::scale(W, vec3(zoom, zoom, zoom));
	W = glm::translate( W, vec3( -(1.0f + (npx+1)/2.0f), 0, -(1.0f + (npy+1)/2.0f) ) );

	m_shader.enable();
		glEnable( GL_DEPTH_TEST );
		glEnable( GL_CULL_FACE );
		glfwSwapInterval(1);

		glUniformMatrix4fv( P_uni, 1, GL_FALSE, value_ptr( proj ) );
		W = view * W;
		glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W ) );
		GLint location = m_shader.getUniformLocation("NormalMatrix");
		mat3 normalMatrix = glm::transpose(glm::inverse(mat3(W)));
		glUniformMatrix3fv(location, 1, GL_FALSE, value_ptr(normalMatrix));
		CHECK_GL_ERRORS;

		// Lighting
		glm::vec3 lpos(0.0f, 10.0f, 0.0f);
		glUniform3fv(light_position, 1, value_ptr(lpos));
		glm::vec3 lcol(0.8f, 0.8f, 0.8f);
		glUniform3fv(light_colour, 1, value_ptr(lcol));
		glm::vec3 lamb(0.1f);
		glUniform3fv(light_ambient, 1, value_ptr(lamb));
		CHECK_GL_ERRORS;

		// Set Material values:
		location = m_shader.getUniformLocation("material.kd");
		glm::vec3 kd(0.5f, 0.7f, 0.5f);
		glUniform3fv(location, 1, value_ptr(kd));
		CHECK_GL_ERRORS;
		location = m_shader.getUniformLocation("material.ks");
		glm::vec3 ks(0.8f, 0.8f, 0.8f);
		glUniform3fv(location, 1, value_ptr(ks));
		CHECK_GL_ERRORS;
		location = m_shader.getUniformLocation("material.shininess");
		float sh = 100.0f;
		glUniform1f(location, sh);
		CHECK_GL_ERRORS;

		// Draw the surface.
		glBindVertexArray( m_surface_vao );
		glDrawArrays( GL_TRIANGLES, 0, surface->get_n_vertices());

	m_shader.disable();

	// Restore defaults
	glBindVertexArray( 0 );

	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
/*
 * Called once, after program is signaled to terminate.
 */
void A1::cleanup()
{}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles cursor entering the window area events.
 */
bool A1::cursorEnterWindowEvent (
		int entered
) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse cursor movement events.
 */
bool A1::mouseMoveEvent(double xPos, double yPos)
{
	bool eventHandled(false);
	if (!ImGui::IsMouseHoveringAnyWindow()) {
		// Rotate the model if mouse is being dragged.
		if (dragging){
			rot_rads = old_rot_rads - 2*3.1415*(old_x - xPos)/m_framebufferWidth;
		} else {
		// Track the position of the mouse in preparation for rotation.
			old_x = xPos;
		}
	}
	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse button events.
 */
bool A1::mouseButtonInputEvent(int button, int actions, int mods) {
	bool eventHandled(false);

	if (!ImGui::IsMouseHoveringAnyWindow()) {
		// The user clicked in the window.  If it's the left
		// mouse button, initiate a rotation.
		if (button == GLFW_MOUSE_BUTTON_LEFT && actions == GLFW_PRESS){
			dragging = true;
		}
		if (button == GLFW_MOUSE_BUTTON_LEFT && actions == GLFW_RELEASE){
			dragging = false;
			rot_rads = fmod(rot_rads, 2*3.1415);
			old_rot_rads = rot_rads;
		}
	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse scroll wheel events.
 */
bool A1::mouseScrollEvent(double xOffSet, double yOffSet) {
	bool eventHandled(false);

	//Scroll down
	if (yOffSet < 0 && zoom >= 0.1){
		zoom -= 0.05;
	}

	//Scroll up
	if (yOffSet > 0 && zoom <= 2.0){
		zoom += 0.05;
	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles window resize events.
 */
bool A1::windowResizeEvent(int width, int height) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles key input events.
 */
bool A1::keyInputEvent(int key, int action, int mods) {
	bool eventHandled(false);

	// Fill in with event handling code...
	if( action == GLFW_PRESS ) {
		//Reset the program.
		if (key == GLFW_KEY_R){
			restart_app();

			eventHandled = true;
		}

		//Quit the program.
		if (key == GLFW_KEY_Q){
			glfwSetWindowShouldClose(m_window, GL_TRUE);

			eventHandled = true;
		}
	}

	return eventHandled;
}
