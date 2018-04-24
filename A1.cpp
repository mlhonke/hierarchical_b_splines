#include "A1.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

#include <iostream>
#include <algorithm>
#include <math.h>

#include <imgui/imgui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "HBSurface.hpp"

using namespace glm;
using namespace std;
int DIM = 2;

//----------------------------------------------------------------------------------------
// Constructor
A1::A1()
{
}

//----------------------------------------------------------------------------------------
// Destructor
A1::~A1()
{
	delete surface;
}

void A1::gen_buffers(GLuint* vaos, GLuint* vbos){
	// Create the vertex array to record buffer assignments.
	glGenVertexArrays( 2, vaos );

	// Create the buffer.
	glGenBuffers( 3, vbos );
	CHECK_GL_ERRORS;
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

	// Build the shader
	b_shader.generateProgramObject();
	b_shader.attachVertexShader( getAssetFilePath( "BasicVertexShader.vs" ).c_str() );
	b_shader.attachFragmentShader( getAssetFilePath( "BasicFragmentShader.fs" ).c_str() );
	b_shader.link();

	//Set up lighting
	light_position = m_shader.getUniformLocation("light.position");
	light_colour = m_shader.getUniformLocation("light.rgbIntensity");
	light_ambient = m_shader.getUniformLocation("ambientIntensity");

	surface = new HBSurface(this, &m_shader, &b_shader, npx, npy, 16);
	//Center the surface.
	W_trans_center = glm::translate(W_trans_center, vec3(- (surface->ncpx-1.0f)/2.0f, 0.0f, - (surface->ncpy-1.0f)/2.0f));
	W_trans = glm::translate(W_trans, vec3(0.0f, 0.0f, - 8.0f));
	update_W();
	// Set up initial view and projection matrices (need to do this here,
	// since it depends on the GLFW window being set up correctly).
	view = glm::lookAt(
		glm::vec3( 0.0f, 0.0f, 0.0f ),
		glm::vec3( 0.0f, 0.0f, -1.0f ),
		glm::vec3( 0.0f, 1.0f, 0.0f ) );

	proj = glm::perspective(
		glm::radians( 45.0f ),
		float( m_framebufferWidth ) / float( m_framebufferHeight ),
		0.1f, 10000.0f );
}

void A1::updateLighting(){
	// Lighting
	m_shader.enable();
		glm::vec3 lpos(0.0f, 10.0f, 10.0f);
		glUniform3fv(light_position, 1, value_ptr(lpos));
		glm::vec3 lcol(0.8f, 0.8f, 0.8f);
		glUniform3fv(light_colour, 1, value_ptr(lcol));
		glm::vec3 lamb(0.1f);
		glUniform3fv(light_ambient, 1, value_ptr(lamb));
		CHECK_GL_ERRORS;
	m_shader.disable();
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
		if (level >= 0){
			ImGui::Text( "Control Point Level: %d", level);
		} else {
			ImGui::Text( "Control Point Level: All");
		}

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

mat4 A1::update_W(){
	//std::cout << "Updating W" << std::endl;
	//std::cout << glm::to_string(W) << std::endl;
	//std::cout << glm::to_string(W_rot) << std::endl;
	W_scale = glm::scale(glm::mat4(), vec3(zoom, zoom, zoom));
	//std::cout << glm::to_string(W) << std::endl;
	W = W_trans*W_rot*W_scale*W_trans_center;

	return W;
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after guiLogic().
 */
void A1::draw()
{
	//std::cout << "New Render Cycle" << std::endl;
	glfwSwapInterval(0);
	// Create a global transformation for the model (centre it).
	//W = glm::translate( W, vec3( 0.0f, 0.0f, 0.0f ) );
	updateLighting();

	surface->render_points(W, proj, view, do_picking, level);
	surface->render_surface(W, proj, view, do_picking);

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

glm::vec3 A1::GetOGLPos(float x, float y, float depth)
{
	//std::cout << x << " " << y << std::endl;
    GLint viewport[4];
	glm::vec4 port(viewport[0], viewport[1], viewport[2], viewport[3]);
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;

	glGetIntegerv( GL_VIEWPORT, viewport );
	//std::cout << viewport[0] << " " << viewport[1] << " " << viewport[2] << " " << viewport[3] << std::endl;
    winX = (float)x;
    winY = (float)viewport[3] - (float)y;
	winZ = depth;
	glm::vec3 win(winX, winY, winZ);
	glm::vec3 pos;
	glm::mat4 modelview = view*W;

    pos = glm::unProject(win, modelview, proj, port);
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;

    return pos;
}

float A1::getDepth(float x, float y){
	GLint viewport[4];
	GLfloat winX, winY, winZ;

	glGetIntegerv( GL_VIEWPORT, viewport );
	winX = (float)x;
	winY = (float)viewport[3] - (float)y;

	glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

	return winZ;
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
		if (dragging && !drag_cp){
			rotate_surface(xPos, yPos);
			//rot_rads = old_rot_rads - 2*3.1415*(old_x - xPos)/m_framebufferWidth;
		} else if (dragging && drag_cp == true){
			//std::cout << xPos << " " << yPos << " " << old_x << " " << old_y << std::endl;
			if (!depth_set){
				depth_set = true;
				depth_val = getDepth(xPos, yPos);
				std::cout << "Depth Read" << std::endl;
			}
			GetOGLPos(old_x, old_y, depth_val);
			glm::vec3 delta_model = GetOGLPos(xPos, yPos, depth_val);// - GetOGLPos(old_x, old_y);

			//std::cout << delta_model.x << " " << delta_model.y << " " << delta_model.z << std::endl;
			surface->move_selected_cp(delta_model);
			old_x = xPos;
			old_y = yPos;
		} else {
		// Track the position of the mouse in preparation for rotation.
			old_x = xPos;
			old_y = yPos;
			W_rot_old = W_rot;
		}
	}
	return eventHandled;
}

int A1::pick_object(){
	do_picking = true;
	double xpos, ypos;
	glfwGetCursorPos(m_window, &xpos, &ypos);
	glClearColor(1.0, 1.0, 1.0, 1.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glClearColor(0.35, 0.35, 0.35, 1.0);

	draw();

	CHECK_GL_ERRORS;

	// Ugly -- FB coordinates might be different than Window coordinates
	// (e.g., on a retina display).  Must compensate.
	xpos *= double(m_framebufferWidth) / double(m_windowWidth);
	// WTF, don't know why I have to measure y relative to the bottom of
	// the window in this case.
	ypos = m_windowHeight - ypos;
	ypos *= double(m_framebufferHeight) / double(m_windowHeight);

	GLubyte buffer[ 4 ] = { 0, 0, 0, 0 };
	// A bit ugly -- don't want to swap the just-drawn false colours
	// to the screen, so read from the back buffer.
	glReadBuffer( GL_BACK );
	// Actually read the pixel at the mouse location.
	glReadPixels( int(xpos), int(ypos), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, buffer );
	CHECK_GL_ERRORS;

	// Reassemble the object ID.
	unsigned int what = buffer[0] + (buffer[1] << 8) + (buffer[2] << 16);

	std::cout << "Selected idx " << what << std::endl;

	do_picking = false;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	return what;
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

			int what = pick_object();
			surface->select_cp(what, lctrl, level);

			if (what == surface->get_selected_cp_idx()){
				drag_cp = true;
			} else {
				drag_cp = false;
			}
		}
		if (button == GLFW_MOUSE_BUTTON_LEFT && actions == GLFW_RELEASE){
			dragging = false;
			depth_set = false;
			rot_rads = fmod(rot_rads, 2*3.1415);
			old_rot_rads = rot_rads;
			drag_x = 0;
			drag_y = 0;
		}
	}

	return eventHandled;
}

void A1::rotate_surface(float xPos, float yPos){
	//std::cout << "Calling rotate surface" << std::endl;
	xPos = 2.0f*xPos/m_framebufferWidth -1.0f;
	yPos = -(2.0f*yPos/m_framebufferHeight -1.0f);
	float xPos_old = 2.0f*old_x/m_framebufferWidth -1.0f;
	float yPos_old = -(2.0f*old_y/m_framebufferHeight -1.0f);
 	float diameter = 2.2f;

	//std::cout << xPos << " " << yPos <<

	float fVecX, fVecY, fVecZ;
	float* pfVecX = &fVecX;
	float* pfVecY = &fVecY;
	float* pfVecZ = &fVecZ;

	vCalcRotVec(xPos, yPos, xPos_old, yPos_old, diameter, pfVecX, pfVecY, pfVecZ);
	//std::cout << glm::to_string(vAxisRotMatrix(fVecX, fVecY, fVecZ)) << std::endl;
	W_rot = vAxisRotMatrix(fVecX, fVecY, fVecZ) * W_rot_old;
	update_W();
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

	update_W();

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

		if (key == GLFW_KEY_S){
			surface->split_selected_cp();
			eventHandled = true;
		}

		if (key == GLFW_KEY_LEFT_CONTROL){
			lctrl = true;
			eventHandled = true;
		}

		if (key == GLFW_KEY_P){
			glm::vec3 coords = surface->get_selected_cp_coords();
			std::cout << glm::to_string(coords) << std::endl;
			eventHandled = true;
		}

		if (key == GLFW_KEY_PAGE_UP){
			level++;
			eventHandled = true;
		}

		if (key == GLFW_KEY_PAGE_DOWN){
			if (level >= 0){
				level--;
			}
			eventHandled = true;
		}

		if (key == GLFW_KEY_E){
			edit_cp = !edit_cp;
			surface->set_mode_edit_cp(edit_cp);
			std::cout << "Setting edit points mode to " << edit_cp << std::endl;

			eventHandled = true;
		}
	}

	if (action == GLFW_RELEASE){
		if (key == GLFW_KEY_LEFT_CONTROL){
			lctrl = false;
			eventHandled = false;
		}
	}

	return eventHandled;
}
