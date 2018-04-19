#include <algorithm>
#include <iostream>

#include "cube.hpp"

using namespace std;

Cube::Cube(float x, float y, float z, int colour):x(x), y(y), z(z), colour(colour){
	/* Cube sides are labelled in the default view. Triangles are
	 * named with respect to the position of their right angle in the
	 * face of the cube
	 */

	//Floor upper-left
	vertices[0] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[1] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[2] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);

	//Floor bottom-right
	vertices[3] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[4] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[5] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);

	//Back bottom-left
	vertices[6] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[7] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[8] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);

	//Back upper-right
	vertices[9] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[10] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[11] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);

	//Left-side bottom-right
	vertices[12] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[13] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[14] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);

	//Left-side upper-left
	vertices[15] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[16] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[17] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);

	//Right-side bottom-right
	vertices[18] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + -0.5f/scale);
	vertices[19] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[20] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);

	//Right-side upper-left
	vertices[21] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[22] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[23] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);

	//Front-side bottom-left
	vertices[24] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);
	vertices[25] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[26] = vec3(x + -0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);

	//Front-side upper-right
	vertices[27] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[28] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[29] = vec3(x + 0.5f/scale, y + -0.5f/scale, z + 0.5f/scale);

	//Top-side upper-left
	vertices[30] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[31] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[32] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);

	//Top-side bottom-right
	vertices[33] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
	vertices[34] = vec3(x + 0.5f/scale, y + 0.5f/scale, z + -0.5f/scale);
	vertices[35] = vec3(x + -0.5f/scale, y + 0.5f/scale, z + 0.5f/scale);
}

Cube::Cube(float x, float y, float z): Cube(x, y, z, default_colour)
{
}

vec3* Cube::get_vertices()
{
	return vertices;
}

bool Cube::draw()
{
	return draw_cube;
}

void Cube::print_properties(){
	cout << x << y << z << endl;
}

bool Cube::operator==(const Cube& test_cube){
	if (x == test_cube.get_x()
	&&  y == test_cube.get_y()
	&&  z == test_cube.get_z())
		return true;
}

/*
 * Use a pointer to allow for easier modification of colours for
 * groups of cubes.
 */
void Cube::set_colour(int new_colour){
	colour = new_colour;
}
int Cube::get_colour(){
	return colour;
}

int Cube::get_x() const{
	return x;
}
int Cube::get_y() const{
	return y;
}
int Cube::get_z() const{
	return z;
}

Cube::Cube(){}

Cube::~Cube()
{
}
