#pragma once

#include <glm/glm.hpp>

using namespace glm;

class Cube
{
public:
	Cube(float x, float y, float z);
	Cube(float x, float y, float z, int colour);
	Cube();
	bool operator==(const Cube& test_cube);

	int get_x() const;
	int get_y() const;
	int get_z() const;

	void set_colour(int new_colour);
	int get_colour();

	~Cube();

	vec3* get_vertices();
	void print_properties();
	bool draw();

private:
	float x;
	float y;
	float z;
	float scale = 10.0f;
	bool draw_cube;
	int default_colour = 0;
	int colour;
	vec3 vertices[36];
};
