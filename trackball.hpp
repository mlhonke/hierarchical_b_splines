#ifndef TRACKBALL_H
#define TRACKBALL_H

#include "events.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>

/* Function prototypes */
void vCalcRotVec(float fNewX, float fNewY,
                 float fOldX, float fOldY,
                 float fDiameter,
                 float *fVecX, float *fVecY, float *fVecZ);
glm::mat4 vAxisRotMatrix(float fVecX, float fVecY, float fVecZ);

#endif
