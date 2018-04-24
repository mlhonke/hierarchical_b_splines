/*
 * Module : events.c
 *
 * Author : Greg Veres
 *
 * Date   : April 10, 1995
 *
 * Purpose: Controls event processing and scene drawing.
 *
 */
#include <stdlib.h>
#include "events.hpp"
#include "trackball.hpp"
#include <X11/Xutil.h>

/* 
 * Local Global Variables
 */
static Matrix mRotations    = {{1.0, 0.0, 0.0, 0.0},
                               {0.0, 0.91, 0.42, 0.0},
                               {0.0, -0.42, 0.91, 0.0},
                               {0.0, 0.0, 0.0, 1.0}};
static Matrix mTranslations = {{1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, -25.0, 1.0}};
Matrix mIdentity     = {{1.0, 0.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0, 0.0},
                        {0.0, 0.0, 1.0, 0.0},
                        {0.0, 0.0, 0.0, 1.0}};

static GLfloat fYellowVec[] = {1.0, 1.0, 0.0};
static GLfloat fWhiteVec[]  = {1.0, 1.0, 1.0};
static GLfloat fBlackVec[]  = {0.0, 0.0, 0.0};
static GLfloat fGreenVec[]  = {0.0, 1.0, 0.0};
static GLfloat fBlueVec[]   = {0.0, 0.0, 1.0};
static GLfloat fGrayVec[]   = {0.5, 0.6, 0.5};
static GLfloat fRedVec[]    = {1.0, 0.0, 0.0};

short sXReference, sYReference;

int nCurrentDir = DIR_NONE;

void vTransposeMatrix(Matrix mSrcDst) {
    GLdouble temp;
    int i,j;

    // Transpose matrix
    for ( i=0; i<4; ++i ) {
        for ( j=i+1; j<4; ++j ) {
            temp = mSrcDst[i][j];
            mSrcDst[i][j] = mSrcDst[j][i];
            mSrcDst[j][i] = temp;
        }
    }
}

/*
 * Name      : void vCopyMatrix(Matrix mSource, Matrix mDestination)
 *
 * Parameters: Matrix mSource      - The source matrix.
 *             Matrix mDestination - The destination matrix.
 *
 * Returns   : void
 *
 * Purpose   : Copies matrix mSource to matrix mDestination.
 *             the result in mDestination.
 */
void vCopyMatrix(Matrix mSource, Matrix mDestination) 
{
    int i, j;

    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++) {
            mDestination[i][j] = mSource[i][j];
        }
    }
}

/*
 * Name      : void vRightMultiply(Matrix mMat1, Matrix mMat2)
 *
 * Parameters: Matrix mMat1 - The first and destination matrix.
 *             Matrix mMat2 - The second matrix.
 *
 * Returns   : void
 *
 * Purpose   : Right multiplies matrix mMat1 by matrix mMat2 and stores
 *             the result in mMat1.
 */
void vRightMultiply(Matrix mMat1, Matrix mMat2) 
{
    int    i, j;
    Matrix mMat3;

    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++) {
            mMat3[i][j] = mMat1[i][0]*mMat2[0][j] + mMat1[i][1]*mMat2[1][j] +
                mMat1[i][2]*mMat2[2][j] + mMat1[i][3]*mMat2[3][j];
        }
    }
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++) {
            mMat1[i][j] = mMat3[i][j];
        }
    }
}

/*
 * Name      : void vTranslate(float fTrans, char cAxis, Matrix mMat)
 *
 * Parameters: float  fAngle - The distance of translation.
 *             char   cAxis  - The axis of rotation.
 *             Matrix mMat   - The matrix to store the result in.
 *
 * Returns   : void
 *
 * Purpose   : Computes the translation along the given axis and stores
 *             the result in mMat.
 */
void vTranslate(float fTrans, char cAxis, Matrix mMat)
{
    vCopyMatrix(mIdentity, mMat);
    switch(cAxis) {
    case 'x':
        mMat[3][0] = fTrans;
        break;

    case 'y':
        mMat[3][1] = fTrans;
        break;

    case 'z':
        mMat[3][2] = fTrans;
        break;
    }
}
