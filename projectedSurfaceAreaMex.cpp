/*
 * projectedSurfaceAreaMex.cpp
 * Usage: s = projectedSurfaceAreaMex(angles, TRIV, X, Y, Z)
 * Computes the projected surface area of s after a rotation.
 * Variables:
 * s - projected surface area.
 * angles - input x, y, z rotation angles.
 * TRIV - triangles of mesh.
 * X - X coordinates of vertices.
 * Y - Y coordinates of vertices.
 * Z - Z coordinates of vertices.
 *
 * David Pickup 2013
 */

#include "mex.h"
#include <math.h>

// Define how to access arrays with two dimentions in matlab.
#define POS(x,y,M) (x + (y*M))

/* Multiply vertices by rotation matrix. */
void matrixMult(double *R, double *X, double *Y, double *Z, double *newX, double *newY, double *newZ, int nVerts)
{
    // Iterate through all vertes.
    for (int i = 0; i < nVerts; i++)
    {
        newX[i] = (R[POS(0,0,3)]*X[i]) + (R[POS(0,1,3)]*Y[i]) + (R[POS(0,2,3)]*Z[i]);
        newY[i] = (R[POS(1,0,3)]*X[i]) + (R[POS(1,1,3)]*Y[i]) + (R[POS(1,2,3)]*Z[i]);
        newZ[i] = (R[POS(2,0,3)]*X[i]) + (R[POS(2,1,3)]*Y[i]) + (R[POS(2,2,3)]*Z[i]);
    }
}

/* Calculate the edge length between two 2D vertices. */
double edgeLength(double X1, double X2, double Y1, double Y2)
{
    return sqrt(((Y1-X1)*(Y1-X1)) + ((Y2-X2)*(Y2-X2)));
}

/* Calculate the area of a 2D polygon from its vertices.  */
double polyArea(double *X, double *Y, int nVerts)
{
    // Initialise area.
    double area = 0;
    // The last vertex is the previous one of the first.
    int j = nVerts - 1;
    
    // Iterate through all vertices.
    for (int i = 0; i < nVerts; i++)
    {
        // Accumulate the area.
        area += (X[j]+X[i]) * (Y[j]-Y[i]);
        
        // Set the previous vertex to the current.
        j = i;
    }
    
    // Return the final area.
    return area / 2;
}

/* Calculate the area of a 2D triangle using Heron's formula. */
double heron(double a, double b, double c)
{
    return ((1/4) * sqrt((((a*a)+(b*b)+(c*c)) * ((a*a)+(b*b)+(c*c))) - (2*((a*a*a*a)+(b*b*b*b)+(c*c*c*c)))));
}

/* The matlab gateway function. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Varify that the correct number of function arguments have been given.
    if (nrhs != 5)
        mexErrMsgTxt("Five input arguments required.");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required.");
    
    // Initialise variables.
    double *TRIV;       // Mesh triangles.
    double *X;          // X coordinates of vertices.
    double *Y;          // Y coordinates of vertices.
    double *Z;          // Z coordinates of vertices.
    double *newX;       // X coordinates of rotated vertices.
    double *newY;       // Y coordinates of rotated vertices.
    double *newZ;       // Z coordinates of rotated vertices.
    double *tmpX;       // temporary X coordinates used during rotation.
    double *tmpY;       // temporary Y coordinates used during rotation.
    double *tmpZ;       // temporary Z coordinates used during rotation.
    double triX[3];     // X coordinates of single triangle.
    double triY[3];     // Y coordinates of single triangle.
    double triZ[3];     // Z coordinates of single triangle.
    double *angles;     // Rotation angles.
    int m,n;            // Size of TRIV.
    double a,b,c;       // Edge lengths.
    
    // Get rotation angles.
    angles = mxGetPr(prhs[0]);
    //mexPrintf("Angles loaded.\n");
    
    // Get triangles and vertex coordinates.
    TRIV = mxGetPr(prhs[1]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    X = mxGetPr(prhs[2]);
    Y = mxGetPr(prhs[3]);
    Z = mxGetPr(prhs[4]);
    int nVerts = mxGetM(prhs[2]);
    newX = (double*) calloc(nVerts,sizeof(double));
    newY = (double*) calloc(nVerts,sizeof(double));
    newZ = (double*) calloc(nVerts,sizeof(double));
    tmpX = (double*) calloc(nVerts,sizeof(double));
    tmpY = (double*) calloc(nVerts,sizeof(double));
    tmpZ = (double*) calloc(nVerts,sizeof(double));
    //mexPrintf("Mesh loaded.\n");
    
    // Create rotation matrices.
    double Rx[] = {1, 0, 0, 0, cos(angles[0]), sin(angles[0]), 0, -sin(angles[0]), cos(angles[0])};
    double Ry[] = {cos(angles[1]), 0, -sin(angles[1]), 0, 1, 0, sin(angles[1]), 0, cos(angles[1])};
    double Rz[] = {cos(angles[2]), sin(angles[2]), 0, -sin(angles[2]), cos(angles[2]), 0, 0, 0, 1};
    //mexPrintf("Created rotation vertices.\n");    

    // Rotate vertices.
    matrixMult(Rx, X, Y, Z, newX, newY, newZ, mxGetM(prhs[2]));
    matrixMult(Ry, newX, newY, newZ, tmpX, tmpY, tmpZ, mxGetM(prhs[2]));
    matrixMult(Rz, tmpX, tmpY, tmpZ, newX, newY, newZ, mxGetM(prhs[2]));
    //mexPrintf("Rotated vertices.\n");
    
    // Initialise output value.
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *r = mxGetPr(plhs[0]);
    r[0] = 0;
    //mexPrintf("Initialised output variable: %f.\n", r[0]);
    
    // Iterate through all triangles.
    for (int i = 0; i < m; i++)
    {
        // Get vertices for current triangle.
        for (int j = 0; j < 3; j++)
        {
            triX[j] = newX[(int)TRIV[POS(i,j,m)]-1];
            triY[j] = newY[(int)TRIV[POS(i,j,m)]-1];
            triZ[j] = newZ[(int)TRIV[POS(i,j,m)]-1];
        }
        //mexPrintf("Got vertices for current triangle.\n");
        
        //mexPrintf("Triangle %d of %d: (%f, %f, %f) (%f, %f, %f) (%f, %f, %f)\n", i, m,
        //          triX[0], triY[0], triZ[0],
        //          triX[1], triY[1], triZ[1],
        //          triX[2], triY[2], triZ[2]);
        
        // Calculate projected area of current triangle and add to total.
        a = polyArea(triX, triY, 3);
        if (a < 0)
            a = a * -1;
        b = polyArea(triX, triZ, 3);
        if (b < 0)
            b = b * -1;
        c = polyArea(triY, triZ, 3);
        if (c < 0)
            c = c * -1;
        r[0] += a + b + c;
    }
    
    free(newX);
    free(newY);
    free(newZ);
    free(tmpX);
    free(tmpY);
    free(tmpZ);
}
