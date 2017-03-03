/*
 * clockmatchModels.cpp
 * Compare all models in database using clock matching. The histograms are
 * compared using DMaxHist measure.
 * Usage: D = clockmatchModels(Hists, permutations, vertexTable, pEdgeTables)
 * Variables:
 * D - distance matrix.
 * Hists - feature histograms for each view of each model.
 * permutations - axis permutations for clockmatching.
 * vertexTable - vertex table of viewpoints.
 * pEdgeTables - edge tables of viewpoints for each permutation.
 *
 * David Pickup 2013
 */

#include "mex.h"
#include <algorithm>

// Define how to access arrays with two dimentions in matlab.
#define POS(x,y,M) ((int)(x + (y*M)))

/* Compares two histogram represented as sparse matrices. */
double DMaxHist(mwSize nzmax1, double *pr1, mwIndex *ir1,
        mwSize nzmax2, double *pr2, mwIndex *ir2)
{
    double numerator, denominator, sum1, sum2;
    int i,j;
    
    // Calculate numerator.
    numerator = 0;
    for (i = 0; i < nzmax1; i++)
        for (j = 0; j < nzmax2; j++)
            if (ir1[i] == ir2[j])
                numerator += std::min(pr1[i],pr2[j]);
    
    // Calculate denominator.
    sum1 = 0;
    for (i = 0; i < nzmax1; i++)
        sum1 += pr1[i];
    
    sum2 = 0;
    for (i = 0; i < nzmax2; i++)
        sum2 += pr2[i];
    
    denominator = std::max(sum1,sum2);
    
//     // Check for zeros.
//     if ((numerator == 0) || (denominator == 0))
//         return 0;
    
    // Return final distance.
    return 1 - (numerator / denominator);
}

/* The matlab gateway function. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Varify that the correct number of function arguments have been given.
    if (nrhs != 4)
        mexErrMsgTxt("Four input arguments required.");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required.");
    
    // Initialise variables.
    char buffer[35];                            // For printing
    const mxArray *Hists, *H1, *H2, *Hs, *Hp;   // histograms.
    const mxArray *perms;                       // permutations.
    const mxArray *vTs;                         // vertex tables.
    const mxArray *pET;                         // edge permutations table.
    double *D;                                  // distance matrix.
    double *edges, *edges1;                     // edges in edge table.
    double *vT, v1, v2;                         // vertex table.
    double *p;                                  // permutation.
    int nModels;                                // number of models.
    int nPerms;                                 // number of permutations.
    mxArray *eT, *eT1;
    double pD, d;
    int i, j, k, l, m, nET, nE, nV, idx1, idx2;
    
    // Get triangles and vertex coordinates.
    Hists = prhs[0];
    perms = prhs[1];
    vTs = prhs[2];
    pET = prhs[3];
    
    // Get number of models/histograms, and permutations.
    nModels = mxGetNumberOfElements(Hists);
    nPerms = mxGetNumberOfElements(perms);
    
    // Initialise output.
    plhs[0] = mxCreateDoubleMatrix(nModels, nModels, mxREAL);
    D = mxGetPr(plhs[0]);
    
    // Get first edge table.
    eT1 = mxGetCell(pET,0);
    
    // Iterate though all pairs of models.
    for (i = 0; i < nModels; i++)
    {
        for(j = i+1; j < nModels; j++)
        {
            // Initialise distance to high value.
            d = 9999;
            
            // Iterate through all axis permutations.
            for (k = 0; k < nPerms; k++)
            {
                // Initialise distance for current permutation.
                pD = 0;

                // Get current permutation and associated edge tables.
                p = mxGetPr(mxGetCell(perms,k));
                eT = mxGetCell(pET,k);
                
                // Get number of edge tables at current permutation.
                nET = mxGetNumberOfElements(eT);
                
                // Get Histograms for models.
                H1 = mxGetCell(Hists,i);
                H2 = mxGetCell(Hists,j);

                // Compute distance for first size vertices.
                for (l = 0; l < 6; l++)
                {
                    // Get stadard and permuted histogram.
                    Hs = mxGetCell(H1,l);
                    Hp = mxGetCell(H2,p[l]-1);
                    pD += DMaxHist(mxGetNzmax(Hs), mxGetPr(Hs),
                            mxGetIr(Hs), mxGetNzmax(Hp), mxGetPr(Hp),
                            mxGetIr(Hp));
                }
                
                // Iterate though all edge tables for current permutation.
                for (l = 0; l < nET-1; l++)
                {
                    // Get edges in current and first edge table.
                    edges = mxGetPr(mxGetCell(eT,l));
                    edges1 = mxGetPr(mxGetCell(eT1,l));
                    
                    // Get current vertex table.
                    vT = mxGetPr(mxGetCell(vTs, l));
                    nV = mxGetM(mxGetCell(vTs,l));
                    
                    // Get number of edges.
                    nE = mxGetM(mxGetCell(eT,l));
                    
                    // Iterate through all edges.
                    for (m = 0; m < nE; m++)
                    {
                        /* Get vertex subdivision for standard and permuted edge tables.*/
                        
                        idx1 = edges1[POS(m,0,nE)]-1;
                        idx2 = edges1[POS(m,1,nE)]-1;
                        
                        v1 = vT[POS(idx1, idx2, nV)];
                        v1 = v1 - 1;
                        
                        idx1 = edges[POS(m,0,nE)]-1;
                        idx2 = edges[POS(m,1,nE)]-1;
                        
                        v2 = vT[POS(idx1, idx2, nV)];
                        v2 = v2 - 1;
                        
                        // Compute distance between histograms.
                        Hs = mxGetCell(H1,v1);
                        Hp = mxGetCell(H2,v2);
                        pD += DMaxHist(mxGetNzmax(Hs), mxGetPr(Hs),
                                mxGetIr(Hs), mxGetNzmax(Hp), mxGetPr(Hp),
                                mxGetIr(Hp));
                    }
                }
                // If permutation distance is lower than stored distance,
                // then store permutation distance.
                if (pD < d)
                        d = pD;
            }
            D[POS(i,j,nModels)] = d;
            D[POS(j,i,nModels)] = d;
        }
        sprintf(buffer,"Finished matching for model %d.\n",i);
        mexPrintf(buffer);
    }
}