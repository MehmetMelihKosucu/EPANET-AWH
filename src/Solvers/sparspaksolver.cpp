/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "sparspaksolver.h"
#include "sparspak.h"

#include <cstring>
#include <limits>
#include <iostream>
#include <ctime>
#include <vector>
#include "network.h"
using namespace std;

// Local module-level functions
//-----------------------------------------------------------------------------
int  compress(
        int n, int nnz, int* xrow, int* xcol, int* xadj, int* adjncy, int* xaij);
void buildAdjncy(
        int n, int nnz, int* xrow, int* xcol, int* xadj,
        int* adjncy, int* adjncy2, int* xaij, int* nz);
void sortAdjncy(
        int n, int nnz, int* xadj, int* adjncy, int* xadj2, int* adjncy2, int* nz);
void transpose(
        int n, int* xadj1, int* adjncy1, int* xadj2, int* adjncy2, int* nz);
int  reorder(
        int n, int* xadj, int* adjncy, int* perm, int* invp, int& nnzl);
int  factorize(
        int n, int& nnzl, int* xadj, int* adjncy, int* perm,
        int* invp, int* xlnz, int* xnzsub, int* nzsub);
void aij2lnz(
        int nnz, int* xrow, int* xcol, int* invp, int* xlnz, int* xnzsub,
        int* nzsub, int* xaij);

//-----------------------------------------------------------------------------

SparspakSolver::SparspakSolver(ostream& logger) :
    nrows(0), nnz(0), nnzl(0), perm(0), invp(0), xlnz(0), xnzsub(0),
    nzsub(0), xaij(0), link(0), first(0), lnz(0), diag(0), rhs(0), temp(0),
    msgLog(logger)
{}

//-----------------------------------------------------------------------------

SparspakSolver::~SparspakSolver()
{
    delete [] perm;
    delete [] invp;
    delete [] xlnz;
    delete [] xnzsub;
    delete [] nzsub;
    delete [] xaij;
    delete [] link;
    delete [] first;
    delete [] lnz;
    delete [] diag;
    delete [] rhs;
    delete [] temp;
}

//-----------------------------------------------------------------------------

int SparspakSolver::init(int nrows_, int nnz_, int* xrow, int* xcol)
{
    // ... save number of equations and number of off-diagonal coeffs.
    nrows = nrows_;
    nnz = nnz_;

    // ... allocate space for pointers from Aij to lnz
    xaij = new int[nnz];
    if ( !xaij ) return 0;
    memset(xaij, 0, nnz*sizeof(int));

    // ... allocate space for row re-ordering
    perm = new int[nrows];
    invp = new int[nrows];
    if ( !perm || !invp ) return 0;


    // ... compress, re-order, and factorize coeff. matrix A
    int* xadj;
    int* adjncy;
    int flag = 0;
    for (;;)
    {
        // ... allocate space for adjacency lists
        xadj = new int[nrows+1];
        adjncy = new int[2*nnz];
        if ( !xadj || !adjncy ) break;

        // ... store matrix A in compressed format
        if ( !compress(nrows, nnz, xrow, xcol, xadj, adjncy, xaij) ) break;

        // ... re-order the rows of A to minimize fill-in
        //clock_t startTime = clock();
        if ( !reorder(nrows, xadj, adjncy, perm, invp, nnzl) ) break;

/************ DEBUG  ******************
    cout << "\n nnzl = " << nnzl;

    for (int i = 0; i < nrows; i++)
    {
        cout << "\n i = " << i << "  perm[i] = " << perm[i] << "  invp[i] = " << invp[i];
    }

//*****************************************/

        // ... allocate space for compressed storage of factorized matrix
        xlnz = new int[nrows+1];
        xnzsub = new int[nrows+1];
        nzsub = new int[nnzl];
        if ( !xlnz || !xnzsub || !nzsub ) break;

        // ... symbolically factorize A to produce L
        if ( !factorize(nrows, nnzl, xadj, adjncy, perm, invp, xlnz,
                        xnzsub, nzsub) ) break;

//************  DEBUG  ********************
        // ... report factorization results

        int nnz0 = xadj[nrows] / 2;
        clock_t startTime = clock();
        double procTime = (double)(clock() - startTime) / //
                          (double)CLOCKS_PER_SEC * 1000.0;

        msgLog << endl;
        msgLog << "  Hydraulic Solution Matrix:" << endl;
        msgLog << "  Number of rows          " << nrows << endl;
        msgLog << "  Off-diagonal non-zeros  " << nnz << endl;
        msgLog << "  Duplicate non-zeros     " << nnz - nnz0 << endl;
        msgLog << "  Amount of fill-in       " << nnzl - nnz0 << endl;
        msgLog << "  Processing time (msec)  " << procTime << endl;
//********************************************/

        // ... all steps were successful
        flag = 1;
        break;
    }

    // ... free memory used for adjacency lists
    delete [] xadj;
    delete [] adjncy;

    // ... return if error condition
    if ( !flag ) return flag;

    // ... map off-diag coeffs. of A to positions in xlnz
    aij2lnz(nnz, xrow, xcol, invp, xlnz, xnzsub, nzsub, xaij);

    // ... allocate space for coeffs. of L and r.h.s vector
    lnz = new double[nnzl];
    diag = new double[nrows];
    rhs = new double[nrows];
    if ( !lnz || ! diag || !rhs ) return 0;

    // ... allocate space for work arrays used by the solve() method
    temp = new double[nrows];
    first = new int[nrows];
    link = new int[nrows];
    if ( !temp || !first || !link ) return 0;
    return 1;
}

//-----------------------------------------------------------------------------

int SparspakSolver::solve(int n, double x[])
{
    // ... call sp_numfct to numerically evaluate the factorized matrix L

/*//*********  DEBUG  ****************************
    --diag;  --rhs; --invp;
    cout << "\n Before call to numfct:";
    for (int i = 1; i <= nrows; i++)
    {
        int j = invp[i] - 1;
        cout << "\n diag[" << j << "] = " << diag[i] << ",  rhs[" << j << "] = " << rhs[i];
    }
    ++diag;  ++rhs;  ++invp;
//*********************************************/

    int flag;
    sp_numfct(nrows, xlnz, lnz, xnzsub, nzsub, diag, link, first, temp, flag);

    // if the matrix was ill-conditioned, return the problematic row
    if ( flag )
    {
        --invp;
        flag = invp[flag] - 1;
        ++invp;
        return flag;
    }

    // call sp_solve() to solve the system LDL'x = b
    sp_solve(nrows, xlnz, lnz, xnzsub, nzsub, diag, rhs);

    // transfer results from rhs to x (recognizing that rhs
    // arrays are offset by 1)
    --x; --rhs; --invp;
    for (int i = 1; i <= nrows; i++)
    {
        x[i] = rhs[invp[i]];
    }
    ++x; ++rhs; ++invp;
    return -1;
}

void SparspakSolver::reset()
{
    memset(diag, 0, (nrows) * sizeof(double));
    memset(lnz, 0, (nnzl) * sizeof(double));
    memset(rhs, 0, (nrows) * sizeof(double));
}

//-----------------------------------------------------------------------------

/*void SparspakSolver::reset()
{
    // First, verify that our arrays are properly allocated
    if (!diag || !lnz || !rhs || nrows <= 0 || nnzl <= 0) {
        std::cerr << "Error: Arrays not properly initialized in reset()\n";
        std::cerr << "diag: " << (diag ? "allocated" : "null") << "\n";
        std::cerr << "lnz: " << (lnz ? "allocated" : "null") << "\n";
        std::cerr << "rhs: " << (rhs ? "allocated" : "null") << "\n";
        std::cerr << "nrows: " << nrows << "\n";
        std::cerr << "nnzl: " << nnzl << "\n";
        return;
    }

    try {
        // Use secure zeroing operations
        std::fill(diag, diag + nrows, 0.0);
        std::fill(lnz, lnz + nnzl, 0.0);
        std::fill(rhs, rhs + nrows, 0.0);

        // Optional: Add verification
#ifdef _DEBUG
        for (int i = 0; i < nnzl; i++) {
            if (lnz[i] != 0.0) {
                std::cerr << "Warning: lnz[" << i << "] not properly zeroed\n";
            }
        }
#endif
    }
    catch (const std::exception& e) {
        std::cerr << "Error during reset: " << e.what() << "\n";
        std::cerr << "Memory locations:\n";
        std::cerr << "diag: " << static_cast<void*>(diag) << "\n";
        std::cerr << "lnz: " << static_cast<void*>(lnz) << "\n";
        std::cerr << "rhs: " << static_cast<void*>(rhs) << "\n";
    }
} // */

//-----------------------------------------------------------------------------


double SparspakSolver::getDiag(int i)
{
    int k = invp[i] - 1;
    return diag[k];
}

//-----------------------------------------------------------------------------

double SparspakSolver::getOffDiag(int i)
{
    int k = xaij[i] - 1;
    return lnz[k];
}

//-----------------------------------------------------------------------------

double SparspakSolver::getRhs(int i)
{
    int k = invp[i] - 1;
    return rhs[k];
}

//-----------------------------------------------------------------------------

void SparspakSolver::setDiag(int i, double value)
{
    int k = invp[i] - 1;
    diag[k] = value;
}

//-----------------------------------------------------------------------------

void SparspakSolver::setRhs(int i, double value)
{
    int k = invp[i] - 1;
    rhs[k] = value;
}

//-----------------------------------------------------------------------------

void SparspakSolver::addToDiag(int i, double value)
{
    int k = invp[i] - 1;
    diag[k] += value;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Adds a value to an off-diagonal coefficient of A.
//
// \param[in]  j  Index of off-diagonal coefficient in A (0 <= j < nnz).
// \param[in]  value  Value to add to the coefficient at position j in A.
//
// The coefficient is stored in the permuted matrix L.
void SparspakSolver::addToOffDiag(int j, double value)
{
    int k = xaij[j] - 1;
    
    // Apply conditioning to the value before adding it
    double conditionedValue = value;
    
    // 1. Apply absolute value thresholding for large values
    const double MAX_COEFFICIENT = 1e6;  // Set based on your system characteristics
    if (std::abs(conditionedValue) > MAX_COEFFICIENT) {
        conditionedValue = MAX_COEFFICIENT * (conditionedValue > 0 ? 1.0 : -1.0);
    }
    
    // 2. Filter out extremely small values that contribute to ill-conditioning
    const double MIN_COEFFICIENT = 1e-10;  // Adjust based on your numerical precision needs
    if (std::abs(conditionedValue) < MIN_COEFFICIENT) {
        conditionedValue = 0.0;
    }
    
    // 3. Apply symmetric scaling to maintain matrix balance
    // This is particularly important for hydraulic network matrices
    if (lnz[k] != 0.0 && (lnz[k] * conditionedValue) < 0) {
        // If changing sign, be more conservative
        conditionedValue *= 0.5;
    }
    
    // Add the conditioned value to the matrix
    lnz[k] += conditionedValue;
}

//-----------------------------------------------------------------------------

void SparspakSolver::addToRhs(int i, double value)
{
    int k = invp[i] - 1;
    rhs[k] += value;
}

bool SparspakSolver::checkDimensions(int expectedSize) {
    // For raw pointer arrays, we use nrows as our size tracker
    // This checks if our current arrays match the expected size
    
    // Basic validation first
    if (expectedSize <= 0) {
        std::cerr << "Error: Invalid expected size (" << expectedSize 
                  << ") in checkDimensions" << std::endl;
        return false;
    }
    
    // Check if arrays are allocated and match expected size
    bool diagOK = (diag != nullptr) && (nrows == expectedSize);
    bool rhsOK = (rhs != nullptr) && (nrows == expectedSize);
    
    // Optional: Add debugging information
    #ifdef DEBUG_SPARSPAK
    std::cout << "Dimension check: expected=" << expectedSize 
              << ", current nrows=" << nrows
              << ", diag allocated=" << (diag != nullptr)
              << ", rhs allocated=" << (rhs != nullptr) << std::endl;
    #endif
    
    return diagOK && rhsOK;
}
    
// Enhanced resize with preservation options
void SparspakSolver::resize(int newSize) {
    // Robust resizing for raw pointer arrays
    // This maintains existing data where possible and safely handles memory
    
    // Input validation
    if (newSize <= 0) {
        std::cerr << "Error: Invalid new size (" << newSize 
                  << ") in resize operation" << std::endl;
        return;
    }
    
    // If we're already the right size, no work needed
    if (nrows == newSize && diag != nullptr && rhs != nullptr) {
        // Just zero out the existing arrays to match std::vector resize behavior
        std::fill(diag, diag + nrows, 0.0);
        std::fill(rhs, rhs + nrows, 0.0);
        return;
    }
    
    // Store old values temporarily if we need to preserve data
    double* oldDiag = diag;
    double* oldRhs = rhs;
    int oldSize = nrows;
    
    // Allocate new arrays
    try {
        diag = new double[newSize];
        rhs = new double[newSize];
        
        // Initialize new arrays to zero (matching std::vector behavior)
        std::fill(diag, diag + newSize, 0.0);
        std::fill(rhs, rhs + newSize, 0.0);
        
        // If we had old data and both old and new arrays are valid,
        // copy over the common elements
        if (oldDiag != nullptr && oldRhs != nullptr && oldSize > 0) {
            int copySize = std::min(oldSize, newSize);
            
            // Copy the diagonal elements
            std::copy(oldDiag, oldDiag + copySize, diag);
            
            // Copy the right-hand-side elements  
            std::copy(oldRhs, oldRhs + copySize, rhs);
            
            #ifdef DEBUG_SPARSPAK
            std::cout << "Resized arrays from " << oldSize << " to " << newSize
                      << ", preserved " << copySize << " elements" << std::endl;
            #endif
        }
        
        // Update our size tracker
        nrows = newSize;
        
        // Clean up old memory
        delete[] oldDiag;
        delete[] oldRhs;
        
    } catch (const std::bad_alloc& e) {
        // If allocation failed, restore old state
        std::cerr << "Error: Memory allocation failed during resize: " 
                  << e.what() << std::endl;
        
        // Restore old pointers if new allocation failed
        diag = oldDiag;
        rhs = oldRhs;
        // nrows stays at oldSize
        
        throw; // Re-throw the exception so caller knows resize failed
    }
}


//=============================================================================

//  Store the matrix non-zero structure in a set of compressed adjacency lists
//  adjncy, where xadj has pointers to the starting index of each list in adjncy.

int compress(int n, int nnz, int* xrow, int* xcol, int* xadj, int* adjncy,
             int* xaij)
{
    int  flag = 0;
    int  *xadj2 = 0;
    int  *adjncy2 = 0;
    int  *nz = 0;

    // ... allocate memory
    xadj2 = new int[n+1];
    adjncy2 = new int[2*nnz];
    nz = new int[n];
    if ( xadj2 && adjncy2 && nz )
    {
        // ... build adjacency lists for the columns of A
        //     (for each column, store the rows indexes with non-zero coeffs.)
        buildAdjncy(n, nnz, xrow, xcol, xadj, adjncy, adjncy2, xaij, nz);

        // ... sort entries stored in each adjacency list
        sortAdjncy(n, nnz, xadj, adjncy, xadj2, adjncy2, nz);

        // ... re-label all row/col indexes to be 1-based
        //     (since the original Sparspak was written in Fortran)
        for (int i = 0; i < 2*nnz; i++) adjncy[i]++;
        for (int i = 0; i <= n; i++) xadj[i]++;
        flag = 1;
    }
    delete [] xadj2;
    delete [] adjncy2;
    delete [] nz;
    return flag;
}

//-----------------------------------------------------------------------------

//  Save the column index of each non-zero coefficient in a list for each row.

void buildAdjncy(
        int n, int nnz, int* xrow, int* xcol, int* xadj,
        int* adjncy, int* adjncy2, int* xaij, int* nz)
{
    int i, j, k, m, dup = 0;
    int adjncy2_size = 2 * nnz;
    // ... use adjncy to temporarily store non-duplicate coeffs.
    //     (parallel links create duplicates)
    int* nondup = adjncy;

    // ... count number of off-diagonal coeffs. in each column of A
    for (i = 0; i < n; i++) nz[i] = 0;
    for (k = 0; k < nnz; k++)
    {
        nz[xrow[k]]++;
        nz[xcol[k]]++;
    }

    // ... initialize adjncy2 to -1 (signifying an empty entry)
    for (i = 0; i < 2*nnz; i++) adjncy2[i] = -1;


    // ... make xadj array point to location in adjncy array where
    //     adjacency list for each column begins
    xadj[0] = 0;
    for (i = 0; i < n; i++)
    {
        xadj[i+1] = xadj[i] + nz[i];
        nz[i] = 0;
    }

    // ... fill adjncy2 array with non-zero row indexes for each column
    for (k = 0; k < nnz; k++)
    {
        i = xrow[k];
        j = xcol[k];

        // ... check for duplicate row/col
        dup = 0;
        for (m = xadj[i]; m < xadj[i]+nz[i]; m++)
        {
            if ( j == adjncy2[m] )
            {
                dup = 1;

                // ... mark xaij with negative of original coeff. index
                xaij[k] = -nondup[m];
                break;
            }
        }

        // ... if not a duplicate, add i and j to adjncy2
        if ( !dup )
        {
            m = xadj[i] + nz[i];
            adjncy2[m] = j;
            nondup[m] = k;
            nz[i]++;
            m = xadj[j] + nz[j];

            if (m >= 0 && m < adjncy2_size) {
                adjncy2[m] = i;
            }
            else {
                std::cerr << "Error: Invalid adjncy2 index: " << m
                    << " (max allowed: " << (adjncy2_size - 1) << ")" << std::endl;
                throw std::runtime_error("Array bounds violation in buildAdjncy");
            }

            adjncy2[m] = i;
            nondup[m] = k;
            nz[j]++;
        }
    }

    // ... re-construct xadj with duplicates removed
    for (i = 0; i < n; i++) xadj[i+1] = xadj[i] + nz[i];

    // ... transfer from adjncy2 to adjncy with duplicates removed
    k = 0;
    for (i = 0; i < 2*nnz; i++)
    {
        if ( adjncy2[i] >= 0 )
        {
            adjncy[k] = adjncy2[i];
            k++;
        }
    }
}

//-----------------------------------------------------------------------------

//  Sort the column indexes stored in each row's adjacency list.

void sortAdjncy(
        int n, int nnz, int* xadj, int* adjncy, int* xadj2,
        int* adjncy2, int* nz)
{
    // ... count number of non-zeros in each row
    //     (xadj[] holds # non-zeros in each column)
    for (int j = 0; j < n; j++) nz[j] = 0;
    for (int i = 0; i < n; i++)
    {
        for (int k = xadj[i]; k < xadj[i+1]; k++)
        {
            int j = adjncy[k];
            nz[j]++;
        }
    }

    // ... fill xadj2 with cumulative # non-zeros in each row
    xadj2[0] = 0;
    for (int i = 0; i < n; i++)
    {
        xadj2[i+1] = xadj2[i] + nz[i];
    }

    // ... transpose adjncy twice to order column indices
    transpose(n, xadj, adjncy, xadj2, adjncy2, nz);
    transpose(n, xadj2, adjncy2, xadj, adjncy, nz);
}

//-----------------------------------------------------------------------------

void transpose(int n, int* xadj1, int* adjncy1, int* xadj2, int* adjncy2, int* nz)
{
     for (int j = 0; j < n; j++) nz[j] = 0;
     for (int i = 0; i < n; i++)
     {
         for (int k = xadj1[i]; k < xadj1[i+1]; k++)
         {
             int j = adjncy1[k];
             int kk = xadj2[j] + nz[j];
             adjncy2[kk] = i;
             nz[j]++;
         }
     }
}

//-----------------------------------------------------------------------------

//  Apply the Multiple Minimum Degree algorithm to re-order the rows of the
//  matrix to minimize the amount of fill-in when the matrix is factorized.

int reorder(int n, int* xadj, int* adjncy, int* perm, int* invp, int& nnzl)
{
    // ... make a copy of the adjacency list
    int nnz2 = xadj[n];
    int* adjncy2 = new int[nnz2];
    if ( ! adjncy2 ) return 0;
    for (int i = 0; i < nnz2; i++)
    {
        adjncy2[i] = adjncy[i];
    }

    // ... create work arrays for row re-ordering
    int flag = 0;
    int *qsize = 0;
    int *llist = 0;
    int *marker = 0;
    int *dhead = 0;
    qsize = new int[n];
    llist = new int[n];
    marker = new int[n];
    dhead = new int[n];
    if ( qsize && llist && marker && dhead )
    {
        // ... call Sparspak sp_genmmd to apply multiple
        //     minimum degree re-ordering to A
        int delta = -1;
        int nofsub = 0;
        int maxint = std::numeric_limits<int>::max();
        sp_genmmd(&n, xadj, adjncy2, invp, perm, &delta, dhead, qsize,
                  llist, marker, &maxint, &nofsub);
        nnzl = nofsub;
        flag = 1;
    }

    // ... delete work arrays
    delete [] adjncy2;
    delete [] qsize;
    delete [] llist;
    delete [] marker;
    delete [] dhead;
    return flag;
}

//-----------------------------------------------------------------------------

//  Symbolically factorize the matrix

int factorize(
        int n, int& nnzl, int* xadj, int* adjncy, int* perm,
        int* invp, int* xlnz, int* xnzsub, int* nzsub)
{
    // ... create work arrays
    int flag = 0;
    int *mrglnk = 0;
    int *rchlnk = 0;
    int *marker = 0;
    mrglnk = new int[n];
    rchlnk = new int[n];
    marker = new int[n];

    // ... call Sparspak sp_smbfct routine
    if ( mrglnk && rchlnk && marker )
    {
        int maxlnz, maxsub = nnzl;
        sp_smbfct(n, xadj, adjncy, perm, invp, xlnz, maxlnz, xnzsub,
                  nzsub, maxsub, mrglnk, rchlnk, marker, flag);

        // ... update nnzl with size needed for lnz
        nnzl = maxlnz;

        // ... a return flag > 0 indicates insufficient memory;
        //     convert it to an error flag
        if ( flag > 0 ) flag = 0;
        else flag = 1;
    }
    delete [] mrglnk;
    delete [] rchlnk;
    delete [] marker;
    return flag;
}

//-----------------------------------------------------------------------------

//  Map the original off-diagonal coeffs. of the matrix to its factorized form.

void aij2lnz(
        int nnz, int* xrow, int* xcol, int* invp, int* xlnz, int* xnzsub,
        int* nzsub, int* xaij)
{
    int i, j, ksub;

    // ... adjust arrays for non-zero offset
    --xlnz; --xnzsub; --nzsub;

    // ... examine each non-zero coefficient
    for (int m = 0; m < nnz; m++)
    {
        // ... skip coeff. if it is marked as being a duplicate
        if ( xaij[m] < 0 ) continue;

        // ... determine its offset row & column below the diagonal
        //     (j is a column index and i is a row index with i > j)
        i = invp[xrow[m]];   // these return indexes starting from 1
        j = invp[xcol[m]];
        if ( i < j )
        {
            ksub = j;
            j = i;
            i = ksub;
        }

        // ... search for row index in nzsub
        ksub = xnzsub[j];
        for (int k = xlnz[j]; k < xlnz[j+1]; k++)
        {
            if ( nzsub[ksub] == i )
            {
                xaij[m] = k;
                break;
            }
            ksub++;
        }
    }

    // ... map any duplicate coeffs. (marked by the negative
    //     of the coeff. index they duplicate)
    for (int m = 0; m < nnz; m++)
    {
        if ( xaij[m] < 0 ) xaij[m] = xaij[-xaij[m]];
    }

    // ... reset arrays for zero offset
    ++xlnz; ++xnzsub; ++nzsub;
}

























int SparspakSolver::initMatrix(int nrows_, int nnz_, int* rowIndices, int* colIndices)
{
    // Basic validation
    if (nrows_ <= 0 || nnz_ <= 0 || !rowIndices || !colIndices) {
        std::cerr << "Invalid input parameters for matrix initialization" << std::endl;
        return 0;
    }

    // Store dimensions
    nrows = nrows_;
    nnz = nnz_;

    try {
        // Allocate core arrays with RAII
        std::vector<int> xadj(nrows + 1);
        std::vector<int> adjncy(2 * nnz);

        // Create compressed sparse matrix structure
        // Note: 'compress' produces 1-based indices.
        if (!compress(nrows, nnz, rowIndices, colIndices,
            xadj.data(), adjncy.data(), xaij)) {
            return 0;
        }

        // Perform matrix ordering (sp_genmmd expects 1-based indexing)
        perm = new int[nrows];
        invp = new int[nrows];
        if (!reorder(nrows, xadj.data(), adjncy.data(), perm, invp, nnzl)) {
            return 0;
        }

        // === Conversion from 1-based to 0-based indexing ===

        // Convert xadj (size: nrows+1)
        for (int i = 0; i < nrows + 1; i++) {
            xadj[i]--;  // Now xadj[i] is 0-based.
        }

        // Convert adjncy (size: 2*nnz)
        for (int i = 0; i < 2 * nnz; i++) {
            adjncy[i]--;  // Now adjncy[i] is 0-based.
        }

        // Convert the permutation arrays (each of size nrows)
        for (int i = 0; i < nrows; i++) {
            perm[i]--;  // Convert permutation from 1-based to 0-based.
            invp[i]--;  // Likewise for inverse permutation.
        }

        // Allocate and initialize solution arrays
        diag = new double[nrows]();
        rhs = new double[nrows];
        if (!rhs) return 0;
        std::fill(rhs, rhs + nrows, 0.0);

        // Calculate fill-in and allocate factorization arrays
        nnzl = calculateFillIn(nrows, xadj.data(), adjncy.data(), perm, invp);
        if (nnzl <= 0) return 0;

        lnz = new double[nnzl]();

        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Matrix initialization failed: " << e.what() << std::endl;
        return 0;
    }
}

int SparspakSolver::initWithInterior(int nrows_, int nnz_, int* xrow, int* xcol, double* segmentLengths, int* segmentTypes)
{
    
    // Add validation at the start
    std::cout << "Starting initialization with:"
        << "\nRequested rows: " << nrows_
        << "\nNon-zero elements: " << nnz_ << std::endl;

    // Verify input parameters
    if (nrows_ <= 0 || nnz_ <= 0 || !xrow || !xcol || !segmentLengths) {
        std::cout << "Invalid input parameters detected" << std::endl;
        return 0;
    }

    // Then, save the dimensions including interior points
    // nrows_ now includes both regular nodes and interior points
    nrows = nrows_;
    nnz = nnz_;   // Number of connections including interior segments

    // Allocate and explicitly initialize invp array
    invp = new int[nrows];
    for (int i = 0; i < nrows; i++) {
        invp[i] = i + 1;  // Proper 1-based indexing
        std::cout << "invp[" << i << "] = " << invp[i] << std::endl;
    }

    // Similar explicit initialization for xaij
    xaij = new int[nnz];
    for (int i = 0; i < nnz; i++) {
        xaij[i] = i + 1;  // Proper initialization
        std::cout << "xaij[" << i << "] = " << xaij[i] << std::endl;
    }
    memset(xaij, 0, nnz * sizeof(int));


    // Allocate additional arrays for interior point properties
    double* segmentWeights = new double[nnz];  // Weights based on segment lengths
    int* interiorFlags = new int[nrows];       // Flags to identify interior points
    if (!segmentWeights || !interiorFlags) return 0;

    // Process segment properties to create weights for matrix ordering
    for (int i = 0; i < nnz; i++) {
        // Shorter segments get higher weights to keep them together in ordering
        segmentWeights[i] = 1.0 / (segmentLengths[i] + 1e-10);
    }

    // Begin the matrix processing loop
    int* xadj;
    int* adjncy;
    double* adjwgt;  // New array for weighted adjacency
    int flag = 0;

    for (;;) {
        // Allocate space for weighted adjacency lists
        xadj = new int[nrows + 1];
        adjncy = new int[2 * nnz];
        adjwgt = new double[2 * nnz];  // Weights for each adjacency
        if (!xadj || !adjncy || !adjwgt) break;

        std::cout << "\nInitializing matrix solver..."
            << "\nTotal nodes: " << nrows
            << "\nTotal non-zero elements: " << nnz
            << "\nChecking segment lengths...\n";

        for (int i = 0; i < nnz; i++) {
            std::cout << "Segment " << i << " length: " << segmentLengths[i] << "\n";
        }

        // Ensure adjacency list arrays are initialized properly
        memset(adjncy, -1, 2 * nnz * sizeof(int));

        // Store matrix A in compressed format with weights
        if (!compressWithWeights(nrows, nnz, xrow, xcol, segmentWeights,
            xadj, adjncy, adjwgt, xaij)) break;

        // Reorder rows with consideration for interior points
        if (!reorderWithInterior(nrows, xadj, adjncy, adjwgt,
            perm, invp, nnzl)) break;

        // Allocate space for compressed storage of factorized matrix
        xlnz = new int[nrows + 1];
        xnzsub = new int[nrows + 1];
        nzsub = new int[nnzl];
        if (!xlnz || !xnzsub || !nzsub) break;

        // Symbolically factorize A to produce L, now handling interior points
        try {
            if (!factorizeWithInterior(nrows, nnzl, xadj, adjncy, adjwgt, perm, invp, xlnz, xnzsub, nzsub)) {
                std::cerr << "Error: Factorization failed due to invalid adjacency matrix!\n";
                return 0;
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Error: Exception in factorization: " << e.what() << "\n";
            return 0;
        }

        // All steps successful
        flag = 1;
        break;
    }

    // Clean up temporary arrays
    delete[] xadj;
    delete[] adjncy;
    delete[] adjwgt;
    delete[] segmentWeights;
    delete[] interiorFlags;

    if (!flag) return flag;

    // Map coefficients with interior point considerations
    aij2lnzWithInterior(nnz, xrow, xcol, invp, xlnz, xnzsub, nzsub, xaij);

    // Allocate solution arrays
    try {
        diag = new double[nrows]();  // () initializes to zero
        rhs = new double[nrows]();

        // Calculate nnzl before allocating lnz
        nnzl = calculateFillIn(nrows, xadj, adjncy, perm, invp);
        std::cout << "Computed nnzl: " << nnzl << "\n";
        if (nnzl <= 0) {
            std::cerr << "Error: Invalid nnzl calculated! Aborting initialization.\n";
            return 0;
        }

        lnz = new double[nnzl]();

        std::cout << "Arrays successfully allocated:"
            << "\ndiag size: " << nrows
            << "\nrhs size: " << nrows
            << "\nlnz size: " << nnzl << std::endl;

    }
    catch (const std::bad_alloc& e) {
        std::cout << "Memory allocation failed: " << e.what() << std::endl;
        return 0;
    }

    // Allocate work arrays for the enhanced solve method
    temp = new double[nrows];
    first = new int[nrows];
    link = new int[nrows];
    if (!temp || !first || !link) return 0;

    // Add verification after matrix processing
    std::cout << "\nPost-processing verification:"
        << "\nFinal nrows: " << nrows
        << "\nFinal nnz: " << nnz
        << "\nFinal nnzl: " << nnzl << std::endl;

    return 1;
}

bool SparspakSolver::compressWithWeights(int nrows, int nnz, int* xrow, int* xcol, double* weights, int* xadj, int* adjncy, double* adjwgt, int* xaij)
{
    // Initialize adjacency list pointers
    memset(xadj, 0, (nrows + 1) * sizeof(int));

    // Ensure adjacency matrix has valid connections
    for (int i = 0; i < nrows; i++) {
        if (xadj[i] < 0 || xadj[i] >= nnz) {
            std::cerr << "Error: xadj[" << i << "] = " << xadj[i]
                << " is out of bounds!\n";
        }
    }
    for (int i = 0; i < 2 * nnz; i++) {
        if (adjncy[i] < 0 || adjncy[i] >= nrows) {
            std::cerr << "Error: adjncy[" << i << "] = " << adjncy[i]
                << " is invalid before factorization!" << std::endl;
        }
    }

    // First pass: count entries in each row
    for (int k = 0; k < nnz; k++) {
        int i = xrow[k];
        int j = xcol[k];
        if (i < 0 || i >= nrows || j < 0 || j >= nrows) return false;
        xadj[i + 1]++;
        xadj[j + 1]++;
    }

    // Create cumulative sum
    for (int i = 1; i <= nrows; i++) {
        xadj[i] += xadj[i - 1];
    }

    // Second pass: fill adjacency lists with weights
    vector<int> pos(nrows, 0);
    for (int k = 0; k < nnz; k++) {
        int i = xrow[k];
        int j = xcol[k];

        // For pipes with segments, we use the segment-based weights
        double weight = weights[k];

        // Add entry (i,j)
        int ii = xadj[i] + pos[i];
        adjncy[ii] = j;
        adjwgt[ii] = weight;
        xaij[k] = ii;
        pos[i]++;

        // Add entry (j,i)
        int jj = xadj[j] + pos[j];
        adjncy[jj] = i;
        adjwgt[jj] = weight;
        pos[j]++;
    }

    std::cout << "\nAdjacency list:";
    for (int i = 0; i < nrows; i++) {
        std::cout << " xadj[" << i << "] = " << xadj[i] << "\n";
    }
    for (int i = 0; i < nnz; i++) {
        std::cout << " adjncy[" << i << "] = " << adjncy[i] << "\n";
    }

    return true;
}

bool SparspakSolver::reorderWithInterior(int nrows, int* xadj, int* adjncy, double* adjwgt, int* perm, int* invp, int& nnzl)
{
    // Initialize arrays for the modified minimum degree algorithm
    vector<int> degree(nrows);
    vector<bool> isInterior(nrows, false);
    vector<int> supernode(nrows);

    // Calculate initial degrees and identify interior nodes
    for (int i = 0; i < nrows; i++) {
        degree[i] = xadj[i + 1] - xadj[i];
        supernode[i] = i;

        // Check if this is an interior node (for pipes only)
        bool allPipeConnections = true;
        for (int j = xadj[i]; j < xadj[i + 1]; j++) {
            if (adjwgt[j] != 1.0) {  // Weight 1.0 indicates pipe segment
                allPipeConnections = false;
                break;
            }
        }
        isInterior[i] = allPipeConnections && degree[i] == 2;
    }

    // Modified minimum degree ordering that keeps pipe segments together
    for (int k = 0; k < nrows; k++) {
        // Find next node with minimum degree, preferring to keep pipe segments together
        int minDeg = INT_MAX;
        int minNode = -1;

        for (int i = 0; i < nrows; i++) {
            if (degree[i] > 0 && degree[i] < minDeg) {
                // For interior pipe nodes, check if adjacent nodes are already ordered
                if (isInterior[i]) {
                    bool adjacentOrdered = false;
                    for (int j = xadj[i]; j < xadj[i + 1]; j++) {
                        if (degree[adjncy[j]] == 0) {
                            adjacentOrdered = true;
                            break;
                        }
                    }
                    if (adjacentOrdered) {
                        minDeg = degree[i];
                        minNode = i;
                        break;
                    }
                }
                minDeg = degree[i];
                minNode = i;
            }
        }

        if (minNode == -1) return false;

        // Update ordering
        perm = new int[nrows];
        perm[k] = minNode;
        invp[minNode] = k;
        degree[minNode] = 0;

        // Update degrees of adjacent nodes
        for (int j = xadj[minNode]; j < xadj[minNode + 1]; j++) {
            int adj = adjncy[j];
            if (degree[adj] > 0) {
                degree[adj]--;
            }
        }
    }

    // Calculate fill-in (nnzl)
    nnzl = calculateFillIn(nrows, xadj, adjncy, perm, invp);
    return true;
}

void SparspakSolver::aij2lnzWithInterior(int nnz, int* xrow, int* xcol, int* invp, int* xlnz, int* xnzsub, int* nzsub, int* xaij)
{
    // Map each off-diagonal coefficient to its position in L
    for (int k = 0; k < nnz; k++) {
        int i = xrow[k];
        int j = xcol[k];
        int r = invp[i];  // Row i in reordered matrix
        int c = invp[j];  // Column j in reordered matrix

        // Ensure we're storing in lower triangular part
        if (r < c) {
            int temp = r;
            r = c;
            c = temp;
        }

        // Find position of (r,c) in L's data structure
        int position = -1;
        for (int m = xlnz[r]; m < xlnz[r + 1]; m++) {
            if (nzsub[m] == c) {
                position = m;
                break;
            }
        }

        if (position >= 0) {
            xaij[k] = position;
        }
    }
}

bool SparspakSolver::factorizeWithInterior(int nrows, int nnzl, int* xadj, int* adjncy, double* adjwgt, int* perm, int* invp, int* xlnz, int* xnzsub, int* nzsub)
{
    // Initialize the structure of L
    for (int i = 0; i <= nrows; i++) {
        xlnz[i] = 0;
        xnzsub[i] = 0;
    }

    // First pass: symbolic factorization
    vector<bool> marker(nrows, false);
    vector<int> temp(nrows);

    for (int k = 0; k < nrows; k++) {
        int pk = perm[k];  // Get original node number

        // Mark nodes in adjacency list
        for (int p = xadj[pk]; p < xadj[pk + 1]; p++) {
            int j = adjncy[p];
            if (invp[j] > k) {  // Only consider nodes not yet eliminated
                marker[invp[j]] = true;
                temp[xlnz[k]++] = invp[j];
            }
        }

        // Add fill-in entries
        for (int p = 0; p < xlnz[k]; p++) {
            int j = temp[p];
            int pj = perm[j];

            for (int q = xadj[pj]; q < xadj[pj + 1]; q++) {
                int i = invp[adjncy[q]];
                if (i > k && !marker[i]) {
                    marker[i] = true;
                    temp[xlnz[k]++] = i;
                }
            }
        }

        // Clear markers
        for (int p = 0; p < xlnz[k]; p++) {
            marker[temp[p]] = false;
        }
    }

    // Convert counts to pointers
    int nz = 0;
    for (int i = 0; i < nrows; i++) {
        int t = xlnz[i];
        xlnz[i] = nz;
        nz += t;
    }
    xlnz[nrows] = nz;

    if (nz != nnzl) return false;

    return true;
}

int SparspakSolver::calculateFillIn(int nrows, int* xadj, int* adjncy, int* perm, int* invp)
{
    // Create temporary storage for the elimination process
    vector<bool> marker(nrows, false);    // Tracks nodes in current elimination step
    vector<int> degree(nrows, 0);         // Tracks remaining connections for each node
    int fillCount = 0;                    // Counts the total nonzeros including fill-in

    // Initialize the degree array with current adjacency counts
    for (int i = 0; i < nrows; i++) {
        degree[i] = xadj[i + 1] - xadj[i];
        fillCount += degree[i];           // Add original nonzeros to count
    }
    fillCount /= 2;                       // Divide by 2 since we counted each edge twice

    // Process each node in the elimination order
    for (int k = 0; k < nrows; k++) {
        int elimNode = perm[k];           // Get the node being eliminated

        // Mark all neighbors of the eliminated node
        for (int p = xadj[elimNode]; p < xadj[elimNode + 1]; p++) {
            int neighbor = adjncy[p];
            if (invp[neighbor] > k) {     // Only consider nodes not yet eliminated
                marker[invp[neighbor]] = true;
            }
        }

        // For each marked neighbor
        for (int p = xadj[elimNode]; p < xadj[elimNode + 1]; p++) {
            int node1 = adjncy[p];
            int node1Pos = invp[node1];

            // Skip if this node is already eliminated
            if (node1Pos <= k) continue;

            // Look at its neighbors to find potential fill-ins
            for (int q = xadj[node1]; q < xadj[node1 + 1]; q++) {
                int node2 = adjncy[q];
                int node2Pos = invp[node2];

                // If this neighbor is not yet eliminated and not marked
                if (node2Pos > k && !marker[node2Pos]) {
                    marker[node2Pos] = true;   // Mark it
                    fillCount++;               // Found a fill-in entry
                }
            }
        }

        // Clear markers for next iteration
        for (int p = xadj[elimNode]; p < xadj[elimNode + 1]; p++) {
            int neighbor = adjncy[p];
            if (invp[neighbor] > k) {
                marker[invp[neighbor]] = false;
            }
        }
    }

    return fillCount;
}
