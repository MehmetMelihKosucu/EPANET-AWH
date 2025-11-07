/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file sparspaksolver.h
//! \brief Description of the SparspakSolver class.

#ifndef SPARSPAKSOLVER_H_
#define SPARSPAKSOLVER_H_

#include "matrixsolver.h"

//! \class SparspakSolver
//! \brief Solves Ax = b using the SPARSPAK routines.
//!
//! This class is derived from the MatrixSolver class and provides an
//! interface to the SPARSPAK routines, originally developed by George
//! and Liu, for re-ordering, factorizing, and solving via Cholesky
//! decomposition a sparse, symmetric, positive definite set of linear
//! equations Ax = b.

class SparspakSolver: public MatrixSolver
{
  public:

    // Constructor/Destructor

    SparspakSolver(std::ostream& logger);
    ~SparspakSolver();

    // Methods

    int    init(int nrows, int nnz, int* xrow, int* xcol);
    void   reset();

    double getDiag(int i);
    double getOffDiag(int i);
    double getRhs(int i);

    void   setDiag(int i, double a);
    void   setRhs(int i, double b);
    void   addToDiag(int i, double a);
    void   addToOffDiag(int j, double a);
    void   addToRhs(int i, double b);
    int    solve(int n, double x[]);
    bool checkDimensions(int expectedSize);
    void resize(int newSize);

    int initMatrix(int nrows_, int nnz_, int* rowIndices, int* colIndices);

    int initWithInterior(int nrows_, int nnz_, int* xrow, int* xcol, double* segmentLengths, int* segmentTypes);                                              // These functions are useless now.
                                                                                                                                                              
    bool compressWithWeights(int nrows, int nnz, int* xrow, int* xcol, double* weights, int* xadj, int* adjncy, double* adjwgt, int* xaij);                   // These functions are useless now.
                                                                                                                                                             
    bool reorderWithInterior(int nrows, int* xadj, int* adjncy, double* adjwgt, int* perm, int* invp, int& nnzl);                                             // These functions are useless now.
                                                                                                                                                             
    bool factorizeWithInterior(int nrows, int nnzl, int* xadj, int* adjncy, double* adjwgt, int* perm, int* invp, int* xlnz, int* xnzsub, int* nzsub);        // These functions are useless now.
                                                                                                                                                             
    void aij2lnzWithInterior(int nnz, int* xrow, int* xcol, int* invp, int* xlnz, int* xnzsub, int* nzsub, int* xaij);                                        // These functions are useless now.

    int calculateFillIn(int nrows, int* xadj, int* adjncy, int* perm, int* invp);

    int getConnectionCount() const;

  private:

    int     nrows;    // number of rows in system Ax = b
    int     nnz;      // number of non-zero off-diag. coeffs. in A
    int     nnzl;     // number of non-zero off-diag. coeffs. in factorized matrix L
    int*    perm;     // permutation of rows in A
    int*    invp;     // inverse row permutation
    int*    xlnz;     // index vector for non-zero entries in L
    int*    xnzsub;   // index vector for entries of nzsub
    int*    nzsub;    // column indexes for non-zero entries in each row of L
    int*    xaij;     // maps off-diag. coeffs. of A to lnz
    int*    link;     // work array
    int*    first;    // work array
    double* lnz;      // off-diag. coeffs. of factorized matrix L
    double* diag;     // diagonal coeffs. of A
    double* rhs;      // right hand side vector
    double* temp;     // work array
    std::ostream& msgLog;

    // Add size tracking members
    int arraySize;  // Size of our arrays

    bool isIndexValid(int index) const {return (index >= 0 && index < arraySize);
    }
};

#endif
