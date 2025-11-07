/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file matrixsolver.h
//! \brief Description of the MatrixSolver class.

#ifndef MATRIXSOLVER_H_
#define MATRIXSOLVER_H_

#include <string>
#include <ostream>
#include <vector>
#include "node.h"

//! \class MatrixSolver
//! \brief Abstract class for solving a set of linear equations.
//!
//! This class defines an interface for solving a set of linear
//! equations generated from a network of links and nodes.
//!
//! The system of equations is expressed as Ax = b where A is a square
//! symmetric coefficient matrix, b is a right hand side vector, and
//! x is a vector of unknowns.

class MatrixSolver
{
  public:

    MatrixSolver();
    virtual ~MatrixSolver();
    static  MatrixSolver* factory(const std::string solver, std::ostream& logger);

    virtual int    init(int nRows, int nOffDiags, int offDiagRow[], int offDiagCol[])= 0;
	  virtual int    initWithInterior(int nRows, int nOffDiags, int offDiagRow[], int offDiagCol[], double segmentLengths[], int segmentTypes[]) = 0;
	  virtual int    initMatrix(int nrows_, int nnz_, int* rowIndices, int* colIndices) = 0;
    virtual void   reset() = 0;

    virtual double getDiag(int i)    {return 0.0;}
    virtual double getOffDiag(int i) {return 0.0;}
    virtual double getRhs(int i)     {return 0.0;}

    virtual void   setDiag(int row, double a) = 0;
    virtual void   setRhs(int row, double b) = 0;
    virtual void   addToDiag(int row, double a) = 0;
    virtual void   addToOffDiag(int offDiag, double a) = 0;
    virtual void   addToRhs(int row, double b) = 0;
    virtual bool   checkDimensions(int expectedSize) { return true; }
    virtual void resize(int newSize) = 0;
    virtual int    solve(int nRows, double x[]) = 0;
    virtual void setIT1Coupling(int offDiagIndex, Node* it1Node, double coefficient);
    
    std::vector<double> offDiag;  // Already defined but showing for clarity
    std::vector<Node*> it1CouplingNodes;  // Vector to store IT1 nodes
    int offDiagIndex; // Current off-diagonal index
    int offDiagSize;  // Size of off-diagonal array
};

#endif
