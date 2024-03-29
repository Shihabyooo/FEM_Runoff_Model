#pragma once
#include "Globals.hpp"

#define DEFAULT_MAX_ITERATION 10000
#define DEFAULT_CONVERGENCE_THRESHOLD -1.0
#define MIN_CONVERGENCE_THRESHOLD 0.000001
#define INITIAL_X_VALUE 1.0 //Note! If initial value = 0.0, some solvers may break

//Note: Initial version of this program opted to give the user the choice for solvers. Later the scope was limited only to\
Successive Over-Relaxtion solver, making most of this module useless. Leaving as is for now.

//These functions assume the system is in the form Ax = b

//Solver entry points, in general, require (for Ax=b system):
	//ref const to A (matrix),
	//ref const to b (column vector),
	//pointer to x (column vector),
	//pointer to a "residual" (column vector), default value set to NULL (i.e. don't compute the residual, if optional).
	//absolute value of maximum accepable residual magnitude (double), default set to -1.0F, indicating it should be chosen\
		 internally from the values of b (e.g. 0.1% of smallest, non zero value in b?)
	//max number of iterations (size_t), default value set to 1000.
//Each function returns a bool, false if the supplied input does not constitute an acceptable problem system (e.g. aMatrix is not square\
or the length does not equal b and x lengths, etc). Returns true otherwise (regardless of whether the solver converged on a solution or not).
//Note: Not all issues with the system are treated as errors. e.g. for Jacobi and SOR, if the supplied bMatrix is not diagonal dominant\
log a warning with the fact but still compute it.

//TODO consider having a generalized struct for solver log (i.e. has statistics about solution process).

bool Solve(	Matrix_f64 const & aMatrix,
			Vector_f64 const & bVector,
			Vector_f64 & outXVector,
			Vector_f64 & outResiduals,
			ModelParameters const & params);

bool SolverSimple(	Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,
					Vector_f64 & outResiduals);

bool Gaussian(	Matrix_f64 const & aMatrix,
				Vector_f64 const & bVector,
				Vector_f64 & outXVector,
				Vector_f64 & outResiduals);

bool SolverJacobi(	Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,					
					Vector_f64 & outResiduals,
					double weight = 1.0F, //Weight must be > 0 and <= 1.0F.
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverSOR(		Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,
					Vector_f64 & outResiduals,
					double weight = 1.0F, //Weight must be from 0.0F to 2.0, clamped to range. 1.0F indicates Gauss-Seidel solution.
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverPCG(		Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,
					Vector_f64 & outResiduals, 
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverBiCG(	Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,
					Vector_f64 & outResiduals,
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverCGS(	Matrix_f64 const & aMatrix,
					Vector_f64 const & bVector,
					Vector_f64 & outXVector,
					Vector_f64 & outResiduals,
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

void ComputeResiduals(	Matrix_f64 const & aMatrix,
						Vector_f64 const & bVector,
						Vector_f64 const & xVector,
						Vector_f64 & outResiduals);