#pragma once
#include <MatricesPP.hpp>
#include "Globals.hpp"

#define DEFAULT_MAX_ITERATION 1000
#define DEFAULT_CONVERGENCE_THRESHOLD -1.0F
#define INITIAL_X_VALUE 1.0F

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

bool SolverSimple(	Matrix_f32 const & aMatrix,
					Vector_f32 const & bVector,
					Vector_f32 & outXVector,
					Vector_f32 & outResiduals);


//bool SolverGaussJordan(Matrix_f32 const & aMatrix,
//						Vector_f32 const & bVector,
//						Vector_f32 & outXVector,
//						Vector_f32 & outResiduals);

bool SolverJacobi(	Matrix_f32 const & aMatrix,
					Vector_f32 const & bVector,
					Vector_f32 & outXVector,					
					Vector_f32 & outResiduals,
					double weight = 1.0F, //Weight must be > 0 and <= 1.0F.
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverSOR(		Matrix_f32 const & aMatrix,
					Vector_f32 const & bVector,
					Vector_f32 & outXVector,
					Vector_f32 & outResiduals,
					double weight = 1.0F, //Weight must be from 0.0F to 2.0, clamped to range. 1.0F indicates Gauss-Seidel solution.
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

bool SolverPCG(		Matrix_f32 const & aMatrix,
					Vector_f32 const & bVector,
					Vector_f32 & outXVector,
					Vector_f32 & outResiduals, 
					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
					size_t maxIterations = DEFAULT_MAX_ITERATION);

//bool SolverBiCG(	Matrix_f32 const & aMatrix,
//					Vector_f32 const & bVector,
//					Vector_f32 & outXVector,
//					Vector_f32 & outResiduals,
//					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
//					size_t maxIterations = DEFAULT_MAX_ITERATION);

//bool SolverGMRES(	Matrix_f32 const & aMatrix,
//					Vector_f32 const & bVector,
//					Vector_f32 & outXVector,
//					Vector_f32 & outResiduals,
//					double threshold = DEFAULT_CONVERGENCE_THRESHOLD, //negative value indicates residual threshold should be chosen automatically.
//					size_t maxIterations = DEFAULT_MAX_ITERATION);

void ComputeResiduals(	Matrix_f32 const & aMatrix,
						Vector_f32 const & bVector,
						Vector_f32 const & xVector,
						Vector_f32 & outResiduals);