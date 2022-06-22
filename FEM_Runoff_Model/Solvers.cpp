#include "Solvers.hpp"
#include "LogManager.hpp"


double ComputeThreshold(Vector_f32 const & bVector)
{
	float minValue = FLT_MAX;
	for (int i = 0; i < bVector.Rows(); i++)
		minValue = Min(abs(bVector.GetValue(i)), minValue);

	float threshold = Max(0.001F * minValue, MIN_CONVERGENCE_THRESHOLD);
	LogMan::Log("Using automatic threshold of: " + std::to_string(threshold));
	return threshold;
}

bool CheckSystem(Matrix_f32 const & aMatrix, Vector_f32 const & bVector)
{
	if (aMatrix.Columns() != bVector.Rows())
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}
	if (!Matrix_f32::IsSquared(aMatrix))
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}
	return true;
}

bool SolverSimple(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals)
{
	LogMan::Log("Using Simple solver.");

	if (!CheckSystem(aMatrix, bVector))
		return false;

	else if( !Matrix_f32::IsInvertible(aMatrix, true))
	{
		LogMan::Log("Error! Supplied factors matrix is not invertible. ", LOG_ERROR);
		return false;
	}

	outXVector = static_cast<Vector_f32>(aMatrix.Invert() * bVector);

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);

	return true;
}

bool Gaussian(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals)
{
	if (!CheckSystem(aMatrix, bVector))
		return false;

	//create augMat = [A|b]
	//create a float ** = float * [systemSize], and then have each element point to address on xVector content. We need this because\
	the augMat may very likely get rearanged (row swapping), this means the order of x elements must change too, but we can't, so we\
	change the float ** instead, and then when assigning x values, we dereference *float[i].
	//loop over first column elements, find row with largest element, swap it with the first row.
		//If all first column elements are zero (unlikely here, but still), return false.
	//typical elementary row ops to reduce matrix to upper tri echelon.
	//backwards substitution -> from systemSize - 1 to 0

	return true;
}

bool SolverJacobi(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double weight,  double threshold, size_t maxIterations)
{
	LogMan::Log("Using Jacobi solver.");

	if (!CheckSystem(aMatrix, bVector))
		return false;

	if (threshold < MIN_CONVERGENCE_THRESHOLD)
		threshold = ComputeThreshold(bVector);

	size_t systemSize = bVector.Rows();
	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);

	//create diagonal (inverted), and lower + upper matrices, so that A = diag + (L + U)
	Matrix_f32 diagInv(systemSize, systemSize), LU(systemSize, systemSize);
	for (int i = 0; i < systemSize; i++)
	{
		for (int j = 0; j < systemSize; j++)
		{
			if (i == j)
				diagInv[i][j] = 1.0F / aMatrix.GetValue(i, j); //invert of a diagonal matrix is reciprocal of each -diagonal- element
			else
				LU[i][j] = aMatrix.GetValue(i, j);
		}
	}

	Vector_f32 tempX = outXVector;
	for (int i = 0; i < maxIterations; i++)
	{
		outXVector = static_cast<Vector_f32>((diagInv * (bVector - static_cast<Vector_f32>(LU * tempX))) * weight) + (tempX * (1.0F - weight));

		ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
		if (outResiduals.Magnitude() <= threshold)
		{
			LogMan::Log("Reached acceptable residual in " + std::to_string(i) + " iterations.", LOG_SUCCESS);
			return true;
		}
		tempX = outXVector;
	}
	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
	
	return true; //while we didn't technically reach a good solution, the routine did execute to the end, so...
}

bool SolverSOR(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double weight, double threshold, size_t maxIterations)
{
	LogMan::Log("Using Succesive Over-Relaxation solver.");

	if (!CheckSystem(aMatrix, bVector))
		return false;
	
	if (threshold < MIN_CONVERGENCE_THRESHOLD)
		threshold = ComputeThreshold(bVector);
	
	size_t systemSize = bVector.Rows();
	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);

	//clamp weight to 0.0F-2.0F range
	weight = weight < 0.0F ? 0.0F : (weight > 2.0F ? 2.0F : weight);

	for (int i = 0; i < maxIterations; i++)
	{
		for (int j = 0; j < systemSize; j++)
		{
			double internalSum = 0.0F;
			for (int k = 0; k <= j - 1; k++)
				internalSum += aMatrix.GetValue(j, k) * outXVector[k];

			for (int k = j + 1; k < systemSize; k++)
				internalSum += aMatrix.GetValue(j, k) * outXVector[k];

			internalSum = (bVector.GetValue(j) - internalSum) / aMatrix.GetValue(j, j);
			outXVector[j] = outXVector[j] + weight * (internalSum - outXVector[j]);
		}

		ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
		if (outResiduals.Magnitude() <= threshold)
		{
			LogMan::Log("Reached acceptable residual in " + std::to_string(i) + " iterations.", LOG_SUCCESS);
			return true;
		}
	}
	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
	
	return true;
}

//TODO double, triple and quadruple check this method. It was translated from text around coup time, and I don't really remember how the hell\
it worked...
bool SolverPCG(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double threshold, size_t maxIterations)
{
	LogMan::Log("Using Preconditioned Conjugate Gradients solver.");
	
	if (!CheckSystem(aMatrix, bVector))
		return false;

	if (!aMatrix.IsSymmetric(0.0001F)) //TODO test/research whether PCG requires exact symmetry
	{
		LogMan::Log("Warning! Supplied factor matrix is assymetric. PCG solvers are designed for symmetric matrices.", LOG_WARN);
	}

	if (threshold < MIN_CONVERGENCE_THRESHOLD)
		threshold = ComputeThreshold(bVector);

	size_t systemSize = bVector.Rows();
	//init xVector
	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);
	/*outXVector[0] = bVector.GetValue(0);
	outXVector[systemSize - 1] = bVector.GetValue(systemSize - 1);*/

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);

	Matrix_f32 conditioner(systemSize, systemSize); //init zeroed conditioner matrix.
	//Invert of a diagonal matrix = replace diagonal elements with reciprocal. No need to use expensive methods of Matrix_f32
	/*for (int i = 0; i < systemSize; i++)
		conditioner[i][i] = 1.0F / aMatrix.GetValue(i, i);*/
	
	conditioner = Matrix_f32::Identity(systemSize); //test

	Vector_f32 dVector = conditioner * outResiduals;

	double delta = (static_cast<Matrix_f32>(outResiduals.Transpose()) * dVector).GetValue(0, 0);
	double allowableTolerance = pow(threshold, 2) * delta; //source paper's algorithm
	allowableTolerance = Min(allowableTolerance, threshold); //this is probably not necessary...

	//while (abs(delta) > allowableTolerance)
	for (int i = 0; i < maxIterations; i++)
	{
		
		Vector_f32 qVector = aMatrix * dVector;
		double alpha = delta / ((static_cast<Matrix_f32>(dVector.Transpose()) * qVector).GetValue(0, 0));
		outXVector = outXVector + dVector * alpha;
		outResiduals = outResiduals - qVector * alpha;
		
		Vector_f32 sVector = conditioner * outResiduals;
		
		double deltaOld = delta;
		delta = (static_cast<Matrix_f32>(outResiduals.Transpose()) * sVector).GetValue(0, 0);
		double beta = delta / deltaOld;

		dVector = sVector + dVector * beta;

		if (abs(delta) < allowableTolerance)
		{
			LogMan::Log("Reached acceptable residual in " + std::to_string(i) + " iterations.", LOG_SUCCESS);
			ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
			return true;
		}
	}

	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
	return true;
}

bool SolverBiCG(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double threshold, size_t maxIterations)
{
	LogMan::Log("Using Bi-Conjugate Gradients solver.");

	if (!CheckSystem(aMatrix, bVector))
		return false;
	
	if (threshold < MIN_CONVERGENCE_THRESHOLD)
		threshold = ComputeThreshold(bVector);

	size_t systemSize = bVector.Rows();
	//init xVector
	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);
	/*outXVector[0] = bVector.GetValue(0);
	outXVector[systemSize - 1] = bVector.GetValue(systemSize - 1);*/

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
	Vector_f32 residuals2 = outResiduals;
	
	//Matrix_f32 conditioner(systemSize, systemSize); //init zeroed conditioner matrix.
	////Invert of a diagonal matrix = replace diagonal elements with reciprocal. No need to use expensive methods of Matrix_f32
	//for (int i = 0; i < systemSize; i++)
	//	conditioner[i][i] = 1.0F / aMatrix.GetValue(i, i);
	////conditioner = Matrix_f32::Identity(systemSize); //test

	Vector_f32 d = outResiduals;
	Vector_f32 d2 = residuals2;

	double delta = (static_cast<Matrix_f32>(residuals2.Transpose()) * outResiduals)[0][0];
	double deltaOld = delta;
	double allowableTolerance = pow(threshold, 2) * delta; //source paper's algorithm
	allowableTolerance = Min(allowableTolerance, threshold); //this is probably not necessary...
	
	for (int i = 0; i < maxIterations; i++)
	{
		Vector_f32 q = aMatrix * d;
		Vector_f32 q2 = static_cast<Matrix_f32>(aMatrix.Transpose()) * d;

		double alpha = delta / (static_cast<Matrix_f32>(d2.Transpose()) * q)[0][0];
		
		if (alpha == 0.0F || isnan(alpha) || isinf(alpha))
		{
			LogMan::Log("Algorithm broke at iteration " + std::to_string(i) + ". Alpha is 0, NaN or Inf.", LOG_ERROR);
			return false;
		}

		outXVector = outXVector + d * alpha;
		outResiduals = outResiduals - q * alpha;
		residuals2 = residuals2 - q2 * alpha;

		deltaOld = delta;
		delta = (static_cast<Matrix_f32>(residuals2.Transpose()) * outResiduals)[0][0];
		
		double beta = delta / deltaOld;
		d = outResiduals + d * beta;
		d2 = residuals2 + d2 * beta;

		//test orthogonality relations:
		//std::cout << "r'T * r: " << (static_cast<Matrix_f32>(residuals2.Transpose()) * outResiduals)[0][0] << std::endl; //test
		//std::cout << "p'T A p: " << (static_cast<Matrix_f32>(d2.Transpose()) * q)[0][0] << std::endl; //test

		if (abs(delta) < allowableTolerance)
		{
			LogMan::Log("Reached acceptable residual in " + std::to_string(i) + " iterations.", LOG_SUCCESS);
			ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
			return true;
		}
	}

	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
	return true;
}

//bool SolverGMRES(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double threshold, size_t maxIterations)
//{
//	LogMan::Log("Using Generalized Minimal Residual solver.");
//
//	if (!CheckSystem(aMatrix, bVector))
//		return false;
//
//	if (threshold < MIN_CONVERGENCE_THRESHOLD)
//		threshold = ComputeThreshold(bVector);
//
//	size_t systemSize = bVector.Rows();
//	//init xVector
//	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);
//
//	for (int i = 0; i < maxIterations; i++)
//	{
//		// r = b - Ax
//		ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
//
//		//v = r / ||r||2
//		//s = ||r||2 e1
//
//		//loop
//			//w = Av
//			//loop
//				//h =
//	}
//
//	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
//	return true;
//}

bool SolverCGS(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double threshold, size_t maxIterations)
{
	LogMan::Log("Using Conjugate Gradients Squared solver.");

	if (!CheckSystem(aMatrix, bVector))
		return false;

	if (threshold < MIN_CONVERGENCE_THRESHOLD)
		threshold = ComputeThreshold(bVector);

	size_t systemSize = bVector.Rows();
	//init xVector
	outXVector = Vector_f32(systemSize, INITIAL_X_VALUE);

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);
	
	Matrix_f32 conditioner(systemSize, systemSize); //init zeroed conditioner matrix.
	//Invert of a diagonal matrix = replace diagonal elements with reciprocal. No need to use expensive methods of Matrix_f32
	/*for (int i = 0; i < systemSize; i++)
		conditioner[i][i] = 1.0F / aMatrix.GetValue(i, i);*/

	conditioner = Matrix_f32::Identity(systemSize); //test

	Vector_f32 residuals2 = outResiduals;

	//u = r
	//p = u
	Vector_f32 u = outResiduals;
	Vector_f32 p = u;
	double delta = (static_cast<Matrix_f32>(residuals2.Transpose()) * outResiduals)[0][0];
	for (int i = 0; i < maxIterations; i++)
	{
		//delta = r'Transpose * r
		double deltaOld = delta;
		delta = (static_cast<Matrix_f32>(residuals2.Transpose()) * outResiduals)[0][0];
		//if delta == zero -> fail
		if (delta == 0.0F || isnan(delta) || isinf(delta))
		{
			LogMan::Log("Algorithm broke at iteration " + std::to_string(i) + ". Delta is 0, NaN or Inf.", LOG_ERROR);
			return false;
		}

		//p' = Cp
		//v' = Ap'
		Vector_f32 p2 = conditioner * p;
		Vector_f32 v = aMatrix * p2;

		//alpha = delta / (r'Transpose * v')
		double alpha = delta / (static_cast<Matrix_f32>(residuals2.Transpose()) * v)[0][0];

		//q = u - alpha * v'
		Vector_f32 q = u - v * alpha;
		//u' = u + q
		Vector_f32 u2 = u + q;
		//x = x + alpha * u'
		outXVector = outXVector + u2 * alpha;
		//q' = Au'
		Vector_f32 q2 = aMatrix * u2;
		//r = r - alpha*q'
		outResiduals = outResiduals - q2 * alpha;
		//beta = alpha / alpha_old
		double beta = delta / deltaOld;
		//u = r + beta * q
		u = outResiduals + q * beta;
		//p = u + beta * (q + beta * p)
		p = u + (q + p * beta) * beta;

		if (outResiduals.Magnitude() <= threshold)
		{
			LogMan::Log("Reached acceptable residual in " + std::to_string(i) + " iterations.", LOG_SUCCESS);
			return true;
		}
	}

	LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
	return true;
}

void ComputeResiduals(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 const & xVector, Vector_f32 & outResiduals)
{
	//LogMan::Log("Residual computation not implemented yet.", LOG_WARN);
	outResiduals = bVector - aMatrix * xVector;
}
