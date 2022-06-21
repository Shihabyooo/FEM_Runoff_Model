#include "Solvers.hpp"
#include "LogManager.hpp"



bool SolverSimple(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals)
{
	if (aMatrix.Columns() != bVector.Rows())
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}
	else if( !Matrix_f32::IsInvertible(aMatrix, true))
	{
		LogMan::Log("Error! Supplied factors matrix is not invertible. ", LOG_ERROR);
		return false;
	}

	outXVector = static_cast<Vector_f32>(aMatrix.Invert() * bVector);

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);

	return true;
}

bool SolverJacobi(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 & outXVector, Vector_f32 & outResiduals, double weight,  double threshold, size_t maxIterations)
{
	LogMan::Log("Using Jacobi solver.");

	if (aMatrix.Columns() != bVector.Rows())
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}

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

	if (aMatrix.Columns() != bVector.Rows())
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}
	
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
	
	if (aMatrix.Columns() != bVector.Rows())
	{
		LogMan::Log("Error! Supplied factors matrix and RHS vector are not of the same size.", LOG_ERROR);
		return false;
	}
	if (!aMatrix.IsSymmetric(0.0001F)) //TODO test/research whether PCG requires exact symmetry
	{
		LogMan::Log("Warning! Supplied factor matrix is highly assymetric. PCG solvers are designed for symmetric matrices.", LOG_WARN);
	}

	size_t systemSize = bVector.Rows();
	outXVector = Vector_f32(systemSize); //TODO check whether the PCG needs the middle elements to be zero. If not, use INITIAL_X_VALUE
	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);

	Matrix_f32 conditioner(aMatrix.Rows(), aMatrix.Columns()); //init zeroed conditioner matrix.

	for (int i = 0; i < systemSize; i++)
		conditioner[i][i] = aMatrix.GetValue(i, i);

	conditioner = conditioner.Invert();

	Vector_f32 dVector = static_cast<Vector_f32>(conditioner * outResiduals);

	double delta = (static_cast<Matrix_f32>(outResiduals.Transpose()) * dVector).GetValue(0, 0);
	double allowableTolerance = pow(threshold, 2.0F) * delta;

	size_t counter = 0;
	while (abs(delta) > allowableTolerance)
	{
		Vector_f32 qVector = static_cast<Vector_f32>(aMatrix * dVector);
		double alpha = delta / ((static_cast<Matrix_f32>(dVector.Transpose()) * qVector).GetValue(0, 0));
		outXVector = outXVector + dVector * alpha;
		outResiduals = outResiduals - qVector * alpha;
		
		Vector_f32 sVector = static_cast<Vector_f32>(conditioner * outResiduals);
		
		double deltaOld = delta;
		delta = (static_cast<Matrix_f32>(outResiduals.Transpose()) * sVector).GetValue(0, 0);
		double beta = delta / deltaOld;

		dVector = sVector + dVector * beta;

		counter++;
		if (counter > maxIterations)
		{
			LogMan::Log("Reached final iteration without converging on an acceptible solution.", LOG_WARN);
			return true;
		}
	}

	LogMan::Log("Reached acceptable residual in " + std::to_string(counter) + " iterations.", LOG_SUCCESS);

	ComputeResiduals(aMatrix, bVector, outXVector, outResiduals);

	return true;
}

void ComputeResiduals(Matrix_f32 const & aMatrix, Vector_f32 const & bVector, Vector_f32 const & xVector, Vector_f32 & outResiduals)
{
	//LogMan::Log("Residual computation not implemented yet.", LOG_WARN);
	outResiduals = bVector - static_cast<Vector_f32>(aMatrix * xVector);
}
