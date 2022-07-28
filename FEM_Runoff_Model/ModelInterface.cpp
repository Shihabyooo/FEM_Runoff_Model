#pragma once
#include "ModelInterface.hpp"
#include "ModelGlobals.hpp"
#include "ModelImplementation.hpp"
#include "SpatialDataModule.hpp"


void ComputeBoundingBox(std::vector<Vector2D> const & points, Vector2D & min, Vector2D & max)
{
	min = points[0];
	max = points[0];

	for (auto it = points.begin(); it < points.end(); it++)
	{
		min.x = Min(it->x, min.x);
		min.y = Min(it->y, min.y);
		max.x = Max(it->x, max.x);
		max.y = Max(it->y, max.y);
	}
}

void ResetMeshes()
{
	nodes.clear();
	triangles.clear();
	boundaryNodes.clear();
}

bool LoadWatershedBoundary(std::string const & boundaryPath)
{
	if (!FileIO::LoadVectorPath(boundaryPath, shedBoundary))
	{
		return false;
	}

	ComputeBoundingBox(shedBoundary, shedSW, shedNE);

	return true;
}

bool GenerateMesh(MeshGeneratorParameters & params)
{
	LogMan::Log("Attempting to generate FEM mesh");
	ResetMeshes();

	if (!MeshGen::ValidateParameters(params))
	{
		LogMan::Log("ERROR! Invalid meshing parameters!", LOG_ERROR);
		return false;
	}

	params.boundary = &shedBoundary;

	void * elementsContainerPtr = NULL;
	
	bool status = MeshGen::GenerateMesh(params, triangles, &nodes, &boundaryNodes);

	if (status)
		ComputeBoundingBox(nodes, nodesSW, nodesNE);
	else
		ResetMeshes();

	return status;
}

Triangle const * GetElementContainingPoint(Vector2D const & pos)
{
	if (triangles.size() < 1)
		return NULL;

	for (auto it = triangles.begin(); it != triangles.end(); ++it)
		if (it->second.ContainsPoint(pos))
			return &it->second;

	return NULL;
}

bool UpdateNode(size_t id, Vector2D const & newPos)
{
	if (id >= nodes.size())
		return false;

	nodes[id] = newPos;

	//look for triangles using this node and log warning if the change made it invalid
	for (auto it = triangles.begin(); it != triangles.end(); ++it)
	{
		if (it->second.ContainsVertex(id))
		{
			if (!it->second.Validate())
				LogMan::Log("Warning! Updated triangle " + std::to_string(it->second.id) + " is invalid (Area: " + std::to_string(it->second.Area()) + ").", LOG_WARN);
		}
	}

	return true;
}

bool IsBoundaryNode(size_t nodeID)
{
	for (auto it = boundaryNodes.begin(); it != boundaryNodes.end(); ++it)
		if (*it == nodeID)
			return true;
	return false;
}

bool LoadTimeSeries(std::string const & path, TimeSeries & ts)
{
	return FileIO::LoadTimeSeries(path, ts);
}

bool CheckParameters(ModelParameters const & params)
{
	LogMan::Log("Checking parameters");
	bool status = true;

	//Check fails
	//TODO for raster file checks, check also that they are georeffed rasters and inclusive of the mesh boundary.
	
	if (triangles.size() < 1) //existence of tris implies existence of nodes (though I probably need to protect both from change outside this file)
	{
		LogMan::Log("ERROR! No loaded mesh.", LOG_ERROR);
		status = false;
	}

	if (!FileIO::FileExists(params.demPath))
	{
		LogMan::Log("ERROR! Must supply a terrain DEM for the region.", LOG_ERROR);
		status = false;
	}

	if (!FileIO::FileExists(params.slopesPath))
	{
		LogMan::Log("ERROR! Must supply a terrain Slopes for the region.", LOG_ERROR);
		status = false;
	}

	if (!FileIO::FileExists(params.fdrPath))
	{
		LogMan::Log("ERROR! Must supply a flow direction map (AGNPS format) for the region.", LOG_ERROR);
		status = false;
	}

	if (!params.variablePrecipitation && !params.unitTimeSeries.IsValid())
	{
		LogMan::Log("ERROR! Must provide a valid Time-series", LOG_ERROR);
		status = false;
	}

	if (params.variableManningCoefficients && !FileIO::FileExists(params.manningCoefficientRasterPath))
	{
		LogMan::Log("ERROR! Must set a Manning coefficients raster when using variable manning coefficients", LOG_ERROR);
		status = false;
	}
	else if (!params.variableManningCoefficients && params.fixedManningCoeffient <= 0.0)
	{
		LogMan::Log("ERROR! Must set a Manning coefficient  when using fixed manning coefficients", LOG_ERROR);
		status = false;
	}

	if (params.timeStep <= 0.0)
	{
		LogMan::Log("ERROR! Time step must be a positive, real number", LOG_ERROR);
		status = false;
	}

	if (params.endTime <= params.startTime || (params.startTime + params.timeStep) > params.endTime)
	{
		LogMan::Log("ERROR! End time must be greater than Start time with at least Time Step", LOG_ERROR);
		status = false;
	}

	if (params.femOmega < 0.0 || params.femOmega > 1.0)
	{
		LogMan::Log("ERROR! Omega must be from 0.0 to 1.0", LOG_ERROR);
		status = false;
	}

	//Check warn
	//TODO add warning here about low iteration time and such.

	return status;
}

bool Simulate(ModelParameters const & params)
{
	return RunSimulation(params);
}

std::vector<Vector2D> const & GetNodes()
{
	return nodes;
}

std::unordered_map<size_t, Triangle> const & GetTriangles()
{
	return triangles;
}

std::vector<size_t> const & GetBoundaryNodes()
{
	return boundaryNodes;
}

std::pair<Vector2D const&, Vector2D const&> GetNodesBoundingBox()
{
	return std::pair<Vector2D const&, Vector2D const&>(nodesSW, nodesNE);
}

std::pair<Vector2D const&, Vector2D const&> GetWatershedBoundingBox()
{
	return std::pair<Vector2D const&, Vector2D const&>(shedSW, shedNE);
}

std::vector<Vector2D> const & GetWatershedBoundary()
{
	return shedBoundary;
}
