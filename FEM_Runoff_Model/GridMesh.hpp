#pragma once
#include <unordered_map>
#include "Globals.hpp"
#include "Rectangle.hpp"
#include "LogManager.hpp"

#define MIN_RECTANGLE_WIDTH 0.01
#define MIN_RECTANGLE_HEIGHT 0.01
#define MIN_BOUNDARY_TO_DISCRETIZE 3	
#define MIN_RESOLUTION 1

bool GenerateGrid(std::vector<Vector2D> const & boundary,
	std::vector<Vector2D> & outNodes,
	std::unordered_map<size_t, Rectangle> & outRectList,
	std::vector<size_t> & outBoundaryNodes,
	size_t resolution = 100,
	double internalPadding = 0.001,
	double rayCastPadding = 1.0);
