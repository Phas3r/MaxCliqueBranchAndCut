#pragma once

#include "graph.h"

namespace CliquerWrapper
{
	Vertices weightedFindSingleFromGraphCompliment(const Graph& graph, const std::vector<double>& vertexWeights, const double minWeight = 0.0);
	Vertices unweightedFindMax(const Graph& graph);
}
