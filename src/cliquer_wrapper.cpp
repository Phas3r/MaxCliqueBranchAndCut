#include "cliquer_wrapper.h"

#include <algorithm>
#include <cliquer.h>

namespace CliquerWrapper
{
	class graph_t_wrapper
	{
	public:
		graph_t_wrapper(const Graph& graph, const std::vector<double>& vertexWeights = std::vector<double>(), bool compliment = false)
		{
			mpGraph = graph_new((int)graph.GetVertexCount());
			for (size_t v1 = 0; v1 < graph.GetVertexCount(); ++v1) {
				for (size_t v2 = v1; v2 < graph.GetVertexCount(); ++v2) {
					if (v1 == v2) {
						continue;
					}
					bool isEdge = graph.HasEdge(v1, v2);
					if (isEdge != compliment) { // != is used as XOR
						GRAPH_ADD_EDGE(mpGraph, v1, v2);
					}
				}
			}
			for (size_t i = 0; i < vertexWeights.size(); ++i) {
				mpGraph->weights[i] = (float)vertexWeights[i];
			}
		}
		~graph_t_wrapper()
		{
			if (mpGraph) {
				graph_free(mpGraph);
				mpGraph = nullptr;
			}
		}
		graph_t* operator&()
		{
			return mpGraph;
		}
		graph_t* operator->()
		{
			return mpGraph;
		}
	private:
		graph_t* mpGraph = nullptr;
	};

	Vertices getVerticesAndFreeSet(set_t cliquerSet)
	{
		Vertices result;
		if (cliquerSet) {
			result.reserve(set_size(cliquerSet));
			int i = -1;
			while ((i = set_return_next(cliquerSet, i)) >= 0) {
				result.push_back(Vertex(i));
			}
			set_free(cliquerSet);
		}
		return result;
	}

	Vertices weightedFindSingleFromGraphCompliment(const Graph& graph, const std::vector<double>& vertexWeights, const double minWeight)
	{
		graph_t_wrapper cliquer_graph(graph, vertexWeights, true);
		set_t resultSet = clique_find_single(&cliquer_graph, (float)(minWeight), 0.0, true, nullptr);
		return getVerticesAndFreeSet(resultSet);
	}

	Vertices unweightedFindMax(const Graph& graph)
	{
		graph_t_wrapper cliquer_graph(graph);
		set_t resultSet = clique_find_single(&cliquer_graph, 0.0, 0.0, true, nullptr);
		return getVerticesAndFreeSet(resultSet);
	}

}
