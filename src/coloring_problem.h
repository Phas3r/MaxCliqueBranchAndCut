#pragma once
#include <vector>

#include "graph.h"

typedef size_t Color;
typedef std::vector<Color> Colors;

class GraphColoring
{
public:
    GraphColoring(size_t vertexCount);

    void colorVertex(Vertex v, Color c);
    Color getVertexColor(Vertex v) const;
    Vertices& getVerticesByColor(Color c);
    size_t getColorCount() const;
    void clear();

private:
    size_t mColorCount = 0;
    Colors mVertexToColor;
    std::vector<Vertices> mColorToVertices;
};

class ColoringSolver
{
public:
    size_t solve(const Graph& graph, std::vector<GraphColoring>& out_solutions);
    size_t solveIter(const Graph& graph, const Vertices& walkOrder, GraphColoring& out_colors, size_t upperBound);
    size_t solveRlf(const Graph& graph, GraphColoring & out_solution, const Vertices& walkOrder);
private:
    std::vector<std::vector<uint8_t>> mAvailableColors;
};