#include "coloring_problem.h"
#include <random>
#include <iostream>
#include <map>
#include <cassert>

GraphColoring::GraphColoring(size_t vertexCount)
    : mVertexToColor(vertexCount, vertexCount)
    , mColorToVertices(vertexCount)
{
    for (Vertices& vertices : mColorToVertices) {
        vertices.reserve(vertexCount);
    }
}

void GraphColoring::colorVertex(Vertex v, Color c)
{
    mVertexToColor[v] = c;
    if (mColorCount <= c) {
        mColorCount = c + 1;
    }
    mColorToVertices[c].push_back(v);
}

Color GraphColoring::getVertexColor(Vertex v) const
{
    return mVertexToColor[v];
}

Vertices& GraphColoring::getVerticesByColor(Color c)
{
    return mColorToVertices[c];
}

size_t GraphColoring::getColorCount() const
{
    return mColorCount;
}

void GraphColoring::clear()
{
    std::fill(mVertexToColor.begin(), mVertexToColor.end(), mVertexToColor.size());
    for (Vertices& vertices : mColorToVertices) {
        vertices.clear();
    }
    mColorCount = 0;
}

const size_t iterationsNeeded = 30;
const size_t numSoultionsNeeded = 25;

size_t ColoringSolver::solve(const Graph& graph, std::vector<GraphColoring>& out_solutions)
{
    out_solutions.clear();
    out_solutions.reserve(numSoultionsNeeded+1);

    GraphColoring rlfSolution(graph.GetVertexCount());
    size_t rlfColorsNum = solveRlf(graph, rlfSolution, graph.GetVertexSet());
    mAvailableColors.resize(graph.GetVertexCount());
    for (Vertex v = 0; v < graph.GetVertexCount(); ++v) {
        mAvailableColors[v].resize(graph.GetVertexCount(), true);
    }

    std::random_device rd;
    std::default_random_engine randEng(rd());
    Vertices walkOrder(graph.GetVertexCount());
    size_t n = 0;
    std::generate(walkOrder.begin(), walkOrder.end(), [&n]() {return n++; });

    std::multimap<size_t, size_t> colorCountToSolutionIdx;
    GraphColoring colors(graph.GetVertexCount());
    for (size_t iterIdx = 0; iterIdx < iterationsNeeded; ++iterIdx) {
        std::shuffle(walkOrder.begin(), walkOrder.end(), randEng);
        colors.clear();

        size_t upperBound = out_solutions.size() < numSoultionsNeeded ? 0 : colorCountToSolutionIdx.rbegin()->first;
        if (!solveIter(graph, walkOrder, colors, upperBound)) {
            continue;
        }

        if (out_solutions.size() < numSoultionsNeeded) {
            out_solutions.push_back(colors);
            colorCountToSolutionIdx.emplace(colors.getColorCount(), out_solutions.size() - 1);
        }
        else {
            auto worstSolutionPairIt = colorCountToSolutionIdx.rbegin();
            size_t worstSolutionSize = worstSolutionPairIt->first;
            size_t worstSolutionIdx = worstSolutionPairIt->second;
            if (colors.getColorCount() < worstSolutionSize) {
                out_solutions[worstSolutionIdx] = colors;
                ++worstSolutionPairIt;
                colorCountToSolutionIdx.erase(worstSolutionPairIt.base());
                colorCountToSolutionIdx.emplace(out_solutions[worstSolutionIdx].getColorCount(), worstSolutionIdx);
            }
        }
    }

    out_solutions.push_back(rlfSolution);
    return colorCountToSolutionIdx.begin()->first;
}

size_t ColoringSolver::solveIter(const Graph& graph, const Vertices& walkOrder, GraphColoring& out_colors, size_t upperBound)
{
    size_t maxColorUsed = 0;
    size_t vertexCount = graph.GetVertexCount();

    for (Vertex v = 0; v < graph.GetVertexCount(); ++v) {
        std::fill(mAvailableColors[v].begin(), mAvailableColors[v].end(), true);
    }

    for (size_t vertexIdx = 0; vertexIdx < walkOrder.size(); ++vertexIdx) {
        Vertex v = walkOrder[vertexIdx];

        Color minMissingColor = std::distance(mAvailableColors[v].begin(), std::find(mAvailableColors[v].begin(), mAvailableColors[v].end(), true));

        out_colors.colorVertex(v, minMissingColor);
        if (minMissingColor > maxColorUsed) {
            maxColorUsed = minMissingColor;
            if (upperBound && (minMissingColor >= upperBound)) {
                return 0;
            }
        }

        for (Vertex neighbor = 0; neighbor < vertexCount; ++neighbor) {
            if (graph.HasEdge(v, neighbor)) {
                mAvailableColors[neighbor][minMissingColor] = false;
            }
        }
    }

    return maxColorUsed + 1;
}

template <typename T, typename Iter>
void fastRemove(std::vector<T>& vec, Iter it)
{
    *it = vec.back();
    vec.pop_back();
}

size_t ColoringSolver::solveRlf(const Graph& graph, GraphColoring & out_solution, const Vertices& sub_set)
{
    Color newColor = 0;
    size_t coloredCount = 0;
    size_t vertexCount = sub_set.size();
    std::vector<Vertex> uncolored(vertexCount);
    {
        size_t n = 0;
        std::generate(uncolored.begin(), uncolored.end(), [&n]() {return n++; });
    }
    std::vector<Vertex> colorable(vertexCount);

    std::vector<size_t> uncoloredNeighborsCount(vertexCount, 0);
    while (coloredCount < vertexCount) {
        std::fill(uncoloredNeighborsCount.begin(), uncoloredNeighborsCount.end(), 0);
        colorable = uncolored;
        while (!colorable.empty()) 
        {
            for (size_t v1Idx = 0; v1Idx < colorable.size(); ++v1Idx) {
                Vertex v1 = colorable[v1Idx];
                const std::vector<uint8_t>& adj = graph.GetAdjacencyMatrixRow(sub_set[v1]);
                for (size_t v2Idx = v1Idx + 1; v2Idx < colorable.size(); ++v2Idx) {
                    Vertex v2 = colorable[v2Idx];
                    if (adj[sub_set[v2]]) {
                        ++uncoloredNeighborsCount[v1];
                        ++uncoloredNeighborsCount[v2];
                    }
                }
            }
            auto vIt = std::max_element(colorable.begin(), colorable.end(), [&uncoloredNeighborsCount](Vertex v1, Vertex v2) { return uncoloredNeighborsCount[v1] < uncoloredNeighborsCount[v2]; });
            Vertex v = *vIt;
            out_solution.colorVertex(v, newColor);
            for (size_t v2Idx = colorable.size() - 1; v2Idx >= 0; --v2Idx) {
                if (graph.HasEdge(sub_set[v], sub_set[colorable[v2Idx]]) || v == colorable[v2Idx]) {
                    fastRemove(colorable, colorable.begin() + v2Idx);
                }
                if (!v2Idx)
                    break;
            }
            fastRemove(uncolored, std::find(uncolored.begin(), uncolored.end(), v));
            ++coloredCount;
        }
        ++newColor;
    }
    return out_solution.getColorCount();
}


