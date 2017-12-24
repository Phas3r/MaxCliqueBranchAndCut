#include "graph.h"

#include <sstream>
#include <algorithm>

size_t Graph::GetVertexCount() const
{
    return m_adjacency_matrix.size();
}

void Graph::SetVertexSet(Vertices arr)
{
    m_vertices = arr;
}

const Vertices& Graph::GetVerticesOrderedByDegree() const
{
    if (m_vertices_by_degree.empty()) {
        m_vertices_by_degree.resize(GetVertexCount());
        size_t n = 0;
        std::generate(m_vertices_by_degree.begin(), m_vertices_by_degree.end(), [&n]() {return n++; });
        std::vector<size_t> vertexDegrees(GetVertexCount());
        for (size_t v1 = 0; v1 < GetVertexCount(); ++v1) {
            for (size_t v2 = 0; v2 < GetVertexCount(); ++v2) {
                if (v1 == v2 || !HasEdge(v1, v2)) {
                    continue;
                }
                ++vertexDegrees[v1];
            }
        }
        std::sort(m_vertices_by_degree.begin(), m_vertices_by_degree.end(),
            [&vertexDegrees](const Vertex& first, const Vertex& second) { return vertexDegrees[first] > vertexDegrees[second]; });
    }
    return m_vertices_by_degree;
}

void Graph::AddEdge(Vertex v1, Vertex v2)
{
    if (HasEdge(v1, v2) || (v1 == v2))
        return;
    m_vertices_by_degree.clear();
    m_adjacency_matrix[v1][v2] = true;
    m_adjacency_matrix[v2][v1] = true;
}

void Graph::EraseEdge(size_t v1, size_t v2)
{
    if (!HasEdge(v1, v2))
        return;
    m_vertices_by_degree.clear();
    m_adjacency_matrix[v1][v2] = false;
    m_adjacency_matrix[v2][v1] = false;
}

bool Graph::HasEdge(size_t v1, size_t v2) const
{
    return !!m_adjacency_matrix[v1][v2];
}

const Graph::AdjacencyMatrixRow& Graph::GetAdjacencyMatrixRow(Vertex v) const
{
    return m_adjacency_matrix[v];
}

const std::set<Vertex>& Graph::GetNotNeighbours(Vertex v) const
{
    return not_neighbours[v];
}

void Graph::SetAdjacencyMatrixRow(Vertex v, const AdjacencyMatrixRow& row)
{
    m_adjacency_matrix[v] = row;
    for (size_t v2 = 0; v2 < GetVertexCount(); ++v2) {
        if (v2 == v) {
            continue;
        }
        m_adjacency_matrix[v2][v] = row[v2];
    }
}

void Graph::InsertVertex(Vertex v, AdjacencyMatrixRow& adjacency)
{
    m_vertices_by_degree.clear();
    m_adjacency_matrix.insert(m_adjacency_matrix.begin() + v, adjacency);
    for (size_t v2 = 0; v2 < GetVertexCount(); ++v2) {
        if (v2 == v) {
            continue;
        }
        m_adjacency_matrix[v2].insert(m_adjacency_matrix[v2].begin() + v, adjacency[v2]);
    }
}

void Graph::EraseVertex(Vertex v)
{
    m_vertices_by_degree.clear();
    for (AdjacencyMatrixRow& row : m_adjacency_matrix) {
        row.erase(row.begin() + v);
    }
    m_adjacency_matrix.erase(m_adjacency_matrix.begin() + v);
}

void Graph::GenerateNotNeighbours()
{
    not_neighbours.clear();
    not_neighbours.resize(m_vertices.size());
    for (size_t i = 0; i < m_vertices.size(); ++i)
    {
        std::set<Vertex>& not_neighnours_set = not_neighbours[i];
        for (size_t j = 0; j < m_vertices.size(); ++j)
        {
            if (i == j)
                continue;
            if (!HasEdge(i, j))
                not_neighnours_set.insert(j);
        }
    }
    
}

const std::string Graph::ToString()
{
    std::stringstream sstream;
    size_t edgeCount = 0;
    for (const AdjacencyMatrixRow& row : m_adjacency_matrix) {
        size_t rowC = 0;
        for (uint8_t hasEdge : row) {
            if (hasEdge) {
                ++rowC;
                ++edgeCount;
            }
        }
        sstream << rowC << std::endl;
    }
    // Edge count is divided by 2 because each edge (v1, v2) is encountered 2 times: in neighbors of v1 and in neighbors of v2, but we only output it once
    sstream << "p " << GetVertexCount() << ' ' << edgeCount / 2 << std::endl;
    for (size_t v1 = 0; v1 < GetVertexCount(); ++v1) {
        Vertices adjacentVertices;
        for (size_t v2 = v1; v2 < GetVertexCount(); ++v2) {
            if (!HasEdge(v1, v2)) {
                continue;
            }
            adjacentVertices.push_back(Vertex(v2));
        }
        std::sort(adjacentVertices.begin(), adjacentVertices.end());
        for (Vertex v2 : adjacentVertices) {
            // Vertex labels start at 1
            sstream << "e " << v1 + 1 << " " << v2 + 1 << std::endl;
        }
    }
    return sstream.str();
}
