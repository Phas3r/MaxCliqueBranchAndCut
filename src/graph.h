#pragma once

#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

typedef size_t Vertex;

typedef std::vector<Vertex> Vertices;


class Graph
{
public:
	typedef std::vector<uint8_t> AdjacencyMatrixRow;

	Graph() = default;
	Graph(size_t num_vertices)
		: m_adjacency_matrix(num_vertices)
	{
        m_vertices.resize(num_vertices);
        {
            size_t n = 0;
            std::generate(m_vertices.begin(), m_vertices.end(), [&n]() {return n++; });
        }
		for (auto& row : m_adjacency_matrix) {
			row.resize(num_vertices);
		}
	};

	size_t GetVertexCount() const;
    const Vertices GetVertexSet() const { return m_vertices; }
    void SetVertexSet(Vertices arr);
	const Vertices& GetVerticesOrderedByDegree() const;

	void AddEdge(Vertex v1, Vertex v2);
	void EraseEdge(size_t v1, size_t v2);
	bool HasEdge(size_t v1, size_t v2) const;

	const AdjacencyMatrixRow& GetAdjacencyMatrixRow(Vertex v) const;
    const std::set<Vertex>& GetNotNeighbours(Vertex v) const;
	void SetAdjacencyMatrixRow(Vertex v, const AdjacencyMatrixRow& row);
	
	void InsertVertex(Vertex v, AdjacencyMatrixRow& adjacency);
	void EraseVertex(Vertex v);
    void GenerateNotNeighbours();
	// Outputs graph in DIMACS format
	const std::string ToString();
private:
	std::vector<AdjacencyMatrixRow> m_adjacency_matrix;
	mutable Vertices m_vertices_by_degree;
    Vertices m_vertices;
    std::vector<std::set<Vertex> > not_neighbours;
};
