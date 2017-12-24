#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <map>
#include "graph.h"
#include "windows.h"
using namespace std;

Graph readDimacs(const char* filename)
{
    Graph result;

    ifstream ifs(filename);
    if (!ifs.is_open()) {
        cout << "0" << endl;
        cout<< "WARNING: Could not open the file"<<endl;
        exit(1);
    }
    string line;
    size_t line_number = 0;
    size_t edgeLineCount = 0;
    size_t edgeLineNum = 0;
    size_t vertexCount = 0;
    while (getline(ifs, line)) {
        ++line_number;
        switch (line[0]) {
        case 'p':
        {
            stringstream liness(line);
            // ignore 'p'
            liness.ignore(1, ' ');
            string format;
            liness >> format >> vertexCount >> edgeLineCount;
            if (format != "edge" && format != "col") {
#ifdef DEBUG
                cout << "Header on line " << line_number << " uses \"" << format << "\" format, but this solution supports only \"edge\" or \"col\" format is supported";
#endif
                return result;
            }
            result = Graph(vertexCount);
            break;
        }
        case 'e':
        {
            ++edgeLineNum;
            if (edgeLineNum > edgeLineCount) {
                continue;
            }
            std::stringstream liness(line);
            // ignore 'e'
            liness.ignore(1, ' ');
            size_t v1;
            size_t v2;
            liness >> v1 >> v2;
            // Vertex labels start at 1, but we need them to start from 0
            --v1;
            --v2;
            if (v1 >= vertexCount || v2 >= vertexCount) {
#ifdef DEBUG
                cout << "Edge on line " << line_number << " has invalid vertex number";
#endif
                return result;
            }
            if (v1 != v2)
                result.AddEdge(Vertex(v1), Vertex(v2));
            break;
        }

        case 'c':
            // Ignore comment lines
            continue;
        default:
        {
#ifdef DEBUG
            cout << "Failed to parse line " << line_number;
#endif
            return result;
            break;
        }
        }
    }
    if (edgeLineNum != edgeLineCount) {

#ifdef DEBUG
        cout << "Invalid edge count. Got " << edgeLineNum << ", expected " << edgeLineCount;
#endif
        return result;
    }
    return result;
}
