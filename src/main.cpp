#include <cstdio>
#include <iostream>
#include <string>
#include <map>
#include <thread>
#include <chrono>
#include <algorithm>
#include "graph.h"
#include "simple_solver.h"
#include "cliquer_wrapper.h"
#include "utils.h"


void StartCheck(Graph graph)
{
        auto t1 = std::chrono::high_resolution_clock::now();

        auto result = CliquerWrapper::unweightedFindMax(graph);

        auto t2 = std::chrono::high_resolution_clock::now();

        // integral duration: requires duration_cast
        auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        

        // fractional duration: no duration_cast needed

        std::cout << "Cliquer took " << int_ms.count() << " whole milliseconds\n";
        cout << "\nRight Answer:" << endl;
        cout << "Max Clique Size:"<< result.size()<<endl;
        cout << "Max Clique:";
        for (auto vertex : result)
            cout << vertex << "\t";
        cout << "\n";
}

void StartTimeoutWatchdog(uint32_t timeout_secs)
{
    std::thread([timeout_secs]() {
        std::this_thread::sleep_for(std::chrono::seconds(timeout_secs));
        std::cout << "0" << std::endl;
        if (!timeout_secs)
            cout << "WARNING: Timeout limit was zero or was not provided!" << endl;
        exit(1);
    }).detach();
}

int main(int argc, char* argv[])
{
    uint32_t timeout = 0;
    try
    {
        timeout = std::stoul(argv[2]);
    }
    catch (std::invalid_argument e) {
        cout << "WARNING: Could not get timeout limit!" << endl;
        return 1;
    }
    StartTimeoutWatchdog(timeout);;
    Graph in_graph = readDimacs(argv[1]);
    Graph out_graph(in_graph.GetVertexCount());
    // only use this for debug
    //StartCheck(in_graph);
    
    SimpleSolver solver;

    auto t1 = std::chrono::high_resolution_clock::now();
    solver.solve(in_graph, out_graph);
    Vertices clique = out_graph.GetVertexSet();

    auto t2 = std::chrono::high_resolution_clock::now();

    // integral duration: requires duration_cast
    auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    
    //std::cout << "Program took " << int_ms.count() << " whole milliseconds\n";
    cout << clique.size() << endl;
    std::sort(clique.begin(), clique.end());
    for (auto vertex : clique)
    {
        cout << vertex + 1 << " ";
    }
    return 0;
}
