set(SOURCES
    main.cpp
    graph.h
    graph.cpp
    utils.h
    simple_solver.h
    simple_solver.cpp
    cliquer_wrapper.h
    cliquer_wrapper.cpp
    coloring_problem.h
    coloring_problem.cpp
)

find_package(Threads REQUIRED)

include_directories(${CMAKE_3RD_PARTY_DIR}/cplex/include)

link_directories(${CMAKE_3RD_PARTY_DIR}/cplex/lib/)

add_executable(SlashchininMaxClique ${SOURCES})

target_link_libraries(SlashchininMaxClique cliquer ${CMAKE_THREAD_LIBS_INIT} cplex1271.lib)
