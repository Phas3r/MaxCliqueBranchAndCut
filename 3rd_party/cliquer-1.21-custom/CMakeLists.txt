project (cliquer C)

set(SOURCES
	cliquer.h
	cliquer.c
	cliquerconf.h
	graph.c
	graph.h
	misc.h
	reorder.c
	reorder.h
	set.h
)

if(UNIX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
elseif(WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MT")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
	set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} /MT")
	set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} /MT")
endif()

add_library(cliquer ${SOURCES})