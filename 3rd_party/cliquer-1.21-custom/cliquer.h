#ifndef CLIQUER_H
#define CLIQUER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>

#include "set.h"
#include "graph.h"
#include "reorder.h"

typedef struct _clique_options clique_options;
struct _clique_options {
	int *(*reorder_function)(graph_t *, boolean);
	int *reorder_map;

	FILE *output;

	boolean (*user_function)(set_t,graph_t *,clique_options *);
	void *user_data;
	set_t *clique_list;
	int clique_list_length;
};

extern clique_options *clique_default_options;

/* Weighted clique functions */
extern float clique_max_weight(graph_t *g,clique_options *opts);
extern set_t clique_find_single(graph_t *g,float min_weight, float max_weight,
				boolean maximal, clique_options *opts);
extern int clique_find_all(graph_t *g, float min_weight, float max_weight,
				boolean maximal, clique_options *opts);

/* Unweighted clique functions */
#define clique_unweighted_max_size clique_unweighted_max_weight
extern float clique_unweighted_max_weight(graph_t *g, clique_options *opts);
extern set_t clique_unweighted_find_single(graph_t *g, float min_size,
					   float max_size,boolean maximal,
					   clique_options *opts);
extern int clique_unweighted_find_all(graph_t *g, float min_size, float max_size,
					boolean maximal, clique_options *opts);

/* Alternate spelling (let's be a little forgiving): */
#define cliquer_options clique_options
#define cliquer_default_options clique_default_options

#ifdef __cplusplus
}
#endif

#endif /* !CLIQUER_H */