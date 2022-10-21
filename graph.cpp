/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc graph.c -O9 -o graph".

To execute:
	- Generate attribute_gen.txt by run attribute_gen.py
	- ./graph k edgelist.txt attribute.txt
	"edgelist.txt" should contain the graph: one edge on each line separated by a space.
	"attribute.txt" should contain the attribute along with the vertice
*/
#include <cstdio>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <set>

// maximum number of edges for memory allocation, will increase if needed
#define NLINKS 100000000

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned n; // # of nodes
	unsigned e; // # of edges
	edge* edges;

    unsigned* attributes; // the attribute array for nodes -> size of n, attributes[i] = A
    unsigned attr_dimension; // attribute set dimension in graph, mainly 2D

	unsigned* ns; // ns[l]: number of nodes in G_l -> size of k+1
	unsigned** sub; //sub[l]: nodes in G_l	-> size of k+1
	unsigned** d; // d[l]: degrees of G_l
	unsigned* cd; // cumulative degree: (starts with 0) length = n+1
	unsigned* adj; // truncated list of neighbors
	unsigned* rank; // ranking of the nodes according to degeneracy ordering -> size of n
	//unsigned *map;// oldID newID correspondance

	unsigned char* lab; // lab[i]: label of node i
	std::vector<std::set<unsigned>> res;	// store the all the clique results
} graph;

void free_graph(graph* g, unsigned char k){
	unsigned char i;
	free(g->ns);
	for (i=2;i<k+1;i++){
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->cd);
	free(g->adj);
	free(g);
}

// Compute the maximum of three unsigned integers.
unsigned max3(unsigned int a,unsigned int b,unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

// read attributes of vertices from file
void read_attribute_file(graph* g, char* attr_file) {
    g->attributes = (unsigned*)malloc(g->n * sizeof(unsigned));
    FILE* f = fopen(attr_file, "r");
    if (f == NULL) {
        printf("Cannot open the attribute file ! \n");
        exit(1);
    }

    int curr_node = 0;
    int curr_attr = 0;
    g->attr_dimension = 0;
	printf("Totol vertices count: %d\n", g->n);
    while (fscanf(f, "%d %d", &curr_node, &curr_attr) == 2) {
        if (curr_node >=  g->n) {
            printf("Erroneous attribute file with vertice %d\n", curr_node);
            exit(1);
        }
        g->attributes[curr_node] = curr_attr;
        if (curr_attr > g->attr_dimension) {
            g->attr_dimension = curr_attr;
        }
    }
    g->attr_dimension++;
    fclose(f);
    printf("Attribute dimension: %d, Finished reading attribute file\n", g->attr_dimension);
    return;
}

// read edgeList from file
graph* read_edgelist(char* edgelist, char* attr_file){
	unsigned e1 = NLINKS;	// maximum number of edges, will increase if needed
	graph* g = (graph*)malloc(sizeof(graph));
	FILE* file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist,"r");
	g->edges = (edge*)malloc(e1 * sizeof(edge));

	while (fscanf(file,"%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) { //Add each edge
		g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
		g->e++;
		if (g->e == e1) {
			e1 += NLINKS;
			g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	g->n++;	// node index starts from 0, thus +1 needed
	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));	// realloc the space to accommodate g->e(updated) # of edges !
    printf("In reading edgelist, total vertices count: %d\n", g->n);
	read_attribute_file(g, attr_file); // fill in attributes array in graph
	return g;
}

void relabel(graph* g){
	unsigned src_rank;
	unsigned target_rank;
	unsigned tmp;
	
	for (unsigned i = 0; i < g->e; i++) {
		src_rank = g->rank[g->edges[i].s];
		target_rank = g->rank[g->edges[i].t];
		if (src_rank < target_rank){	// make sure src_rank always greater than target rank
			tmp = src_rank;
			src_rank = target_rank;
			target_rank = tmp;
		}
		// relabel vertice by their rank
		g->edges[i].s = src_rank;
		g->edges[i].t = target_rank;
	}
}

///// CORE ordering /////////////////////
typedef struct {
	unsigned id;
	unsigned degree;
} vertice;

// Building the heap structure with (id, degree) = (node_index, degree) for each node
typedef struct {
	unsigned n_max;	// max number of nodes
	unsigned n;	// number of nodes
	unsigned* pt;	// pointers to nodes
	vertice* v; // nodes
} bheap;


bheap* construct(unsigned n_max){
	bheap* heap = (bheap*)malloc(sizeof(bheap));	// size of 24B, 4+4+8+8
	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned int*)malloc(n_max * sizeof(unsigned));
	for (unsigned i = 0; i < n_max; i++) {
		(heap->pt)[i] = -1;
	}
	heap->v = (vertice*)malloc(n_max * sizeof(vertice));
	return heap;
}

void swap(bheap* heap, unsigned i, unsigned j) {
	vertice v_tmp = heap->v[i];
	unsigned pt_tmp = heap->pt[v_tmp.id];
	heap->pt[heap->v[i].id] = heap->pt[heap->v[j].id];
	heap->v[i] = heap->v[j];
	heap->pt[heap->v[j].id] = pt_tmp;
	heap->v[j] = v_tmp;
}

void bubble_up(bheap* heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->v[j].degree > heap->v[i].degree) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

void bubble_down(bheap* heap) {
	unsigned i = 0, j1 = 1, j2 = 2;
	unsigned j;
	while (j1 < heap->n) {
		j=( (j2 < heap->n) && (heap->v[j2].degree < heap->v[j1].degree) ) ? j2 : j1 ;
		if (heap->v[j].degree < heap->v[i].degree) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap* heap, vertice v){
	(heap->n)++;
	(heap->pt)[v.id] = heap->n;
	heap->v[heap->n - 1] = v;
	bubble_up(heap, heap->n - 1);
}

void update(bheap* heap, unsigned id){
	unsigned i = (heap->pt)[id];
	if (i != -1) {
		((heap->v[i]).degree)--;
		bubble_up(heap, i);
	}
}

vertice pop_min(bheap* heap){
	vertice min = heap->v[0];
	heap->pt[min.id] = -1;
	heap->v[0] = heap->v[--(heap->n)];
	heap->pt[heap->v[0].id] = 0;
	bubble_down(heap);
	return min;
}

// Building the heap for each node
bheap* mk_heap(unsigned n, unsigned* degrees){
	vertice v;
	bheap* heap = construct(n);
	for (unsigned i = 0; i < n; i++){
		v.id = i;
		v.degree = degrees[i];
		insert(heap, v);
	}
	return heap;
}

void freeheap(bheap* heap){
	free(heap->pt);
	free(heap->v);
	free(heap);
}

// computing degeneracy ordering and core degree for each vertice
void ord_core(graph* g) {
	unsigned* d0 = (unsigned int*)calloc(g->n, sizeof(unsigned));	// record degree for each vertice
	unsigned* cd0 = (unsigned int*)malloc((g->n + 1) * sizeof(unsigned));	// cumulative degrees of all vertices
	unsigned* adj0 = (unsigned int*)malloc(2 * (g->e) * sizeof(unsigned));	// ??

	// update the degree of each vertice in undirected graph
	for (unsigned i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	
	// add up all the degrees 
	cd0[0] = 0;
	for (unsigned i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i-1] + d0[i-1];
	}

	for (unsigned i = 0; i < g->e; i++) {
		unsigned src = g->edges[i].s;
		adj0[ cd0[src] ] = src;
		unsigned dst = g->edges[i].t;
		adj0[ cd0[dst] ] = dst;
	}

	bheap* heap = mk_heap(g->n, d0);	// construct a min heap ordered by the degree
	g->rank = (unsigned int*)malloc(g->n * sizeof(unsigned));	// initialise rank array
	for (unsigned i = 0, r = 0; i < g->n; i++){
		vertice v = pop_min(heap); // pop out the vertice with min degree
		(g->rank)[v.id] = g->n - (++r);	// update new ranking for each vertice in descending order
		for (unsigned j = cd0[v.id]; j < cd0[v.id + 1]; j++){
			update(heap, adj0[j]);
		}
	}

	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
// Building the special graph structure
void mkspecial(graph* g, unsigned char k) {
	unsigned* d = (unsigned int*)calloc(g->n, sizeof(unsigned));

	for (unsigned i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
	}

	g->cd = (unsigned int*)malloc((g->n + 1) * sizeof(unsigned));
	unsigned ns = 0;
	g->cd[0] = 0;
	unsigned max = 0;
	unsigned* sub = (unsigned int*)malloc(g->n * sizeof(unsigned));
	unsigned char* lab = (unsigned char*)malloc(g->n * sizeof(unsigned char));
	for (unsigned i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}
	printf("max degree = %u\n",max);
	fflush(stdout);

	g->adj = (unsigned int*)malloc(g->e*sizeof(unsigned));

	for (unsigned i = 0; i < g->e; i++) {
		g->adj[ g->cd[g->edges[i].s] + d[ g->edges[i].s ]++ ]=g->edges[i].t;
	}
	free(g->edges);

	g->ns = (unsigned int*)malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	g->d = (unsigned int**)malloc((k + 1) * sizeof(unsigned*));
	g->sub = (unsigned int**)malloc((k + 1) * sizeof(unsigned*));

	for (unsigned i = 2; i < k; i++){
		g->d[i] = (unsigned int*)malloc(g->n * sizeof(unsigned));
		g->sub[i] = (unsigned int*)malloc(max * sizeof(unsigned));
	}

	g->d[k] = d;
	g->sub[k] = sub;
	g->lab = lab;
}


// n stored the number of k-cliques
void kclique(unsigned l, graph* g, unsigned long long* n, std::set<unsigned>& R) {
	unsigned end, u, v, w;
	if (l == 2){
		for (unsigned i = 0; i < g->ns[2]; i++) { // loop through all the nodes in the subgraph
			u = g->sub[2][i];
			end = g->cd[u] + g->d[2][u];
			// printf("u is %d, end is %d\n", u, end);
			for (unsigned j = g->cd[u]; j < end; j++) {
				// printf("u is %d, j is %d, g->adj[j] is %d\n", u, j, g->adj[j]);
				std::set<unsigned> temp = {u, g->adj[j]};
				std::set<unsigned> target_clique;
				std::set_union(R.begin(), R.end(), temp.begin(), temp.end(),
					std::inserter(target_clique, target_clique.begin()));
				g->res.push_back(target_clique);
				// NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only) !!
				(*n)++; 
			}
		}
		return;
	}

	// loop over the vertices of G_l
	for (unsigned i = 0; i < g->ns[l]; i++){
		u = g->sub[l][i]; // get the current vertice

		std::set<unsigned> tmp = {u};
		std::set<unsigned> current_clique;
		std::set_union(R.begin(), R.end(),
				tmp.begin(), tmp.end(),                  
				std::inserter(current_clique, current_clique.begin()));
		// printf("current vertice is %u\n", u);

		g->ns[l - 1] = 0;
		
		end = g->cd[u] + g->d[l][u];	// whole cumulative degree count from vertice 0 to u
		for (unsigned j = g->cd[u]; j < end; j++){ // relabeling nodes and forming U
			v = g->adj[j];
			//if (g->lab[v]==l){
			g->lab[v] = l - 1;	// mark v (neighbor of u) to l-1
			g->sub[l - 1][g->ns[l - 1]++] = v;
			g->d[l-1][v] = 0; //new degrees
			//}
		}

		for (unsigned j = 0; j < g->ns[l-1]; j++){ // reodering adjacency list and computing new degrees
			v = g->sub[l-1][j];
			end = g->cd[v] + g->d[l][v];
			for (unsigned k = g->cd[v]; k < end; k++){
				w = g->adj[k];
				if (g->lab[w] == l-1){
					g->d[l-1][v]++;
				} else{
					g->adj[k--] = g->adj[--end];
					g->adj[end] = w;
				}
			}
		}

		kclique(l - 1, g, n, current_clique);
		
		//restoring labels
		for (unsigned j = 0; j < g->ns[l-1]; j++){
			v = g->sub[l - 1][j];
			g->lab[v] = l;
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 4) {
		printf("Not enough args: ./graph k edgelist.txt attributes.txt \n");
		return 0;
	}

	graph* g; // declare a graph g
	unsigned char k = atoi(argv[1]);
	unsigned long long n;

	time_t t0,t1,t2;
	t1 = time(NULL);
	t0 = t1;
	
	printf("Reading graph from file %s in edgelist format\n", argv[2]);
	fflush(stdout);

	printf("Reading attribute file from file %s\n", argv[3]);
	fflush(stdout);

	g = read_edgelist(argv[2], argv[3]);
	fflush(stdout);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
	t1 = t2;

	printf("Building the graph structure\n");
	fflush(stdout);

	//compute core number of each vertice in G, 𝑘(𝑣) = the largest k such that the k-core contains v
	ord_core(g);
	
	// relabled as DAG according to k-core
	// relabel(g);

	mkspecial(g, k);

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
	t1 = t2;

	printf("Iterate over all cliques\n");
	fflush(stdout);

	n = 0;
	std::set<unsigned> R;
	kclique(k, g, &n, R);	 // list all k-cliques
	std::cout << "result size is " << g->res.size() << std::endl;

	for (auto& i : g->res) {
		for (auto& j : i) {
			std::cout << j <<" ";
		}
		std::cout << "\n";
	}

	printf("Number of %u-cliques: %llu\n\n", k, n);

	printf("Start to Filter out unfair clique...\n");

	int threshold = 2;
	int count0 = 0;
	int count1 = 1;
	std::vector<std::set<unsigned>> finalRes;
	for (auto& clique : g->res) {
		count0 = 0;
		count1 = 0;
		for (auto& i : clique) {
			if (g->attributes[i] == 0) {
				count0++;
			} else {
				count1++;
			}
		}
		if (count0 >= threshold && count1 >= threshold) {
			finalRes.push_back(clique);
		}
	}

	for (auto& i : finalRes) {
		for (auto& j : i) {
			std::cout << j <<" ";
		}
		std::cout << "\n";
	}

	printf("Number of Weak fair %u-cliques: %llu\n\n", k, finalRes.size());

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
	fflush(stdout);
	t1 = t2;

	free_graph(g, k);

	printf("- Overall time = %ldh%ldm%lds\n", (t2-t0)/3600, ((t2-t0)%3600)/60, ((t2-t0)%60));
	fflush(stdout);
	
	return 0;
}
