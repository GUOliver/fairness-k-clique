/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc graph.c -O9 -o graph".

To execute:
"./graph k edgelist.txt attribute.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
"attribute.txt" should contain the attribute along with the vertice
Will print the number of k-cliques.
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

// maximum number of edges for memory allocation, will increase if needed
#define NLINKS 100000000

typedef struct {
	unsigned s;
	unsigned t;
} edge;

// typedef struct {
// 	unsigned node;
// 	unsigned deg;
// } nodedeg;

typedef struct {
	unsigned n; // # of nodes
	unsigned e; // # of edges
	edge* edges; // array of edges: size of e

    unsigned* attributes; // the attribute dynamic array for nodes: size of n
    unsigned attr_dimension; // attribute set dimension in graph, mainly will be 2D

	unsigned* ns; // ns[l]: number of nodes in G_l: sie
	unsigned** d; // d[l]: degrees of G_l
	unsigned* cd; // cumulative degree: (starts with 0) length = n+1
	unsigned* adj; // truncated list of neighbors
	unsigned* rank; // ranking of the nodes according to degeneracy ordering
	//unsigned *map;// oldID newID correspondance

	unsigned char* lab; // lab[i] label of node i
	unsigned** sub; //sub[l]: nodes in G_l

	int** resultSet; // results containing the cliques
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
    g->attributes = malloc(g->n * sizeof(int));
    FILE* f = fopen(attr_file, "r");
    if (f == NULL) {
        printf("Cannot open the attribute file ! \n");
        exit(1);
    }

    int curr_node = 0;
    int curr_attr = 0;
    g->attr_dimension = 0;
    while (fscanf(f, "%d %d", &curr_node, &curr_attr) == 2) {
        if (curr_node >=  g->n) {
            printf("Erroneous attribute file with node %d\n", curr_node);
            exit(1);
        }
        g->attributes[curr_node] = curr_attr;
        if (curr_attr > g->attr_dimension) {
            g->attr_dimension = curr_attr;
        }
    }
    g->attr_dimension++;
    fclose(f);
    printf("Finished reading attribute file\n");
    return;
}

// read edgeList from file
graph* read_edgelist(char* edgelist, char* attr_file){
	unsigned e1 = NLINKS;	// maximum number of edges, will increase if needed
	graph* g = malloc(sizeof(graph));
	FILE *file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist,"r");
	g->edges = malloc(e1 * sizeof(edge));

	while (fscanf(file,"%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) { //Add each edge
		g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
		g->e++;
		if (g->e == e1) {
			e1 += NLINKS;
			g->edges=realloc(g->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	g->n++;	// node index starts from 0, thus +1 needed
	g->edges = realloc(g->edges, g->e * sizeof(edge));	// realloc the space to accommodate g->e(updated) # of edges !
    
	read_attribute_file(g, attr_file); // fill in attributes array in graph
	return g;
}

// relabled as DAG 
void relabel(graph* g){
	unsigned source, target, tmp;
	for (unsigned i = 0; i < g->e; i++) {
		source = g->rank[g->edges[i].s];
		target = g->rank[g->edges[i].t];
		if (source < target){
			tmp = source;
			source = target;
			target = tmp;
		}
		g->edges[i].s = source;
		g->edges[i].t = target;
	}
}

///// CORE ordering /////////////////////
typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

// Building the heap structure with (key, value) = (node, degree) for each node
typedef struct {
	unsigned n_max;	// max number of nodes
	unsigned n;	// number of nodes
	unsigned* pt;	// pointers to nodes
	keyvalue* kv; // nodes
} bheap;


bheap* construct(unsigned n_max){
	bheap* heap=malloc(sizeof(bheap));
	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = malloc(n_max * sizeof(unsigned));
	for (unsigned i = 0; i < n_max; i++) {
		(heap->pt)[i] = -1;
	}
	heap->kv = malloc(n_max * sizeof(keyvalue));
	return heap;
}

void swap(bheap* heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

void bubble_up(bheap* heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
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
		j=( (j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap* heap, keyvalue kv){
	(heap->pt)[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap,heap->n - 1);
}

void update(bheap* heap, unsigned key){
	unsigned i = heap->pt[key];
	if (i != -1){
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

keyvalue pop_min(bheap* heap){
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

// Building the heap structure with (key, value) = (node, degree) for each node
bheap* mk_heap(unsigned n, unsigned* degrees){
	keyvalue kv;
	bheap* heap = construct(n);
	for (unsigned i = 0; i < n; i++){
		kv.key = i;
		kv.value = degrees[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap* heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

// computing degeneracy ordering and core value for each vertice
void ord_core(graph* g){
	unsigned i, j, r = 0;
	unsigned n = g->n;
	keyvalue kv;
	bheap* heap;

	unsigned* d0 = calloc(g->n, sizeof(unsigned));	// allocate space and initialize all to 0
	unsigned* cd0 = malloc((g->n + 1) * sizeof(unsigned));	// ???
	unsigned* adj0 = malloc(2 * (g->e) * sizeof(unsigned));	// store the neighbors of each vertices, thus need 2E space

	// fill in the degrees of each vertice in undirected graph
	for (i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}
	
	cd0[0] = 0;
	for (i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i-1] + d0[i-1];
		d0[i-1] = 0;
	}

	for (i = 0;i < g->e; i++) {
		adj0[ cd0[g->edges[i].s] + d0[ g->edges[i].s ]++ ] = g->edges[i].t;
		adj0[ cd0[g->edges[i].t] + d0[ g->edges[i].t ]++ ] = g->edges[i].s;
	}

	heap = mk_heap(n, d0);

	g->rank = malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++){
		kv = pop_min(heap);
		g->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++){
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
	
	unsigned* d = calloc(g->n, sizeof(unsigned));

	for (unsigned i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
	}

	g->cd = malloc((g->n + 1) * sizeof(unsigned));
	unsigned ns = 0;
	g->cd[0] = 0;
	unsigned max = 0;
	unsigned* sub = malloc(g->n * sizeof(unsigned));
	unsigned char* lab = malloc(g->n * sizeof(unsigned char));
	for (unsigned i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}
	printf("max degree = %u\n",max);
	fflush(stdout);

	g->adj=malloc(g->e*sizeof(unsigned));

	for (unsigned i = 0; i < g->e; i++) {
		g->adj[ g->cd[g->edges[i].s] + d[ g->edges[i].s ]++ ]=g->edges[i].t;
	}
	free(g->edges);

	g->ns=malloc((k+1)*sizeof(unsigned));
	g->ns[k]=ns;

	g->d=malloc((k+1)*sizeof(unsigned*));
	g->sub=malloc((k+1)*sizeof(unsigned*));
	for (unsigned i = 2; i < k; i++){
		g->d[i] = malloc(g->n*sizeof(unsigned));
		g->sub[i] = malloc(max*sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
}


// n stored the number of k-cliques
void kclique(unsigned l, graph* g, unsigned long long* n) {
	unsigned end, u, v, w;


	if (l == 2){
		for (unsigned i = 0; i < g->ns[2]; i++) { // loop through all the nodes in the subgraph
			u = g->sub[2][i];
			printf("u is %d\n", u);
			// printf("cd[u] is %d\n", g->cd[u]);
			// printf("d[2][u] is %d\n", g->d[2][u]);
			//(*n)+=g->d[2][u];
			end = g->cd[u] + g->d[2][u];
			// printf("end is %d\n", end);
			for (unsigned j = g->cd[u]; j < end; j++) {
				//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
				(*n)++; 
			}
		}
		return;
	}

	for (unsigned i = 0; i < g->ns[l]; i++){
		u = g->sub[l][i];
		printf("%u %u\n",i ,u);
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];

		for (unsigned j = g->cd[u]; j < end; j++){ // relabeling nodes and forming U
			v = g->adj[j];
			//if (g->lab[v]==l){
			g->lab[v] = l - 1;
			g->sub[l - 1][g->ns[l - 1]++] = v;
			g->d[l-1][v] = 0;//new degrees
			//}
		}

		for (unsigned j = 0; j < g->ns[l-1]; j++){ //reodering adjacency list and computing new degrees
			v=g->sub[l-1][j];
			end=g->cd[v]+g->d[l][v];
			for (unsigned k = g->cd[v]; k < end; k++){
				w = g->adj[k];
				if (g->lab[w] == l-1){
					g->d[l-1][v]++;
				}
				else{
					g->adj[k--] = g->adj[--end];
					g->adj[end] = w;
				}
			}
		}

		kclique(l-1, g, n);
		
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

	// compute core number of each vertice in G, 𝑘(𝑣) = the largest k such that the k-core contains v
	ord_core(g);	
	// relabled as DAG according to k-core
	relabel(g);

	mkspecial(g, k);

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
	t1 = t2;

	printf("Iterate over all cliques\n");
	fflush(stdout);

	n = 0;
	kclique(k, g, &n);	 // list all k-cliques 

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2-t1)/3600, ((t2-t1)%3600)/60, ((t2-t1)%60));
	fflush(stdout);
	t1 = t2;

	free_graph(g, k);

	printf("- Overall time = %ldh%ldm%lds\n", (t2-t0)/3600, ((t2-t0)%3600)/60, ((t2-t0)%60));
	fflush(stdout);
	
	return 0;
}
