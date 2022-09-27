#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

class graph {
private:
    int n;
    int m;
    int attr_size; // number of attributes
    int threshold; // given
    int max_color;
    std::string dir;
    int* offset;    // index of nodes in edge_list
    int* edge_list; // store all neighbors
    int* attribute; // the attribute of node
    int* pend;
    int* nei_vis;
    int* fairness_d;

    long long now_mem = 0;
    long long max_mem = 0;
    long long max_local_mem = 0;
    long long now_local_mem = 0;
    std::vector<int> left;
    std::vector<int> component;
    std::vector<int> union_X;
    std::vector<std::vector<int>> ResultSet;

    int* peeling_idx;
    int* color;

    int* X_vis;
    int* nvis;
    int* rvis;

    bool debug_is_on = false;

public:
    graph(const char* dir);
    ~graph();
};