#pragma once

#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include <memory>

using namespace std;

// ============ Tree nodes ============
struct UPGMA_Node {
    string name;
    int    level;   // merge level, leaf=-1, incremented each round
};

// ============ Edges ============
struct UPGMA_Edge {
    UPGMA_Node child;          // child node
    UPGMA_Node parent;         // parent node
    double     branch_length;  // branch length from child to parent = parent_height - child_accumulated_height
};

// ============ Tree building ============
struct UPGMA_Tree {
    vector<UPGMA_Node>              nodes;             // all nodes, last one is root
    vector<UPGMA_Edge>              edges;             // all merge edges
    vector<vector<double>>          distance_matrix;   // working distance matrix (shrunk during tree build)
    vector<string>                  cluster_names;     // current cluster names (updated during merge)
};

vector<string> clean_chain_names(const vector<string>& names);
vector<vector<double>> build_distance_matrix(const vector<vector<double>>& tm_mat, int chain_num);
unique_ptr<UPGMA_Tree> build_upgma_tree(const vector<string>& names, const vector<vector<double>>& dist_mat);
string get_tree_string(const UPGMA_Tree* tree, int precision = 5);
void print_tree(const UPGMA_Tree* tree, ostream& os = cout);
bool save_svg(const UPGMA_Tree* tree, const string& filename);
void save_distance_matrix(const vector<string>& names, const vector<vector<double>>& dist_mat, int chain_num);
void output_upgma_tree(const vector<string>& names, const vector<vector<double>>& tm_mat, int chain_num);
