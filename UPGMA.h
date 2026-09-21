#pragma once

#include <string>
#include <vector>
#include <ostream>
#include <iostream>
#include <memory>

// ============ Tree nodes ============
struct UPGMA_Node {
    std::string name;
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
    std::vector<UPGMA_Node>              nodes;             // all nodes, last one is root
    std::vector<UPGMA_Edge>              edges;             // all merge edges
    std::vector<std::vector<double>>          distance_matrix;   // working distance matrix (shrunk during tree build)
    std::vector<std::string>                  cluster_names;     // current cluster names (updated during merge)
};

std::vector<std::string> clean_chain_names(const std::vector<std::string>& names);
std::vector<std::vector<double>> build_distance_matrix(const std::vector<std::vector<double>>& tm_mat, int chain_num);
std::unique_ptr<UPGMA_Tree> build_upgma_tree(const std::vector<std::string>& names, const std::vector<std::vector<double>>& dist_mat);
std::string get_tree_string(const UPGMA_Tree* tree, int precision = 5);
void print_tree(const UPGMA_Tree* tree, std::ostream& os = std::cout);
bool save_svg(const UPGMA_Tree* tree, const std::string& filename);
void save_distance_matrix(const std::vector<std::string>& names, const std::vector<std::vector<double>>& dist_mat, int chain_num);
void output_upgma_tree(const std::vector<std::string>& names, const std::vector<std::vector<double>>& tm_mat, int chain_num);
