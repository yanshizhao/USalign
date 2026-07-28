#include "UPGMA.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <map>
#include <functional>

static vector<string> split_name(const string& name)
{
    vector<string> parts;
    size_t j = 0;
    for (size_t i = 0; i <= name.size(); i++) {
        if (i == name.size() || name[i] == '+') {
            parts.push_back(name.substr(j, i - j));
            j = i + 1;
        }
    }
    return parts;
}

static pair<int,int> find_min_distance(
    const vector<vector<double>>& mat, int len)
{
    double min_dis = 1e10;
    int col = -1, row = -1;
    for (int i = 0; i < len; i++)
        for (int j = 0; j < i; j++)
            if (mat[i][j] < min_dis) {
                min_dis = mat[i][j];
                col = i; row = j;
            }
    return {col, row};
}

// @input  nodes     - all nodes in tree
// @input  edges     - all edges in tree
// @input  name      - node name to update
// @input  new_level - merge round number to assign
// @output updates level field of matching nodes and edges in-place
static void change_node_level(
    vector<UPGMA_Node>& nodes, vector<UPGMA_Edge>& edges,
    const string& name, int new_level)
{
    for (auto& n : nodes)
    {
        if (n.name == name) n.level = new_level;
    }

    for (auto& e : edges) 
    {
        if (e.child.name == name) 
        {
            e.child.level = new_level;
        } 
        else if (e.parent.name == name) 
        {
            e.parent.level = new_level;
        }
    }
}

// UPGMA average linkage: compute avg distance between two groups
// @input  group_a   - group A members (e.g. "1e94e+1g41a")
// @input  group_b   - group B members (e.g. "1e32a")
// @input  orig_names - original structure names
// @input  orig_dist  - original distance matrix
// @output average of all cross-group pair distances
static double get_avg_distance_between(
    const string& group_a, const string& group_b,
    const vector<string>& orig_names,
    const vector<vector<double>>& orig_dist)
{
    vector<string> members_a = split_name(group_a);
    vector<string> members_b = split_name(group_b);
    double sum = 0;
    int cnt = 0;
    for (size_t i = 0; i < orig_names.size(); i++)
        for (size_t j = 0; j < orig_names.size(); j++) {
            bool in_a = find(members_a.begin(), members_a.end(), orig_names[i]) != members_a.end();
            bool in_b = find(members_b.begin(), members_b.end(), orig_names[j]) != members_b.end();
            if (in_a && in_b) { sum += orig_dist[i][j]; cnt++; }
        }
    return (cnt > 0) ? sum / cnt : 0;
}

// shrink distance matrix by one row/col after merging two clusters
// @input  mat       - current distance matrix (will be shrunk)
// @input  names     - current cluster names (will be shrunk)
// @input  merged_a  - index of merged cluster A
// @input  merged_b  - index of merged cluster B
// @input  orig_names - original structure names (for backtracking)
// @input  orig_dist  - original distance matrix (for backtracking)
// @output mat   updated to (n-1)x(n-1)
// @output names updated to n-1 names (merged cluster added)
static void renew_matrix(
    vector<vector<double>>& mat, vector<string>& names,
    int merged_a, int merged_b,
    const vector<string>& orig_names,
    const vector<vector<double>>& orig_dist)
{
    int max_idx = max(merged_a, merged_b);
    int min_idx = min(merged_a, merged_b);
    int old_len = (int)names.size();
    int new_len = old_len - 1;

    // build new name list: keep unmerged clusters, merge two into "A+B" at end
    vector<string> new_names;
    for (int i = 0; i < old_len; i++)
        if ((i != min_idx) && (i != max_idx))
            new_names.push_back(names[i]);
    new_names.push_back(names[max_idx] + "+" + names[min_idx]);
    vector<vector<double>> new_mat(new_len, vector<double>(new_len, 0.0));

    // build old name->index map for O(1) lookup
    map<string, int> name_to_idx;
    for (int t = 0; t < old_len; t++)
        name_to_idx[names[t]] = t;

    for (int i = 0; i < new_len; i++) {
        for (int k = 0; k < i; k++) {
            map<string, int>::iterator iter_i = name_to_idx.find(new_names[i]);
            map<string, int>::iterator iter_k = name_to_idx.find(new_names[k]);
            if (iter_i != name_to_idx.end() && iter_k != name_to_idx.end()) {
                int idx_i = iter_i->second;
                int idx_k = iter_k->second;
                new_mat[i][k] = mat[max(idx_i, idx_k)][min(idx_i, idx_k)];
            } else {
                new_mat[i][k] = get_avg_distance_between(new_names[i], new_names[k], orig_names, orig_dist);
            }
        }
    }
    mat = move(new_mat);
    names = move(new_names);
}


// ============ UPGMA clustering ============
// @input  names    - structure names (cleaned, no chain ID)
// @input  dist_mat - distance matrix
// @output UPGMA_Tree* - UPGMA tree (nodes + edges)
// 1) Initialize each structure as its own cluster
// 2) Merge the closest pair each iteration
// 3) Update distance matrix (average linkage)
// 4) Repeat until one cluster remains

unique_ptr<UPGMA_Tree> build_upgma_tree(
    const vector<string>& names,
    const vector<vector<double>>& dist_mat)
{
    auto tree = unique_ptr<UPGMA_Tree>(new UPGMA_Tree);
    int structure_num = (int)names.size();
    tree->cluster_names = names;
    tree->distance_matrix = dist_mat;
    tree->nodes.reserve(structure_num * 2 - 1);

    for (int i = 0; i < structure_num; i++)
        tree->nodes.push_back({names[i], -1});

    // initialize leaf heights to 0 for correct branch length calculation
    map<string, double> node_height;
    for (int i = 0; i < structure_num; i++)
        node_height[names[i]] = 0.0;

    for (int iter = 0; iter < structure_num - 1; iter++) {

        pair<int,int> min_pair = find_min_distance(tree->distance_matrix, (int)tree->cluster_names.size());
        int closest_idx_a = min_pair.first;
        int closest_idx_b = min_pair.second;

        string cluster_a_name = tree->cluster_names[closest_idx_a];
        string cluster_b_name = tree->cluster_names[closest_idx_b];
        double merge_dist = tree->distance_matrix[closest_idx_a][closest_idx_b];
        double parent_height = merge_dist / 2.0;
        string parent_name = cluster_a_name + "+" + cluster_b_name;

        tree->nodes.push_back({parent_name, -1});
        for (auto& node : tree->nodes) {
            if (node.name == cluster_a_name || node.name == cluster_b_name)
            {
                change_node_level(tree->nodes, tree->edges, node.name, iter + 1);
                // branch length from parent to child = parent height - child accumulated height
                // child height was accumulated in previous merges; only add remaining increment to new parent
                double branch_len = parent_height - node_height[node.name];
                tree->edges.push_back({node, {parent_name, -1}, branch_len});
            }
        }
        node_height[parent_name] = parent_height;

        renew_matrix(tree->distance_matrix, tree->cluster_names, closest_idx_a, closest_idx_b, names, dist_mat);
    }
    return tree;
}



// ============ Newick output ============
// recursive Newick builder (helper)
static void build_newick_node(const string& node,
    const map<string, vector<pair<string, double>>>& children,
    stringstream& ss)
{
    auto it = children.find(node);
    if (it == children.end()) 
    { 
        ss << node; 
        return; 
    }
    ss << "(";
    for (size_t i = 0; i < it->second.size(); i++) 
    {
        if (i > 0) ss << ",";
        build_newick_node(it->second[i].first, children, ss);
        ss << ":" << it->second[i].second;
    }
    ss << ")";
}

// Serialize tree to Newick format string
string get_tree_string(const UPGMA_Tree* tree, int precision)
{
    if (tree->nodes.empty() || tree->edges.empty()) return "";
    map<string, vector<pair<string, double>>> children;
    for (auto& e : tree->edges) {
        string child  = e.child.name;
        string parent = e.parent.name;
        children[parent].push_back({child, e.branch_length});
    }
    stringstream ss;
    ss << fixed << setprecision(precision);
    build_newick_node(tree->nodes.back().name, children, ss);
    ss << ";";
    return ss.str();
}


// build parent->children name map from edges (no branch lengths)
static void build_child_tree_map(const UPGMA_Tree* tree,
    map<string, vector<string>>& children)
{
    for (auto& e : tree->edges)
        children[e.parent.name].push_back(e.child.name);
}


// ============ ASCII tree output ============
// Print a vertical ASCII tree to terminal
// Uses only ASCII characters for compatibility

// Recursively print a node and its subtree as ASCII tree
// @input  node     - current node name to print
// @input  children - map: parent -> list of children
// @input  os       - output stream
// @input  prefix   - indentation string for tree connectors
// For each child:
//   - if leaf:   print "prefix + connector + node_name"
//   - if branch: print "prefix + connector + *", then recurse with deeper prefix
static void print_ascii_node(
    const string& node,
    const map<string, vector<string>>& children,
    ostream& os,
    const string& prefix)
{
    auto it = children.find(node);
    bool is_leaf = (it == children.end() || it->second.empty());

    if (is_leaf) {
        os << prefix << node << "\n";
        return;
    }

    for (size_t i = 0; i < it->second.size(); i++) {
        bool last = (i == it->second.size() - 1);
        string conn = last ? "`-- " : "|-- ";
        string next = prefix + (last ? "    " : "|   ");
        auto child_it = children.find(it->second[i]);
        bool child_is_leaf = (child_it == children.end() || child_it->second.empty());

        if (child_is_leaf) {
            os << prefix + conn + it->second[i] << "\n";
        } else {
            os << prefix + conn + "*\n";
            print_ascii_node(it->second[i], children, os, next);
        }
    }
}

// Print entire tree as ASCII diagram
// @input tree - UPGMA tree with nodes and edges
// @input os   - output stream
void print_tree(const UPGMA_Tree* tree, ostream& os)
{
    if (tree->nodes.empty() || tree->edges.empty()) return;

    // Build parent-child map and find root
    map<string, vector<string>> children;
    build_child_tree_map(tree, children);
    string root = tree->nodes.back().name;
    os << "# UPGMA Tree:" << "\n";

    // Print root ".", then recursively print its children
    os << ".\n";
    auto root_it = children.find(root);
    if (root_it != children.end()) {
        for (size_t i = 0; i < root_it->second.size(); i++) {
            bool last = (i == root_it->second.size() - 1);
            string conn = last ? "`-- " : "|-- ";
            bool child_is_leaf = (children.find(root_it->second[i]) == children.end() ||
                                  children.at(root_it->second[i]).empty());
            if (child_is_leaf) {
                os << conn + root_it->second[i] << "\n";
            } else {
                os << conn + "*\n";
                print_ascii_node(root_it->second[i], children, os, "    ");
            }
        }
    }
}



// ============ SVG vector output ============
// Generate horizontal SVG tree diagram
// Root on left, leaves on right, monospace font
// Zero external dependencies, viewable in browser

// calculate SVG layout: depth and Y coordinates
static void calc_svg_layout(
    const map<string, vector<string>>& children, const string& root,
    map<string, int>& depth, map<string, double>& ypos,
    int& leaf_count, int& svg_height)
{
    // calculate node depths
    function<void(const string&, int)> set_depth = [&](const string& node, int d) {
        depth[node] = d;
        auto it = children.find(node);
        if (it != children.end())
            for (auto& c : it->second) set_depth(c, d + 1);
    };
    set_depth(root, 0);

    // assign Y coordinates
    int leaf_idx = 0;
    const double top_margin = 30;
    const double leaf_margin = 24;

    function<void(const string&)> assign_y = [&](const string& node) {
        auto it = children.find(node);
        if (it == children.end() || it->second.empty()) {
            ypos[node] = top_margin + leaf_idx * leaf_margin;
            leaf_idx++;
            return;
        }
        for (auto& c : it->second) assign_y(c);
        double sum = 0;
        for (auto& c : it->second) sum += ypos[c];
        ypos[node] = sum / (double)it->second.size();
    };
    assign_y(root);
    leaf_count = leaf_idx;
    svg_height = max(200, (int)(top_margin * 2 + leaf_count * leaf_margin));
}

// recursive SVG tree drawer
static void draw_svg_tree(
    ofstream& ofs,
    const map<string, vector<string>>& children,
    const map<string, int>& depth,
    const map<string, double>& ypos,
    const string& node)
{
    const double left_margin = 50;
    const double x_step = 80;
    double nx = left_margin + depth.at(node) * x_step;
    double ny = ypos.at(node);
    auto it = children.find(node);

    if (it != children.end() && !it->second.empty()) {
        double first_y = ypos.at(it->second.front());
        double last_y  = ypos.at(it->second.back());

        if (first_y != last_y)
            ofs << "<line x1=\"" << nx << "\" y1=\"" << first_y
                << "\" x2=\"" << nx << "\" y2=\"" << last_y
                << "\" stroke=\"black\" stroke-width=\"1\"/>\n";

        for (auto& child : it->second) {
            double cy = ypos.at(child);
            double cx = left_margin + depth.at(child) * x_step;

            ofs << "<line x1=\"" << nx << "\" y1=\"" << cy
                << "\" x2=\"" << cx << "\" y2=\"" << cy
                << "\" stroke=\"black\" stroke-width=\"1\"/>\n";

            draw_svg_tree(ofs, children, depth, ypos, child);
        }
    }

    if (it == children.end() || it->second.empty())
        ofs << "<text x=\"" << (nx + 5) << "\" y=\"" << (ny + 4)
            << "\">" << node << "</text>\n";
}

bool save_svg(const UPGMA_Tree* tree, const string& filename)
{
    if (tree->nodes.empty() || tree->edges.empty()) return false;

    map<string, vector<string>> children;
    // build parent-child map
    build_child_tree_map(tree, children);
    string root_name = tree->nodes.back().name;

    map<string, int> depth;
    map<string, double> ypos;
    int leaf_count = 0, svg_height = 200;
    // calculate depth and Y positions
    calc_svg_layout(children, root_name, depth, ypos, leaf_count, svg_height);

    if (leaf_count == 0) return false;
    int svg_width = 600;

    ofstream ofs(filename);
    if (!ofs) return false;

    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << svg_width << "\" height=\"" << svg_height << "\">\n";
    ofs << "<style>text{font-family:monospace;font-size:12px;}</style>\n";
    // recursively draw lines and leaf labels
    draw_svg_tree(ofs, children, depth, ypos, root_name);
    ofs << "</svg>\n";
    return true;
}


// ============ Distance matrix ============
// @input  tm_mat    - TM-score matrix (chain_num x chain_num)
// @input  chain_num - number of chains
// @output dist_mat  - distance matrix, formula: dist = 1 / (1 + TM)
//         0 < dist <= 1, smaller = more similar

vector<vector<double>> build_distance_matrix(
    const vector<vector<double>>& tm_mat,
    int chain_num)
{
    vector<vector<double>> dist_mat(chain_num, vector<double>(chain_num, 0.0));
    for (int i = 0; i < chain_num; i++)
        for (int j = i + 1; j < chain_num; j++) {
            double d = 1.0 / (1.0 + tm_mat[i][j]);
            dist_mat[i][j] = d;
            dist_mat[j][i] = d;
        }
    return dist_mat;
}


// ============ Clean chain names ============
// @input  names - chain names with "filename:chainID" format (e.g. "1d2na.atm:D")
// @output vector<string> - cleaned filenames only (e.g. "1d2na.atm")
vector<string> clean_chain_names(const vector<string>& names)
{
    vector<string> result = names;
    for (size_t i = 0; i < result.size(); i++) {
        size_t pos = result[i].rfind(':');
        if (pos != string::npos)
            result[i] = result[i].substr(0, pos);
    }
    return result;
}


// ============ Distance matrix output ============
// Save distance matrix in PHYLIP format
void save_distance_matrix(const vector<string>& names,
    const vector<vector<double>>& dist_mat, int chain_num)
{
    ofstream dmf("upgma_tree.dist");
    dmf << chain_num << "\n";
    for (int i = 0; i < chain_num; i++) {
        string name = names[i].substr(0, 10);
        if (name.size() < 10) name.append(10 - name.size(), ' ');
        dmf << name;
        for (int j = 0; j < chain_num; j++) {
            double d = (i == j) ? 0.0 : dist_mat[i][j];
            dmf << "  " << fixed << setprecision(5) << d;
        }
        dmf << "\n";
    }
    dmf.close();
}


// ============ One-click build + output ============
// Build UPGMA tree from TM-score matrix
// Output: upgma_tree.txt, upgma_tree.dist,
//         upgma_tree.svg, ASCII tree
void output_upgma_tree(
    const vector<string>& names,
    const vector<vector<double>>& tm_mat,
    int chain_num)
{
    if (chain_num < 3) return;

    // 1. Build UPGMA tree
    vector<string> clean_names = clean_chain_names(names);
    vector<vector<double>> dist_mat = build_distance_matrix(tm_mat, chain_num);
    auto tree = build_upgma_tree(clean_names, dist_mat);

    // 2. Save Newick format tree
    string newick = get_tree_string(tree.get(), 5);
    ofstream nf("upgma_tree.txt");
    nf << newick << "\n";
    nf.close();

    // 3. Save SVG and distance matrix
    save_svg(tree.get(), "upgma_tree.svg");
    save_distance_matrix(clean_names, dist_mat, chain_num);

    // 4. Print tree to terminal
    cout << "\n# UPGMA Phylogenetic Tree:\n";
    cout << "# Newick: " << newick << "\n";
    print_tree(tree.get(), cout);
    cout << "\n";
}

