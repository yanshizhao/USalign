#pragma once

#include "common.h"

struct BiounitasymOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 0;                // all models
    int infmt_opt = -1;             // PDB or PDBx/mmCIF format
    int split_opt = 1;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

void print_help();
bool parse_numeric_option(int argc, char* argv[], int& i, BiounitasymOptions& biounitasym_opts);
bool parse_string_option(int argc, char* argv[], int& i, BiounitasymOptions& biounitasym_opts);
void parse_arguments(int argc, char* argv[], BiounitasymOptions& biounitasym_opts);
void postprocess_arguments(BiounitasymOptions& biounitasym_opts);
double count_inter_chain_clashes(const CoordArray& xa, const CoordArray& ya);
void compute_pair_clash_ratios(const std::string& filename, const std::vector<std::vector<std::string> >& PDB_lines, const int chainnum, std::vector<double>& clashratio_vec);
void report_biounit_verdict(const std::string& filename, const std::vector<double>& clashratio_vec);
void report_multimodel_clash(const std::string& filename, const BiounitasymOptions& biounitasym_opts);
