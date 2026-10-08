#pragma once

#include "common.h"
#include "TMalign.h"

struct Pdb2ssOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 3;                // TER, END, or different chainID
    int infmt_opt = -1;             // PDB format
    int split_opt = 0;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

void print_help();
bool parse_numeric_option(int argc, char* argv[], int& i, Pdb2ssOptions& pdb2ss_opts);
bool parse_string_option(int argc, char* argv[], int& i, Pdb2ssOptions& pdb2ss_opts);
void parse_arguments(int argc, char* argv[], Pdb2ssOptions& pdb2ss_opts);
void postprocess_arguments(Pdb2ssOptions& pdb2ss_opts);
void output_ss_chain(const std::string& filename, const Pdb2ssOptions& pdb2ss_opts, const std::vector<std::string>& chain_lines, const int mol_type, const std::string& chain_id);
void convert_file_to_ss(const std::string& filename, const Pdb2ssOptions& pdb2ss_opts);
