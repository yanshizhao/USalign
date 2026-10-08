#pragma once

#include "common.h"

struct Pdb2xyzOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    int ter_opt = 3;                // TER, END, or different chainID
    int infmt_opt = -1;             // PDB or PDBx/mmCIF format
    int split_opt = 0;              // do not split chain
    int het_opt = 0;                // do not read HETATM residues
    bool normalize_atom_name = false;
};

void print_help();
bool parse_numeric_option(int argc, char* argv[], int& i, Pdb2xyzOptions& pdb2xyz_opts);
bool parse_string_option(int argc, char* argv[], int& i, Pdb2xyzOptions& pdb2xyz_opts);
void parse_arguments(int argc, char* argv[], Pdb2xyzOptions& pdb2xyz_opts);
void check_split_conflicts(Pdb2xyzOptions& pdb2xyz_opts);
void postprocess_arguments(Pdb2xyzOptions& pdb2xyz_opts);
void output_xyz_chain(const std::string& filename, const Pdb2xyzOptions& pdb2xyz_opts, const std::vector<std::string>& chain_lines, const std::string& chain_id);
void convert_file_to_xyz(const std::string& filename, const Pdb2xyzOptions& pdb2xyz_opts);
