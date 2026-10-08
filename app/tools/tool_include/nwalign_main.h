#pragma once

#include "common.h"

struct NWalignOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string yname = "";
    std::string atom_opt = "auto";  // use C alpha atom for protein and C3' for RNA
    std::string mol_opt = "auto";   // auto-detect the molecule type as protein/RNA
    std::string suffix_opt = "";    // set -suffix to empty
    std::string dir_opt = "";       // set -dir to empty
    std::string dir1_opt = "";      // set -dir1 to empty
    std::string dir2_opt = "";      // set -dir2 to empty
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    int infmt1_opt = -1;            // FASTA sequence
    int infmt2_opt = -1;            // FASTA sequence
    int ter_opt = 3;                // TER, END, or different chainID
    int split_opt = 0;              // do not split chain
    int outfmt_opt = 0;             // set -outfmt to full output
    int glocal = 0;
    int het_opt = 0;                // do not read HETATM residues
    bool h_opt = false;             // print full help message
    bool normalize_atom_name = false;
};

void print_extra_help();
void print_help(bool h_opt=false);
bool parse_numeric_option(int argc, char* argv[], int& i, NWalignOptions& nwalign_opts);
bool parse_string_option(int argc, char* argv[], int& i, NWalignOptions& nwalign_opts);
void parse_arguments(int argc, char* argv[], NWalignOptions& nwalign_opts);
void check_input_presence(NWalignOptions& nwalign_opts);
void postprocess_arguments(NWalignOptions& nwalign_opts);
void extract_chain_sequence(const std::vector<std::string>& chain_lines, const int infmt_opt, const int xlen, std::string& seq);
void align_all_chain_pairs(std::vector<std::string>& chain1_list, std::vector<std::string>& chain2_list, const NWalignOptions& nwalign_opts);
