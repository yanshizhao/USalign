#pragma once

#include "common.h"
#include "se.h"
#include "output.h"

struct SeOptions
{
    std::string outdir_opt = "";
    std::string xname       = "";
    std::string yname       = "";
    std::string fname_lign  = ""; // file name for user alignment
    std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    std::string suffix_opt="";    // set -suffix to empty
    std::string dir_opt   ="";    // set -dir to empty
    std::string dir1_opt  ="";    // set -dir1 to empty
    std::string dir2_opt  ="";    // set -dir2 to empty
    std::vector<std::string> sequence; // get value from alignment file
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double Lnorm_ass;
    double d0_scale;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    het_opt=0;        // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    bool h_opt = false; // print full help message
    bool i_opt = false; // flag for -i, stick to user given alignment
    bool a_opt = false; // flag for -a, normalized by average length
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool do_opt= false; // flag for -do, output distance of i-th aligned pair
};

void print_extra_help();
void print_help(bool h_opt=false);
bool parse_file_option(int argc, char* argv[], int& i, SeOptions& se_opts);
bool parse_numeric_option(int argc, char* argv[], int& i, SeOptions& se_opts);
bool parse_string_option(int argc, char* argv[], int& i, SeOptions& se_opts);
void parse_arguments(int argc, char* argv[], SeOptions& se_opts);
void check_byresi_split_conflicts(SeOptions& se_opts);
void apply_derived_options(SeOptions& se_opts);
void postprocess_arguments(SeOptions& se_opts);
void read_chain(const std::vector<std::string>& chain_lines, const std::string& filename, const std::string& chain_id, const SeOptions& se_opts, const int mol_type, ParsedChain& chain);
void align_one_chain_pair(ParsedChain& chain1, ParsedChain& chain2, const std::vector<std::string>& pdb_lines1, const std::vector<std::string>& pdb_lines2, const SeOptions& se_opts);
void align_all_chain_pairs(std::vector<std::string>& chain1_list, std::vector<std::string>& chain2_list, SeOptions& se_opts);
