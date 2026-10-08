#pragma once

#include "common.h"
#include "MMalign.h"

struct MMalignOptions
{
    std::string outdir_opt = "";
    std::string xname       = "";
    std::string yname       = "";
    std::string fname_super = ""; // file name for superposed structure
    std::string fname_lign  = ""; // file name for user alignment
    std::string fname_matrix= ""; // file name for output matrix
    std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    std::string suffix_opt="";    // set -suffix to empty
    std::string dir1_opt  ="";    // set -dir1 to empty
    std::string dir2_opt  ="";    // set -dir2 to empty
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double d0_scale    =0;
    double TMcut     =-1;
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =1;     // ENDMDL or END
    int    split_opt =2;     // split by chain
    int    outfmt_opt=0;     // set -outfmt to full output
    int    mirror_opt=0;     // do not align mirror
    int    het_opt   =0;     // do not read HETATM residues
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool o_opt = false; // flag for -o, output superposed structure
    bool d_opt = false; // flag for -d, user specified d0
    bool   full_opt  = false;// do not show chain level alignment
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
};

void print_extra_help();
void print_help(bool h_opt=false);
bool parse_file_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts);
bool parse_numeric_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts);
bool parse_string_option(int argc, char* argv[], int& i, MMalignOptions& mmalign_opts);
void parse_arguments(int argc, char* argv[], MMalignOptions& mmalign_opts);
void check_input_presence(MMalignOptions& mmalign_opts);
void check_dir_conflicts(MMalignOptions& mmalign_opts);
void check_numeric_conflicts(MMalignOptions& mmalign_opts);
void check_ter_split_conflicts(MMalignOptions& mmalign_opts);
void postprocess_arguments(MMalignOptions& mmalign_opts);
void build_file_lists(const MMalignOptions& mmalign_opts, std::vector<std::string>& chain1_list, std::vector<std::string>& chain2_list);
void parse_complexes(const MMalignOptions& mmalign_opts, const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list, ComplexData& complex1, ComplexData& complex2);
void align_monomer_pair(const ComplexData& complex1, const ComplexData& complex2, const MMalignOptions& mmalign_opts, std::vector<std::string>& sequence);
void align_oligomer_pairs(const ComplexData& complex1, const ComplexData& complex2, const MMalignOptions& mmalign_opts, std::vector<std::string>& sequence);
