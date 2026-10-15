#pragma once

#include "common.h"
#include "flexalign.h"

struct USFlexalignOptions
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
    std::string dir_opt   ="";    // set -dir to empty
    std::string dir1_opt  ="";    // set -dir1 to empty
    std::string dir2_opt  ="";    // set -dir2 to empty
    std::string dirpair_opt="";   // set -dirpair to empty
    std::vector<std::string> sequence; // get value from alignment file
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    double Lnorm_ass;
    double d0_scale;
    double TMcut     =-1;
    int    i_opt = 0;     // 1 for -i, 3 for -I
    int    o_opt = 0;     // 1 for -o, 2 for -rasmol
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
    int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
    int    ter_opt   =-1;    // default change to 2 (END, or different chainID)
    int    split_opt =-1;    // default change to 2 (split each chain)
    int    outfmt_opt=0;     // set -outfmt to full output
    int    mirror_opt=0;     // do not align mirror
    int    het_opt   =0;     // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    int    hinge_opt =9;     // maximum number of hinge allowed
    int    read_resi =0;     // whether to read residue index
    double TMpass_opt=0.85;  // early stopping threshold for AFP
    bool hinge_set = false;  // whether -hinge is set by user
    bool h_opt = false; // print full help message
    bool v_opt = false; // print version
    bool m_opt = false; // flag for -m, output rotation matrix
    bool u_opt = false; // flag for -u, normalized by user specified length
    bool d_opt = false; // flag for -d, user specified d0
    bool usbcat_opt=false;   // flag for -afp, AFP-enhancement mechanism
    bool fast_opt  =false; // flags for -fast, fTM-align algorithm
};

void print_extra_help();
void print_help(bool h_opt=false);
inline void set_output_format(int& o_opt, std::string& fname_super, const int incoming, const std::string& value);
bool parse_file_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts);
bool parse_numeric_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts);
bool parse_string_option(int argc, char* argv[], int& i, USFlexalignOptions& usflexalign_opts);
void parse_arguments(int argc, char* argv[], USFlexalignOptions& usflexalign_opts);
void check_dir_conflicts(USFlexalignOptions& usflexalign_opts);
void check_hinge(USFlexalignOptions& usflexalign_opts);
void apply_derived_options(USFlexalignOptions& usflexalign_opts);
void check_flexalign_input_presence(USFlexalignOptions& usflexalign_opts);
void postprocess_arguments(USFlexalignOptions& usflexalign_opts);
void fill_flexalign_params(FlexalignParams& flex_params, const USFlexalignOptions& usflexalign_opts);
void fill_common_inputs(AlignCommonInput& common_inputs, const USFlexalignOptions& usflexalign_opts, const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list);
void align_all_chain_pairs(AlignCommonInput& common_inputs, const FlexalignParams& flex_params);
