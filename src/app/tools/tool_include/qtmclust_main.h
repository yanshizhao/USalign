#pragma once

#include "common.h"
#include "HwRMSD.h"
#include "TMalign.h"

struct QTMclustOptions
{
    std::string outdir_opt = "";
    std::string xname = "";
    std::string fname_clust = ""; // file name for output cluster result
    std::string fname_lign  = ""; // file name for user alignment
    std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
    std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
    std::string suffix_opt="";    // set -suffix to empty
    std::string dir_opt   ="";    // set -dir to empty
    std::vector<std::string> sequence; // get value from alignment file
    std::vector<std::string> chain2parse;
    std::vector<std::string> model2parse;
    double TMcut       = 0.5;
    double Lnorm_ass = 0;
    double d0_scale = 0;
    int    i_opt = 0;     // 3 for -I, stick to user given alignment
    int    a_opt = 0;     // flag for -a, do not normalized by average length
    int    s_opt = 2;     // flag for -s, normalized by longer length
    int    infmt_opt =-1;    // PDB or PDBx/mmCIF format
    int    ter_opt   =3;     // TER, END, or different chainID
    int    split_opt =0;     // do not split chain
    int    het_opt   =0;     // do not read HETATM residues
    int    byresi_opt=0;     // set -byresi to 0
    bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
    bool   u_opt = false; // flag for -u, normalized by user specified length
    bool   d_opt = false; // flag for -d, user specified d0
    bool   h_opt = false; // print full help message
};

struct QTMclustRun
{
    const QTMclustOptions& qtmclust_opts;
    std::vector<int>& mol_vec;
    std::vector<std::string>& chainID_list;
    CharMatrix& seq_vec;
    std::vector<std::string>& sec_vec;
    std::vector<std::vector<std::vector<float> > >& xyz_vec;
    const std::vector<size_t>& clust_repr_vec;
    std::vector<size_t>& clust_mem_vec;
    std::map<size_t,size_t>& clust_repr_map;
    std::vector<size_t>& index_vec;
    std::map<std::string, std::map<std::string,bool> >& init_cluster;
    CoordArray& xa;
    size_t chain_i;
    int xlen;
    double TMcut;
    double fast_lb;
    double fast_ub;
    double lb_HwRMSD;
    double lb_TMfast;
    double ub_HwRMSD;
    double ub_TMfast;
    int glocal;
    int iter_opt;
    int min_repr_num;
    int max_repr_num;
};

void print_extra_help();
void print_help(bool h_opt=false);
void filter_lower_bound(double &lb_HwRMSD, double &lb_TMfast, const double TMcut, const int s_opt,const int mol_type);
void read_init_cluster(const std::string&filename, std::map<std::string, std::map<std::string,bool> > &init_cluster);
bool parse_file_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts, std::map<std::string, std::map<std::string,bool> >& init_cluster);
bool parse_numeric_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts);
bool parse_string_option(int argc, char* argv[], int& i, QTMclustOptions& qtmclust_opts);
void parse_arguments(int argc, char* argv[], QTMclustOptions& qtmclust_opts, std::map<std::string, std::map<std::string,bool> >& init_cluster);
void check_input_presence(QTMclustOptions& qtmclust_opts);
void check_dir_conflicts(QTMclustOptions& qtmclust_opts);
void check_numeric_and_split(QTMclustOptions& qtmclust_opts);
void apply_derived_options(QTMclustOptions& qtmclust_opts);
void postprocess_arguments(QTMclustOptions& qtmclust_opts);
void build_file_lists(const QTMclustOptions& qtmclust_opts, std::vector<std::string>& chain_list);
void read_all_chains(const std::vector<std::string>& chain_list, const QTMclustOptions& qtmclust_opts, std::vector<int>& mol_vec, std::vector<std::string>& chainID_list, std::vector<std::pair<int,size_t> >& chainLen_list, CharMatrix& seq_vec, std::vector<std::string>& sec_vec, std::vector<std::vector<std::vector<float> > >& xyz_vec);
void fill_chain_candidates(QTMclustRun& run);
void run_hwrmsd_prefilter_stage(QTMclustRun& run);
void fill_qtmclust_align_options(ChainPairAlignOptions& align_opts, const QTMclustOptions& qtmclust_opts, const bool fast_opt, const int mol_type);
bool run_tmalign_refine_stage(QTMclustRun& run);
void report_clusters(const std::vector<size_t>& clust_repr_vec, const std::vector<size_t>& clust_mem_vec, const std::vector<std::string>& chainID_list, const std::string& fname_clust);
