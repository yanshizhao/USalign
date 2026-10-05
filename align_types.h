#pragma once

#include "basic_fun.h"

#include <string>
#include <vector>

// ---- single-chain parse result (used by mm0/mm5-6/mm7) ----
struct ParsedChain
{
    // --- 24-byte heavy objects (vector / string) ---
    CoordArray     chain_coords;            // 3D coordinates
    std::string    chain_seq;              // sequence
    std::string    chain_sec;              // secondary structure
    std::vector<std::string> resi_vec;      // residue index (for -do output)
    std::string    chain_id;               // chain ID
    std::string    filename;               // source filename (for output)
    std::vector<std::string> pdb_lines;    // raw PDB lines (for -do output)
    // --- 4-byte scalars packed together ---
    int            chain_len;              // length
    int            cur_complex_mol_list;   // molecule type (-1=protein, 1=RNA)
};

// ---- multi-chain parse result (used by mm1/mm2/mm4) ----
struct ComplexData
{
    DoubleCube coords;           // coordinates of each chain (was xa_vec)
    CharMatrix seqs;             // sequence of each chain (was seqx_vec)
    CharMatrix secs;             // secondary structure of each chain (was secx_vec)
    std::vector<int> mol_types;  // molecule type of each chain (was mol_vec)
    std::vector<int> lengths;    // length of each chain (was xlen_vec)
    std::vector<std::string> chain_ids;  // chain ID (was chainID_list)
    std::vector<std::string> resi;       // residue indices (was resi_vec)
    int total_len_aa;
    int total_len_na;
};

struct TrimmedComplex
{
    DoubleCube coords;
    CharMatrix seqs;
    CharMatrix secs;
    std::vector<int> lengths;
    int chain_count;
    int max_aa_len;
    int max_na_len;
};

struct SoiAlignContext
{
    CoordArray k_nearest1;
    CoordArray k_nearest2;
    IntPairArray chain1_bonds;
    IntPairArray chain2_bonds;
    std::vector<double> pair_distances;
};

struct ChainPairDesc
{
    std::string name1;
    std::string name2;
    std::string chain_id1;
    std::string chain_id2;
    int len1;
    int len2;
};

// ---- per chain-pair structure alignment result (shared by every pair in mm0-4) ----
struct ChainPairAlignResult
{
    Vec3 t0;
    RotMat u0;
    double TM1;
    double TM2;
    double TM3;
    double TM4;
    double TM5;
    double d0_0;
    double TM_0;
    double d0A;
    double d0B;
    double d0u;
    double d0a;
    double d0_out;
    std::string seqM;
    std::string seqxA;
    std::string seqyA;
    double rmsd0;
    double Liden;
    double TM_ali;
    double rmsd_ali;
    int L_ali;
    int n_ali;
    int n_ali8;
    std::vector<double> do_vec;
    std::vector<int> invmap;
};

inline void init_rt_mat(Vec3& t0, RotMat& u0)
{
    u0[0][0]=u0[1][1]=u0[2][2]=1;
    u0[0][1]=         u0[0][2]=
    u0[1][0]=         u0[1][2]=
    u0[2][0]=         u0[2][1]=
    t0[0]   =t0[1]   =t0[2]   =0;
}

struct MstaIterationState
{
    std::vector<std::string> msa;
    std::vector<std::pair<double,int> > TM_pair_vec;
    int repr_idx;
    int compare_num;
    int max_iter;
    double TM4_total_max;
    std::vector<double> TM_vec;
    std::vector<double> d0_vec;
    std::vector<double> seqID_vec;
    DoubleMatrix TM_mat;
    DoubleMatrix d0_mat;
    DoubleMatrix seqID_mat;
    std::vector<int> assign_list;
    std::vector<std::string> xname_vec;
    std::vector<std::string> yname_vec;
    ChainPairAlignResult totals;
    int xlen_total;
    int ylen_total;
};

// ---- command-line options container (shared by all mm0-7 modes) ----
struct UserOptions
{
    // file paths: structure files and output files
    std::string xname;
    std::string yname;
    std::string fname_super;
    std::string fname_lign;
    std::string fname_matrix;

    // directory options: directory arguments given on the command line
    std::string dir_opt;
    std::string dirpair_opt;
    std::string dir1_opt;
    std::string dir2_opt;

    // command-line options: format, filtering and other configuration
    int infmt1_opt;
    int infmt2_opt;
    int ter_opt;
    int split_opt;
    int het_opt;
    std::string atom_opt;
    std::string mol_opt;
    int mirror_opt;
    std::vector<std::string> chain2parse1;
    std::vector<std::string> chain2parse2;
    std::vector<std::string> model2parse1;
    std::vector<std::string> model2parse2;
    int byresi_opt;
    bool fast_opt;
    int i_opt;
    int o_opt;
    int a_opt;
    bool m_opt;
    bool u_opt;
    bool d_opt;
    int outfmt_opt;

    // count / ratio thresholds
    double Lnorm_ass;
    double d0_scale;
    double TMcut;

    UserOptions() : infmt1_opt(-1), infmt2_opt(-1), ter_opt(-1), split_opt(-1),
                    het_opt(0), atom_opt("auto"), mol_opt("auto"), mirror_opt(0),
                    byresi_opt(0), fast_opt(false), i_opt(0), o_opt(0), a_opt(0),
                    m_opt(false), u_opt(false), d_opt(false), outfmt_opt(0),
                    Lnorm_ass(0.0), d0_scale(0.0), TMcut(-1) {}
};

struct ParsedInput
{
    // parse result: chain lists, sequence and derived config obtained from a PDB file
    std::vector<std::string> chain1_list;
    std::vector<std::string> chain2_list;
    std::vector<std::string> sequence;
    bool normalize_atom_name;

    ParsedInput() : normalize_atom_name(false) {}
};

struct ControlOptions
{
    int         mm_opt;
    int         cp_opt;
    bool        full_opt;
    int         closeK_opt;
    std::string chainmapfile;
    bool        se_opt;
    bool        do_opt;
    int         parallel_threads;
    bool        h_opt;
    bool        v_opt;
    std::string suffix_opt;
    bool        usbcat_opt;
    int         hinge_opt;
    bool        hinge_set;
    double      TMpass_opt;

    ControlOptions() : mm_opt(0), cp_opt(0), full_opt(false), closeK_opt(-1),
                       se_opt(false), do_opt(false), parallel_threads(0),
                       h_opt(false), v_opt(false), chainmapfile(""),
                       suffix_opt(""), usbcat_opt(false), hinge_opt(9),
                       hinge_set(false), TMpass_opt(0.85) {}
};

struct AlignCommonInput
{
    UserOptions user_options;
    ParsedInput parsed_input;
    ControlOptions control_options;
};


struct MMalignParams
{
    std::string dir1_opt;              // normalized dir1
    std::string dir2_opt;              // normalized dir2
    std::vector<std::string> chain1_list;  // chain file list of structure 1 for this call
    std::vector<std::string> chain2_list;  // chain file list of structure 2 for this call
};


struct AllChainPairsResult
{
    DoubleMatrix tm_matrix;                 // TM-score matrix of chain pairs
    RotArray rotations;                     // rotation of each chain pair
    std::vector<std::vector<std::string> > aligned_seq1;   // chain-pair residue alignment 1 
    std::vector<std::vector<std::string> > aligned_seq2;   // chain-pair residue alignment 2 
    std::vector<std::vector<std::string> > aligned_consensus;  // consensus sequence
    double best_pair_tm;                 // best monomer chain-pair TM 
    int best_pair_chain1_idx;                     // was maxTMmono_i
    int best_pair_chain2_idx;                     // was maxTMmono_j
};

struct ChainAssignResult
{
    std::vector<int> chain2_of_chain1;   // index in structure 2 for each chain of structure 1
    std::vector<int> chain1_of_chain2;   // index in structure 1 for each chain of structure 2
};


struct MMalignFinalParams
{
    std::string xname;
    std::string yname;
    std::string fname_super;
    std::string fname_lign;
    std::string fname_matrix;
    std::string seqx_arg;
    std::string seqy_arg;
    int    len_aa;
    int    len_na;
    int    chain1_num;
    int    chain2_num;
    double d0_scale;
    bool   m_opt;
    int    o_opt;
    int    outfmt_opt;
    int    ter_opt;
    int    split_opt;
    bool   a_opt;
    bool   d_opt;
    bool   fast_opt;
    bool   full_opt;
    int    mirror_opt;
};


struct ChainPairAlignOptions
{
    int    i_opt;              // alignment mode: 0=auto 1=strategy+user-aln 2=user-aln only 3=-I
    int    a_opt;              // whether to output TM3
    int    u_opt;              // 0/1/2: 1=TM4 by user Lnorm; 2=also override search d0 with Lnorm
    bool   d_opt;              // whether to output TM5
    bool   fast_opt;           // fast mode
    bool   se_opt;             // use se_main instead of TMalign_main
    bool   cp_opt;             // use CPalign_main
    double Lnorm;              // normalization length
    double d0_scale;           // d0 scale
    double TMcut;              // early-termination threshold
    int    parallel_threads;   // intra-engine parallelism
    int    ss_opt;             // 1=skip the secondary-structure-based initial strategy
    int    mol_type;           // sum of the molecule types of the two chains
};

struct OutputOptions
{
    std::string fname_matrix;
    std::string fname_super;
    double Lnorm_ass;
    double d0_scale;
    int    outfmt_opt;
    int    ter_opt;
    int    mm_opt;
    int    split_opt;
    int    o_opt;
    int    i_opt;
    int    a_opt;
    int    mirror_opt;
    bool   u_opt;
    bool   d_opt;
};

inline void fill_chain_pair_align_options(ChainPairAlignOptions& align_opts,
    const UserOptions& user_opts, const ControlOptions& ctrl_opts)
{
    align_opts.i_opt = user_opts.i_opt;
    align_opts.a_opt = user_opts.a_opt;
    align_opts.u_opt = user_opts.u_opt;
    align_opts.d_opt = user_opts.d_opt;
    align_opts.fast_opt = user_opts.fast_opt;
    align_opts.se_opt = ctrl_opts.se_opt;
    align_opts.cp_opt = ctrl_opts.cp_opt;
    align_opts.Lnorm = user_opts.Lnorm_ass;
    align_opts.d0_scale = user_opts.d0_scale;
    align_opts.TMcut = user_opts.TMcut;
    align_opts.parallel_threads = ctrl_opts.parallel_threads;
    align_opts.ss_opt = 0;
}

inline void fill_out_opts_from_user_opts(OutputOptions& out_opts, const UserOptions& user_opts)
{
    out_opts.fname_matrix = user_opts.m_opt ? user_opts.fname_matrix : "";
    out_opts.fname_super = user_opts.fname_super;
    out_opts.Lnorm_ass = user_opts.Lnorm_ass;
    out_opts.d0_scale = user_opts.d0_scale;
    out_opts.outfmt_opt = user_opts.outfmt_opt;
    out_opts.ter_opt = user_opts.ter_opt;
    out_opts.mm_opt = 0;
    out_opts.split_opt = user_opts.split_opt;
    out_opts.o_opt = user_opts.o_opt;
    out_opts.i_opt = user_opts.i_opt;
    out_opts.a_opt = user_opts.a_opt;
    out_opts.mirror_opt = user_opts.mirror_opt;
    out_opts.u_opt = user_opts.u_opt;
    out_opts.d_opt = user_opts.d_opt;
}

inline void fill_out_opts_mmalign_final(OutputOptions& out_opts, const MMalignFinalParams& params)
{
    out_opts.fname_matrix = params.m_opt ? params.fname_matrix : "";
    out_opts.fname_super = params.fname_super;
    out_opts.Lnorm_ass = 0;
    out_opts.d0_scale = params.d0_scale;
    out_opts.outfmt_opt = params.outfmt_opt;
    out_opts.ter_opt = params.ter_opt;
    out_opts.mm_opt = 0;
    out_opts.split_opt = params.split_opt;
    out_opts.o_opt = params.o_opt;
    out_opts.i_opt = 0;
    out_opts.a_opt = params.a_opt;
    out_opts.mirror_opt = params.mirror_opt;
    out_opts.u_opt = false;
    out_opts.d_opt = params.d_opt;
}
