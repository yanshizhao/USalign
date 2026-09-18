#pragma once

// 常设公共层共享类型：跨模式共享的数据结构。
// 仅依赖 basic_fun.h 的基础类型（CoordArray/Vec3/RotMat/DoubleMatrix 等），
// 不引用任何引擎头（TMalign.h/MMalign.h/flexalign.h/…），避免 include 环。
// 全部使用 std:: 前缀，不在共享头引入 using namespace std。
#include "basic_fun.h"

#include <string>
#include <vector>

// ---- 单链解析结果（mm0/mm5-6/mm7 用）----
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

// ---- 多链解析结果（mm1/mm2/mm4 用）----
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

// ---- 单条链对的结构比对结果（mm0-4 的每对链共用）----
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

// ---- 命令行选项容器（mm0-7 全部模式共用）----
struct UserOptions
{
    // 文件路径：结构文件和输出文件
    std::string xname;
    std::string yname;
    std::string fname_super;
    std::string fname_lign;
    std::string fname_matrix;

    // 目录选项：命令行提供的目录参数
    std::string dir_opt;
    std::string dirpair_opt;
    std::string dir1_opt;
    std::string dir2_opt;

    // 命令行选项：格式、过滤等配置
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

    // 计数/比例阈值
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
    // 解析结果：由 PDB 文件解析得到的链列表、序列及派生配置
    std::vector<std::string> chain1_list;
    std::vector<std::string> chain2_list;
    std::vector<std::string> sequence;
    bool autojustify;

    ParsedInput() : autojustify(false) {}
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

// ---- MMalign 每次调用的 per-call 差异参数（dir 与链列表随调用分支变化）----
struct MMalignParams
{
    std::string dir1_opt;              // 归一化后的 dir1（main 计算）
    std::string dir2_opt;              // 归一化后的 dir2（main 计算）
    std::vector<std::string> chain1_list;  // 本次调用结构1 的链文件列表
    std::vector<std::string> chain2_list;  // 本次调用结构2 的链文件列表
};

// ---- mm1 全链对打分结果（stage B）、链分配结果（stage C） ----
struct AllChainPairsResult
{
    DoubleMatrix tm_matrix;                 // TM-score matrix of chain pairs (was TMave_mat)
    RotArray rotations;                     // rotation of each chain pair (was ut_mat)
    std::vector<std::vector<std::string> > aligned_seq1;   // chain-pair residue alignment 1 (was seqxA_mat)
    std::vector<std::vector<std::string> > aligned_seq2;   // chain-pair residue alignment 2 (was seqyA_mat)
    std::vector<std::vector<std::string> > aligned_consensus;  // consensus sequence (was seqM_mat)
    double best_pair_tm;                 // best monomer chain-pair TM (was maxTMmono)
    int best_pair_chain1_idx;                     // was maxTMmono_i
    int best_pair_chain2_idx;                     // was maxTMmono_j
};

struct ChainAssignResult
{
    std::vector<int> chain2_of_chain1;   // index in structure 2 for each chain of structure 1 (was assign1_list)
    std::vector<int> chain1_of_chain2;   // index in structure 1 for each chain of structure 2 (was assign2_list)
};

// ---- mm1 最终输出（MMalign_final）的逐调用命名/长度/选项 ----
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

// ---- 单条链对结构比对的引擎选项（per-pair，跨模式共享）----
struct ChainPairAlignOptions
{
    int    i_opt;              // 比对模式：0=自动 1=策略+用户比对 2=仅用户比对 3=-I
    int    a_opt;              // 是否产出 TM3（平均长度归一）
    int    u_opt;              // 0/1/2：1=按用户 Lnorm 出 TM4；2=额外用 Lnorm 覆盖搜索 d0
    bool   d_opt;              // 是否产出 TM5（缩放 d0）
    bool   fast_opt;           // 快速模式
    bool   se_opt;             // 走 se_main 而非 TMalign_main
    bool   cp_opt;             // 走 CPalign_main（仅 mm0/sn 置位）
    double Lnorm;              // 归一化长度
    double d0_scale;           // d0 缩放
    double TMcut;              // 早退阈值
    int    parallel_threads;   // 引擎内并行度
    int    ss_opt;             // 1=跳过基于二级结构的初始策略
    int    mol_type;           // 两链分子类型之和（-1/0/1，逐对覆盖）
};