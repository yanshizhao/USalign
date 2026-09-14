#pragma once

// 常设公共层共享类型：跨模式共享的数据结构。
// 仅依赖 basic_fun.h 的基础类型（CoordArray/Vec3/RotMat/DoubleMatrix 等），
// 不引用任何引擎头（TMalign.h/MMalign.h/flexalign.h/…），避免 include 环。
#include "basic_fun.h"

#include <string>
#include <vector>

using namespace std;

// ---- 单链解析结果（mm0/mm5-6/mm7 用）----
struct ParsedChain
{
    // --- 24-byte heavy objects (vector / string) ---
    CoordArray     chain_coords;            // 3D coordinates
    string         chain_seq;              // sequence
    string         chain_sec;              // secondary structure
    vector<string> resi_vec;                 // residue index (for -do output)
    string         chain_id;                 // chain ID
    string         filename;                 // source filename (for output)
    vector<string> pdb_lines;                // raw PDB lines (for -do output)
    // --- 4-byte scalars packed together ---
    int            chain_len;                // length
    int            cur_complex_mol_list;     // molecule type (-1=protein, 1=RNA)
};

// ---- 多链解析结果（mm1/mm2/mm4 用）----
struct ComplexData
{
    DoubleCube coords;           // coordinates of each chain (was xa_vec)
    CharMatrix seqs;             // sequence of each chain (was seqx_vec)
    CharMatrix secs;             // secondary structure of each chain (was secx_vec)
    vector<int> mol_types;       // molecule type of each chain (was mol_vec)
    vector<int> lengths;         // length of each chain (was xlen_vec)
    vector<string> chain_ids;    // chain ID (was chainID_list)
    vector<string> resi;         // residue indices (was resi_vec)
    int total_len_aa;
    int total_len_na;
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
    string seqM;
    string seqxA;
    string seqyA;
    double rmsd0;
    double Liden;
    double TM_ali;
    double rmsd_ali;
    int L_ali;
    int n_ali;
    int n_ali8;
    vector<double> do_vec;
    vector<int> invmap;
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