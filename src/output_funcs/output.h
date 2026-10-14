#pragma once
#include "common.h"

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

extern const char USALIGN_VERSION[];

void output_pair_alignment(const AlignCommonInput& common_inputs,
    int xlen, int ylen,
    const ChainPairAlignResult& chain_align_res,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg = NULL,
    std::ostream& os = std::cout);

void output_do_distance_table(std::ostream& os,
    const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::vector<double>& do_vec,
    size_t right_num);

void output_soi_dist_table(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2);

void copy_structure2_to_outdir(const std::string& yname,
    const std::string& fname_super, std::string& struct2_basename);

// script format: 0 - no script; 1 - pymol; 3 - chimerax
void output_pymol(const std::string xname, const std::string yname,
    const std::string fname_super, const Vec3& t, const RotMat& u, const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const std::string chainID1, const std::string chainID2, const int o_opt=1);

void read_mTMalign_chain_pdb(const std::string& name, const std::string& fname_super,
    const int chain_idx, const RotArray& ut_mat, int& compress_type, std::string& filename);

// pymol/chimerax color cycle for multi-chain (mTMalign) output
const std::vector<std::string> mTMalign_color_list = {
    "red", "green", "blue", "yellow", "violet", "cyan", "salmon", "lime",
    "pink", "slate", "magenta", "orange", "marine", "olive", "purple",
    "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey",
};

void output_mTMalign_pymol(const std::vector<std::string>& chain_list,
    const RotArray& ut_mat, const std::string& fname_super, const int o_opt=1);

struct TMalignRasmolState
{
    std::stringstream buf;
    std::stringstream buf_all;
    std::stringstream buf_atm;
    std::stringstream buf_all_atm;
    std::stringstream buf_all_atm_lig;
    std::stringstream buf_tm;
    std::string line;
    Vec3 x;  // before transform
    Vec3 x1; // after transform
    bool after_ter; // true if passed the "TER" line in PDB
    std::string asym_id; // chain ID
    std::map<std::string,int> _atom_site;
    int atom_site_pos;
    std::vector<std::string> line_vec;
    std::string atom; // 4-character atom name
    std::string AA;   // 3-character residue name
    std::string resi; // 4-character residue sequence number
    std::string inscode; // 1-character insertion code
    std::string model_index; // model index
    bool is_mmcif;
    int ca_idx1; // all CA atoms
    int lig_idx1; // all atoms
    std::vector <int> idx_vec;
    int ca_idx2; // all CA atoms
    int lig_idx2; // all atoms
    std::vector<std::string> resi_aln1;
    std::vector<std::string> resi_aln2;
    int i1;
    int i2;
    std::string chain1_sele;
    std::string chain2_sele;
};

void normalize_mmcif_atom_fields(
    std::map<std::string,int>& atom_site,
    const std::vector<std::string>& line_vec,
    std::string& atom, std::string& AA, std::string& resi, std::string& inscode);

inline void write_rasmol_mmcif_chain_a(std::ostream& os, TMalignRasmolState& rasmol_state, const int serial)
{
    os<<"ATOM  "<<std::setw(6)<<std::setw(5)<<serial%100000<<' '<<rasmol_state.atom<<' '
        <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
        <<setiosflags(std::ios::fixed)<<std::setprecision(3)
        <<std::setw(8)<<rasmol_state.x1[0]
        <<std::setw(8)<<rasmol_state.x1[1]
        <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
}

void read_rasmol_structure1(TMalignRasmolState& rasmol_state,
    const std::string& xname, const Vec3& t, const RotMat& u,
    const int ter_opt, const int mm_opt, const int mirror_opt);

inline void write_rasmol_mmcif_chain_b(std::ostream& os, TMalignRasmolState& rasmol_state, const int serial)
{
    os<<"ATOM  "<<std::setw(6)<<std::setw(5)<<serial%100000<<' '
        <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_x"]]
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_y"]]
        <<std::setw(8)<<rasmol_state.line_vec[rasmol_state._atom_site["Cartn_z"]]
        <<'\n';
}

void read_rasmol_structure2(TMalignRasmolState& rasmol_state,
    const std::string& yname, const int ter_opt, const int mm_opt);

inline void write_rasmol_to_file(const TMalignRasmolState& rasmol_state,
    const std::string& fname_super, const int mm_opt)
{
    std::ofstream fp;
    if (!mm_opt)
    {
        fp.open((fname_super).c_str());
        fp<<rasmol_state.buf.str();
        fp.close();
    }
    fp.open((fname_super+"_all").c_str());
    fp<<rasmol_state.buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open((fname_super+"_atm").c_str());
        fp<<rasmol_state.buf_atm.str();
        fp.close();
    }
    fp.open((fname_super+"_all_atm").c_str());
    fp<<rasmol_state.buf_all_atm.str();
    fp.close();
    fp.open((fname_super+"_all_atm_lig").c_str());
    fp<<rasmol_state.buf_all_atm_lig.str();
    fp.close();
}

void output_rasmol(const std::string xname, const std::string yname, const std::string fname_super, const Vec3& t,
    const RotMat& u, const int ter_opt, const int mm_opt, const int split_opt,
    const int mirror_opt, const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2, const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden);

// extract rotation matrix based on TMscore8
void output_rotation_matrix(const std::string& fname_matrix,
    const Vec3& t, const RotMat& u,
    std::ostream& os = std::cout);

void output_results(const ChainPairDesc& chain_pair,
    const ChainPairAlignResult& chain_align_res,
    const OutputOptions& out_opts,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    std::ostream& os = std::cout);

void output_mTMalign_results(const std::string xname, const std::string yname,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const ChainPairAlignResult& chain_align_res,
    const double Lnorm_ass, const double d0_scale,
    const std::string& fname_matrix, const int outfmt_opt, const int ter_opt,
    const int mm_opt, const int split_opt, const int o_opt,
    const std::string fname_super, const int i_opt, const int a_opt,
    const bool u_opt, const bool d_opt, const int mirror_opt,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const double ccTM_score,
    std::ostream& os = std::cout);

void output_cp(const std::string&xname, const std::string&yname,
    const std::string &seqxA, const int outfmt_opt, int &right_num,
    std::ostream& os = std::cout);
