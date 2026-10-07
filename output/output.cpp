#include "output.h"
#include "tmalign_output.h"

#include <iomanip>
#include <iostream>

const char USALIGN_BANNER[] =
"\n"
" ********************************************************************\n"
" * US-align (Version 20260819)                                      *\n"
" * Universal Structure Alignment of Proteins and Nucleic Acids      *\n"
" * Reference: C Zhang, L Freddolino, Y Zhang. (2026) Nat Protoc     *\n"
" *            C Zhang, M Shine, AM Pyle, Y Zhang. (2022) Nat Methods*\n"
" *            C Zhang, AM Pyle (2022) iScience.                     *\n"
" * Please email comments and suggestions to zhang@zhanggroup.org    *\n"
" ********************************************************************";

// ---- Output one chain-pair alignment result (version / CP / results / do-block) ----
void output_pair_alignment(const AlignCommonInput& common_inputs,
    int xlen, int ylen,
    const ChainPairAlignResult& chain_align_res,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::string& chain_id1, const std::string& chain_id2,
    const std::vector<std::string>& resi_vec1,
    const std::vector<std::string>& resi_vec2,
    const std::vector<double>* dist_arg)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    if (user_opts.outfmt_opt == 0) print_banner(USALIGN_BANNER);
    int right_num=0;
    if (ctrl_opts.cp_opt) output_cp(
        user_opts.xname.substr(user_opts.dir1_opt.size() + user_opts.dir_opt.size()),
        user_opts.yname.substr(user_opts.dir2_opt.size() + user_opts.dir_opt.size()),
        chain_align_res.seqxA, user_opts.outfmt_opt, right_num);
    OutputOptions out_opts;
    fill_out_opts_from_user_opts(out_opts, user_opts);
    const std::string structure1_display = user_opts.xname.substr(
        user_opts.dir1_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size());
    const std::string structure2_display = user_opts.yname.substr(
        user_opts.dir2_opt.size() + user_opts.dir_opt.size() + user_opts.dirpair_opt.size());
    output_results(
        ChainPairDesc{structure1_display, structure2_display, chain_id1, chain_id2, xlen, ylen},
        chain_align_res, out_opts, resi_vec1, resi_vec2);
    if (dist_arg) output_soi_dist_table(std::cout, ylen,
        chain_align_res.invmap, *dist_arg, pdb_lines1, pdb_lines2);
    else if (ctrl_opts.do_opt || (ctrl_opts.cp_opt && user_opts.outfmt_opt <= 0))
    {
        output_do_distance_table(std::cout, chain_align_res.seqxA, chain_align_res.seqyA,
            pdb_lines1, pdb_lines2, chain_align_res.do_vec, right_num);
    }
}

void output_do_distance_table(std::ostream& os,
    const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2,
    const std::vector<double>& do_vec,
    size_t right_num)
{
    os << "###############\t###############\t#########" << std::endl;
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#" << std::endl;
    size_t r1 = right_num;
    size_t r2 = 0;
    int    postcp = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        r1 += seqxA[r] != '-';
        r2 += seqyA[r] != '-';
        if (seqxA[r] == '*')
        {
            os << "###### Circular\tPermutation ###\t#########\n";
            r1 = 0;
            postcp = 1;
        }
        else if (seqxA[r] != '-' && seqyA[r] != '-')
        {
            os << pdb_lines1[r1 - 1].substr(12, 15) << '\t'
               << pdb_lines2[r2 - 1].substr(12, 15) << '\t'
               << std::setw(9) << std::setiosflags(std::ios::fixed) << std::setprecision(3)
               << do_vec[r - postcp] << '\n';
        }
    }
    os << "###############\t###############\t#########" << std::endl;
}

void output_soi_dist_table(std::ostream& os, int ylen,
    const std::vector<int>& invmap, const std::vector<double>& dist_list,
    const std::vector<std::string>& pdb_lines1,
    const std::vector<std::string>& pdb_lines2)
{
    os << "###############\t###############\t#########\n";
    os << "#Aligned atom 1\tAligned atom 2 \tDistance#\n";
    int r1;
    int r2;
    for (r2=0;r2<ylen;r2++)
    {
        r1=invmap[r2];
        if (r1<0) continue;
        os<<pdb_lines1[r1].substr(12,15)<<'\t'
            <<pdb_lines2[r2].substr(12,15)<<'\t'
            <<std::setw(9)<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)
            <<dist_list[r2]<<'\n';
    }
    os << "###############\t###############\t#########\n";
}
