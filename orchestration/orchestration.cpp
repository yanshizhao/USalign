#include "orchestration.h"

#include "banner.h"
#include "output_do.h"
#include "se.h"

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

// ---- Single chain-pair structure alignment dispatcher (cp / se / TMalign_main) ----
void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (align_opts.cp_opt) CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
    else if (align_opts.se_opt)
    {
        chain_align_res.invmap.assign(ylen + 1, 0);
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);
        se_main(xa, ya, seqx, seqy,
            chain_align_res,
            xlen, ylen, sequence, align_opts, outfmt_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0; r2<ylen; r2++)
            {
                r1 = chain_align_res.invmap[r2];
                if (r1 < 0) continue;
                chain_align_res.L_ali += 1;
                chain_align_res.Liden += (seqx[r1] == seqy[r2]);
            }
        }
    }
    else TMalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
}

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
