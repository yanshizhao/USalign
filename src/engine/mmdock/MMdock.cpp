#include "MMdock.h"
#include "MMalign.h"
#include "common.h"
#include "output.h"
#include "se.h"

void mmdock_align_trimmed(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const TrimmedComplex& trimmed, int trim_idx,
    int mol_type, double Lnorm_pair, double d0_scale, double TMcut,
    bool fast_opt, int parallel_threads,
    std::vector<std::string>& sequence)
{
    int ylen_trim = trimmed.lengths[trim_idx];
    CoordArray ya_trim(ylen_trim);
    std::string seqy_trim;
    std::string secy_trim;
    secy_trim.resize(ylen_trim + 1);
    copy_chain_data(trimmed.coords[trim_idx], trimmed.seqs[trim_idx], trimmed.secs[trim_idx],
        ylen_trim, ya_trim, seqy_trim, secy_trim);

    ChainPairAlignOptions trim_opts;
    trim_opts.i_opt = 0;
    trim_opts.a_opt = false;
    trim_opts.u_opt = true;
    trim_opts.d_opt = false;
    trim_opts.fast_opt = fast_opt;
    trim_opts.se_opt = false;
    trim_opts.cp_opt = false;
    trim_opts.Lnorm = Lnorm_pair;
    trim_opts.d0_scale = d0_scale;
    trim_opts.TMcut = TMcut;
    trim_opts.parallel_threads = parallel_threads;
    trim_opts.ss_opt = 0;
    trim_opts.mol_type = mol_type;

    TMalign_main(xa, ya_trim, seqx, seqy_trim, secx, secy_trim,
        chain_align_res, xlen, ylen_trim, sequence, trim_opts);
    chain_align_res.seqxA.clear();
    chain_align_res.seqyA.clear();

    CoordArray xt(xlen);
    do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
    chain_align_res.invmap.assign(ylen + 1, 0);
    trim_opts.i_opt = 0;
    trim_opts.u_opt = 2;
    trim_opts.parallel_threads = parallel_threads;
    se_main(xt, ya, seqx, seqy,
        chain_align_res, xlen, ylen, sequence, trim_opts, 1);

    if (sequence.size() < 2) sequence.push_back("");
    if (sequence.size() < 2) sequence.push_back("");
    sequence[0] = chain_align_res.seqxA;
    sequence[1] = chain_align_res.seqyA;

    trim_opts.i_opt = 2;
    trim_opts.u_opt = true;
    trim_opts.parallel_threads = 1;
    TMalign_main(xt, ya, seqx, seqy, secx.c_str(), secy.c_str(),
        chain_align_res, xlen, ylen, sequence, trim_opts);
}

void run_mmdock_parallel(
    const DoubleCube& xa_vec, const DoubleCube& ya_vec,
    const CharMatrix& seqx_vec, const CharMatrix& seqy_vec,
    const CharMatrix& secx_vec, const CharMatrix& secy_vec,
    const std::vector<int>& xlen_vec, const std::vector<int>& ylen_vec,
    const std::vector<int>& mol_vec1, const std::vector<int>& mol_vec2,
    DoubleMatrix& TMave_mat, std::vector<std::vector<std::string>>& seqxA_mat,
    std::vector<std::vector<std::string>>& seqyA_mat, int chain1_num, int chain2_num,
    int len_aa, int len_na, int outfmt_opt, double TMcut, double d0_scale,
    bool fast_opt, const TrimmedComplex& trimmed, int parallel_threads = 1)
{
    #pragma omp parallel for schedule(dynamic, 1) num_threads(parallel_threads)
    for (int i = 0; i < chain1_num; i++)
    {
        int xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (int j = 0; j < chain2_num; j++)
                TMave_mat[i][j] = -1;   // no symmetric write
            continue;
        }

        std::vector<std::string> pair_sequence(2, "");

        CoordArray xa(xlen);
        std::string seqx, secx;
        secx.resize(xlen + 1);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i],
            xlen, xa, seqx, secx);

        for (int j = 0; j < chain2_num; j++)
        {
            // Skip protein-RNA cross-type alignment
            if (mol_vec1[i] * mol_vec2[j] < 0)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            int ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            CoordArray ya(ylen);
            std::string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j],
                ylen, ya, seqy, secy);

            int Lnorm_pair = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0) Lnorm_pair = len_na;

            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;

            if (trimmed.chain_count && trimmed.lengths[j] < ylen)
            {
                mmdock_align_trimmed(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, trimmed, j, mol_vec1[i] + mol_vec2[j],
                    Lnorm_pair, d0_scale, TMcut, fast_opt, 1, pair_sequence);
            }
            else
            {    // ---- no trimComplex branch ----
                ChainPairAlignOptions align_opts;
                align_opts.i_opt = 0;
                align_opts.a_opt = 0;
                align_opts.u_opt = 1;
                align_opts.d_opt = false;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = false;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_pair;
                align_opts.d0_scale = d0_scale;
                align_opts.TMcut = TMcut;
                align_opts.parallel_threads = 1;
                align_opts.ss_opt = 0;
                align_opts.mol_type = mol_vec1[i] + mol_vec2[j];

                align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy, xlen, ylen, align_opts, pair_sequence, outfmt_opt);
            }

            // Store chain_align_res (rectangular matrix, no symmetric write)
            seqxA_mat[i][j] = chain_align_res.seqxA;
            seqyA_mat[i][j] = chain_align_res.seqyA;
            TMave_mat[i][j] = chain_align_res.TM4 * Lnorm_pair;
        }
    }
}

void mmdock_assign_and_output(AlignCommonInput& common_inputs,
    const ComplexData& complex1, const ComplexData& complex2,
    DoubleMatrix& TMave_mat, std::vector<std::vector<std::string> >& seqxA_mat,
    std::vector<std::vector<std::string> >& seqyA_mat, int chain1_num, int chain2_num, bool fast_opt)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    int    i,j;                    // chain index
    int    xlen, ylen;             // chain length
    std::string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    std::string secx;                   // for the secondary structure
    std::string secy;

    std::vector<int> assign1_list(chain1_num);
    std::vector<int> assign2_list(chain2_num);
    enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);

    if (user_opts.outfmt_opt==0) print_version(USALIGN_VERSION);
    RotArray ut_mat; // rotation matrices for all-against-all alignment
    ut_mat.resize(chain1_num);
    int ui;
    int uj;
    std::vector<std::string>xname_vec;
    std::vector<std::string>yname_vec;
    std::vector<double>TM_vec;
    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        xname_vec.push_back(user_opts.xname+complex1.chain_ids[i]);
        if (j<0)
        {
            std::cerr<<"Warning! "<<complex1.chain_ids[i]<<" cannot be alighed"<<std::endl;
            for (ui=0;ui<3;ui++)
            {
                for (uj=0;uj<4;uj++) ut_mat[i][ui*3+uj]=0;
                ut_mat[i][ui*3+ui]=1;
            }
            yname_vec.push_back(user_opts.yname);
            continue;
        }
        yname_vec.push_back(user_opts.yname+complex2.chain_ids[j]);

        xlen =complex1.lengths[i];
        secx.resize(xlen+1);
        xa.clear();
        xa.reserve(xlen);
        copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i], xlen,xa,seqx,secx);

        ylen =complex2.lengths[j];
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(complex2.coords[j],complex2.seqs[j],complex2.secs[j], ylen,ya,seqy,secy);

        ChainPairAlignResult chain_align_res = { 0};
        chain_align_res.d0_out = 5.0;
        ChainPairAlignOptions align_opts;
        fill_chain_pair_align_options(align_opts, user_opts, common_inputs.control_options);
        align_opts.i_opt = 3;
        align_opts.fast_opt = fast_opt;
        align_opts.se_opt = false;
        align_opts.cp_opt = false;
        align_opts.TMcut = -1;
        align_opts.parallel_threads = 1;
        align_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];

        int c;
        for (c=0; c<parsed_input.sequence.size(); c++) parsed_input.sequence[c].clear();
        parsed_input.sequence.clear();
        parsed_input.sequence.push_back(seqxA_mat[i][j]);
        parsed_input.sequence.push_back(seqyA_mat[i][j]);
            
        align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
            xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);
        
        for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=chain_align_res.u0[ui][uj];
        for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=chain_align_res.t0[uj];

        TM_vec.push_back(chain_align_res.TM1);
        TM_vec.push_back(chain_align_res.TM2);

        OutputOptions out_opts;
        fill_out_opts_from_user_opts(out_opts, user_opts);
        out_opts.fname_matrix = "";
        out_opts.fname_super = "";
        out_opts.o_opt = 0;
        const std::string structure1_display = user_opts.xname;
        const std::string structure2_display = user_opts.yname;
        const ChainPairDesc pair_desc{structure1_display, structure2_display, complex1.chain_ids[i], complex2.chain_ids[j], xlen, ylen};
        if (user_opts.outfmt_opt<2) output_results(
            pair_desc, chain_align_res, out_opts, complex1.resi, complex2.resi);
        
        }
    if (user_opts.outfmt_opt==2)
    {
        double TM=0;
        for (i=0;i<TM_vec.size();i++) TM+=TM_vec[i]*TM_vec[i];
        TM=sqrt(TM/TM_vec.size());
        std::string query_name=user_opts.xname;
        std::string template_name=user_opts.yname;

        for (i=0;i<chain1_num;i++)
        {
            j=assign1_list[i];
            if (j<0) continue;
            query_name   +=complex1.chain_ids[i];
            template_name+=complex2.chain_ids[j];
        }
        fcout("%s\t%s\t%.4f\n", query_name, template_name, TM);
    }

    if (user_opts.m_opt) output_dock_rotation_matrix(user_opts.fname_matrix,
        xname_vec,yname_vec, ut_mat, assign1_list);

    if (user_opts.o_opt) output_superposed_chains(parsed_input.chain1_list, user_opts.ter_opt, user_opts.split_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, user_opts.mirror_opt, ut_mat, user_opts.fname_super);
}

void mmdock_pairwise_alignment(AlignCommonInput& common_inputs, ComplexData& complex1, ComplexData& complex2,
    const TrimmedComplex& trimmed, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >& seqxA_mat, std::vector<std::vector<std::string> >& seqyA_mat,
    const int chain1_num, const int chain2_num, const int len_aa, const int len_na, const bool fast_opt)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    int i;
    int j;
    int xlen;
    int ylen;
    CoordArray xa;
    CoordArray ya;
    std::string seqx;
    std::string seqy;
    std::string secx;
    std::string secy;

    bool mmdock_parallel_done = false;
#ifdef _OPENMP
    if (ctrl_opts.parallel_threads > 1 && chain1_num > 1) {
        run_mmdock_parallel(
            complex1.coords, complex2.coords, complex1.seqs, complex2.seqs,
            complex1.secs, complex2.secs, complex1.lengths, complex2.lengths,
            complex1.mol_types, complex2.mol_types, TMave_mat,
            seqxA_mat, seqyA_mat, chain1_num, chain2_num, len_aa, len_na,
            user_opts.outfmt_opt, user_opts.TMcut, user_opts.d0_scale, fast_opt,
            trimmed, ctrl_opts.parallel_threads);
        mmdock_parallel_done = true;
    }
#endif
    if (!mmdock_parallel_done)
    {
        for (i=0;i<chain1_num;i++)
        {
            xlen=complex1.lengths[i];
            if (xlen<3)
            {
                for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
                continue;
            }
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(complex1.coords[i],complex1.seqs[i],complex1.secs[i],
                xlen,xa,seqx,secx);

            for (j=0;j<chain2_num;j++)
            {
                if (complex1.mol_types[i]*complex2.mol_types[j]<0) //no protein-RNA alignment
                {
                    TMave_mat[i][j]=-1;
                    continue;
                }

                ylen=complex2.lengths[j];
                if (ylen<3)
                {
                    TMave_mat[i][j]=-1;
                    continue;
                }
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(complex2.coords[j],complex2.seqs[j],complex2.secs[j],
                    ylen,ya,seqy,secy);

                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                int Lnorm_pair=len_aa;
                if (complex1.mol_types[i]+complex2.mol_types[j]>0) Lnorm_pair=len_na;

                if (trimmed.chain_count && trimmed.lengths[j]<ylen)
                {
                    mmdock_align_trimmed(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, trimmed, j, complex1.mol_types[i]+complex2.mol_types[j],
                        Lnorm_pair, user_opts.d0_scale, user_opts.TMcut, fast_opt,
                        ctrl_opts.parallel_threads, parsed_input.sequence);
                }
                else
                {
                    ChainPairAlignOptions align_opts;
                    align_opts.i_opt = 0;
                    align_opts.a_opt = 0;
                    align_opts.u_opt = 1;
                    align_opts.d_opt = false;
                    align_opts.fast_opt = fast_opt;
                    align_opts.se_opt = false;
                    align_opts.cp_opt = false;
                    align_opts.Lnorm = Lnorm_pair;
                    align_opts.d0_scale = user_opts.d0_scale;
                    align_opts.TMcut = user_opts.TMcut;
                    align_opts.parallel_threads = 1;
                    align_opts.ss_opt = 0;
                    align_opts.mol_type = complex1.mol_types[i]+complex2.mol_types[j];
                    align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);
                }

                seqxA_mat[i][j]=chain_align_res.seqxA;
                seqyA_mat[i][j]=chain_align_res.seqyA;
                TMave_mat[i][j]=chain_align_res.TM4*Lnorm_pair;
            }
        }
    }
}

int MMdock(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    bool fast_opt = user_opts.fast_opt;

    ComplexData complex1;
    ComplexData complex2;
    int    i;                    // chain index
    int    xlen, ylen;             // chain length
    std::string seqx, seqy;             // for the protein sequence
    CoordArray xa;                     // structure of single chain
    CoordArray ya;
    std::string secx;                   // for the secondary structure
    std::string secy;

    // parse complex
    parse_chain_list(parsed_input.chain1_list, complex1,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, parsed_input.normalize_atom_name, user_opts.mirror_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    if (complex1.coords.size()==0) PrintErrorAndQuit("ERROR! 0 individual chain");

    parse_chain_list(parsed_input.chain2_list, complex2,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt2_opt,
        user_opts.atom_opt, parsed_input.normalize_atom_name, 0, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2);
    if (complex1.coords.size()>complex2.coords.size())  PrintErrorAndQuit("ERROR! more individual chains to align than number of chains in complex template");

    int len_aa=getmin(complex1.total_len_aa,complex2.total_len_aa);
    int len_na=getmin(complex1.total_len_na,complex2.total_len_na);
    if (user_opts.a_opt)
    {
        len_aa=(complex1.total_len_aa+complex2.total_len_aa)/2;
        len_na=(complex1.total_len_na+complex2.total_len_na)/2;
    }

    // perform monomer alignment if there is only one chain
    if (complex1.coords.size()==1 && complex2.coords.size()==1)
    {
        xlen = complex1.lengths[0];
        ylen = complex2.lengths[0];
        secx.resize(xlen+1);
        secy.resize(ylen+1);
        xa.clear();
        xa.reserve(xlen);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(complex1.coords[0],complex1.seqs[0],complex1.secs[0], xlen,xa,seqx,secx);
        copy_chain_data(complex2.coords[0],complex2.seqs[0],complex2.secs[0], ylen,ya,seqy,secy);

        ChainPairAlignResult chain_align_res = { 0};
        chain_align_res.d0_out = 5.0;
        ChainPairAlignOptions align_opts;
        fill_chain_pair_align_options(align_opts, user_opts, common_inputs.control_options);
        align_opts.i_opt = 0;
        align_opts.fast_opt = fast_opt;
        align_opts.se_opt = false;
        align_opts.cp_opt = false;
        align_opts.parallel_threads = 1;
        align_opts.mol_type = complex1.mol_types[0]+complex2.mol_types[0];

        align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy, xlen, ylen, align_opts, parsed_input.sequence, user_opts.outfmt_opt);

        OutputOptions out_opts;
        fill_out_opts_from_user_opts(out_opts, user_opts);
        out_opts.outfmt_opt = (user_opts.outfmt_opt==2?user_opts.outfmt_opt:3);
        out_opts.mm_opt = true;

        const std::string structure1_display = user_opts.xname.substr(user_opts.dir1_opt.size());
        const std::string structure2_display = user_opts.yname.substr(user_opts.dir2_opt.size());
        const ChainPairDesc pair_desc{structure1_display, structure2_display, complex1.chain_ids[0], complex2.chain_ids[0], xlen, ylen};
        output_results(
            pair_desc, chain_align_res, out_opts, complex1.resi, complex2.resi);
        if (user_opts.outfmt_opt==2) fcout("%s%s\t%s%s\t%.4f\n",
            user_opts.xname.substr(user_opts.dir1_opt.size()), complex1.chain_ids[0],
            user_opts.yname.substr(user_opts.dir2_opt.size()), complex2.chain_ids[0],
            sqrt((chain_align_res.TM1*chain_align_res.TM1+chain_align_res.TM2*chain_align_res.TM2)/2));

        return 0;
    }

    // declare TM-score tables
    int chain1_num=complex1.coords.size();
    int chain2_num=complex2.coords.size();
    std::vector<std::string> tmp_str_vec(chain2_num,"");
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain1_num,std::vector<double>(chain2_num));
    std::vector<std::vector<std::string> >seqxA_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> > seqM_mat(chain1_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_mat(chain1_num,tmp_str_vec);

    // trimComplex
    TrimmedComplex trimmed;
    trimmed.max_aa_len=0;
    trimmed.max_na_len=0;
    for (i=0;i<chain1_num;i++)
    {
        xlen=complex1.lengths[i];
        if      (complex1.mol_types[i]>0  && xlen>trimmed.max_na_len) trimmed.max_na_len=xlen;
        else if (complex1.mol_types[i]<=0 && xlen>trimmed.max_aa_len) trimmed.max_aa_len=xlen;
    }
    trimComplex(trimmed.coords,trimmed.seqs,
        trimmed.secs,trimmed.lengths,complex2.coords,complex2.seqs,complex2.secs,complex2.lengths,
        complex2.mol_types,trimmed.max_aa_len,trimmed.max_na_len,trimmed.chain_count);

    // Auto-enable fast mode BEFORE parallel entry to keep both paths consistent
    if (len_aa + len_na > 500) fast_opt = true;

    mmdock_pairwise_alignment(common_inputs, complex1, complex2, trimmed, TMave_mat,
        seqxA_mat, seqyA_mat, chain1_num, chain2_num, len_aa, len_na, fast_opt);
    DoubleCube().swap(trimmed.coords);
    CharMatrix().swap(trimmed.seqs);
    CharMatrix().swap(trimmed.secs);
    std::vector<int> ().swap(trimmed.lengths);

    mmdock_assign_and_output(common_inputs, complex1, complex2, TMave_mat,
        seqxA_mat, seqyA_mat, chain1_num, chain2_num, fast_opt);

    std::vector<std::vector<std::string> >().swap(seqxA_mat);
    std::vector<std::vector<std::string> >().swap(seqM_mat);
    std::vector<std::vector<std::string> >().swap(seqyA_mat);
    return 1;
}

