#include "mTMalign.h"
#include "MMalign.h"
#include "common.h"
#include "output.h"
#include "UPGMA.h"
#include "se.h"

// ============ Helper functions for ccTM-score ============
//Build msa_to_res mapping: msa_to_res[i][l] = original residue index at MSA column l (-1=gap)
static void build_msa_to_res(const std::vector<std::string>& msa_seqs, int chain_num, int msa_len, IntMatrix& msa_to_res) 
{
    msa_to_res.assign(chain_num, std::vector<int>(msa_len, -1));
    for (int i = 0; i < chain_num; i++) 
    {
        int res_idx = 0;
        for (int l = 0; l < msa_len; l++) {
            if (msa_seqs[i][l] != '-') {
                msa_to_res[i][l] = res_idx;
                res_idx++;
            }
        }
    }
}

//Find candidate columns: all structures have residues
static void find_no_gap_cols(const IntMatrix& msa_to_res, int chain_num, int msa_len, std::vector<int>& no_gap_cols) 
{
    for (int l = 0; l < msa_len; l++) {
        int count = 0;
        for (int i = 0; i < chain_num; i++)
            if (msa_to_res[i][l] != -1) count++;
        if (count == chain_num)
            no_gap_cols.push_back(l);
    }
}

// Select reference structure (most shared columns)
static int select_ref_struct(const IntMatrix& msa_to_res, int msa_len, int chain_num) 
{
    int ref_idx = 0, max_cnt = 0;
    for (int i = 0; i < chain_num; i++) 
    {
        int cnt = 0;
        for (int l = 0; l < msa_len; l++) 
        {
            if (msa_to_res[i][l] == -1) continue;
            for (int k = 0; k < chain_num; k++)
            {
                if (k != i && msa_to_res[k][l] != -1) 
                { 
                    cnt++; 
                    break; 
                }
            }
        }
        if (cnt > max_cnt) 
        { 
            max_cnt = cnt; 
            ref_idx = i; 
        }
    }
    return ref_idx;
}

// Coordinate re-alignment: from ua_vec to eliminate chain error
static void align_to_ref(const DoubleCube& orig_coords, const IntMatrix& msa_to_res, const std::vector<int>& no_gap_cols, int chain_num, int ref_idx, std::vector<CoordArray>& aligned_coord) 
{
    aligned_coord.resize(chain_num);

    // Pre-compute reference coordinates for all candidate columns
    CoordArray ref_core_coords;
    for (int k = 0; k < (int)no_gap_cols.size(); k++) 
    {
        int col_idx = no_gap_cols[k];
        int ref_res_idx = msa_to_res[ref_idx][col_idx];
        ref_core_coords.push_back({ { orig_coords[ref_idx][ref_res_idx][0], orig_coords[ref_idx][ref_res_idx][1], orig_coords[ref_idx][ref_res_idx][2] } });
    }

    for (int i = 0; i < chain_num; i++)
    {
        CoordArray cur_struct_coords;
        for (int r = 0; r < (int)orig_coords[i].size(); r++)
            cur_struct_coords.push_back({ { orig_coords[i][r][0], orig_coords[i][r][1], orig_coords[i][r][2] } });

        if (i == ref_idx)
        {
            aligned_coord[i] = cur_struct_coords;
            continue;
        }

        CoordArray cur_core_coords;
        for (int k = 0; k < (int)no_gap_cols.size(); k++) {
            int col_idx = no_gap_cols[k];
            int cur_res_idx = msa_to_res[i][col_idx];
            cur_core_coords.push_back({ { orig_coords[i][cur_res_idx][0], orig_coords[i][cur_res_idx][1], orig_coords[i][cur_res_idx][2] } });
        }

        if (cur_core_coords.size() < 4)
        {
            aligned_coord[i] = cur_struct_coords;
            continue;
        }
        Vec3 t; RotMat u; double rms;
        Kabsch(cur_core_coords, ref_core_coords, (int)cur_core_coords.size(), 1, rms, t, u);
        CoordArray rotated_coords(cur_struct_coords.size());
        do_rotation(cur_struct_coords, rotated_coords, (int)cur_struct_coords.size(), t, u);
        aligned_coord[i] = rotated_coords;
    }
}

//Calculate pairwise CA distances at each candidate column
static void calc_pairwise_distances(const std::vector<CoordArray>& aligned_coord, const IntMatrix& msa_to_res, const std::vector<int>& no_gap_cols, int chain_num, int tot_num_pair, DoubleMatrix& ca_dist_matrix) 
{
    ca_dist_matrix.assign(tot_num_pair, std::vector<double>(no_gap_cols.size(), -1));
    for (int p = 0, i = 0; i < chain_num; i++) 
    {
        for (int j = i + 1; j < chain_num; j++, p++) 
        {
            for (int k = 0; k < (int)no_gap_cols.size(); k++) 
            {
                int l = no_gap_cols[k];
                if (msa_to_res[i][l] != -1 && msa_to_res[j][l] != -1) 
                {
                    double dx = aligned_coord[i][msa_to_res[i][l]][0] - aligned_coord[j][msa_to_res[j][l]][0];
                    double dy = aligned_coord[i][msa_to_res[i][l]][1] - aligned_coord[j][msa_to_res[j][l]][1];
                    double dz = aligned_coord[i][msa_to_res[i][l]][2] - aligned_coord[j][msa_to_res[j][l]][2];
                    ca_dist_matrix[p][k] = sqrt(dx*dx + dy*dy + dz*dz);
                }
            }
        }
    }
}

//Filter Common Core columns: all pair distances <= 4.0 Å
static void select_common_core_cols(const std::vector<std::vector<double>>& ca_dist_matrix, const std::vector<int>& no_gap_cols, int tot_num_pair, std::vector<int>& core_cols) 
{
    for (int k = 0; k < (int)no_gap_cols.size(); k++) 
    {
        bool col_pass = true;
        for (int p = 0; p < tot_num_pair; p++) 
        {
            if (ca_dist_matrix[p][k] > 4.0) 
            { 
                col_pass = false; 
                break; 
            }
        }
        if (col_pass) core_cols.push_back(no_gap_cols[k]);
    }
}

//Extract CC residue indices: core_res_idx[i][j] = msa_to_res[i][core_cols[j]]
static void extract_common_core_res_idx(const IntMatrix& msa_to_res, const std::vector<int>& core_cols, int chain_num, IntMatrix& core_res_idx) 
{
    core_res_idx.assign(chain_num, std::vector<int>((int)core_cols.size(), -1));
    for (int i = 0; i < chain_num; i++)
    {
        for (int j = 0; j < (int)core_cols.size(); j++)
        {
            core_res_idx[i][j] = msa_to_res[i][core_cols[j]];
        }
    }
}

//Calculate pairwise TM-score using only CC residues
static double calc_common_core_TM_sum(const DoubleCube& orig_coords, const IntMatrix& core_res_idx, const std::vector<int>& core_cols, int chain_num, const std::vector<int>& len_vec, int cur_complex_mol_list) 
{
    double TM_sum = 0;
    int simplify_step = 1;
    int score_sum_method = 0;
    double score_d8 = 0;
    int core_col_count = (int)core_cols.size();

    for (int i = 0; i < chain_num; i++) 
    {
        for (int j = 0; j < i; j++) 
        {
            CoordArray r1(core_col_count), r2(core_col_count), xtm(core_col_count), ytm(core_col_count), xt(core_col_count);
            for (int k = 0; k < core_col_count; k++)
            {
                int ri = core_res_idx[i][k]; 
                int rj = core_res_idx[j][k];
                xtm[k][0] = orig_coords[i][ri][0]; 
                xtm[k][1] = orig_coords[i][ri][1]; 
                xtm[k][2] = orig_coords[i][ri][2];
                ytm[k][0] = orig_coords[j][rj][0]; 
                ytm[k][1] = orig_coords[j][rj][1]; 
                ytm[k][2] = orig_coords[j][rj][2];
            }

            double pair_Lnorm = std::min(len_vec[i], len_vec[j]);
            double D0_MIN, Lnorm_out, d0, d0_search;
            parameter_set4final(pair_Lnorm, D0_MIN, Lnorm_out, d0, d0_search, cur_complex_mol_list);
            Vec3 t0; 
            RotMat u0; 
            double rmsd;

            double TM = TMscore8_search(r1, r2, xtm, ytm, xt, core_col_count,
                t0, u0, simplify_step, score_sum_method, rmsd,
                d0_search, pair_Lnorm, score_d8, d0);
            TM_sum += TM;
        }
    }
    return TM_sum;
}

// ============ ccTM-score calculation ============
// Calculate average pairwise TM-score based on Common Core columns
// @param  orig_coords     - original coordinates (for TM-score calculation)
// @param  seqxA_mat  - pairwise alignment matrix (diagonal = MSA sequences)
// @param  chain_num  - number of structures
// @param  len_vec    - length of each structure
// @param  cur_complex_mol_list   - molecule type (0=protein, >0=RNA)
// @return ccTM-score (0 if no Common Core found)
double calc_ccTM_score(
    const DoubleCube& orig_coords,
    const std::vector<std::vector<std::string>>& seqxA_mat,
    int chain_num,
    const std::vector<int>& len_vec,
    int cur_complex_mol_list)
{
    // Extract MSA sequences from seqxA_mat diagonal
    std::vector<std::string> msa_seqs(chain_num);
    for (int i = 0; i < chain_num; i++) {
        msa_seqs[i] = seqxA_mat[i][i];
    }

    // 1. Build msa_to_res: msa_to_res[struct][col] = original residue index (-1 = gap)
    int msa_len = (int)msa_seqs[0].size();
    IntMatrix msa_to_res;
    build_msa_to_res(msa_seqs, chain_num, msa_len, msa_to_res);

    // 2. Find no-gap columns: all structures have residues at these columns
    std::vector<int> no_gap_cols;
    find_no_gap_cols(msa_to_res, chain_num, msa_len, no_gap_cols);
    if (no_gap_cols.empty())
    {
        fcout(std::cerr, "Warning! No no-gap MSA columns in ccTM-score; "
                         "the ccTM-score is 0.\n");
        return 0.0;
    }
    // Kabsch re-alignment needs >= 4 non-colinear points; with fewer no-gap
    // columns the alignment is degenerate, so the ccTM-score is set to 0.
    if ((int)no_gap_cols.size() < 4)
    {
        fcout(std::cerr, "Warning! Only %d no-gap MSA column(s) (< 4) in "
                         "ccTM-score: Kabsch re-alignment is degenerate; the "
                         "ccTM-score is set to 0.\n",
              (int)no_gap_cols.size());
        return 0.0;
    }

    // 3. Re-align all structures from original coords to reference (eliminate chain error)
    //    select_ref_struct => align_to_ref (Kabsch + do_rotation)
    int ref_idx = select_ref_struct(msa_to_res, msa_len, chain_num);
    std::vector<CoordArray> aligned_coord;
    align_to_ref(orig_coords, msa_to_res, no_gap_cols, chain_num, ref_idx, aligned_coord);

    // 4. Compute pairwise CA distances between every pair at no-gap columns
    int tot_num_pair = chain_num * (chain_num - 1) / 2;
    if(tot_num_pair == 0) return 0.0;
    DoubleMatrix ca_dist_matrix;
    calc_pairwise_distances(aligned_coord, msa_to_res, no_gap_cols, chain_num, tot_num_pair, ca_dist_matrix);

    // 5. Select Common Core columns: all pairs have CA distance <= 4A
    std::vector<int> core_cols;
    select_common_core_cols(ca_dist_matrix, no_gap_cols, tot_num_pair, core_cols);
    if (core_cols.empty())
    {
        fcout(std::cerr, "Warning! No common core columns in ccTM-score (all "
                         "pairs have CA distance > 4.0 A at every no-gap MSA "
                         "column, usually caused by a distant structure); the "
                         "ccTM-score is 0.\n");
        return 0.0;
    }

    // 6. Extract original residue indices at Common Core columns
    IntMatrix core_res_idx;
    extract_common_core_res_idx(msa_to_res, core_cols, chain_num, core_res_idx);

    // 7. Compute all pairwise TM-scores using ONLY Common Core residues, sum them
    double TM_sum = calc_common_core_TM_sum(orig_coords, core_res_idx, core_cols, chain_num, len_vec, cur_complex_mol_list);
    return TM_sum / tot_num_pair;
}


// Parallel all-against-all pairwise alignment for MSTA (mTMalign)
static void run_mTMalign_pairwise_parallel(
    const DoubleCube& a_vec, const CharMatrix& seq_vec, const CharMatrix& sec_vec,
    const std::vector<int>& len_vec,
    const std::vector<std::string>& chain_list, const std::vector<std::string>& chainID_list,
    std::vector<std::vector<std::string>>& seqxA_mat, std::vector<std::vector<std::string>>& seqyA_mat,
    DoubleMatrix& TMave_mat, const std::vector<std::string>& resi_vec,
    int chain_num, double Lnorm_ass, double d0_scale,
    bool u_opt, int cur_complex_mol_list, int outfmt_opt, bool fast_opt, double TMcut,
    bool full_opt, bool se_opt, int ter_opt, int split_opt, int o_opt, int a_opt, bool d_opt,
    int parallel_threads)
{

    std::vector<std::string> full_out(chain_num * chain_num);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(parallel_threads)
    for (int chain_i = 0; chain_i < chain_num; chain_i++)
    {
        int xlen = len_vec[chain_i];
        if (xlen < 3) continue;
        CoordArray xa(xlen);
        std::string seqx, secx;
        secx.resize(xlen + 1);
        copy_chain_data(a_vec[chain_i], seq_vec[chain_i], sec_vec[chain_i], xlen, xa, seqx, secx);
        seqxA_mat[chain_i][chain_i] = seqyA_mat[chain_i][chain_i] = seqx;
        for (int chain_j = chain_i + 1; chain_j < chain_num; chain_j++)
        {
            int ylen = len_vec[chain_j];
            if (ylen < 3) continue;
            CoordArray ya(ylen);
            std::string seqy, secy;
            secy.resize(ylen + 1);
            copy_chain_data(a_vec[chain_j], seq_vec[chain_j], sec_vec[chain_j], ylen, ya, seqy, secy);

            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            std::vector<std::string> local_seq(2);
            local_seq[0] = seqxA_mat[chain_i][chain_j];
            local_seq[1] = seqyA_mat[chain_i][chain_j];
            ChainPairAlignOptions align_opts;
            align_opts.i_opt = 0;
            align_opts.a_opt = 0;
            align_opts.u_opt = u_opt;
            align_opts.d_opt = false;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = se_opt;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.d0_scale = d0_scale;
            align_opts.TMcut = TMcut;
            align_opts.parallel_threads = 1;
            align_opts.ss_opt = 0;
            align_opts.mol_type = cur_complex_mol_list;
            align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, local_seq, outfmt_opt);

            TMave_mat[chain_i][chain_j] = chain_align_res.TM4; TMave_mat[chain_j][chain_i] = chain_align_res.TM4;
            seqxA_mat[chain_i][chain_j] = seqyA_mat[chain_j][chain_i] = chain_align_res.seqxA;
            seqyA_mat[chain_i][chain_j] = seqxA_mat[chain_j][chain_i] = chain_align_res.seqyA;
            if (full_opt)
            {
                std::stringstream ss;
                OutputOptions out_opts;
                out_opts.fname_matrix = "";
                out_opts.fname_super = "";
                out_opts.Lnorm_ass = Lnorm_ass;
                out_opts.d0_scale = d0_scale;
                out_opts.outfmt_opt = outfmt_opt;
                out_opts.ter_opt = ter_opt;
                out_opts.mm_opt = true;
                out_opts.split_opt = split_opt;
                out_opts.o_opt = o_opt;
                out_opts.i_opt = 0;
                out_opts.a_opt = a_opt;
                out_opts.mirror_opt = 0;
                out_opts.u_opt = false;
                out_opts.d_opt = d_opt;
                const std::string structure1_display = chain_list[chain_i];
                const std::string structure2_display = chain_list[chain_j];
                output_results(
                    ChainPairDesc{structure1_display, structure2_display, chainID_list[chain_i], chainID_list[chain_j], xlen, ylen},
                    chain_align_res, out_opts, resi_vec, resi_vec, ss);
                full_out[chain_i * chain_num + chain_j] = ss.str();
            }
        }
    }

    // Serial output in (chain_i, chain_j) order, matching the serial path
    if (full_opt)
    {
        for (int chain_i = 0; chain_i < chain_num; chain_i++)
            for (int chain_j = chain_i + 1; chain_j < chain_num; chain_j++)
                if (!full_out[chain_i * chain_num + chain_j].empty())
                    std::cout << full_out[chain_i * chain_num + chain_j];
    }
}

// ---- Pick the chain with the largest column sum of the TM-score matrix (mm4) ----
int select_representative(const DoubleMatrix& TMave_mat, int chain_num)
{
    std::vector<double> TMave_list(chain_num, 0.0);
    for (int j=0; j<chain_num; j++) TMave_list[j]=0;
    for (int i=0; i<chain_num; i++)
        for (int j=0; j<chain_num; j++)
            TMave_list[j]+=TMave_mat[i][j];
    int repr_idx=0;
    double repr_TM=0;
    for (int j=0; j<chain_num; j++)
    {
        if (TMave_list[j]<repr_TM) continue;
        repr_TM=TMave_list[j];
        repr_idx=j;
    }
    return repr_idx;
}

struct MstaIterBuffers
{
    CoordArray member_coords;
    CoordArray partner_coords;
    std::string member_sec;
    std::string partner_sec;
    int member_len=0;
    int partner_len=0;
    int member_chain_idx=0;
    int partner_chain_idx=0;
    int residue_idx=0;
    int tm_order_idx=0;
};

struct MstaIterParams
{
    UserOptions& user_opts;
    ParsedInput& parsed_input;
    ControlOptions& ctrl_opts;
    const ComplexData& complex;
    DoubleCube& a_vec;
    const DoubleMatrix& TMave_mat;
    std::vector<std::vector<std::string> >& seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat;
    const std::vector<std::string>& chainID_list;
    const std::vector<std::string>& resi_vec;
    int chain_num;
    double Lnorm_ass;
    bool u_opt;
    bool fast_opt;
    int mol_type_sum;
    MstaIterBuffers buffers;

    MstaIterParams(AlignCommonInput& common_inputs, ComplexData& complex_in,
        DoubleCube& a_vec_in, const DoubleMatrix& TMave_mat_in,
        std::vector<std::vector<std::string> >& seqxA_mat_in, std::vector<std::vector<std::string> >& seqyA_mat_in,
        int chain_num_in, double Lnorm_ass_in, bool u_opt_in, bool fast_opt_in,
        int cur_complex_mol_list_in)
        : user_opts(common_inputs.user_options),
          parsed_input(common_inputs.parsed_input),
          ctrl_opts(common_inputs.control_options),
          complex(complex_in),
          a_vec(a_vec_in),
          TMave_mat(TMave_mat_in),
          seqxA_mat(seqxA_mat_in),
          seqyA_mat(seqyA_mat_in),
          chainID_list(complex_in.chain_ids),
          resi_vec(complex_in.resi),
          chain_num(chain_num_in),
          Lnorm_ass(Lnorm_ass_in),
          u_opt(u_opt_in),
          fast_opt(fast_opt_in),
          mol_type_sum(cur_complex_mol_list_in) {}
};


void msta_superpose_to_rep(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    ControlOptions& ctrl_opts = msta_iter_params.ctrl_opts;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    const std::vector<std::string>& chainID_list = msta_iter_params.chainID_list;
    const std::vector<std::string>& resi_vec = msta_iter_params.resi_vec;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    const DoubleMatrix& TMave_mat = msta_iter_params.TMave_mat;
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;

        user_opts.yname=parsed_input.chain1_list[msta_iter_state.repr_idx].substr(user_opts.dir_opt.size())+chainID_list[msta_iter_state.repr_idx];
        CoordArray xt;
        msta_iter_state.TM_pair_vec.clear();

        for (i=0; i<chain_num; i++) msta_iter_state.assign_list[i]=-1;
        msta_iter_state.assign_list[msta_iter_state.repr_idx]=msta_iter_state.repr_idx;
        for (r=0;r<parsed_input.sequence.size();r++) parsed_input.sequence[r].clear(); parsed_input.sequence.clear();
        parsed_input.sequence.push_back("");
        parsed_input.sequence.push_back("");
        for (i=0;i<chain_num;i++)
        {
            msta_iter_state.yname_vec.push_back(user_opts.yname);
            xlen = len_vec[i];
            if (i==msta_iter_state.repr_idx || xlen<3) continue;
            msta_iter_state.TM_pair_vec.push_back(std::make_pair(-TMave_mat[i][msta_iter_state.repr_idx],i));
        }
        sort(msta_iter_state.TM_pair_vec.begin(),msta_iter_state.TM_pair_vec.end());
    
        int tm_idx;
        if (user_opts.outfmt_opt<0) std::cout<<"#PDBchain1\tPDBchain2\tTM1\tTM2\t"
                               <<"RMSD\tID1\tID2\tIDali\tL1\tL2\tLali"<<std::endl;
        for (tm_idx=0; tm_idx<msta_iter_state.TM_pair_vec.size(); tm_idx++)
        {
            i=msta_iter_state.TM_pair_vec[tm_idx].second;
            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);

            double maxTM=TMave_mat[i][msta_iter_state.repr_idx];
            int maxj=msta_iter_state.repr_idx;
            for (j=0;j<chain_num;j++)
            {
                if (i==j || msta_iter_state.assign_list[j]<0 || TMave_mat[i][j]<=maxTM) continue;
                maxj=j;
                maxTM=TMave_mat[i][j];
            }
            j=maxj;
            msta_iter_state.assign_list[i]=j;
            ylen = len_vec[j];
            std::string seqy;
            secy.resize(ylen+1);
            ya.clear();
            ya.reserve(ylen);
            copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j], ylen,ya,seqy,secy);

            parsed_input.sequence[0]=seqxA_mat[i][j];
            parsed_input.sequence[1]=seqyA_mat[i][j];

            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            ChainPairAlignOptions align_opts;
            fill_chain_pair_align_options(align_opts, user_opts, ctrl_opts);
            align_opts.i_opt = 2;
            align_opts.u_opt = u_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.mol_type = cur_complex_mol_list;
            // entry function for structure alignment
            align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                xlen, ylen, align_opts, parsed_input.sequence,
                user_opts.outfmt_opt);

            OutputOptions out_opts;
            out_opts.fname_matrix = "";
            out_opts.fname_super = "";
            out_opts.Lnorm_ass = Lnorm_ass;
            out_opts.d0_scale = user_opts.d0_scale;
            out_opts.outfmt_opt = 2;
            out_opts.ter_opt = user_opts.ter_opt;
            out_opts.mm_opt = 0;
            out_opts.split_opt = user_opts.split_opt;
            out_opts.o_opt = 0;
            out_opts.i_opt = 0;
            out_opts.a_opt = user_opts.a_opt;
            out_opts.mirror_opt = 0;
            out_opts.u_opt = u_opt;
            out_opts.d_opt = user_opts.d_opt;
            const std::string structure1_display = msta_iter_state.xname_vec[i];
            const std::string structure2_display = msta_iter_state.xname_vec[j];
            if (user_opts.outfmt_opt<0) output_results(
                ChainPairDesc{structure1_display, structure2_display, "", "", xlen, ylen},
                chain_align_res, out_opts, resi_vec, resi_vec);
         
            xt.resize(xlen);
            do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
            for (r=0;r<xlen;r++)
            {
                a_vec[i][r][0]=xt[r][0];
                a_vec[i][r][1]=xt[r][1];
                a_vec[i][r][2]=xt[r][2];
            }
        
            parsed_input.sequence[0].clear();
            parsed_input.sequence[1].clear();
        }
}


void msta_build_msa(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        ylen = len_vec[msta_iter_state.repr_idx];
        std::string seqy;
        secy.resize(ylen+1);
        ya.clear();
        ya.reserve(ylen);
        copy_chain_data(a_vec[msta_iter_state.repr_idx],seq_vec[msta_iter_state.repr_idx],sec_vec[msta_iter_state.repr_idx], ylen,ya,seqy,secy);

        // recover alignment
        int    ylen_ext=ylen;        // chain length
        CoordArray ya_ext;               // structure of single chain
        std::string seqy_ext;            // for the protein sequence
        std::string secy_ext;            // for the secondary structure
        for (r=0;r<msta_iter_state.msa.size();r++) msta_iter_state.msa[r].clear(); msta_iter_state.msa.clear();
        msta_iter_state.msa.assign(ylen,""); // row is position along msta_iter_state.msa; column is sequence
        std::vector<std::string> msa_ext;      // row is position along msta_iter_state.msa; column is sequence
        for (r=0;r<ylen;r++) msta_iter_state.msa[r]=seqy[r];
        msta_iter_state.assign_list[msta_iter_state.repr_idx]=0;
        for (tm_idx=0; tm_idx<msta_iter_state.TM_pair_vec.size(); tm_idx++)
        {
            i=msta_iter_state.TM_pair_vec[tm_idx].second;
            msta_iter_state.assign_list[i]=tm_idx+1;

            xlen = len_vec[i];
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
        
            ChainPairAlignResult chain_align_res = { 0};
            chain_align_res.d0_out = 5.0;
            chain_align_res.invmap.assign(ylen + 1, 0);
            std::string& seqxA = chain_align_res.seqxA;
            std::string& seqyA = chain_align_res.seqyA;
            ChainPairAlignOptions align_opts;
            fill_chain_pair_align_options(align_opts, user_opts, msta_iter_params.ctrl_opts);
            align_opts.i_opt = 0;
            align_opts.u_opt = u_opt;
            align_opts.fast_opt = fast_opt;
            align_opts.se_opt = true;
            align_opts.cp_opt = false;
            align_opts.Lnorm = Lnorm_ass;
            align_opts.parallel_threads = 1;
            align_opts.mol_type = cur_complex_mol_list;
            se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen,
                parsed_input.sequence, align_opts, 1);

            int rx=0;
            int ry=0;
            ylen_ext=seqxA.size();
            ya_ext.resize(ylen_ext);                    // structure of single chain
            seqy_ext.resize(ylen_ext+1);            // for the protein sequence
            secy_ext.resize(ylen_ext+1);            // for the secondary structure
            std::string tmp_gap="";
            for (r=0;r<msta_iter_state.msa[0].size();r++) tmp_gap+='-';
            for (r=msa_ext.size();r<ylen_ext;r++) msa_ext.push_back("");
            for (r=0;r<ylen_ext;r++)
            {
                if (seqyA[r]=='-')
                {
                    msa_ext[r]=tmp_gap+seqxA[r];
                    ya_ext[r][0]=xa[rx][0];
                    ya_ext[r][1]=xa[rx][1];
                    ya_ext[r][2]=xa[rx][2];
                    seqy_ext[r]=seqx[rx];
                    secy_ext[r]=secx[rx];
                }
                else
                {
                    msa_ext[r]=msta_iter_state.msa[ry]+seqxA[r];
                    ya_ext[r][0]=ya[ry][0];
                    ya_ext[r][1]=ya[ry][1];
                    ya_ext[r][2]=ya[ry][2];
                    seqy_ext[r]=seqy[ry];
                    secy_ext[r]=secy[ry];
                }
                rx+=(seqxA[r]!='-');
                ry+=(seqyA[r]!='-');
            }

            // copy ya_ext to ya


            ylen=ylen_ext;
            ya.resize(ylen);
            seqy.assign(seqy_ext, 0, ylen);
            secy.resize(ylen+1);
            for (r=0;r<ylen;r++)
            {
                ya[r][0]=ya_ext[r][0];
                ya[r][1]=ya_ext[r][1];
                ya[r][2]=ya_ext[r][2];
                secy[r]=secy_ext[r];
            }
            for (r=0;r<ylen;r++)
            {
                if (r<msta_iter_state.msa.size()) msta_iter_state.msa[r]=msa_ext[r];
                else msta_iter_state.msa.push_back(msa_ext[r]);
            }

            // clean up
            tmp_gap.clear();
        }
        std::vector<std::string>().swap(msa_ext);
        std::vector<std::pair<double,int> >().swap(msta_iter_state.TM_pair_vec);
}

void msta_rebuild_pair_matrices(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& ylen = buffers.partner_len;
    int& r = buffers.residue_idx;
    int& tm_idx = buffers.tm_order_idx;

        for (i=0; i<chain_num; i++)
        {
            tm_idx=msta_iter_state.assign_list[i];
            if (tm_idx<0) continue;
            seqyA_mat[i][i]="";
            for (r=0 ;r<ylen ; r++) seqyA_mat[i][i]+=msta_iter_state.msa[r][tm_idx];
            seqxA_mat[i][i]=seqyA_mat[i][i];
        }
        for (i=0;i<chain_num; i++)
        {
            if (msta_iter_state.assign_list[i]<0) continue;
            std::string seqxA=seqxA_mat[i][i];
            for (j=0; j<chain_num; j++)
            {
                if (i==j || msta_iter_state.assign_list[j]<0) continue;
                std::string seqyA=seqyA_mat[j][j];
                seqxA_mat[i][j]=seqyA_mat[i][j]="";
                for (r=0;r<ylen;r++)
                {
                    if (seqxA[r]=='-' && seqyA[r]=='-') continue;
                    seqxA_mat[i][j]+=seqxA[r];
                    seqyA_mat[i][j]+=seqyA[r];
                }
                seqyA.clear();
            }
            seqxA.clear();
        }
}

bool msta_accumulate_stats(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    UserOptions& user_opts = msta_iter_params.user_opts;
    ParsedInput& parsed_input = msta_iter_params.parsed_input;
    const CharMatrix& seq_vec = msta_iter_params.complex.seqs;
    const CharMatrix& sec_vec = msta_iter_params.complex.secs;
    const std::vector<int>& len_vec = msta_iter_params.complex.lengths;
    DoubleCube& a_vec = msta_iter_params.a_vec;
    std::vector<std::vector<std::string> >& seqxA_mat = msta_iter_params.seqxA_mat;
    std::vector<std::vector<std::string> >& seqyA_mat = msta_iter_params.seqyA_mat;
    const int chain_num = msta_iter_params.chain_num;
    const double Lnorm_ass = msta_iter_params.Lnorm_ass;
    const bool u_opt = msta_iter_params.u_opt;
    const bool fast_opt = msta_iter_params.fast_opt;
    const int cur_complex_mol_list = msta_iter_params.mol_type_sum;
    MstaIterBuffers& buffers = msta_iter_params.buffers;
    int& i = buffers.member_chain_idx;
    int& j = buffers.partner_chain_idx;
    int& xlen = buffers.member_len;
    int& ylen = buffers.partner_len;
    CoordArray& xa = buffers.member_coords;
    CoordArray& ya = buffers.partner_coords;
    std::string& secx = buffers.member_sec;
    std::string& secy = buffers.partner_sec;

        // recover statistics such as TM-score
        msta_iter_state.compare_num=0;
        msta_iter_state.totals.TM1=0, msta_iter_state.totals.TM2=0;
        msta_iter_state.totals.TM3=0, msta_iter_state.totals.TM4=0, msta_iter_state.totals.TM5=0;
        msta_iter_state.totals.d0_0=0, msta_iter_state.totals.TM_0=0;
        msta_iter_state.totals.d0A=0, msta_iter_state.totals.d0B=0, msta_iter_state.totals.d0u=0, msta_iter_state.totals.d0a=0;
        msta_iter_state.totals.d0_out=0;
        msta_iter_state.totals.rmsd0 = 0.0;
        msta_iter_state.totals.L_ali=0;
        msta_iter_state.totals.Liden=0;
        msta_iter_state.totals.TM_ali=0, msta_iter_state.totals.rmsd_ali=0;
        msta_iter_state.totals.n_ali=0;
        msta_iter_state.totals.n_ali8=0;
        msta_iter_state.xlen_total=0, msta_iter_state.ylen_total=0;
        for (i=0; i< chain_num; i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i], xlen,xa,seqx,secx);
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                msta_iter_state.compare_num++;
                std::string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);
                parsed_input.sequence[0]=seqxA_mat[i][j];
                parsed_input.sequence[1]=seqyA_mat[i][j];
            
                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                chain_align_res.invmap.assign(ylen + 1, 0);
                ChainPairAlignOptions align_opts;
                fill_chain_pair_align_options(align_opts, user_opts, msta_iter_params.ctrl_opts);
                align_opts.i_opt = 1;
                align_opts.u_opt = u_opt;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = true;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_ass;
                align_opts.parallel_threads = 1;
                align_opts.mol_type = cur_complex_mol_list;
                se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen,
                    parsed_input.sequence, align_opts, 1);

                if (xlen<=ylen)
                {
                    msta_iter_state.xlen_total+=xlen;
                    msta_iter_state.ylen_total+=ylen;
                    msta_iter_state.totals.TM1+=chain_align_res.TM1;
                    msta_iter_state.totals.TM2+=chain_align_res.TM2;
                    msta_iter_state.totals.d0A+=chain_align_res.d0A;
                    msta_iter_state.totals.d0B+=chain_align_res.d0B;
                }
                else
                {
                    msta_iter_state.xlen_total+=ylen;
                    msta_iter_state.ylen_total+=xlen;
                    msta_iter_state.totals.TM1+=chain_align_res.TM2;
                    msta_iter_state.totals.TM2+=chain_align_res.TM1;
                    msta_iter_state.totals.d0A+=chain_align_res.d0B;
                    msta_iter_state.totals.d0B+=chain_align_res.d0A;
                }
                msta_iter_state.TM_mat[i][j]=chain_align_res.TM2;
                msta_iter_state.TM_mat[j][i]=chain_align_res.TM1;
                msta_iter_state.d0_mat[i][j]=chain_align_res.d0B;
                msta_iter_state.d0_mat[j][i]=chain_align_res.d0A;
                msta_iter_state.seqID_mat[i][j]=1.*chain_align_res.Liden/xlen;
                msta_iter_state.seqID_mat[j][i]=1.*chain_align_res.Liden/ylen;

                msta_iter_state.totals.TM3+=chain_align_res.TM3;
                msta_iter_state.totals.TM4+=chain_align_res.TM4;
                msta_iter_state.totals.TM5+=chain_align_res.TM5;
                msta_iter_state.totals.d0_0+=chain_align_res.d0_0;
                msta_iter_state.totals.TM_0+=chain_align_res.TM_0;
                msta_iter_state.totals.d0u+=chain_align_res.d0u;
                msta_iter_state.totals.d0_out+=chain_align_res.d0_out;
                msta_iter_state.totals.rmsd0+=chain_align_res.rmsd0;
                msta_iter_state.totals.L_ali+=chain_align_res.L_ali;        // Aligned length in standard_TMscore
                msta_iter_state.totals.Liden+=chain_align_res.Liden;
                msta_iter_state.totals.TM_ali+=chain_align_res.TM_ali;
                msta_iter_state.totals.rmsd_ali+=chain_align_res.rmsd_ali;  // TMscore and rmsd in standard_TMscore
                msta_iter_state.totals.n_ali+=chain_align_res.n_ali;
                msta_iter_state.totals.n_ali8+=chain_align_res.n_ali8;
            }
            
        }
        if (msta_iter_state.totals.TM4<=msta_iter_state.TM4_total_max) return false;
        msta_iter_state.TM4_total_max=msta_iter_state.totals.TM4;
        return true;
}


void msta_iterate(MstaIterationState& msta_iter_state, MstaIterParams& msta_iter_params)
{
    const DoubleMatrix& TMave_mat = msta_iter_params.TMave_mat;
    const int chain_num = msta_iter_params.chain_num;
    int iter=0;

    for (iter=0; iter<msta_iter_state.max_iter; iter++)
    {
        msta_iter_state.repr_idx=select_representative(TMave_mat, chain_num);

        msta_superpose_to_rep(msta_iter_state, msta_iter_params);
        msta_build_msa(msta_iter_state, msta_iter_params);
        msta_rebuild_pair_matrices(msta_iter_state, msta_iter_params);

        if (!msta_accumulate_stats(msta_iter_state, msta_iter_params)) break;
    }
}

void mTMalign_pairwise_alignment(AlignCommonInput& common_inputs, const ComplexData& complex,
    const DoubleCube& a_vec, const int chain_num, const double Lnorm_ass, const bool u_opt,
    const bool fast_opt, const int cur_complex_mol_list, DoubleMatrix& TMave_mat,
    std::vector<std::vector<std::string> >& seqxA_mat,
    std::vector<std::vector<std::string> >& seqyA_mat)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;
    const CharMatrix& seq_vec = complex.seqs;
    const CharMatrix& sec_vec = complex.secs;
    const std::vector<int>& len_vec = complex.lengths;
    const std::vector<std::string>& chainID_list = complex.chain_ids;
    const std::vector<std::string>& resi_vec = complex.resi;
    int i;
    int j;
    int xlen;
    int ylen;
    CoordArray xa;
    CoordArray ya;
    std::string secx;
    std::string secy;
    bool pair_parallel_done = false;
#ifdef _OPENMP
    if (ctrl_opts.parallel_threads > 1) {
        run_mTMalign_pairwise_parallel(
            a_vec, seq_vec, sec_vec, len_vec,
            parsed_input.chain1_list, chainID_list, seqxA_mat, seqyA_mat,
            TMave_mat, resi_vec, chain_num, Lnorm_ass, user_opts.d0_scale,
            u_opt, cur_complex_mol_list, user_opts.outfmt_opt, fast_opt, user_opts.TMcut,
            ctrl_opts.full_opt, ctrl_opts.se_opt, user_opts.ter_opt, user_opts.split_opt, user_opts.o_opt, user_opts.a_opt, user_opts.d_opt,
            ctrl_opts.parallel_threads);
        pair_parallel_done = true;
    }
#endif
    if (!pair_parallel_done)
    {
        // serial fallback
        for (i=0;i<chain_num;i++)
        {
            xlen=len_vec[i];
            if (xlen<3) continue;
            std::string seqx;
            secx.resize(xlen+1);
            xa.clear();
            xa.reserve(xlen);
            copy_chain_data(a_vec[i],seq_vec[i],sec_vec[i],xlen,xa,seqx,secx);
            seqxA_mat[i][i]=seqyA_mat[i][i]=seqx;
            for (j=i+1;j<chain_num;j++)
            {
                ylen=len_vec[j];
                if (ylen<3) continue;
                std::string seqy;
                secy.resize(ylen+1);
                ya.clear();
                ya.reserve(ylen);
                copy_chain_data(a_vec[j],seq_vec[j],sec_vec[j],ylen,ya,seqy,secy);

                ChainPairAlignResult chain_align_res = { 0};
                chain_align_res.d0_out = 5.0;
                ChainPairAlignOptions align_opts;
                align_opts.i_opt = 0;
                align_opts.a_opt = 0;
                align_opts.u_opt = u_opt;
                align_opts.d_opt = false;
                align_opts.fast_opt = fast_opt;
                align_opts.se_opt = ctrl_opts.se_opt;
                align_opts.cp_opt = false;
                align_opts.Lnorm = Lnorm_ass;
                align_opts.d0_scale = user_opts.d0_scale;
                align_opts.TMcut = user_opts.TMcut;
                align_opts.parallel_threads = 1;
                align_opts.ss_opt = 0;
                align_opts.mol_type = cur_complex_mol_list;
                // entry function for structure alignment
                align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                    xlen, ylen, align_opts, parsed_input.sequence,
                    user_opts.outfmt_opt);

                // store chain_align_res
                TMave_mat[i][j]=chain_align_res.TM4; TMave_mat[j][i]=chain_align_res.TM4;
                seqxA_mat[i][j]=seqyA_mat[j][i]=chain_align_res.seqxA;
                seqyA_mat[i][j]=seqxA_mat[j][i]=chain_align_res.seqyA;
                OutputOptions out_opts;
                out_opts.fname_matrix = "";
                out_opts.fname_super = "";
                out_opts.Lnorm_ass = Lnorm_ass;
                out_opts.d0_scale = user_opts.d0_scale;
                out_opts.outfmt_opt = user_opts.outfmt_opt;
                out_opts.ter_opt = user_opts.ter_opt;
                out_opts.mm_opt = true;
                out_opts.split_opt = user_opts.split_opt;
                out_opts.o_opt = user_opts.o_opt;
                out_opts.i_opt = 0;
                out_opts.a_opt = user_opts.a_opt;
                out_opts.mirror_opt = 0;
                out_opts.u_opt = false;
                out_opts.d_opt = user_opts.d_opt;
                const std::string structure1_display = parsed_input.chain1_list[i];
                const std::string structure2_display = parsed_input.chain1_list[j];
                if (ctrl_opts.full_opt) output_results(
                    ChainPairDesc{structure1_display, structure2_display, chainID_list[i], chainID_list[j], xlen, ylen},
                    chain_align_res, out_opts, resi_vec, resi_vec);
            }
        }
    }
}

int mTMalign(AlignCommonInput& common_inputs)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    double Lnorm_ass = user_opts.Lnorm_ass;
    bool u_opt = user_opts.u_opt;
    bool fast_opt = user_opts.fast_opt;

    // declare previously global variables
    int    i,j;                    // chain index
    int    xlen=0, ylen=0;         // chain length (serial path updates them in pair loop)
    CoordArray xa;                     // structure of single chain
    CoordArray ya;

    // parse chain list
    ComplexData complex;
    parse_chain_list(parsed_input.chain1_list, complex,
        user_opts.ter_opt, user_opts.split_opt, user_opts.mol_opt, user_opts.infmt1_opt,
        user_opts.atom_opt, parsed_input.normalize_atom_name, false, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    int chain_num=complex.coords.size();
    if (chain_num<=1) PrintErrorAndQuit("ERROR! <2 chains for multiple alignment");
    DoubleCube a_vec = complex.coords;
    DoubleCube ua_vec = complex.coords;
    const std::vector<int>& mol_vec = complex.mol_types;
    const std::vector<std::string>& chainID_list = complex.chain_ids;
    const std::vector<int>& len_vec = complex.lengths;
    const std::vector<std::string>& resi_vec = complex.resi;
    int    cur_complex_mol_list=0;
    int    total_len=0;
    xlen=0;
    for (i=0; i<chain_num; i++)
    {
        if (len_vec[i]>xlen) xlen=len_vec[i];
        total_len+=len_vec[i];
        cur_complex_mol_list+=mol_vec[i];
    }
    if (!u_opt) Lnorm_ass=total_len/chain_num;
    u_opt=true;
    total_len-=xlen;
    if (total_len>750) fast_opt=true;

    // get all-against-all alignment
    DoubleMatrix TMave_mat;
    TMave_mat.assign(chain_num,std::vector<double>(chain_num));
    std::vector<std::string> tmp_str_vec(chain_num,"");
    std::vector<std::vector<std::string> >seqxA_mat(chain_num,tmp_str_vec);
    std::vector<std::vector<std::string> >seqyA_mat(chain_num,tmp_str_vec);
    for (i=0;i<chain_num;i++) for (j=0;j<chain_num;j++) TMave_mat[i][j]=0;
    mTMalign_pairwise_alignment(common_inputs, complex, a_vec, chain_num,
        Lnorm_ass, u_opt, fast_opt, cur_complex_mol_list,
        TMave_mat, seqxA_mat, seqyA_mat);
    // representative related variables
    int r;
    MstaIterationState state;
    // Empty until rebuilt inside the iteration loop (msa.assign(ylen,"")),
    // so it does not depend on stale xlen/ylen left by the pairwise loop.
    for (i=0;i<chain_num;i++) state.xname_vec.push_back(
        parsed_input.chain1_list[i].substr(user_opts.dir_opt.size())+chainID_list[i]);

    // build and output UPGMA phylogenetic tree
    output_upgma_tree(state.xname_vec, TMave_mat, chain_num, user_opts.outdir_opt);

    state.assign_list.assign(chain_num, 0);
    state.compare_num = 0;
    state.TM4_total_max=0;

    state.max_iter=5-static_cast<int>(total_len/200);
    if (state.max_iter<2) state.max_iter=2;
    state.TM_vec.assign(chain_num,0);
    state.d0_vec.assign(chain_num,0);
    state.seqID_vec.assign(chain_num,0);
    state.TM_mat.assign(chain_num,state.TM_vec);
    state.d0_mat.assign(chain_num,state.d0_vec);
    state.seqID_mat.assign(chain_num,state.seqID_vec);
    MstaIterParams msta_iter_params(common_inputs, complex, a_vec, TMave_mat, seqxA_mat,
        seqyA_mat, chain_num, Lnorm_ass, u_opt, fast_opt, cur_complex_mol_list);
    msta_iterate(state, msta_iter_params);
    for (i=0;i<chain_num;i++)
    {
        for (j=0;j<chain_num;j++)
        {
            if (i==j) continue;
            state.TM_vec[i]+=state.TM_mat[i][j];
            state.d0_vec[i]+=state.d0_mat[i][j];
            state.seqID_vec[i]+=state.seqID_mat[i][j];
        }
        state.TM_vec[i]/=(chain_num-1);
        state.d0_vec[i]/=(chain_num-1);
        state.seqID_vec[i]/=(chain_num-1);
    }
    state.xlen_total    /=state.compare_num;
    state.ylen_total    /=state.compare_num;
    state.totals.TM1     /=state.compare_num;
    state.totals.TM2     /=state.compare_num;
    state.totals.d0A     /=state.compare_num;
    state.totals.d0B     /=state.compare_num;
    state.totals.TM3     /=state.compare_num;
    state.totals.TM4     /=state.compare_num;
    state.totals.TM5     /=state.compare_num;
    state.totals.d0_0    /=state.compare_num;
    state.totals.TM_0    /=state.compare_num;
    state.totals.d0u     /=state.compare_num;
    state.totals.d0_out  /=state.compare_num;
    state.totals.rmsd0   /=state.compare_num;
    state.totals.L_ali   /=state.compare_num;
    state.totals.Liden   /=state.compare_num;
    state.totals.TM_ali  /=state.compare_num;
    state.totals.rmsd_ali/=state.compare_num;
    state.totals.n_ali   /=state.compare_num;
    state.totals.n_ali8  /=state.compare_num;
    user_opts.xname="shorter";
    user_opts.yname="longer";
    std::string seqM="";
    std::string seqxA="";
    std::string seqyA="";
    std::stringstream buf;
    for (i=0; i<chain_num; i++)
    {
        if (state.assign_list[i]<0) continue;
        buf <<">"<<state.xname_vec[i]<<"\tL="<<len_vec[i]
            <<"\td0="<<std::setiosflags(std::ios::fixed)<<std::setprecision(2)<<state.d0_vec[i]
            <<"\tseqID="<<setiosflags(std::ios::fixed)<<std::setprecision(3)<<state.seqID_vec[i]
            <<"\tTM-score="<<setiosflags(std::ios::fixed)<<std::setprecision(5)<<state.TM_vec[i];
        if (i==state.repr_idx) buf<<"\t*";
        buf<<'\n'<<seqxA_mat[i][i]<<std::endl;
    }
    seqM=buf.str();
    seqM=seqM.substr(0,seqM.size()-1);
    buf.str(std::string());
    if (user_opts.outfmt_opt==0) print_version(USALIGN_VERSION);
    // calculate ccTM-score
    double ccTM_score = calc_ccTM_score(ua_vec, seqxA_mat, chain_num, len_vec, cur_complex_mol_list);
    state.totals.seqM=seqM;
    state.totals.seqxA=seqxA;
    state.totals.seqyA=seqyA;
    output_mTMalign_results( user_opts.xname,user_opts.yname, "","",
        state.xlen_total, state.ylen_total, state.totals,
        Lnorm_ass, user_opts.d0_scale,
        "", user_opts.outfmt_opt, user_opts.ter_opt, 0, user_opts.split_opt, false,
        "", false, user_opts.a_opt, u_opt, user_opts.d_opt, false,
        resi_vec, resi_vec, ccTM_score );

    if (user_opts.m_opt || user_opts.o_opt)
    {
        RotArray ut_mat; // rotation matrices for all-against-all alignment
        ut_mat.resize(chain_num);
        int ui;
        int uj;
        Vec3 t;
        RotMat u;
        double rmsd;
        for (i=0;i<chain_num;i++)
        {
            xlen=ylen=a_vec[i].size();
            xa.clear();
            xa.resize(xlen);
            ya.clear();
            ya.resize(ylen);
            for (r=0;r<xlen;r++)
            {
                xa[r][0]=ua_vec[i][r][0];
                xa[r][1]=ua_vec[i][r][1];
                xa[r][2]=ua_vec[i][r][2];
                ya[r][0]= a_vec[i][r][0];
                ya[r][1]= a_vec[i][r][1];
                ya[r][2]= a_vec[i][r][2];
            }
            Kabsch(xa, ya, xlen, 1, rmsd, t, u);
            for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) ut_mat[i][ui*3+uj]=u[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[i][9+uj]=t[uj];
        }
        DoubleCube().swap(ua_vec);

        if (user_opts.m_opt)
        {
            state.assign_list[state.repr_idx]=-1;
            output_dock_rotation_matrix(user_opts.fname_matrix,
                state.xname_vec,state.yname_vec, ut_mat, state.assign_list);
        }

        
        
        if (user_opts.o_opt) output_mTMalign_pymol(parsed_input.chain1_list,
            ut_mat, user_opts.fname_super, user_opts.o_opt);

    }

    // clean up
    return 1;
}

