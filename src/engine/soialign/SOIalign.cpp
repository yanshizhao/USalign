#include "SOIalign.h"

#ifdef _OPENMP
struct SOIEgsMaxCell { double score; int i; int j; };
#pragma omp declare reduction(soiegsmax : SOIEgsMaxCell : \
    omp_out = (omp_in.score > omp_out.score || \
        (omp_in.score == omp_out.score && (omp_in.i < omp_out.i || \
            (omp_in.i == omp_out.i && omp_in.j < omp_out.j)))) ? omp_in : omp_out) \
    initializer(omp_priv = SOIEgsMaxCell{0, -1, -1})
#endif

void soi_egs(DoubleMatrix& score, const int xlen, const int ylen, std::vector<int>& invmap,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads){

    std::vector<int> fwdmap(xlen, -1); // j=fwdmap[i];
    for (int i=0; i<xlen; i++) fwdmap[i]=-1;
    for (int j=0; j<ylen; j++)
    {
        int i=invmap[j];
        if (i>=0) fwdmap[i]=j;
    }

    // stage 1 - make initial assignment, starting from the highest score pair
    bool omp_active = false;
#ifdef _OPENMP
    omp_active = (parallel_threads>1 && omp_get_active_level()==0);
#endif
    while(1)
    {
        double max_score=0;
        int maxi=-1;
        int maxj=-1;
        if (omp_active)
        {
#ifdef _OPENMP
            SOIEgsMaxCell best = {0, -1, -1};
            #pragma omp parallel for schedule(static) num_threads(parallel_threads) reduction(soiegsmax:best)
            for (int i=0;i<xlen;i++)
            {
                if (fwdmap[i]>=0) continue;
                for (int j=0;j<ylen;j++)
                {
                    if (invmap[j]>=0 || score[i+1][j+1]<=best.score) continue;
                    if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond, fwdmap,invmap)) continue;
                    best.score=score[i+1][j+1];
                    best.i=i;
                    best.j=j;
                }
            }
            max_score=best.score;
            maxi=best.i;
            maxj=best.j;
#endif
        }
        else
        {
            for (int i=0;i<xlen;i++)
            {
                if (fwdmap[i]>=0) continue;
                for (int j=0;j<ylen;j++)
                {
                    if (invmap[j]>=0 || score[i+1][j+1]<=max_score) continue;
                    if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond, fwdmap,invmap)) continue;
                    maxi=i;
                    maxj=j;
                    max_score=score[i+1][j+1];
                }
            }
        }
        if (maxi<0) break; // no assignment;
        invmap[maxj]=maxi;
        fwdmap[maxi]=maxj;
    }

    double total_score=0;
    for (int j=0;j<ylen;j++)
    {
        int i=invmap[j];
        if (i>=0) total_score+=score[i+1][j+1];
    }

    // stage 2 - swap assignment until total score cannot be improved
    for (int iter=0; iter<getmin(xlen,ylen)*5; iter++)
    {
        double delta_score=-1;
        for (int i=0;i<xlen;i++)
        {
            int oldj=fwdmap[i];
            for (int j=0;j<ylen;j++)
            {
                int oldi=invmap[j];
                if (score[i+1][j+1]<=0 || oldi==i) continue;
                if (mm_opt==6 && (!sec2sq(i,j,secx_bond,secy_bond,fwdmap,invmap) ||
                            !sec2sq(oldi,oldj,secx_bond,secy_bond,fwdmap,invmap)))
                    continue;
                delta_score=score[i+1][j+1];
                if (oldi>=0 && oldj>=0) delta_score+=score[oldi+1][oldj+1];
                if (oldi>=0) delta_score-=score[oldi+1][j+1];
                if (oldj>=0) delta_score-=score[i+1][oldj+1];

                if (delta_score>0) // successful swap
                {
                    fwdmap[i]=j;
                    if (oldi>=0) fwdmap[oldi]=oldj;
                    invmap[j]=i;
                    if (oldj>=0) invmap[oldj]=oldi;
                    total_score+=delta_score;
                    break;
                }
            }
        }
        if (delta_score<=0) break; // cannot make further swap
    }
}

void fill_soi_params(SoiAlignParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ControlOptions& ctrl_opts = common_inputs.control_options;

    params.read_resi  = (user_opts.o_opt != 0) ? 2 : 0;
    params.closeK_opt = ctrl_opts.closeK_opt;
    params.mm_opt     = ctrl_opts.mm_opt;
    fill_chain_pair_align_options(params.align, user_opts, ctrl_opts);
}


int parse_chain_for_soi(const std::vector<std::string>& pdb_lines, const std::string& chain_name,
    int& mol, const UserOptions& user_opts, const SoiAlignParams& soi_params,
    const bool mirror_opt, CoordArray& coords, std::string& seq, std::string& sec,
    std::vector<std::string>& resi_vec, CoordArray& k_nearest, IntPairArray& sec_bond)
{
    int len=pdb_lines.size();
    if (user_opts.mol_opt=="RNA") mol=1;
    else if (user_opts.mol_opt=="protein") mol=-1;
    if (!len)
    {
        std::cerr<<"Warning! Cannot parse file: "<<chain_name
            <<". Chain length 0."<<std::endl;
        return len;
    }
    else if (len<3)
    {
        std::cerr<<"Sequence is too short <3!: "<<chain_name<<std::endl;
        return len;
    }
    coords.clear();
    coords.reserve(len);
    if (soi_params.closeK_opt>=3) k_nearest.resize(len*soi_params.closeK_opt);
    sec.resize(len + 1);
    len = read_PDB(pdb_lines, coords, seq, resi_vec, soi_params.read_resi);
    if (mirror_opt) for (int r=0;r<len;r++) coords[r][2]=-coords[r][2];
    if (mol>0) make_sec(seq, coords, len, sec, user_opts.atom_opt);
    else make_sec(coords, len, sec); // secondary structure assignment
    if (soi_params.closeK_opt>=3) getCloseK(coords, len, soi_params.closeK_opt, k_nearest);
    if (soi_params.mm_opt==6)
    {
        sec_bond.resize(len);
        assign_sec_bond(sec_bond, sec, len);
    }
    return len;
}

struct SoiAlignFileData
{
    std::vector<std::vector<std::string> > pdb_lines;
    std::vector<std::string> chain_ids;
    std::vector<int> mol_types;
    int chain_num;
};

struct SoiAlignPairTask
{
    int file1_idx;
    int chain1_idx;
    int file2_idx;
    int chain2_idx;
    int order;
};

const SoiAlignFileData* load_soi_file_cache(std::map<std::string, SoiAlignFileData>& cache,
    const std::string& name, int infmt_opt, const PdbParseOptions& pdb_parse_opts)
{
    std::map<std::string, SoiAlignFileData>::iterator cached = cache.find(name);
    if (cached != cache.end()) return &cached->second;
    SoiAlignFileData& file_data = cache[name];
    file_data.chain_num = (int)get_PDB_lines(name, file_data.pdb_lines, file_data.chain_ids,
        file_data.mol_types, infmt_opt, pdb_parse_opts);
    return &file_data;
}

void run_soi_task(const SoiAlignPairTask& task,
    const AlignCommonInput& common_inputs, const SoiAlignParams& soi_params,
    const SoiAlignFileData& file1, const SoiAlignFileData& file2,
    const std::string& xname, const std::string& yname,
    std::ostream& out_stream)
{
    UserOptions task_user_opts = common_inputs.user_options;
    task_user_opts.xname = xname;
    task_user_opts.yname = yname;

    SoiAlignContext soi_context;
    std::string secx, secy, seqx, seqy;
    CoordArray xa, ya;
    std::vector<std::string> resi_vec1, resi_vec2;

    int mol1 = file1.mol_types[task.chain1_idx];
    int mol2 = file2.mol_types[task.chain2_idx];
    int xlen = parse_chain_for_soi(file1.pdb_lines[task.chain1_idx], xname,
        mol1, task_user_opts, soi_params, task_user_opts.mirror_opt,
        xa, seqx, secx, resi_vec1, soi_context.k_nearest1, soi_context.chain1_bonds);
    if (xlen<3) return;

    int ylen = parse_chain_for_soi(file2.pdb_lines[task.chain2_idx], yname,
        mol2, task_user_opts, soi_params, false, ya, seqy, secy, resi_vec2,
        soi_context.k_nearest2, soi_context.chain2_bonds);
    if (ylen<3) return;

    ChainPairAlignResult chain_align_res = { 0};
    chain_align_res.d0_out = 5.0;
    ChainPairAlignOptions soi_align_opts = soi_params.align;
    soi_align_opts.fast_opt = (getmin(xlen,ylen)>1500)?true:task_user_opts.fast_opt;
    soi_align_opts.mol_type = mol1+mol2;
#ifdef _OPENMP
    soi_align_opts.parallel_threads = 1;
#endif
    soi_context.pair_distances.assign(ylen+1, 0.0);

    soi_align_pair(chain_align_res, xa, ya, soi_params.closeK_opt,
        seqx, seqy, secx, secy, xlen, ylen, soi_align_opts,
        common_inputs.parsed_input.sequence, soi_context,
        soi_params.mm_opt, task_user_opts.outfmt_opt, out_stream);

    AlignCommonInput task_input = common_inputs;
    task_input.user_options = task_user_opts;
    task_input.control_options.parallel_threads = 1;
    output_pair_alignment(task_input, xlen, ylen, chain_align_res,
        file1.pdb_lines[task.chain1_idx], file2.pdb_lines[task.chain2_idx],
        file1.chain_ids[task.chain1_idx], file2.chain_ids[task.chain2_idx],
        resi_vec1, resi_vec2, &soi_context.pair_distances, out_stream);
}

void build_soi_pair_tasks(const std::vector<std::string>& chain1_list,
    const std::vector<std::string>& chain2_list,
    const std::vector<const SoiAlignFileData*>& files1,
    const std::vector<const SoiAlignFileData*>& files2,
    bool upper_triangular, bool dirpair, std::vector<SoiAlignPairTask>& tasks)
{
    for (int file1_idx = 0; file1_idx < (int)chain1_list.size(); file1_idx++)
    {
        const SoiAlignFileData& file1 = *files1[file1_idx];
        if (!file1.chain_num)
        {
            std::cerr<<"Warning! Cannot parse file: "<<chain1_list[file1_idx]
                <<". Chain number 0."<<std::endl;
            continue;
        }

        for (int chain1_idx = 0; chain1_idx < file1.chain_num; chain1_idx++)
        {
            int file2_start = upper_triangular ? (file1_idx + 1) : 0;
            for (int file2_idx = file2_start; file2_idx < (int)chain2_list.size(); file2_idx++)
            {
                if (dirpair && file1_idx!=file2_idx) continue;
                const SoiAlignFileData& file2 = *files2[file2_idx];
                if (!file2.chain_num)
                {
                    std::cerr<<"Warning! Cannot parse file: "<<chain2_list[file2_idx]
                        <<". Chain number 0."<<std::endl;
                    continue;
                }
                for (int chain2_idx = 0; chain2_idx < file2.chain_num; chain2_idx++)
                    tasks.push_back({file1_idx, chain1_idx, file2_idx, chain2_idx, (int)tasks.size()});
            }
        }
    }
}

// sequence order independent alignment
int SOIalign(AlignCommonInput& common_inputs, const SoiAlignParams& soi_params)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ParsedInput& parsed_input = common_inputs.parsed_input;

    const PdbParseOptions pdb_parse_opts1 = make_pdb_parse_opts(
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    const PdbParseOptions pdb_parse_opts2 = make_pdb_parse_opts(
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2);

    std::map<std::string, SoiAlignFileData> file1_cache;
    std::map<std::string, SoiAlignFileData> file2_cache;
    std::vector<const SoiAlignFileData*> files1(parsed_input.chain1_list.size());
    std::vector<const SoiAlignFileData*> files2(parsed_input.chain2_list.size());

    for (int i=0;i<(int)parsed_input.chain1_list.size();i++)
        files1[i] = load_soi_file_cache(file1_cache, parsed_input.chain1_list[i], user_opts.infmt1_opt, pdb_parse_opts1);

    for (int j=0;j<(int)parsed_input.chain2_list.size();j++)
        files2[j] = load_soi_file_cache(file2_cache, parsed_input.chain2_list[j], user_opts.infmt2_opt, pdb_parse_opts2);

    std::vector<SoiAlignPairTask> tasks;
    build_soi_pair_tasks(parsed_input.chain1_list, parsed_input.chain2_list, files1, files2, user_opts.dir_opt.size() > 0, user_opts.dirpair_opt.size() != 0, tasks);

    int task_num = (int)tasks.size();
    std::vector<std::string> out_lines(task_num);

#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads) if(common_inputs.control_options.parallel_threads>1 && task_num>1)
    for (int task_idx = 0; task_idx < task_num; task_idx++)
    {
        const SoiAlignPairTask& task = tasks[task_idx];
        std::stringstream out_stream;
        run_soi_task(task, common_inputs, soi_params, *files1[task.file1_idx], *files2[task.file2_idx],
            parsed_input.chain1_list[task.file1_idx], parsed_input.chain2_list[task.file2_idx], out_stream);
        out_lines[task.order] = out_stream.str();
    }

    for (int task_idx = 0; task_idx < task_num; task_idx++)
        std::cout << out_lines[task_idx];
    return 0;
}

int SOIalign_main(CoordArray& xa_c, CoordArray& ya_c,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    const std::string &seqx, const std::string &seqy, const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    std::ostream& os)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    std::vector<int>& invmap = chain_align_res.invmap;
    const bool fast_opt = align_opts.fast_opt;
    int i;
    int j;
    SoiAlignSearchState soi_search_state;

    /***********************/
    // allocate memory
    /***********************/
    int minlen = std::min(xlen, ylen);
    int maxlen = (xlen>ylen)?xlen:ylen;
    soi_search_state.score.assign( xlen+1, std::vector<double>(ylen+1));
    soi_search_state.scoret.assign(ylen+1, std::vector<double>(xlen+1));
    soi_search_state.path.assign( maxlen+1, std::vector<char>(maxlen+1));
    soi_search_state.val.assign(  maxlen+1, std::vector<double>(maxlen+1));
    soi_search_state.xtm.resize(minlen);
    soi_search_state.ytm.resize(minlen);
    soi_search_state.xt.resize(xlen);
    soi_search_state.yt.resize(ylen);
    soi_search_state.r1.resize(minlen);
    soi_search_state.r2.resize(minlen);

    /***********************/
    //    parameter set
    /***********************/
    parameter_set4search(xlen, ylen, soi_search_state.D0_MIN, soi_search_state.Lnorm,
        soi_search_state.score_d8, soi_search_state.d0, soi_search_state.d0_search, soi_search_state.dcu0);
    soi_search_state.simplify_step    = 40; //for simplified search engine
    soi_search_state.score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    soi_search_state.fwdmap0.assign(xlen+1, -1);
    soi_search_state.invmap0.assign(ylen+1, -1);
    soi_search_state.TMmax=-1;
    soi_search_state.TM=-1;
    soi_search_state.local_d0_search = soi_search_state.d0_search;
    soi_search_state.iteration_max=(fast_opt)?2:30;

    soi_align_initial_seq(xa_c, ya_c, seqx, seqy, secx, secy,
        xlen, ylen, sequence, align_opts, soi_search_state, chain_align_res, secx_bond, secy_bond, mm_opt,
        align_opts.parallel_threads, os);
    soi_align_initial_closek(xa_c, ya_c, xk, yk, closeK_opt,
        xlen, ylen, soi_search_state, chain_align_res, secx_bond, secy_bond, mm_opt,
        align_opts.parallel_threads);

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for (i=0; i<xlen; i++) soi_search_state.fwdmap0[i]=-1;
    for (j=0; j<ylen; j++)
    {
        i=soi_search_state.invmap0[j];
        invmap[j]=i;
        if (i>=0)
        {
            soi_search_state.fwdmap0[i]=j;
            flag=true;
        }
    }
    if(!flag)
    {
        os << "There is no alignment between the two structures! "
             << "Program stop with no result!" << std::endl;
        TM1=TM2=TM3=TM4=TM5=0;
        init_rt_mat(t0, u0);
        return 1;
    }

    soi_align_final_tmscore(xa_c, ya_c, xlen, ylen, align_opts, soi_search_state, chain_align_res);
    soi_derive_alignment_strings(xa_c, ya_c, seqx, seqy, xlen, ylen, soi_search_state, chain_align_res, dist_list);

    return 0;
}

double SOI_iter(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads, const bool init_invmap)
{
    double rmsd;
    std::vector<int> invmap(ylen+1);

    int iteration;
    int i;
    int j;
    int k;
    double tmscore;
    double tmscore_max;
    double tmscore_old=0;
    tmscore_max=-1;

    for (iteration=0; iteration<iteration_max; iteration++)
    {
        if (iteration==0 && init_invmap)
            for (j=0;j<ylen;j++) invmap[j]=invmap0[j];
        else
        {
            for (j=0; j<ylen; j++) invmap[j]=-1;
            if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
        }
        soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt, parallel_threads);

        k=0;
        for (j=0; j<ylen; j++)
        {
            i=invmap[j];
            if (i<0) continue;

            xtm[k][0]=xa[i][0];
            xtm[k][1]=xa[i][1];
            xtm[k][2]=xa[i][2];

            ytm[k][0]=ya[j][0];
            ytm[k][1]=ya[j][1];
            ytm[k][2]=ya[j][2];
            k++;
        }

        tmscore = TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
            40, 8, rmsd, local_d0_search, Lnorm, score_d8, d0);

        if (tmscore>tmscore_max)
        {
            tmscore_max=tmscore;
            for (j=0; j<ylen; j++) invmap0[j]=invmap[j];
        }

        if (iteration>0 && fabs(tmscore_old-tmscore)<0.000001) break;
        tmscore_old=tmscore;
        do_rotation(xa, xt, xlen, t, u);
        SOI_super2score(xt, ya, xlen, ylen, score, d0, score_d8);
    }

    return tmscore_max;
}
void get_SOI_initial_assign(CoordArray& xk, CoordArray& yk, const int closeK_opt,
    DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, const int xlen, const int ylen,
    Vec3& t, RotMat& u, std::vector<int>& invmap,
    double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads)
{
    int i;
    int j;
    int k;
    CoordArray xfrag, xtran, yfrag;
    xfrag.resize(closeK_opt);
    xtran.resize(closeK_opt);
    yfrag.resize(closeK_opt);
    double rmsd;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;

    for (i=0;i<xlen;i++)
    {
        for (k=0;k<closeK_opt;k++)
        {
            xfrag[k][0]=xk[i*closeK_opt+k][0];
            xfrag[k][1]=xk[i*closeK_opt+k][1];
            xfrag[k][2]=xk[i*closeK_opt+k][2];
        }

        for (j=0;j<ylen;j++)
        {
            for (k=0;k<closeK_opt;k++)
            {
                yfrag[k][0]=yk[j*closeK_opt+k][0];
                yfrag[k][1]=yk[j*closeK_opt+k][1];
                yfrag[k][2]=yk[j*closeK_opt+k][2];
            }

            Kabsch(xfrag, yfrag, closeK_opt, 1, rmsd, t, u);

            do_rotation(xfrag, xtran, closeK_opt, t, u);

            k=closeK_opt-1;
            d2=dist(xtran[k], yfrag[k]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+d2/d02);
        }
    }

    for (j=0;j<ylen;j++) invmap[j]=-1;
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);
    for (j=0; j<ylen;j++) i=invmap[j];
    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt, parallel_threads);

}

void soi_align_initial_seq(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads,
    std::ostream& os)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &local_d0_search = soi_search_state.local_d0_search;
    DoubleMatrix &score = soi_search_state.score;
    DoubleMatrix &scoret = soi_search_state.scoret;
    CharMatrix &path = soi_search_state.path;
    DoubleMatrix &val = soi_search_state.val;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &yt = soi_search_state.yt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double &TMmax = soi_search_state.TMmax;
    double &TM = soi_search_state.TM;
    int &iteration_max = soi_search_state.iteration_max;
    int i;
    int j;
    CPalign_main(xa_c, ya_c, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts, os);
    chain_align_res.do_vec.clear();
    if (mm_opt==6)
    {
        i=0;
        j=0;
        for (int r=0;r<seqxA.size();r++)
        {
            if (seqxA[r]=='*') // circular permutation point
            {
                for (int jj=0;jj<j;jj++) if (invmap0[jj]>=0)
                    invmap0[jj]+=xlen - i;
                i=0;
                continue;
            }
            if (seqyA[r]!='-')
            {
                if (seqxA[r]!='-') invmap0[j]=i;
                j++;
            }
            if (seqxA[r]!='-') i++;
        }
        for (j=0;j<ylen;j++)
        {
            i=invmap0[j];
            if (i>=0) fwdmap0[i]=j;
        }
    }
    do_rotation(xa_c, xt, xlen, t0, u0);
    SOI_super2score(xt, ya_c, xlen, ylen, score, d0, score_d8);
    for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];
    TMmax=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa_c, ya_c,
        xlen, ylen, t0, u0, invmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt, parallel_threads, true);
    TM   =SOI_iter(r2, r1, ytm, xtm, yt, scoret, path, val, ya_c, xa_c,
        ylen, xlen, t0, u0, fwdmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secy_bond, secx_bond, mm_opt, parallel_threads, true);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (j=0; j<ylen; j++) invmap0[j]=-1;
        for (i=0; i<xlen; i++) 
        {
            j=fwdmap0[i];
            if (j>=0) invmap0[j]=i;
        }
    }
}
void soi_align_initial_closek(CoordArray& xa_c, CoordArray& ya_c,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    int xlen, int ylen, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    const int parallel_threads)
{
    std::vector<int> &invmap = chain_align_res.invmap;
    Vec3 &t = soi_search_state.t;
    RotMat &u = soi_search_state.u;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &local_d0_search = soi_search_state.local_d0_search;
    DoubleMatrix &score = soi_search_state.score;
    DoubleMatrix &scoret = soi_search_state.scoret;
    CharMatrix &path = soi_search_state.path;
    DoubleMatrix &val = soi_search_state.val;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &yt = soi_search_state.yt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double &TMmax = soi_search_state.TMmax;
    double &TM = soi_search_state.TM;
    int &iteration_max = soi_search_state.iteration_max;
    int i;
    int j;
    /***************************************************************/
    // initial alignment with sequence order independent alignment
    /***************************************************************/
    if (closeK_opt>=3)
    {
        get_SOI_initial_assign(xk, yk, closeK_opt, score, path, val,
            xlen, ylen, t, u, invmap, d0, score_d8,
            secx_bond, secy_bond, mm_opt, parallel_threads);
        for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];

        SOI_assign2super(r1, r2, xtm, ytm, xt, xa_c, ya_c,
            xlen, ylen, t, u, invmap, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa_c, ya_c,
            xlen, ylen, t, u, invmap, iteration_max,
            local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt, parallel_threads);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j = 0; j<ylen; j++) invmap0[j] = invmap[j];
        }

        for (i=0;i<xlen;i++) fwdmap0[i]=-1;
        if (mm_opt==6) NWDP_TM(scoret, path, val, ylen, xlen, -0.6, fwdmap0);
        soi_egs(scoret, ylen, xlen, fwdmap0, secy_bond, secx_bond, mm_opt, parallel_threads);
        SOI_assign2super(r2, r1, ytm, xtm, yt, ya_c, xa_c,
            ylen, xlen, t, u, fwdmap0, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r2, r1, ytm, xtm, yt, scoret, path, val, ya_c, xa_c, ylen, xlen, t, u,
            fwdmap0, iteration_max, local_d0_search, Lnorm, d0, score_d8,secy_bond, secx_bond, mm_opt, parallel_threads);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j=0; j<ylen; j++) invmap0[j]=-1;
            for (i=0; i<xlen; i++) 
            {
                j=fwdmap0[i];
                if (j>=0) invmap0[j]=i;
            }
        }
    }
}
void soi_align_final_tmscore(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM2 = chain_align_res.TM2;
    double &d0B = chain_align_res.d0B;
    double &rmsd0 = chain_align_res.rmsd0;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    const int mol_type = align_opts.mol_type;
    const bool fast_opt = align_opts.fast_opt;
    double &D0_MIN = soi_search_state.D0_MIN;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &d0_search = soi_search_state.d0_search;
    double &local_d0_search = soi_search_state.local_d0_search;
    int &simplify_step = soi_search_state.simplify_step;
    int &score_sum_method = soi_search_state.score_sum_method;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    Vec3 &t = soi_search_state.t;
    RotMat &u = soi_search_state.u;
    double &TM = soi_search_state.TM;
    int i;
    int j;

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alignment
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
        invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
        false, Lnorm, score_d8, d0);
    
    double rmsd;
    simplify_step=1;
    score_sum_method=0;

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    double d;
    std::vector<int> m1(xlen); //alignd index in x
    std::vector<int> m2(ylen); //alignd index in y
    copy_t_u(t, u, t0, u0);
    
    //****************************************//
    //              Final TMscore 1           //
    //****************************************//

    do_rotation(xa_c, xt, xlen, t, u);
    k=0;
    n_ali=0;
    for (i=0; i<xlen; i++)
    {
        j=fwdmap0[i];
        if(j>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xt[i], ya_c[j]));
            if (d <= score_d8)
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa_c[i][0];
                xtm[k][1]=xa_c[i][1];
                xtm[k][2]=xa_c[i][2];

                ytm[k][0]=ya_c[j][0];
                ytm[k][1]=ya_c[j][1];
                ytm[k][2]=ya_c[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya_c[j][0];
                r2[k][1] = ya_c[j][1];
                r2[k][2] = ya_c[j][2];

                k++;
            }
            else fwdmap0[i]=-1;
        }
    }
    n_ali8=k;

    Kabsch(r1, r2, n_ali8, 0, rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);

    //normalized by length of structure A
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    soi_align_final_tmscore_inv(xa_c, ya_c, xlen, ylen, align_opts, soi_search_state, chain_align_res);
}
void soi_align_final_tmscore_inv(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM1 = chain_align_res.TM1;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    double &d0_0 = chain_align_res.d0_0;
    double &TM_0 = chain_align_res.TM_0;
    double &d0A = chain_align_res.d0A;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    double &d0_out = chain_align_res.d0_out;
    int &n_ali8 = chain_align_res.n_ali8;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;
    const int a_opt = align_opts.a_opt;
    const bool u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    const int mol_type = align_opts.mol_type;
    double &D0_MIN = soi_search_state.D0_MIN;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &d0_search = soi_search_state.d0_search;
    double &local_d0_search = soi_search_state.local_d0_search;
    int &simplify_step = soi_search_state.simplify_step;
    int &score_sum_method = soi_search_state.score_sum_method;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double rmsd;
    double Lnorm_0=ylen;
    int i;
    int j;
    int k;
    double d;
    std::vector<int> m1(xlen);
    std::vector<int> m2(ylen);
    //****************************************//
    //              Final TMscore 2           //
    //****************************************//

    do_rotation(xa_c, xt, xlen, t0, u0);
    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            d=sqrt(dist(xt[i], ya_c[j]));
            if (d <= score_d8)
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa_c[i][0];
                xtm[k][1]=xa_c[i][1];
                xtm[k][2]=xa_c[i][2];

                ytm[k][0]=ya_c[j][0];
                ytm[k][1]=ya_c[j][1];
                ytm[k][2]=ya_c[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya_c[j][0];
                r2[k][1] = ya_c[j][1];
                r2[k][2] = ya_c[j][2];

                k++;
            }
            else invmap[j]=invmap0[j]=-1;
        }
    }

    //normalized by length of structure B
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    TM_0 = TM1;

    if (a_opt>0)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0, d0_search, mol_type);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm, d0, d0_search);
        d0_out=d0_scale;
        d0_0=d0_scale;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM5;
    }
}
void soi_derive_alignment_strings(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy, int xlen, int ylen,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Liden = chain_align_res.Liden;
    const double d0_out = chain_align_res.d0_out;
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    CoordArray &xt = soi_search_state.xt;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    int i;
    int j;
    int k;
    double d;
    // derive alignment from superposition
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap0[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    do_rotation(xa_c, xt, xlen, t0, u0);

    Liden=0;
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap0[j];
        dist_list[j]=-1;
        if (i<0) continue;
        d=sqrt(dist(xt[i], ya_c[j]));
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        dist_list[j]=d;
        seqxA[j]=seqx[i];
        Liden+=(seqx[i]==seqy[j]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap0[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }

}

void soi_build_score_seed_invmap(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    ChainPairAlignResult& chain_align_res, IntPairArray& secx_bond, IntPairArray& secy_bond,
    const int mm_opt, DoubleMatrix& score, const int parallel_threads)
{
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const int mol_type = align_opts.mol_type;
    double D0_MIN;
    double Lnorm;
    double score_d8,d0,d0_search,dcu0;
    CharMatrix  path;
    DoubleMatrix val;
    int i;
    int j;

    /***********************/
    // allocate memory
    /***********************/
    score.assign(xlen+1, std::vector<double>(ylen+1));
    path.assign( xlen+1, std::vector<char>(ylen+1));
    val.assign(  xlen+1, std::vector<double>(ylen+1));

    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
        score_d8, d0, d0_search, dcu0); // set score_d8
    parameter_set4final(xlen, D0_MIN, Lnorm,
        d0B, d0_search, mol_type); // set d0B
    parameter_set4final(ylen, D0_MIN, Lnorm,
        d0A, d0_search, mol_type); // set d0A
    if (a_opt)
        parameter_set4final((xlen+ylen)*0.5, D0_MIN, Lnorm,
            d0a, d0_search, mol_type); // set d0a
    if (u_opt)
    {
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0u, d0_search, mol_type); // set d0u
        if (u_opt==2)
        {
            parameter_set4search(Lnorm_ass, Lnorm_ass, D0_MIN, Lnorm,
                score_d8, d0, d0_search, dcu0); // set score_d8
        }
    }

    for(j=0; j<ylen; j++) invmap[j]=-1;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for(i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xa[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);

    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt, parallel_threads);
}
bool soi_se_score_alignment(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    const DoubleMatrix& score, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list)
{
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    double &rmsd0 = chain_align_res.rmsd0;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    std::vector<int> m1;
    std::vector<int> m2;
    int i;
    int j;
    double d;
    if (outfmt_opt<2)
    {
        m1.resize(xlen);
        m2.resize(ylen);
    }

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        dist_list[j]=-1;
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xa[i], ya[j]));
            dist_list[j]=d;
            if (score[i+1][j+1]>0)
            {
                if (outfmt_opt<2)
                {
                    m1[k]=i;
                    m2[k]=j;
                }
                k++;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);
    if (outfmt_opt>=2)
    {
        return true;
    }
    return false;
}
void soi_se_extract_strings(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    ChainPairAlignResult& chain_align_res)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Liden = chain_align_res.Liden;
    const double d0_out = chain_align_res.d0_out;
    std::vector<int> &invmap = chain_align_res.invmap;
    int i;
    int j;
    int k = chain_align_res.n_ali8;
    double d;

    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    std::vector<int> fwdmap(xlen+1, -1);
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap[j];
        if (i<0) continue;
        d=sqrt(dist(xa[i], ya[j]));
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        fwdmap[i]=j;
        seqxA[j]=seqx[i];
        Liden+=(seqxA[k]==seqyA[k]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }
}
