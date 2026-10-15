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
