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
    int i;
    int j;
    std::vector<int> fwdmap(xlen, -1); // j=fwdmap[i];
    for (i=0; i<xlen; i++) fwdmap[i]=-1;
    for (j=0; j<ylen; j++)
    {
        i=invmap[j];
        if (i>=0) fwdmap[i]=j;
    }

    // stage 1 - make initial assignment, starting from the highest score pair
    double max_score;
    int maxi;
    int maxj;
    bool omp_active = false;
#ifdef _OPENMP
    omp_active = (parallel_threads>1 && omp_get_active_level()==0);
#endif
    while(1)
    {
        max_score=0;
        maxi=maxj=-1;
        if (omp_active)
        {
#ifdef _OPENMP
            SOIEgsMaxCell best = {0, -1, -1};
            #pragma omp parallel for schedule(static) num_threads(parallel_threads) \
                reduction(soiegsmax:best) private(j)
            for (i=0;i<xlen;i++)
            {
                if (fwdmap[i]>=0) continue;
                for (j=0;j<ylen;j++)
                {
                    if (invmap[j]>=0 || score[i+1][j+1]<=best.score) continue;
                    if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond,
                        fwdmap,invmap)) continue;
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
            for (i=0;i<xlen;i++)
            {
                if (fwdmap[i]>=0) continue;
                for (j=0;j<ylen;j++)
                {
                    if (invmap[j]>=0 || score[i+1][j+1]<=max_score) continue;
                    if (mm_opt==6 && !sec2sq(i,j,secx_bond,secy_bond,
                        fwdmap,invmap)) continue;
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
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i>=0) total_score+=score[i+1][j+1];
    }

    // stage 2 - swap assignment until total score cannot be improved
    int iter;
    int oldi;
    int oldj;
    double delta_score;
    for (iter=0; iter<getmin(xlen,ylen)*5; iter++)
    {
        delta_score=-1;
        for (i=0;i<xlen;i++)
        {
            oldj=fwdmap[i];
            for (j=0;j<ylen;j++)
            {
                oldi=invmap[j];
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

    // clean up

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
        mol2, task_user_opts, soi_params, false,
        ya, seqy, secy, resi_vec2, soi_context.k_nearest2, soi_context.chain2_bonds);
    if (ylen<3) return;

    ChainPairAlignResult chain_align_res = { 0};
    chain_align_res.d0_out = 5.0;
    ChainPairAlignOptions soi_align_opts = soi_params.align;
    soi_align_opts.fast_opt = (getmin(xlen,ylen)>1500)?true:task_user_opts.fast_opt;
    soi_align_opts.mol_type = mol1+mol2;
#ifdef _OPENMP
    if (omp_get_active_level()>0) soi_align_opts.parallel_threads = 1;
#endif
    soi_context.pair_distances.assign(ylen+1, 0.0);

    soi_align_pair(chain_align_res, xa, ya, soi_context.k_nearest1, soi_context.k_nearest2, soi_params.closeK_opt,
        seqx, seqy, secx, secy, xlen, ylen, soi_align_opts,
        common_inputs.parsed_input.sequence, soi_context.pair_distances, soi_context.chain1_bonds, soi_context.chain2_bonds,
        soi_params.mm_opt, task_user_opts.outfmt_opt, out_stream);

    AlignCommonInput task_input = common_inputs;
    task_input.user_options = task_user_opts;
    output_pair_alignment(task_input,
        xlen, ylen, chain_align_res,
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

    const PdbParseOptions pdb_parse_opts1{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2};

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

