#include "flexalign.h"
#include "output.h"

void output_flexalign_pymol(const std::string xname, const std::string yname,
                            const std::string fname_super, const std::vector<std::vector<double> > &tu_vec,
                            const int ter_opt,
                            const int mm_opt, const int split_opt, const int mirror_opt,
                            const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
                            const std::vector<std::string> &resi_vec1, const std::vector<std::string> &resi_vec2,
                            const std::string chainID1, const std::string chainID2){
    FlexalignPymolState pymol_state;
    build_pymol_hinge_map(pymol_state, seqM, seqxA, resi_vec1);
    read_pymol_structure(pymol_state, xname, resi_vec1, tu_vec, mirror_opt, ter_opt);
    build_pymol_chain_suffix(pymol_state, chainID1, chainID2, mm_opt, split_opt, ter_opt);
    build_pymol_resi_suffix(pymol_state, seqM, seqxA, seqyA, resi_vec1, resi_vec2, mm_opt);
    write_pymol_to_file(pymol_state, fname_super, yname, mm_opt, ter_opt);
}
void run_flexalign_main(
    CoordArray& xa, CoordArray& ya, const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen,
    const AlignCommonInput& common_inputs,
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &flexalign_res,
    std::ostream& os){
    flexalign_res.hingeNum = flexalign_main(
        xa, ya, seqx, seqy, secx, secy, xlen, ylen, common_inputs, mol_type, hinge_opt, ss_opt, flexalign_res, os);

    bool refine_fallback = (hinge_opt) && (flexalign_res.hingeNum <= 1) && (flexalign_res.n_ali8 < 0.6 * getmin(xlen, ylen));
    if (refine_fallback)
    {
        FlexAlignResult refine_res;
        refine_res.tu_vec.push_back(flexalign_res.tu_vec[0]);
        tu2t_u(flexalign_res.tu_vec[0], refine_res.t0, refine_res.u0);

        refine_res.hingeNum = flexalign_main(
            xa, ya, seqx, seqy, secx, secy, xlen, ylen, common_inputs, mol_type, hinge_opt, ss_opt, refine_res, os);

        double TM = (flexalign_res.TM1 > flexalign_res.TM2) ? flexalign_res.TM1 : flexalign_res.TM2;
        double TM_h = (refine_res.TM1 > refine_res.TM2) ? refine_res.TM1 : refine_res.TM2;
        if (TM_h > TM)
        {
            flexalign_res = refine_res; 
        }
    }
}
void flexalign_usbcat_main(
    CoordArray& xa,
    CoordArray& ya,
    const std::string &seqx,
    const std::string &seqy,
    const std::string &secx,
    const std::string &secy,
    int xlen,
    int ylen,
    const AlignCommonInput& common_inputs,
    const FlexAlignResult& flexalign_main_res,
    int mol_type, double best_global_max_TM,
    int sparse_val,
    const FlexalignParams& flex_params,
    FlexAlignResult& flexalign_res,
    std::ostream& os){
    FlexAlignResult global_best_align = flexalign_main_res;

    USBCATParams usb_cat_para;
    fill_usbcat_params(usb_cat_para, flex_params.hinge_opt);

    LocalDistTables dist_tables;
    build_local_dist_tables(xa, ya, usb_cat_para, dist_tables);

    double min_resid_num = best_global_max_TM * std::min(xlen, ylen);

    RegionBoundsAllChain bounds_default;
    generate_bounds(xa, ya, dist_tables, usb_cat_para, xlen, ylen, min_resid_num,
                    3.0, 4.0, 4.0, sparse_val, bounds_default);
    RegionBoundsAllChain bounds_strict;
    generate_bounds(xa, ya, dist_tables, usb_cat_para, xlen, ylen, min_resid_num,
                    2.0, 3.0, 2.0, sparse_val, bounds_strict);

    RegionBoundPool region_bound_pool;
    fill_region_bound_pools(bounds_default, bounds_strict, region_bound_pool);

    update_global_best_align(
        region_bound_pool, flex_params.hinge_set, 
        flex_params.hinge_opt, usb_cat_para, xa, 
        ya, seqx, seqy, secx, secy, common_inputs, 
        mol_type, xlen, ylen, flexalign_res.d0_out,
        global_best_align, best_global_max_TM, os);

    if (best_global_max_TM < 0)
    {
        flexalign_res.hingeNum = 0;
        return;
    }

    flexalign_res = global_best_align;
    flexalign_res.hingeNum = flexalign_res.tu_vec.size();
}
void fill_flexalign_params(
    FlexalignParams& flex_params,
    const ControlOptions& control)
{
    flex_params.mode = control.usbcat_opt ? FLEX_USBCAT : FLEX_BEST;
    flex_params.hinge_opt = control.hinge_opt;
    flex_params.hinge_set = (flex_params.mode == FLEX_USBCAT) ? control.hinge_set : false;
    flex_params.TMpass = (flex_params.mode == FLEX_USBCAT) ? control.TMpass_opt : 0.85;
    flex_params.sparse_val = 0;
}

struct FlexalignFileData
{
    std::vector<std::vector<std::string> > pdb_lines;
    std::vector<std::string> chain_ids;
    std::vector<int> mol_types;
    int chain_num;
};

struct FlexalignPairTask
{
    int file1_idx;
    int chain1_idx;
    int file2_idx;
    int chain2_idx;
    int order;
};

const FlexalignFileData* load_flexalign_file_cache(std::map<std::string, FlexalignFileData>& cache,
    const std::string& name, int infmt_opt, const PdbParseOptions& pdb_parse_opts)
{
    std::map<std::string, FlexalignFileData>::iterator cached = cache.find(name);
    if (cached != cache.end()) return &cached->second;
    FlexalignFileData& file_data = cache[name];
    file_data.chain_num = (int)get_PDB_lines(name, file_data.pdb_lines, file_data.chain_ids,
        file_data.mol_types, infmt_opt, pdb_parse_opts);
    return &file_data;
}

void run_flexalign_task(const FlexalignPairTask& task,
    const AlignCommonInput& common_inputs, const FlexalignParams& flex_params,
    const FlexalignFileData& file1, const FlexalignFileData& file2,
    const std::string& xname, const std::string& yname,
    int read_resi, std::ostream& os)
{
    UserOptions task_user_opts = common_inputs.user_options;
    task_user_opts.xname = xname;
    task_user_opts.yname = yname;

    ParsedChain chain1_data;
    if (!parse_chain(xname, file1.pdb_lines, file1.chain_ids, file1.mol_types,
            task.chain1_idx, task_user_opts.mol_opt, task_user_opts.mirror_opt, read_resi,
            task_user_opts.atom_opt, chain1_data))
        return;

    ParsedChain chain2_data;
    if (!parse_chain(yname, file2.pdb_lines, file2.chain_ids, file2.mol_types,
            task.chain2_idx, task_user_opts.mol_opt, 0, read_resi,
            task_user_opts.atom_opt, chain2_data))
        return;

    AlignCommonInput task_input = common_inputs;
    task_input.user_options = task_user_opts;
    if (task_user_opts.byresi_opt)
        extract_aln_from_resi(task_input.parsed_input.sequence, chain1_data.chain_seq, chain2_data.chain_seq,
            chain1_data.resi_vec, chain2_data.resi_vec, task_user_opts.byresi_opt);

    FlexAlignResult task_result;
    run_flexalign(chain1_data, chain2_data, task_input, flex_params, task_result, os);

    if (task_user_opts.outfmt_opt==0) print_version(USALIGN_VERSION, os);
    output_flexalign_results(chain1_data, chain2_data, task_result, task_user_opts, false, os);
}

void build_flexalign_pair_tasks(const std::vector<std::string>& chain1_list,
    const std::vector<std::string>& chain2_list,
    const std::vector<const FlexalignFileData*>& files1,
    const std::vector<const FlexalignFileData*>& files2,
    bool upper_triangular, bool dirpair, std::vector<FlexalignPairTask>& tasks)
{
    for (int file1_idx = 0; file1_idx < (int)chain1_list.size(); file1_idx++)
    {
        const FlexalignFileData& file1 = *files1[file1_idx];
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
                const FlexalignFileData& file2 = *files2[file2_idx];
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

int Flexalign(AlignCommonInput& common_inputs, const FlexalignParams& flex_params, FlexAlignResult& flex_result)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ParsedInput& parsed_input = common_inputs.parsed_input;

    int read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;  // whether to read residue index

    const PdbParseOptions pdb_parse_opts1{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2};

    std::map<std::string, FlexalignFileData> file1_cache;
    std::map<std::string, FlexalignFileData> file2_cache;
    std::vector<const FlexalignFileData*> files1(parsed_input.chain1_list.size());
    std::vector<const FlexalignFileData*> files2(parsed_input.chain2_list.size());
    for (int i=0;i<(int)parsed_input.chain1_list.size();i++)
        files1[i] = load_flexalign_file_cache(file1_cache, parsed_input.chain1_list[i],
            user_opts.infmt1_opt, pdb_parse_opts1);
    for (int j=0;j<(int)parsed_input.chain2_list.size();j++)
        files2[j] = load_flexalign_file_cache(file2_cache, parsed_input.chain2_list[j],
            user_opts.infmt2_opt, pdb_parse_opts2);

    std::vector<FlexalignPairTask> tasks;
    build_flexalign_pair_tasks(parsed_input.chain1_list, parsed_input.chain2_list,
        files1, files2, user_opts.dir_opt.size() > 0, user_opts.dirpair_opt.size() != 0, tasks);

    int task_num = (int)tasks.size();
    std::vector<std::string> out_lines(task_num);

#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads)
    for (int task_idx = 0; task_idx < task_num; task_idx++)
    {
        const FlexalignPairTask& task = tasks[task_idx];
        std::stringstream out_stream;
        run_flexalign_task(task, common_inputs, flex_params,
            *files1[task.file1_idx], *files2[task.file2_idx],
            parsed_input.chain1_list[task.file1_idx], parsed_input.chain2_list[task.file2_idx],
            read_resi, out_stream);
        out_lines[task.order] = out_stream.str();
    }

    for (int task_idx = 0; task_idx < task_num; task_idx++)
        std::cout << out_lines[task_idx];

    flex_result = FlexAlignResult();
    return 0;
}

