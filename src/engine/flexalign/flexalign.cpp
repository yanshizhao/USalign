#include "flexalign.h"
#include "output.h"

void output_flexalign_pymol(const std::string xname, const std::string yname,
                            const std::string fname_super, const std::vector<std::vector<double> > &tu_vec,
                            const int ter_opt,
                            const int mm_opt, const int split_opt, const int mirror_opt,
                            const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
                            const std::vector<std::string> &resi_vec1, const std::vector<std::string> &resi_vec2,
                            const std::string chainID1, const std::string chainID2, const int o_opt){
    FlexalignPymolState pymol_state;
    const int num_hinges = static_cast<int>(tu_vec.size());
    build_pymol_hinge_map(pymol_state, seqM, seqxA, resi_vec1);
    read_pymol_structure(pymol_state, xname, resi_vec1, tu_vec, mirror_opt, ter_opt);
    build_pymol_chain_suffix(pymol_state, chainID1, chainID2, mm_opt, split_opt, ter_opt, o_opt);
    build_pymol_resi_suffix(pymol_state, seqM, seqxA, seqyA, resi_vec1, resi_vec2, mm_opt, num_hinges, o_opt);
    write_pymol_to_file(pymol_state, fname_super, yname, mm_opt, ter_opt, o_opt);
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
    if (!parse_chain(xname, file1.pdb_lines, file1.chain_ids, file1.mol_types, task.chain1_idx,
        task_user_opts.mol_opt, task_user_opts.mirror_opt, read_resi, task_user_opts.atom_opt, chain1_data))
        return;

    ParsedChain chain2_data;
    if (!parse_chain(yname, file2.pdb_lines, file2.chain_ids, file2.mol_types, task.chain2_idx,
        task_user_opts.mol_opt, 0, read_resi, task_user_opts.atom_opt, chain2_data))
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

    const PdbParseOptions pdb_parse_opts1 = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse1, user_opts.model2parse1);
    const PdbParseOptions pdb_parse_opts2 = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt, user_opts.chain2parse2, user_opts.model2parse2);

    std::map<std::string, FlexalignFileData> file1_cache;
    std::map<std::string, FlexalignFileData> file2_cache;
    std::vector<const FlexalignFileData*> files1(parsed_input.chain1_list.size());
    std::vector<const FlexalignFileData*> files2(parsed_input.chain2_list.size());

    for (int i=0;i<(int)parsed_input.chain1_list.size();i++)
        files1[i] = load_flexalign_file_cache(file1_cache, parsed_input.chain1_list[i], user_opts.infmt1_opt, pdb_parse_opts1);

    for (int j=0;j<(int)parsed_input.chain2_list.size();j++)
        files2[j] = load_flexalign_file_cache(file2_cache, parsed_input.chain2_list[j], user_opts.infmt2_opt, pdb_parse_opts2);

    std::vector<FlexalignPairTask> tasks;
    build_flexalign_pair_tasks(parsed_input.chain1_list, parsed_input.chain2_list,
        files1, files2, user_opts.dir_opt.size() > 0, user_opts.dirpair_opt.size() != 0, tasks);

    int task_num = (int)tasks.size();
    std::vector<std::string> out_lines(task_num);

#pragma omp parallel for schedule(dynamic, 1) num_threads(common_inputs.control_options.parallel_threads) if(common_inputs.control_options.parallel_threads>1 && task_num>1)
    for (int task_idx = 0; task_idx < task_num; task_idx++)
    {
        const FlexalignPairTask& task = tasks[task_idx];
        std::stringstream out_stream;

        run_flexalign_task(task, common_inputs, flex_params, *files1[task.file1_idx], *files2[task.file2_idx],
            parsed_input.chain1_list[task.file1_idx], parsed_input.chain2_list[task.file2_idx], read_resi, out_stream);
        out_lines[task.order] = out_stream.str();
    }

    for (int task_idx = 0; task_idx < task_num; task_idx++)
        std::cout << out_lines[task_idx];

    flex_result = FlexAlignResult();
    return 0;
}

int flexalign_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy, const std::string &secx, const std::string &secy,
    int xlen, int ylen,
    const AlignCommonInput& common_inputs,
    int mol_type, int hinge_opt, int ss_opt,
    FlexAlignResult& flexalign_res,
    std::ostream& os)
{
    const UserOptions& user_opts = common_inputs.user_options;
    const ParsedInput& parsed_input = common_inputs.parsed_input;
    bool fast_opt = (std::min(xlen, ylen) > 1500) ? true : user_opts.fast_opt;

    ChainPairAlignOptions align_opts;
    fill_chain_pair_align_options(align_opts, user_opts, common_inputs.control_options);
    align_opts.fast_opt = fast_opt;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.TMcut = -1;
    align_opts.parallel_threads = 1;
#ifdef _OPENMP
    if (omp_get_active_level()==0)
        align_opts.parallel_threads = common_inputs.control_options.parallel_threads;
#endif
    align_opts.ss_opt = ss_opt;
    align_opts.mol_type = mol_type;

    std::vector<double> tu_tmp(12,0);
    int round2=flexalign_res.tu_vec.size();
    if (round2==0)
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, flexalign_res,
            xlen, ylen, parsed_input.sequence, align_opts, os);

        t_u2tu(flexalign_res.t0,flexalign_res.u0,tu_tmp);
        flexalign_res.tu_vec.push_back(tu_tmp);
    }

    flexalign_res.invmap.assign(ylen+1, -1);

    CoordArray xt;
    xt.resize(xlen);

    flexalign_res.TM1 = flexalign_res.TM2 = flexalign_res.TM3 = flexalign_res.TM4 = flexalign_res.TM5 = flexalign_res.rmsd0 = 0;
    flexalign_res.seqM = "";
    flexalign_res.seqxA = "";
    flexalign_res.seqyA = "";
    flexalign_res.n_ali = flexalign_res.n_ali8 = 0;

    score_with_cur_transform(
        xa, ya, xt, seqx, seqy,
        xlen, ylen, parsed_input.sequence, align_opts, 1, flexalign_res);
    if (round2) refine_via_split_alignment(
        xa, ya, flexalign_res.tu_vec, seqx, seqy, secx, secy,
        xlen, ylen, parsed_input.sequence, align_opts,
        flexalign_res, os);
    search_hinge_regions(
        xa, ya, xt, flexalign_res.tu_vec, seqx, seqy, secx, secy,
        xlen, ylen, hinge_opt, parsed_input.sequence, align_opts,
        flexalign_res, os);

    if (flexalign_res.tu_vec.size()<=1)
    {
        return flexalign_res.tu_vec.size();
    }
    // re-derive alignment based on tu_vec
    denoise_segments(
        xa, ya, xt, flexalign_res.tu_vec,
        xlen, ylen, align_opts, flexalign_res);

    return flexalign_res.tu_vec.size();
}

void output_flexalign_matrix(const std::string& fname_matrix,
    const DoubleMatrix&tu_vec, std::ostream& os)
{
    Vec3 t; RotMat u;
    std::stringstream ss;
    for (int hinge=0;hinge<tu_vec.size();hinge++)
    {
        tu2t_u(tu_vec[hinge],t,u);
        ss << "------ The rotation matrix to rotate Structure_1 to Structure_2 ------\n";
        ss << strfmt("m %18s %14s %14s %14s\n", "t[m]", "u[m][0]", "u[m][1]", "u[m][2]");
        for (int k = 0; k < 3; k++)
        {
            ss << strfmt("%d %18.10f %14.10f %14.10f %14.10f\n", k, clean_fmt(t[k]), clean_fmt(u[k][0]), clean_fmt(u[k][1]), clean_fmt(u[k][2]));
        }
    }
    ss << "\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):\n"
            "for(i=0; i<L; i++)\n"
            "{\n"
            "   X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i];\n"
            "   Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i];\n"
            "   Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i];\n"
            "}\n";
    if (fname_matrix == "-")
       os<<ss.str();
    else
    {
        std::fstream fout;
        fout.open(fname_matrix, std::ios::out | std::ios::trunc);
        if (fout)
        {
            fout<<ss.str();
            fout.close();
        }
        else os << "Open file to output rotation matrix fail.\n";
    }
    ss.str(std::string());
}

void read_rasmol_structure1(FlexalignRasmolState& rasmol_state,
    const std::string& xname, const DoubleMatrix& tu_vec,
    const int mm_opt, const int ter_opt, const int mirror_opt)
{
    std::ifstream fin;
    int i;
    Vec3 t; RotMat u;
    rasmol_state.after_ter=false;
    rasmol_state.asym_id="";
    fin.open(xname.c_str());
    int hinge=0;
    while (fin.good())
    {
        std::getline(fin, rasmol_state.line);
        if (ter_opt>=3 && rasmol_state.line.compare(0,3,"TER")==0) rasmol_state.after_ter=true;
        if (rasmol_state.is_mmcif==false && rasmol_state.line.size()>=54 &&
           (rasmol_state.line.compare(0, 6, "ATOM  ")==0 ||
            rasmol_state.line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (rasmol_state.line[16]!='A' && rasmol_state.line[16]!=' ') continue;
            rasmol_state.x[0]=safe_stod(rasmol_state.line.substr(30,8).c_str());
            rasmol_state.x[1]=safe_stod(rasmol_state.line.substr(38,8).c_str());
            rasmol_state.x[2]=safe_stod(rasmol_state.line.substr(46,8).c_str());
            if (mirror_opt) rasmol_state.x[2]=-rasmol_state.x[2];
            if (rasmol_state.read_resi==1) rasmol_state.resi=rasmol_state.line.substr(22,5);
            else rasmol_state.resi=rasmol_state.line.substr(22,5)+rasmol_state.line[21];
            hinge=0;
            if (rasmol_state.resi2hinge_dict.count(rasmol_state.resi)) hinge=rasmol_state.resi2hinge_dict[rasmol_state.resi];
            tu2t_u(tu_vec[hinge],t,u);
            transform(t, u, rasmol_state.x, rasmol_state.x1);

            if (rasmol_state.after_ter && rasmol_state.line.compare(0,6,"ATOM  ")==0) continue;
            rasmol_state.lig_idx1++;
            rasmol_state.buf_all_atm_lig<<rasmol_state.line.substr(0,6)<<std::setw(5)<<rasmol_state.lig_idx1
                <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (rasmol_state.chain1_sele.size() && rasmol_state.line[21]!=rasmol_state.chain1_sele[0]) continue;
            if (rasmol_state.after_ter || rasmol_state.line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && rasmol_state.asym_id!=rasmol_state.line.substr(21,1)) 
                {
                    rasmol_state.after_ter=true;
                    continue;
                }
                rasmol_state.asym_id=rasmol_state.line[21];
            }
            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1
                <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),rasmol_state.resi_aln1.end(),
                rasmol_state.line.substr(22,4))!=rasmol_state.resi_aln1.end())
            {
                rasmol_state.buf_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1
                    <<rasmol_state.line.substr(11,9)<<" A"<<rasmol_state.line.substr(22,8)
                    <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                    <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1] <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            }
            if (rasmol_state.line.substr(12,4)!=" CA " && rasmol_state.line.substr(12,4)!=" C3'") continue;
            rasmol_state.ca_idx1++;
            rasmol_state.buf_all<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1]<<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            if (std::find(rasmol_state.resi_aln1.begin(),rasmol_state.resi_aln1.end(),
                rasmol_state.line.substr(22,4))==rasmol_state.resi_aln1.end()) continue;
            if (!mm_opt) rasmol_state.buf<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" A"<<rasmol_state.line.substr(22,8)
                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                <<std::setw(8)<<rasmol_state.x1[0]<<std::setw(8)<<rasmol_state.x1[1]<<std::setw(8)<<rasmol_state.x1[2]<<'\n';
            rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1);
        }
        else if (rasmol_state.line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (rasmol_state.line.size()) break;
            }
            if (rasmol_state.line.compare(0,11,"_atom_site.")) continue;
            rasmol_state.atom_site.clear();
            rasmol_state.atom_site_pos=0;
            rasmol_state.atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=rasmol_state.atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+xname);
                if (rasmol_state.line.size()==0) continue;
                if (rasmol_state.line.compare(0,11,"_atom_site.")) break;
                rasmol_state.atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=++rasmol_state.atom_site_pos;
            }

            if (rasmol_state.is_mmcif==false)
            {
                rasmol_state.is_mmcif=true;
            }

            while(1)
            {
                rasmol_state.line_vec.clear();
                split(rasmol_state.line,rasmol_state.line_vec);
                if (rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]!="ATOM" &&
                    rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]!="HETATM") break;
                if (rasmol_state.atom_site.count("pdbx_PDB_model_num"))
                {
                    if (rasmol_state.model_index.size() && rasmol_state.model_index!=
                        rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_model_num"]])
                        break;
                    rasmol_state.model_index=rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_model_num"]];
                }

                rasmol_state.x[0]=safe_stod(rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]].c_str());
                rasmol_state.x[1]=safe_stod(rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]].c_str());
                rasmol_state.x[2]=safe_stod(rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]].c_str());
                if (mirror_opt) rasmol_state.x[2]=-rasmol_state.x[2];


                if (rasmol_state.atom_site.count("auth_seq_id"))
                    rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["auth_seq_id"]];
                else rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["label_seq_id"]];
                if (rasmol_state.atom_site.count("pdbx_PDB_ins_code") && 
                    rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]]!="?")
                    rasmol_state.resi+=rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]][0];
                else rasmol_state.resi+=" ";
                if (rasmol_state.read_resi>=2)
                {
                    if (rasmol_state.atom_site.count("auth_asym_id"))
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]];
                    else rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]];
                    if (rasmol_state.asym_id==".") rasmol_state.asym_id=" ";
                    rasmol_state.resi+=rasmol_state.asym_id[0];
                }
                hinge=0;
                if (rasmol_state.resi2hinge_dict.count(rasmol_state.resi)) hinge=rasmol_state.resi2hinge_dict[rasmol_state.resi];
                tu2t_u(tu_vec[hinge],t,u);
                transform(t, u, rasmol_state.x, rasmol_state.x1);

                if (rasmol_state.atom_site.count("label_alt_id")==0 || 
                    rasmol_state.line_vec[rasmol_state.atom_site["label_alt_id"]]=="." ||
                    rasmol_state.line_vec[rasmol_state.atom_site["label_alt_id"]]=="A")
                {
                    rasmol_state.atom=rasmol_state.line_vec[rasmol_state.atom_site["label_atom_id"]];
                    if (rasmol_state.atom[0]=='"') rasmol_state.atom=rasmol_state.atom.substr(1);
                    if (rasmol_state.atom.size() && rasmol_state.atom[rasmol_state.atom.size()-1]=='"')
                        rasmol_state.atom=rasmol_state.atom.substr(0,rasmol_state.atom.size()-1);
                    if      (rasmol_state.atom.size()==0) rasmol_state.atom="    ";
                    else if (rasmol_state.atom.size()==1) rasmol_state.atom=" "+rasmol_state.atom+"  ";
                    else if (rasmol_state.atom.size()==2) rasmol_state.atom=" "+rasmol_state.atom+" ";
                    else if (rasmol_state.atom.size()==3) rasmol_state.atom=" "+rasmol_state.atom;
                    else if (rasmol_state.atom.size()>=5) rasmol_state.atom=rasmol_state.atom.substr(0,4);
            
                    rasmol_state.AA=rasmol_state.line_vec[rasmol_state.atom_site["label_comp_id"]]; // residue name
                    if      (rasmol_state.AA.size()==1) rasmol_state.AA="  "+rasmol_state.AA;
                    else if (rasmol_state.AA.size()==2) rasmol_state.AA=" " +rasmol_state.AA;
                    else if (rasmol_state.AA.size()>=4) rasmol_state.AA=rasmol_state.AA.substr(0,3);
                
                    if (rasmol_state.atom_site.count("auth_seq_id"))
                        rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["auth_seq_id"]];
                    else rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["label_seq_id"]];
                    while (rasmol_state.resi.size()<4) rasmol_state.resi=' '+rasmol_state.resi;
                    if (rasmol_state.resi.size()>4) rasmol_state.resi=rasmol_state.resi.substr(0,4);
                
                    rasmol_state.inscode=' ';
                    if (rasmol_state.atom_site.count("pdbx_PDB_ins_code") && 
                        rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]]!="?")
                        rasmol_state.inscode=rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]][0];

                    if (rasmol_state.atom_site.count("auth_asym_id"))
                    {
                        if (rasmol_state.chain1_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]]!=rasmol_state.chain1_sele;
                        else if (ter_opt>=2 && rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]];
                    }
                    else if (rasmol_state.atom_site.count("label_asym_id"))
                    {
                        if (rasmol_state.chain1_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]]!=rasmol_state.chain1_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx1 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]];
                    }

                    if (rasmol_state.after_ter==false ||
                        rasmol_state.line_vec[rasmol_state.atom_site["group_pdb"]]=="HETATM")
                    {
                        rasmol_state.lig_idx1++;
                        rasmol_state.buf_all_atm_lig<<std::left<<std::setw(6)
                            <<rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]<<std::right
                            <<std::setw(5)<<rasmol_state.lig_idx1%100000<<' '<<rasmol_state.atom<<' '
                            <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                            <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                            <<std::setw(8)<<rasmol_state.x1[0]
                            <<std::setw(8)<<rasmol_state.x1[1]
                            <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                        if (rasmol_state.after_ter==false &&
                            rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]=="ATOM")
                        {
                            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(6)
                                <<std::setw(5)<<rasmol_state.lig_idx1%100000<<' '<<rasmol_state.atom<<' '
                                <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                                <<std::setw(8)<<rasmol_state.x1[0]
                                <<std::setw(8)<<rasmol_state.x1[1]
                                <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                            if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),
                                rasmol_state.resi_aln1.end(),rasmol_state.resi)!=rasmol_state.resi_aln1.end())
                            {
                                rasmol_state.buf_atm<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<rasmol_state.lig_idx1%100000<<' '
                                    <<rasmol_state.atom<<' '<<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                                    <<std::setw(8)<<rasmol_state.x1[0]
                                    <<std::setw(8)<<rasmol_state.x1[1]
                                    <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                            }
                            if (rasmol_state.atom==" CA " || rasmol_state.atom==" C3'")
                            {
                                rasmol_state.ca_idx1++;
                                rasmol_state.buf_all<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<rasmol_state.ca_idx1%100000<<' '<<rasmol_state.atom<<' '
                                    <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                                    <<std::setw(8)<<rasmol_state.x1[0]
                                    <<std::setw(8)<<rasmol_state.x1[1]
                                    <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                                if (!mm_opt && std::find(rasmol_state.resi_aln1.begin(),
                                    rasmol_state.resi_aln1.end(),rasmol_state.resi)!=rasmol_state.resi_aln1.end())
                                {
                                    rasmol_state.buf<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<rasmol_state.ca_idx1%100000<<' '<<rasmol_state.atom<<' '
                                    <<rasmol_state.AA<<" A"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<setiosflags(std::ios::fixed)<<std::setprecision(3)
                                    <<std::setw(8)<<rasmol_state.x1[0]
                                    <<std::setw(8)<<rasmol_state.x1[1]
                                    <<std::setw(8)<<rasmol_state.x1[2]<<'\n';
                                    rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1);
                                }
                            }
                        }
                    }
                }

                while(1)
                {
                    if (fin.good()) getline(fin, rasmol_state.line);
                    else break;
                    if (rasmol_state.line.size()) break;
                }
            }
        }
        else if (rasmol_state.line.size() && rasmol_state.is_mmcif==false)
        {
            if (ter_opt>=1 && rasmol_state.line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) rasmol_state.buf<<"TER\n";
    rasmol_state.buf_all<<"TER\n";
    if (!mm_opt) rasmol_state.buf_atm<<"TER\n";
    rasmol_state.buf_all_atm<<"TER\n";
    rasmol_state.buf_all_atm_lig<<"TER\n";
    for (i=1;i<rasmol_state.ca_idx1;i++) rasmol_state.buf_all<<"CONECT"
        <<std::setw(5)<<i%100000<<std::setw(5)<<(i+1)%100000<<'\n';
    if (!mm_opt) for (i=1;i<rasmol_state.idx_vec.size();i++) rasmol_state.buf<<"CONECT"
        <<std::setw(5)<<rasmol_state.idx_vec[i-1]%100000<<std::setw(5)<<rasmol_state.idx_vec[i]%100000<<'\n';
    rasmol_state.idx_vec.clear();
}

void read_rasmol_structure2(FlexalignRasmolState& rasmol_state,
    const std::string& yname, const int mm_opt, const int ter_opt)
{
    std::ifstream fin;
    int i;
    rasmol_state.after_ter=false;
    rasmol_state.asym_id="";
    fin.open(yname.c_str());
    while (fin.good())
    {
        std::getline(fin, rasmol_state.line);
        if (ter_opt>=3 && rasmol_state.line.compare(0,3,"TER")==0) rasmol_state.after_ter=true;
        if (rasmol_state.line.size()>=54 && (rasmol_state.line.compare(0, 6, "ATOM  ")==0 ||
            rasmol_state.line.compare(0, 6, "HETATM")==0)) // PDB format
        {
            if (rasmol_state.line[16]!='A' && rasmol_state.line[16]!=' ') continue;
            if (rasmol_state.after_ter && rasmol_state.line.compare(0,6,"ATOM  ")==0) continue;
            rasmol_state.lig_idx2++;
            rasmol_state.buf_all_atm_lig<<rasmol_state.line.substr(0,6)<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (rasmol_state.chain1_sele.size() && rasmol_state.line[21]!=rasmol_state.chain1_sele[0]) continue;
            if (rasmol_state.after_ter || rasmol_state.line.compare(0,6,"ATOM  ")) continue;
            if (ter_opt>=2)
            {
                if (rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && rasmol_state.asym_id!=rasmol_state.line.substr(21,1))
                {
                    rasmol_state.after_ter=true;
                    continue;
                }
                rasmol_state.asym_id=rasmol_state.line[21];
            }
            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),rasmol_state.resi_aln2.end(),
                rasmol_state.line.substr(22,4))!=rasmol_state.resi_aln2.end())
            {
                rasmol_state.buf_atm<<"ATOM  "<<std::setw(5)<<rasmol_state.lig_idx1+rasmol_state.lig_idx2
                    <<rasmol_state.line.substr(11,9)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            }
            if (rasmol_state.line.substr(12,4)!=" CA " && rasmol_state.line.substr(12,4)!=" C3'") continue;
            rasmol_state.ca_idx2++;
            rasmol_state.buf_all<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1+rasmol_state.ca_idx2<<' '<<rasmol_state.line.substr(12,4)
                <<' '<<rasmol_state.line.substr(17,3)<<" B"<<rasmol_state.line.substr(22,32)<<'\n';
            if (std::find(rasmol_state.resi_aln2.begin(),rasmol_state.resi_aln2.end(),rasmol_state.line.substr(22,4)
                )==rasmol_state.resi_aln2.end()) continue;
            if (!mm_opt) rasmol_state.buf<<"ATOM  "<<std::setw(5)<<rasmol_state.ca_idx1+rasmol_state.ca_idx2<<' '
                <<rasmol_state.line.substr(12,4)<<' '<<rasmol_state.line.substr(17,3)<<" B"
                <<rasmol_state.line.substr(22,32)<<'\n';
            rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1+rasmol_state.ca_idx2);
        }
        else if (rasmol_state.line.compare(0,5,"loop_")==0) // PDBx/mmCIF
        {
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (rasmol_state.line.size()) break;
            }
            if (rasmol_state.line.compare(0,11,"_atom_site.")) continue;
            rasmol_state.atom_site.clear();
            rasmol_state.atom_site_pos=0;
            rasmol_state.atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=rasmol_state.atom_site_pos;
            while(1)
            {
                if (fin.good()) getline(fin, rasmol_state.line);
                else PrintErrorAndQuit("ERROR! Unexpected end of "+yname);
                if (rasmol_state.line.size()==0) continue;
                if (rasmol_state.line.compare(0,11,"_atom_site.")) break;
                rasmol_state.atom_site[rasmol_state.line.substr(11,rasmol_state.line.size()-12)]=++rasmol_state.atom_site_pos;
            }

            while(1)
            {
                rasmol_state.line_vec.clear();
                split(rasmol_state.line,rasmol_state.line_vec);
                if (rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]!="ATOM" &&
                    rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]!="HETATM") break;
                if (rasmol_state.atom_site.count("pdbx_PDB_model_num"))
                {
                    if (rasmol_state.model_index.size() && rasmol_state.model_index!=
                        rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_model_num"]])
                        break;
                    rasmol_state.model_index=rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_model_num"]];
                }

                if (rasmol_state.atom_site.count("label_alt_id")==0 || 
                    rasmol_state.line_vec[rasmol_state.atom_site["label_alt_id"]]=="." ||
                    rasmol_state.line_vec[rasmol_state.atom_site["label_alt_id"]]=="A")
                {
                    rasmol_state.atom=rasmol_state.line_vec[rasmol_state.atom_site["label_atom_id"]];
                    if (rasmol_state.atom[0]=='"') rasmol_state.atom=rasmol_state.atom.substr(1);
                    if (rasmol_state.atom.size() && rasmol_state.atom[rasmol_state.atom.size()-1]=='"')
                        rasmol_state.atom=rasmol_state.atom.substr(0,rasmol_state.atom.size()-1);
                    if      (rasmol_state.atom.size()==0) rasmol_state.atom="    ";
                    else if (rasmol_state.atom.size()==1) rasmol_state.atom=" "+rasmol_state.atom+"  ";
                    else if (rasmol_state.atom.size()==2) rasmol_state.atom=" "+rasmol_state.atom+" ";
                    else if (rasmol_state.atom.size()==3) rasmol_state.atom=" "+rasmol_state.atom;
                    else if (rasmol_state.atom.size()>=5) rasmol_state.atom=rasmol_state.atom.substr(0,4);
            
                    rasmol_state.AA=rasmol_state.line_vec[rasmol_state.atom_site["label_comp_id"]]; // residue name
                    if      (rasmol_state.AA.size()==1) rasmol_state.AA="  "+rasmol_state.AA;
                    else if (rasmol_state.AA.size()==2) rasmol_state.AA=" " +rasmol_state.AA;
                    else if (rasmol_state.AA.size()>=4) rasmol_state.AA=rasmol_state.AA.substr(0,3);
                
                    if (rasmol_state.atom_site.count("auth_seq_id"))
                        rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["auth_seq_id"]];
                    else rasmol_state.resi=rasmol_state.line_vec[rasmol_state.atom_site["label_seq_id"]];
                    while (rasmol_state.resi.size()<4) rasmol_state.resi=' '+rasmol_state.resi;
                    if (rasmol_state.resi.size()>4) rasmol_state.resi=rasmol_state.resi.substr(0,4);
                
                    rasmol_state.inscode=' ';
                    if (rasmol_state.atom_site.count("pdbx_PDB_ins_code") && 
                        rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]]!="?")
                        rasmol_state.inscode=rasmol_state.line_vec[rasmol_state.atom_site["pdbx_PDB_ins_code"]][0];
                    
                    if (rasmol_state.atom_site.count("auth_asym_id"))
                    {
                        if (rasmol_state.chain2_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]]!=rasmol_state.chain2_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["auth_asym_id"]];
                    }
                    else if (rasmol_state.atom_site.count("label_asym_id"))
                    {
                        if (rasmol_state.chain2_sele.size()) rasmol_state.after_ter
                            =rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]]!=rasmol_state.chain2_sele;
                        if (ter_opt>=2 && rasmol_state.ca_idx2 && rasmol_state.asym_id.size() && 
                            rasmol_state.asym_id!=rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]])
                            rasmol_state.after_ter=true;
                        rasmol_state.asym_id=rasmol_state.line_vec[rasmol_state.atom_site["label_asym_id"]];
                    }
                    if (rasmol_state.after_ter==false || 
                        rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]=="HETATM")
                    {
                        rasmol_state.lig_idx2++;
                        rasmol_state.buf_all_atm_lig<<std::left<<std::setw(6)
                            <<rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]<<std::right
                            <<std::setw(5)<<(rasmol_state.lig_idx1+rasmol_state.lig_idx2)%100000<<' '
                            <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]]
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]]
                            <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]]
                            <<'\n';
                        if (rasmol_state.after_ter==false &&
                            rasmol_state.line_vec[rasmol_state.atom_site["group_PDB"]]=="ATOM")
                        {
                            rasmol_state.buf_all_atm<<"ATOM  "<<std::setw(6)
                                <<std::setw(5)<<(rasmol_state.lig_idx1+rasmol_state.lig_idx2)%100000<<' '
                                <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]]
                                <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]]
                                <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]]
                                <<'\n';
                            if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),
                                rasmol_state.resi_aln2.end(),rasmol_state.resi)!=rasmol_state.resi_aln2.end())
                            {
                                rasmol_state.buf_atm<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<(rasmol_state.lig_idx1+rasmol_state.lig_idx2)%100000<<' '
                                    <<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]]
                                    <<'\n';
                            }
                            if (rasmol_state.atom==" CA " || rasmol_state.atom==" C3'")
                            {
                                rasmol_state.ca_idx2++;
                                rasmol_state.buf_all<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<(rasmol_state.ca_idx1+rasmol_state.ca_idx2)%100000
                                    <<' '<<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]]
                                    <<'\n';
                                if (!mm_opt && std::find(rasmol_state.resi_aln2.begin(),
                                    rasmol_state.resi_aln2.end(),rasmol_state.resi)!=rasmol_state.resi_aln2.end())
                                {
                                    rasmol_state.buf<<"ATOM  "<<std::setw(6)
                                    <<std::setw(5)<<(rasmol_state.ca_idx1+rasmol_state.ca_idx2)%100000
                                    <<' '<<rasmol_state.atom<<' '<<rasmol_state.AA<<" B"<<rasmol_state.resi<<rasmol_state.inscode<<"   "
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_x"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_y"]]
                                    <<std::setw(8)<<rasmol_state.line_vec[rasmol_state.atom_site["Cartn_z"]]
                                    <<'\n';
                                    rasmol_state.idx_vec.push_back(rasmol_state.ca_idx1+rasmol_state.ca_idx2);
                                }
                            }
                        }
                    }
                }

                if (fin.good()) getline(fin, rasmol_state.line);
                else break;
            }
        }
        else if (rasmol_state.line.size())
        {
            if (ter_opt>=1 && rasmol_state.line.compare(0,3,"END")==0) break;
        }
    }
    fin.close();
    if (!mm_opt) rasmol_state.buf<<"TER\n";
    rasmol_state.buf_all<<"TER\n";
    if (!mm_opt) rasmol_state.buf_atm<<"TER\n";
    rasmol_state.buf_all_atm<<"TER\n";
    rasmol_state.buf_all_atm_lig<<"TER\n";
    for (i=rasmol_state.ca_idx1+1;i<rasmol_state.ca_idx1+rasmol_state.ca_idx2;i++) rasmol_state.buf_all<<"CONECT"
        <<std::setw(5)<<i%100000<<std::setw(5)<<(i+1)%100000<<'\n';
    for (i=1;i<rasmol_state.idx_vec.size();i++) rasmol_state.buf<<"CONECT"
        <<std::setw(5)<<rasmol_state.idx_vec[i-1]%100000<<std::setw(5)<<rasmol_state.idx_vec[i]%100000<<'\n';
    rasmol_state.idx_vec.clear();

}

void output_flexalign_rasmol(const std::string xname, const std::string yname,
    const std::string fname_super,const DoubleMatrix&tu_vec,
    const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden)
{
    FlexalignRasmolState rasmol_state;
    rasmol_state.ca_idx1=0;
    rasmol_state.lig_idx1=0;
    rasmol_state.ca_idx2=0;
    rasmol_state.lig_idx2=0;

    build_resi_hinge_map(rasmol_state, seqM, seqxA, resi_vec1, xlen);
    rasmol_state.resi=resi_vec1[0];
    rasmol_state.read_resi=rasmol_state.resi.size()-4;

    build_rasmol_headers(rasmol_state, xname, yname, chainID1, chainID2,
        xlen, ylen, d0A, n_ali8, rmsd, TM1, Liden, mm_opt);

    if (!mm_opt)
    {
        if (split_opt==2 && ter_opt>=1) // align one chain from model 1
        {
            rasmol_state.chain1_sele=chainID1.substr(1);
            rasmol_state.chain2_sele=chainID2.substr(1);
        }
        else if (split_opt==2 && ter_opt==0) // align one chain from each model
        {
            int i;
            for (i=1;i<chainID1.size();i++) if (chainID1[i]==',') break;
            rasmol_state.chain1_sele=chainID1.substr(i+1);
            for (i=1;i<chainID2.size();i++) if (chainID2[i]==',') break;
            rasmol_state.chain2_sele=chainID2.substr(i+1);
        }
    }


    rasmol_state.is_mmcif=false;

    extract_aligned_resi(rasmol_state, seqM, seqxA, seqyA, resi_vec1, resi_vec2, mm_opt,
        static_cast<int>(tu_vec.size()));

    read_rasmol_structure1(rasmol_state, xname, tu_vec, mm_opt, ter_opt, mirror_opt);

    read_rasmol_structure2(rasmol_state, yname, mm_opt, ter_opt);

    write_rasmol_to_file(rasmol_state, fname_super, mm_opt);
}

void read_pymol_mmcif_loop(FlexalignPymolReadState& read_state,
    const std::string& xname,
    const DoubleMatrix& tu_vec, const int mirror_opt)
{
    std::stringstream& buf = read_state.buf;
    read_state.infmt = 3;
    buf << read_state.line << '\n';
    read_pymol_line_or_quit(read_state, xname);
    while (!read_state.line.size())
        read_pymol_line_or_quit(read_state, xname);
    buf << read_state.line << '\n';
    if (read_state.line.compare(0, 11, "_atom_site."))
        return;
    std::map<std::string, int> atom_site;
    size_t atom_site_pos = 0;
    atom_site[Trim(read_state.line.substr(11))] = atom_site_pos;
    while (1)
    {
        read_pymol_line_or_quit(read_state, xname);
        while (!read_state.line.size())
            read_pymol_line_or_quit(read_state, xname);
        if (read_state.line.compare(0, 11, "_atom_site."))
            break;
        atom_site[Trim(read_state.line.substr(11))] = ++atom_site_pos;
        buf << read_state.line << '\n';
    }

    if (atom_site.count("group_PDB") *
            atom_site.count("Cartn_x") *
            atom_site.count("Cartn_y") *
            atom_site.count("Cartn_z") ==
        0)
    {
        buf << read_state.line << '\n';
        std::cerr << "Warning! Missing one of the following _atom_site data items: group_PDB, Cartn_x, Cartn_y, Cartn_z" << std::endl;
        return;
    }

    std::vector<std::string> line_vec;
    while (1)
    {
        line_vec.clear();
        split(read_state.line, line_vec);
        if (line_vec[atom_site["group_PDB"]] != "ATOM" &&
            line_vec[atom_site["group_PDB"]] != "HETATM")
            break;

        Vec3 x;
        x[0] = atof(line_vec[atom_site["Cartn_x"]].c_str());
        x[1] = atof(line_vec[atom_site["Cartn_y"]].c_str());
        x[2] = atof(line_vec[atom_site["Cartn_z"]].c_str());
        if (mirror_opt)
            x[2] = -x[2];

        if (atom_site.count("auth_seq_id"))
            read_state.resi = line_vec[atom_site["auth_seq_id"]];
        else
            read_state.resi = line_vec[atom_site["label_seq_id"]];
        if (atom_site.count("pdbx_PDB_ins_code") &&
            line_vec[atom_site["pdbx_PDB_ins_code"]] != "?")
            read_state.resi += line_vec[atom_site["pdbx_PDB_ins_code"]][0];
        else
            read_state.resi += " ";
        if (read_state.read_resi >= 2)
        {
            if (atom_site.count("auth_asym_id"))
                read_state.asym_id = line_vec[atom_site["auth_asym_id"]];
            else
                read_state.asym_id = line_vec[atom_site["label_asym_id"]];
            if (read_state.asym_id == ".")
                read_state.asym_id = " ";
            read_state.resi += read_state.asym_id[0];
        }
        int hinge = 0;
        if (read_state.resi2hinge_dict.count(read_state.resi))
            hinge = read_state.resi2hinge_dict[read_state.resi];
        tu2t_u(tu_vec[hinge], read_state.t, read_state.u);
        Vec3 x1;
        transform(read_state.t, read_state.u, x, x1);

        for (atom_site_pos = 0; atom_site_pos < atom_site.size(); atom_site_pos++)
        {
            if (atom_site_pos == atom_site["Cartn_x"])
                buf << std::setiosflags(std::ios::fixed) << std::setprecision(3)
                    << std::setw(8) << x1[0] << ' ';
            else if (atom_site_pos == atom_site["Cartn_y"])
                buf << std::setiosflags(std::ios::fixed) << std::setprecision(3)
                    << std::setw(8) << x1[1] << ' ';
            else if (atom_site_pos == atom_site["Cartn_z"])
                buf << std::setiosflags(std::ios::fixed) << std::setprecision(3)
                    << std::setw(8) << x1[2] << ' ';
            else
                buf << line_vec[atom_site_pos] << ' ';
        }
        buf << '\n';

        if (!read_pymol_line(read_state))
            break;
    }
    if (read_state.compress_type ? read_state.fin_gz.good() : read_state.fin.good())
        buf << read_state.line << '\n';
}

void read_pymol_structure(FlexalignPymolState& pymol_state,
    const std::string& xname, const std::vector<std::string>& resi_vec1,
    const DoubleMatrix& tu_vec, const int mirror_opt, const int ter_opt)
{
    std::stringstream& buf = pymol_state.buf;
    FlexalignPymolReadState read_state(pymol_state);
    read_state.compress_type = COMPRESS_FILE;
    open_file_for_read(xname, read_state.fin, read_state.fin_gz, read_state.compress_type);
    read_state.resi = resi_vec1[0];
    read_state.read_resi = read_state.resi.size() - 4;
    read_state.infmt = -1; // 0 - PDB, 3 - PDBx/mmCIF
    while (read_pymol_line(read_state))
    {
        const std::string& line = read_state.line;
        if (line.compare(0, 6, "ATOM  ") == 0 ||
            line.compare(0, 6, "HETATM") == 0) // PDB format
        {
            read_pymol_pdb_atom_line(read_state, tu_vec, mirror_opt);
        }
        else if (line.compare(0, 5, "loop_") == 0) // PDBx/mmCIF
        {
            read_pymol_mmcif_loop(read_state, xname, tu_vec, mirror_opt);
        }
        else if (line.size())
        {
            buf << line << '\n';
            if (ter_opt >= 1 && line.compare(0, 3, "END") == 0)
                break;
        }
    }
    if (read_state.compress_type)
        read_state.fin_gz.close();
    else
        read_state.fin.close();
    pymol_state.infmt = read_state.infmt;
}

void build_pymol_resi_suffix(FlexalignPymolState& pymol_state,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2,
    const int mm_opt, const int num_hinges, const int o_opt)
{
    int i1 = -1;
    int i2 = -1;
    // Arrays to store selections and bonds separately for each hinge
    pymol_state.resi1_sele.assign(num_hinges, "");
    pymol_state.resi2_sele.assign(num_hinges, "");
    pymol_state.resi1_bond.assign(num_hinges, "");
    pymol_state.resi2_bond.assign(num_hinges, "");
    std::vector<std::string> prev_resi1(num_hinges, "");
    std::vector<std::string> prev_resi2(num_hinges, "");
    if (!mm_opt)
    {
        for (int i = 0; i < seqM.size(); i++)
        {
            i1 += (seqxA[i] != '-' && seqxA[i] != '*');
            i2 += (seqyA[i] != '-');
            if (seqM[i] == ' ' || seqxA[i] == '*')
                continue;
            std::string curr_resi1 = resi_vec1[i1].substr(0, 4);
            std::string curr_resi2 = resi_vec2[i2].substr(0, 4);
            // Extract hinge index from the alignment mapping (seqM)
            int hinge_idx = hinge_char_to_index(seqM[i]);
            // Safety check to prevent index out of bounds
            if (hinge_idx >= num_hinges)
                hinge_idx = num_hinges - 1;
            if (o_opt == 3)
            {
                append_chimerax_resi_sele(pymol_state.resi1_sele[hinge_idx], curr_resi1);
                append_chimerax_resi_sele(pymol_state.resi2_sele[hinge_idx], curr_resi2);
            }
            else
            {
                if (pymol_state.resi1_sele[hinge_idx].size() == 0)
                    pymol_state.resi1_sele[hinge_idx] = "i. " + curr_resi1;
                else
                {
                    pymol_state.resi1_sele[hinge_idx] += " or i. " + curr_resi1;
                    pymol_state.resi1_bond[hinge_idx] += "bond structure1 and i. " + prev_resi1[hinge_idx] +
                                  ", i. " + curr_resi1 + "\n";
                }
                if (pymol_state.resi2_sele[hinge_idx].size() == 0)
                    pymol_state.resi2_sele[hinge_idx] = "i. " + curr_resi2;
                else
                {
                    pymol_state.resi2_sele[hinge_idx] += " or i. " + curr_resi2;
                    pymol_state.resi2_bond[hinge_idx] += "bond structure2 and i. " + prev_resi2[hinge_idx] +
                                  ", i. " + curr_resi2 + "\n";
                }
            }
            prev_resi1[hinge_idx] = curr_resi1;
            prev_resi2[hinge_idx] = curr_resi2;
        }
        // Wrap final selections in parentheses
        if (o_opt != 3)
        {
            for (int h = 0; h < num_hinges; h++)
            {
                if (pymol_state.resi1_sele[h].size())
                    pymol_state.resi1_sele[h] = " and ( " + pymol_state.resi1_sele[h] + ")";
                if (pymol_state.resi2_sele[h].size())
                    pymol_state.resi2_sele[h] = " and ( " + pymol_state.resi2_sele[h] + ")";
            }
        }
    }
}

void write_pymol_to_file(const FlexalignPymolState& pymol_state,
    const std::string& fname_super,
    const std::string& yname, const int mm_opt, const int ter_opt, const int o_opt)
{
    std::string fname_super_full = fname_super;
    if (pymol_state.infmt == 0)
        fname_super_full += ".pdb";
    else if (pymol_state.infmt == 3)
        fname_super_full += ".cif";
    const std::string struct1_basename = get_basename(fname_super_full);
    std::string struct2_basename;
    copy_structure2_to_outdir(yname, fname_super, struct2_basename);
    std::ofstream fp;
    fp.open(fname_super_full.c_str());
    fp << pymol_state.buf.str();
    fp.close();

    const std::string& chain1_sele = pymol_state.chain1_sele;
    const std::string& chain2_sele = pymol_state.chain2_sele;
    const std::vector<std::string>& resi1_sele = pymol_state.resi1_sele;
    const std::vector<std::string>& resi2_sele = pymol_state.resi2_sele;
    const std::vector<std::string>& resi1_bond = pymol_state.resi1_bond;
    const std::vector<std::string>& resi2_bond = pymol_state.resi2_bond;
    const int num_hinges = static_cast<int>(resi1_sele.size());
    // Define 35 distinct PyMOL colors for Structure 1
    // The first 10 are highly contrasting for standard cases
    const char* color1_list[] = {
        "blue", "red", "green", "yellow", "magenta", "cyan", "orange", "purple", "brown", "pink", // Top 10 highly contrasting
        "marine", "firebrick", "forest", "gold", "violet", "deepteal", "salmon", "aquamarine", "chocolate", "hotpink",
        "lime", "olive", "slate", "paleblue", "palegreen", "paleyellow", "lightpink", "lightblue", "lightmagenta",
        "dirtyviolet", "dash", "tv_blue", "tv_red", "tv_green" // Extended palette to safely cover >31 hinges
    };
    const int color_palette_size = static_cast<int>(sizeof(color1_list) / sizeof(color1_list[0]));
    std::stringstream buf_pymol;
    std::vector<std::string> pml_list;
    pml_list.push_back(fname_super);
    pml_list.push_back(fname_super + "_atm");
    pml_list.push_back(fname_super + "_all");
    pml_list.push_back(fname_super + "_all_atm");
    pml_list.push_back(fname_super + "_all_atm_lig");

    for (int p = 0; p < pml_list.size(); p++)
    {
        if (mm_opt && p <= 1)
            continue;
        if (o_opt == 3)
        {
            buf_pymol
                << "#!/usr/bin/env chimerax --script\n"
                << "open \"" << struct1_basename << "\" name structure1\n"
                << "open \"" << struct2_basename << "\" name structure2\n";
            if (p == 0)
            {
                buf_pymol
                    << "hide all cartoon\n"
                    << "hide all atoms\n"
                    << "hide all bonds\n";
                for (int h = 0; h < num_hinges; h++)
                {
                    if (resi1_sele[h].size())
                        buf_pymol << "show #1" << chain1_sele << resi1_sele[h] << " cartoon\n";
                    if (resi2_sele[h].size())
                        buf_pymol << "show #2" << chain2_sele << resi2_sele[h] << " cartoon\n";
                }
            }
            else if (p == 1)
            {
                buf_pymol
                    << "hide all cartoon\n"
                    << "hide all atoms\n"
                    << "hide all bonds\n";
                for (int h = 0; h < num_hinges; h++)
                {
                    if (resi1_sele[h].size())
                        buf_pymol
                            << "show #1" << chain1_sele << resi1_sele[h] << " cartoon\n"
                            << "show #1" << chain1_sele << resi1_sele[h] << " & nucleic cartoon,bonds\n";
                    if (resi2_sele[h].size())
                        buf_pymol
                            << "show #2" << chain2_sele << resi2_sele[h] << " cartoon\n"
                            << "show #2" << chain2_sele << resi2_sele[h] << " & nucleic cartoon,bonds\n";
                }
            }
            else if (p == 2)
            {
                buf_pymol
                    << "hide all atoms\n"
                    << "hide all bonds\n"
                    << "cartoon style protein modeh default arrows f xsect oval width 1 thick 1\n"
                    << "cartoon style nucleic xsect oval width 1 thick 1\n"
                    << "hide all cartoon\n"
                    << "show #1" << chain1_sele << " cartoon\n"
                    << "show #2" << chain2_sele << " cartoon\n";
            }
            else if (p == 3)
            {
                buf_pymol
                    << "hide all bonds\n"
                    << "hide all atoms\n"
                    << "hide all cartoon\n"
                    << "show #1" << chain1_sele << " cartoon\n"
                    << "show #1" << chain1_sele << " & nucleic cartoon,bonds\n"
                    << "show #2" << chain2_sele << " cartoon\n"
                    << "show #2" << chain2_sele << " & nucleic cartoon,bonds\n";
            }
            else if (p == 4)
            {
                buf_pymol
                    << "hide all bonds\n"
                    << "hide all atoms\n"
                    << "hide all cartoon\n"
                    << "show #1" << chain1_sele << " cartoon\n"
                    << "show #2" << chain2_sele << " cartoon\n"
                    << "show solvent atoms,bonds\n";
            }
            buf_pymol
                << "color #1 grey\n"
                << "color #2 white\n";
            for (int h = 0; h < num_hinges; h++)
            {
                if (resi1_sele[h].size())
                    buf_pymol << "color #1" << chain1_sele << resi1_sele[h] << " "
                        << color1_list[h % color_palette_size] << "\n";
            }
            buf_pymol << "set bgColor white\n" << std::endl;
            fp.open((pml_list[p] + ".cxc").c_str());
            fp << buf_pymol.str();
            fp.close();
            buf_pymol.str(std::string());
            continue;
        }
        buf_pymol
            << "#!/usr/bin/env pymol\n"
            << "cmd.load(\"" << struct1_basename << "\", \"structure1\")\n"
            << "cmd.load(\"" << struct2_basename << "\", \"structure2\")\n"
            << "hide all\n"
            << "set all_states, " << ((ter_opt == 0) ? "on" : "off") << '\n';
        if (p == 0) // .pml
        {
            if (chain1_sele.size())
                buf_pymol
                    << "remove structure1 and not " << chain1_sele.substr(4) << "\n";
            if (chain2_sele.size())
                buf_pymol
                    << "remove structure2 and not " << chain2_sele.substr(4) << "\n";
            buf_pymol
                << "remove not n. CA and not n. C3'\n";
            // Add bonds for each hinge
            for (int h = 0; h < num_hinges; h++)
                buf_pymol << resi1_bond[h] << resi2_bond[h];
            // Show sticks for each aligned hinge region
            for (int h = 0; h < num_hinges; h++)
            {
                if (resi1_sele[h].size())
                    buf_pymol << "show stick, structure1" << chain1_sele << resi1_sele[h] << "\n";
                if (resi2_sele[h].size())
                    buf_pymol << "show stick, structure2" << chain2_sele << resi2_sele[h] << "\n";
            }
        }
        else if (p == 1) // _atm.pml
        {
            // Show cartoons for aligned hinge regions
            for (int h = 0; h < num_hinges; h++)
            {
                if (resi1_sele[h].size())
                    buf_pymol << "show cartoon, structure1" << chain1_sele << resi1_sele[h] << "\n";
                if (resi2_sele[h].size())
                    buf_pymol << "show cartoon, structure2" << chain2_sele << resi2_sele[h] << "\n";
            }
        }
        else if (p == 2) // _all.pml
        {
            buf_pymol
                << "show ribbon, structure1" << chain1_sele << "\n"
                << "show ribbon, structure2" << chain2_sele << "\n";
        }
        else if (p == 3) // _all_atm.pml
        {
            buf_pymol
                << "show cartoon, structure1" << chain1_sele << "\n"
                << "show cartoon, structure2" << chain2_sele << "\n";
        }
        else if (p == 4) // _all_atm_lig.pml
        {
            buf_pymol
                << "show cartoon, structure1\n"
                << "show cartoon, structure2\n"
                << "show stick, not polymer\n"
                << "show sphere, not polymer\n";
        }
        // Create selection groups in PyMOL for easy user interaction
        for (int h = 0; h < num_hinges; h++)
        {
            if (resi1_sele[h].size())
                buf_pymol << "select hinge_" << h << "_str1, structure1" << chain1_sele << resi1_sele[h] << "\n";
            if (resi2_sele[h].size())
                buf_pymol << "select hinge_" << h << "_str2, structure2" << chain2_sele << resi2_sele[h] << "\n";
        }
        // Global color setting: Structure 1 unaligned gets light gray, Structure 2 gets solid darker gray
        buf_pymol
            << "color grey, structure1\n"
            << "color white, structure2\n";
        // Color each hinge in Structure 1 with distinct colors
        for (int h = 0; h < num_hinges; h++)
        {
            if (resi1_sele[h].size())
                buf_pymol << "color " << color1_list[h % color_palette_size] << ", hinge_" << h << "_str1\n";
        }
        buf_pymol
            << "set ribbon_width, 6\n"
            << "set stick_radius, 0.3\n"
            << "set sphere_scale, 0.25\n"
            << "set ray_shadow, 0\n"
            << "bg_color white\n"
            << "set transparency=0.2\n"
            << "zoom polymer and ((structure1" << chain1_sele
            << ") or (structure2" << chain2_sele << "))\n"
            << std::endl;

        fp.open((pml_list[p] + ".pml").c_str());
        fp << buf_pymol.str();
        fp.close();
        buf_pymol.str(std::string());
    }
}

void output_flexalign_results(
    const ParsedChain& chain1_data, const ParsedChain& chain2_data,
    const FlexAlignResult& flexalign_res,
    const UserOptions& opts,
    const int mm_opt, std::ostream& os)
{
    const std::string xname = opts.xname.substr(opts.dir1_opt.size() + opts.dir_opt.size() + opts.dirpair_opt.size());
    const std::string yname = opts.yname.substr(opts.dir2_opt.size() + opts.dir_opt.size() + opts.dirpair_opt.size());

    const std::string& chainID1 = chain1_data.chain_id;
    const std::string& chainID2 = chain2_data.chain_id;
    const int xlen = chain1_data.chain_len;
    const int ylen = chain2_data.chain_len;
    const std::vector<std::string>& resi_vec1 = chain1_data.resi_vec;
    const std::vector<std::string>& resi_vec2 = chain2_data.resi_vec;

    if (opts.outfmt_opt<=0)
    {
        fcout(os, "\nName of Structure_1: %s%s (to be superimposed onto Structure_2)\n",
            xname, chainID1);
        fcout(os, "Name of Structure_2: %s%s\n", yname, chainID2);
        fcout(os, "Length of Structure_1: %d residues\n", xlen);
        fcout(os, "Length of Structure_2: %d residues\n\n", ylen);

        if (opts.i_opt)
            fcout(os, "User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", flexalign_res.TM_ali, flexalign_res.L_ali, flexalign_res.rmsd_ali);

        fcout(os, "Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", flexalign_res.n_ali8, flexalign_res.rmsd0, (flexalign_res.n_ali8>0)?flexalign_res.Liden/flexalign_res.n_ali8:0);
        fcout(os, "TM-score= %6.5f (normalized by length of Structure_1: L=%d, d0=%.2f)\n", flexalign_res.TM2, xlen, flexalign_res.d0B);
        fcout(os, "TM-score= %6.5f (normalized by length of Structure_2: L=%d, d0=%.2f)\n", flexalign_res.TM1, ylen, flexalign_res.d0A);

        if (opts.a_opt==1)
            fcout(os, "TM-score= %6.5f (if normalized by average length of two structures: L=%.1f, d0=%.2f)\n", flexalign_res.TM3, (xlen+ylen)*0.5, flexalign_res.d0a);
        if (opts.u_opt)
            fcout(os, "TM-score= %6.5f (normalized by user-specified L=%.2f and d0=%.2f)\n", flexalign_res.TM4, opts.Lnorm_ass, flexalign_res.d0u);
        if (opts.d_opt)
            fcout(os, "TM-score= %6.5f (scaled by user-specified d0=%.2f, and L=%d)\n", flexalign_res.TM5, opts.d0_scale, ylen);
        os << "(You should use TM-score normalized by length of the reference structure)\n";

        //output alignment
        os << "\n([0-9,a-z,A-Z] denote different aligned fragment pairs separated by different hinges)\n";
        os << flexalign_res.seqxA << "\n";
        os << flexalign_res.seqM << "\n";
        os << flexalign_res.seqyA << "\n";
    }
    else if (opts.outfmt_opt==1)
    {
        fcout(os, ">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname, chainID1, xlen, flexalign_res.d0B, flexalign_res.Liden/xlen, flexalign_res.TM2);
        os << flexalign_res.seqxA << "\n";
        fcout(os, ">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname, chainID2, ylen, flexalign_res.d0A, flexalign_res.Liden/ylen, flexalign_res.TM1);
        os << flexalign_res.seqyA << "\n";

        fcout(os, "# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            flexalign_res.n_ali8, flexalign_res.rmsd0, (flexalign_res.n_ali8>0)?flexalign_res.Liden/flexalign_res.n_ali8:0);

        if (opts.i_opt)
            fcout(os, "# User-specified initial alignment: TM=%.5lf\tLali=%4d\trmsd=%.3lf\n", flexalign_res.TM_ali, flexalign_res.L_ali, flexalign_res.rmsd_ali);

        if(opts.a_opt)
            fcout(os, "# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", flexalign_res.TM3, (xlen+ylen)*0.5, flexalign_res.d0a);

        if(opts.u_opt)
            fcout(os, "# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", flexalign_res.TM4, opts.Lnorm_ass, flexalign_res.d0u);

        if(opts.d_opt)
            fcout(os, "# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", flexalign_res.TM5, opts.d0_scale, ylen);

        os << "$$$$\n";
    }
    else if (opts.outfmt_opt==2)
    {
        fcout(os, "%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            flexalign_res.TM2, flexalign_res.TM1, flexalign_res.rmsd0, flexalign_res.Liden/xlen, flexalign_res.Liden/ylen, (flexalign_res.n_ali8>0)?flexalign_res.Liden/flexalign_res.n_ali8:0,
            xlen, ylen, flexalign_res.n_ali8);
    }
    os << std::endl;

    if (opts.m_opt && !opts.fname_matrix.empty()) output_flexalign_matrix(
            opts.fname_matrix, flexalign_res.tu_vec, os);

    if (opts.o_opt==1 || opts.o_opt==3) output_flexalign_pymol(xname, yname, opts.fname_super, flexalign_res.tu_vec,
            opts.ter_opt, mm_opt, opts.split_opt, opts.mirror_opt, flexalign_res.seqM, flexalign_res.seqxA, flexalign_res.seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2, opts.o_opt);
    else if (opts.o_opt==2)
        output_flexalign_rasmol(xname, yname, opts.fname_super, flexalign_res.tu_vec,
            opts.ter_opt, mm_opt, opts.split_opt, opts.mirror_opt, flexalign_res.seqM, flexalign_res.seqxA, flexalign_res.seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2,
            xlen, ylen, flexalign_res.d0A, flexalign_res.n_ali8, flexalign_res.rmsd0, flexalign_res.TM1, flexalign_res.Liden);
}

void recompute_flexalign_scores(const std::string &seqM, const std::string &seqxA,
    const std::string &seqyA, std::vector<char>& seqM_char, std::vector<double>& di_vec,
    const std::vector<int>& invmap, DoubleMatrix& tu_vec, CoordArray& xa,
    CoordArray& ya, CoordArray& xt, Vec3& t0, RotMat& u0,
    const int xlen, const int ylen,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &rmsd0, double &Liden,
    const double d0A, const double d0B, const double d0a, const double d0u,
    const double d0_scale, const double Lnorm_ass,
    const int a_opt, const bool u_opt, const bool d_opt, const int n_ali8)
{
    for (int hinge=(int)tu_vec.size()-1; hinge>=0; hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        do_rotation(xa, xt, xlen, t0, u0);
        char cur_seg = hinge + '0';
        for (int j=0; j<ylen; j++)
        {
            int i = invmap[j];
            if (i < 0) continue;
            if (seqM_char[j] != cur_seg) continue;
            double d = sqrt(dist(xt[i], ya[j]));
            if (di_vec[j] < 0 || d <= di_vec[j])
                di_vec[j] = d;
        }
    }

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    Liden=0;
    for (int r=0; r<(int)seqM.size(); r++)
    {
        if (seqM[r]!=' ')
            Liden+=(seqxA[r]==seqyA[r]);
    }

    for(int j=0; j<ylen; j++)
    {
        int i=invmap[j];
        if(i<0) continue;
        {
            double d=di_vec[j];
            TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
            TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
            if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
            if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
            if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
            rmsd0+=d*d;
        }
    }
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);
}

void extract_split_substruct(
        const CoordArray& xa, const CoordArray& ya,
        const std::string &seqx, const std::string &seqy,
        const std::string &secx, const std::string &secy,
        const std::string &seqxA, const std::string &seqyA,
        const int xlen, const int ylen, const int n_ali8,
        AlignSubstruct sub[NUM_SPLIT_VARIANTS])
{
    AlignSubstruct &sub0 = sub[0];
    AlignSubstruct &sub1 = sub[1];
    sub0.seqx_hinge_trimmed.clear();
    sub0.seqy_hinge_trimmed.clear();
    sub0.secx_hinge_trimmed.resize(xlen + 1);
    sub0.secy_hinge_trimmed.resize(ylen + 1);
    sub0.xa_hinge_trimmed.resize(xlen);
    sub0.ya_hinge_trimmed.resize(ylen);
    sub0.xlen_hinge_trimmed = n_ali8;
    sub0.ylen_hinge_trimmed = ylen - n_ali8;
    sub1.seqx_hinge_trimmed.clear();
    sub1.seqy_hinge_trimmed.clear();
    sub1.secx_hinge_trimmed.resize(xlen + 1);
    sub1.secy_hinge_trimmed.resize(ylen + 1);
    sub1.xa_hinge_trimmed.resize(xlen);
    sub1.ya_hinge_trimmed.resize(ylen);
    sub1.xlen_hinge_trimmed = xlen - n_ali8;
    sub1.ylen_hinge_trimmed = n_ali8;
    int i = -1;
    int j = -1;
    int a0 = 0;
    int b0 = 0;
    int a1 = 0;
    int b1 = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        i += (seqxA[r] != '-');
        j += (seqyA[r] != '-');
        if (seqxA[r] != '-' && seqyA[r] != '-')
        {
            sub0.seqx_hinge_trimmed += seqx[i];
            sub0.secx_hinge_trimmed[a0] = secx[i];
            sub0.xa_hinge_trimmed[a0][0] = xa[i][0];
            sub0.xa_hinge_trimmed[a0][1] = xa[i][1];
            sub0.xa_hinge_trimmed[a0][2] = xa[i][2];
            a0++;
            sub1.seqy_hinge_trimmed += seqy[j];
            sub1.secy_hinge_trimmed[b1] = secy[j];
            sub1.ya_hinge_trimmed[b1][0] = ya[j][0];
            sub1.ya_hinge_trimmed[b1][1] = ya[j][1];
            sub1.ya_hinge_trimmed[b1][2] = ya[j][2];
            b1++;
        }
        else if (seqxA[r] != '-')
        {
            sub1.seqx_hinge_trimmed += seqx[i];
            sub1.secx_hinge_trimmed[a1] = secx[i];
            sub1.xa_hinge_trimmed[a1][0] = xa[i][0];
            sub1.xa_hinge_trimmed[a1][1] = xa[i][1];
            sub1.xa_hinge_trimmed[a1][2] = xa[i][2];
            a1++;
        }
        else if (seqyA[r] != '-')
        {
            sub0.seqy_hinge_trimmed += seqy[j];
            sub0.secy_hinge_trimmed[b0] = secy[j];
            sub0.ya_hinge_trimmed[b0][0] = ya[j][0];
            sub0.ya_hinge_trimmed[b0][1] = ya[j][1];
            sub0.ya_hinge_trimmed[b0][2] = ya[j][2];
            b0++;
        }
    }
}

void score_split_substruct(
        AlignSubstruct& sub,
        CoordArray& xa, CoordArray& ya,
        const std::string &seqx, const std::string &seqy,
        std::vector<double>& do_vec,
        const int xlen, const int ylen,
        const std::vector<std::string> &sequence,
        const ChainPairAlignOptions& align_opts,
        Vec3& t0, RotMat& u0,
        SubstructAlignRes& sub_align_res,
        std::ostream& os)
{
    ChainPairAlignResult chain_align_res = { 0};
    chain_align_res.d0_out = 5.0;
    chain_align_res.t0 = t0;
    chain_align_res.u0 = u0;
    chain_align_res.invmap.assign(ylen + 1, -1);
    chain_align_res.do_vec = do_vec;
    TMalign_main(sub.xa_hinge_trimmed, sub.ya_hinge_trimmed, sub.seqx_hinge_trimmed, sub.seqy_hinge_trimmed, sub.secx_hinge_trimmed, sub.secy_hinge_trimmed,
        chain_align_res, sub.xlen_hinge_trimmed, sub.ylen_hinge_trimmed, sequence, align_opts, os);
    t0 = chain_align_res.t0;
    u0 = chain_align_res.u0;
    do_vec = chain_align_res.do_vec;
    sub_align_res.t_u_result = chain_align_res.t0;
    sub_align_res.u_u_result = chain_align_res.u0;
    chain_align_res.TM1 = chain_align_res.TM2 = chain_align_res.TM3 = chain_align_res.TM4 = chain_align_res.TM5 = chain_align_res.rmsd0 = 0;
    chain_align_res.seqM = "";
    chain_align_res.seqxA = "";
    chain_align_res.seqyA = "";
    chain_align_res.n_ali = chain_align_res.n_ali8 = 0;
    CoordArray xt;
    xt.resize(xlen);
    score_with_cur_transform(
        xa, ya, xt, seqx, seqy, xlen, ylen,
        sequence, align_opts, 1, chain_align_res);
    do_vec = chain_align_res.do_vec;
    sub_align_res.TM1 = chain_align_res.TM1;
    sub_align_res.TM2 = chain_align_res.TM2;
    sub_align_res.TM3 = chain_align_res.TM3;
    sub_align_res.TM4 = chain_align_res.TM4;
    sub_align_res.TM5 = chain_align_res.TM5;
    sub_align_res.seqM = chain_align_res.seqM;
    sub_align_res.seqxA = chain_align_res.seqxA;
    sub_align_res.seqyA = chain_align_res.seqyA;
    sub_align_res.n_ali = chain_align_res.n_ali;
    sub_align_res.n_ali8 = chain_align_res.n_ali8;
    sub_align_res.rmsd0 = chain_align_res.rmsd0;
    sub_align_res.d0_0 = chain_align_res.d0_0;
    sub_align_res.TM_0 = chain_align_res.TM_0;
    sub_align_res.d0A = chain_align_res.d0A;
    sub_align_res.d0B = chain_align_res.d0B;
    sub_align_res.d0u = chain_align_res.d0u;
    sub_align_res.d0a = chain_align_res.d0a;
    sub_align_res.d0_out = chain_align_res.d0_out;
    sub_align_res.L_ali = chain_align_res.L_ali;
    sub_align_res.Liden = chain_align_res.Liden;
    sub_align_res.TM_ali = chain_align_res.TM_ali;
    sub_align_res.rmsd_ali = chain_align_res.rmsd_ali;
    sub_align_res.invmap = chain_align_res.invmap;
}

void fill_sub_struct_align_result(
        const SubstructAlignRes sub_align_res[NUM_SPLIT_VARIANTS],
        DoubleMatrix& tu_vec,
        ChainPairAlignResult& chain_align_res)
{
    double TM_h=(sub_align_res[0].TM1>sub_align_res[0].TM2)?sub_align_res[0].TM1:sub_align_res[0].TM2;
    double TM  =(sub_align_res[1].TM1>sub_align_res[1].TM2)?sub_align_res[1].TM1:sub_align_res[1].TM2;
    const SubstructAlignRes* best_align;
    if (TM_h>TM)
        best_align=&sub_align_res[0];
    else
        best_align=&sub_align_res[1];
    chain_align_res.t0 = best_align->t_u_result;
    chain_align_res.u0 = best_align->u_u_result;
    t_u2tu(chain_align_res.t0, chain_align_res.u0, tu_vec[0]);
    chain_align_res.TM1=best_align->TM1;
    chain_align_res.TM2=best_align->TM2;
    chain_align_res.TM3=best_align->TM3;
    chain_align_res.TM4=best_align->TM4;
    chain_align_res.TM5=best_align->TM5;
    chain_align_res.seqM=best_align->seqM;
    chain_align_res.seqxA=best_align->seqxA;
    chain_align_res.seqyA=best_align->seqyA;
    chain_align_res.rmsd0=best_align->rmsd0;
    chain_align_res.n_ali=best_align->n_ali;
    chain_align_res.n_ali8=best_align->n_ali8;
    chain_align_res.L_ali=best_align->L_ali;
    chain_align_res.Liden=best_align->Liden;
    chain_align_res.TM_ali=best_align->TM_ali;
    chain_align_res.rmsd_ali=best_align->rmsd_ali;
    chain_align_res.d0_0=best_align->d0_0;
    chain_align_res.TM_0=best_align->TM_0;
    chain_align_res.d0A=best_align->d0A;
    chain_align_res.d0B=best_align->d0B;
    chain_align_res.d0u=best_align->d0u;
    chain_align_res.d0a=best_align->d0a;
    chain_align_res.d0_out=best_align->d0_out;
    chain_align_res.invmap = best_align->invmap;
}

void extract_unalign_subcoords(
        const std::string& seqxA, const std::string& seqyA,
        const std::string& seqx, const std::string& seqy,
        const std::string& secx, const std::string& secy,
        const CoordArray& xa, const CoordArray& ya,
        const int xlen_hinge_trimmed, const int ylen_hinge_trimmed,
        CoordArray& xa_hinge_trimmed, CoordArray& ya_hinge_trimmed,
        std::string& seqx_hinge_trimmed, std::string& seqy_hinge_trimmed,
        std::string& secx_hinge_trimmed, std::string& secy_hinge_trimmed,
        std::vector<int>& r1toi, std::vector<int>& r2toj)
{
    xa_hinge_trimmed.resize(xlen_hinge_trimmed);
    ya_hinge_trimmed.resize(ylen_hinge_trimmed);
    secx_hinge_trimmed.resize(xlen_hinge_trimmed + 1);
    secy_hinge_trimmed.resize(ylen_hinge_trimmed + 1);
    r1toi.assign(xlen_hinge_trimmed, 0);
    r2toj.assign(ylen_hinge_trimmed, 0);
    int i = -1, j = -1, r1 = 0, r2 = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        i += (seqxA[r] != '-');
        j += (seqyA[r] != '-');
        if (seqyA[r] == '-')
        {
            seqx_hinge_trimmed.push_back(seqx[i]);
            secx_hinge_trimmed[r1] = secx[i];
            xa_hinge_trimmed[r1][0] = xa[i][0];
            xa_hinge_trimmed[r1][1] = xa[i][1];
            xa_hinge_trimmed[r1][2] = xa[i][2];
            r1toi[r1] = i;
            r1++;
        }
        if (seqxA[r] == '-')
        {
            seqy_hinge_trimmed.push_back(seqy[j]);
            secy_hinge_trimmed[r2] = secy[j];
            ya_hinge_trimmed[r2][0] = ya[j][0];
            ya_hinge_trimmed[r2][1] = ya[j][1];
            ya_hinge_trimmed[r2][2] = ya[j][2];
            r2toj[r2] = j;
            r2++;
        }
    }
}

void search_hinge_regions(
        CoordArray& xa, CoordArray& ya, CoordArray& xt, DoubleMatrix& tu_vec,
        const std::string &seqx, const std::string &seqy, const std::string &secx,
        const std::string &secy,
        const int xlen, const int ylen, const int hinge_opt, const std::vector<std::string> &sequence,
        const ChainPairAlignOptions& align_opts,
        ChainPairAlignResult& chain_align_res,
        std::ostream& os)
{
    for (int r2=0;r2<chain_align_res.seqM.size();r2++)
    {
        if (chain_align_res.seqM[r2]=='1')
            chain_align_res.seqM[r2]='0';
    }

    int minlen = std::min(xlen, ylen);
    for (int hinge=0;hinge<hinge_opt;hinge++)
    {
        if (minlen-chain_align_res.n_ali8<5) break;

        int xlen_hinge_trimmed=xlen - chain_align_res.n_ali8;
        int ylen_hinge_trimmed=ylen - chain_align_res.n_ali8;
        CoordArray xa_hinge_trimmed;
        CoordArray ya_hinge_trimmed;
        std::string seqx_hinge_trimmed, seqy_hinge_trimmed, secx_hinge_trimmed, secy_hinge_trimmed;
        std::vector<int> r1toi, r2toj;

        extract_unalign_subcoords(chain_align_res.seqxA, chain_align_res.seqyA, seqx, seqy, secx, secy, xa, ya,
            xlen_hinge_trimmed, ylen_hinge_trimmed, xa_hinge_trimmed, ya_hinge_trimmed, seqx_hinge_trimmed, seqy_hinge_trimmed, secx_hinge_trimmed, secy_hinge_trimmed, r1toi, r2toj);

        CurHingeResult cur_hing_res;
        cur_hing_res.tmalign_res.t0 = chain_align_res.t0;
        cur_hing_res.tmalign_res.u0 = chain_align_res.u0;
        cur_hing_res.tmalign_res.do_vec = chain_align_res.do_vec;
        TMalign_main(xa_hinge_trimmed, ya_hinge_trimmed, seqx_hinge_trimmed, seqy_hinge_trimmed,
            secx_hinge_trimmed, secy_hinge_trimmed, cur_hing_res.tmalign_res,
            xlen_hinge_trimmed, ylen_hinge_trimmed, sequence, align_opts, os);
        chain_align_res.t0 = cur_hing_res.tmalign_res.t0;
        chain_align_res.u0 = cur_hing_res.tmalign_res.u0;
        chain_align_res.do_vec = cur_hing_res.tmalign_res.do_vec;

        do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
        fill_hing_res(cur_hing_res, chain_align_res.TM1, chain_align_res.TM2, chain_align_res.TM3, chain_align_res.TM4, chain_align_res.TM5, chain_align_res.rmsd0,
            chain_align_res.seqM, chain_align_res.seqxA, chain_align_res.seqyA, chain_align_res.n_ali, chain_align_res.n_ali8, chain_align_res.invmap, ylen);

        cur_hing_res.tmalign_res.d0_0 = chain_align_res.d0_0;
        cur_hing_res.tmalign_res.TM_0 = chain_align_res.TM_0;
        cur_hing_res.tmalign_res.d0A = chain_align_res.d0A;
        cur_hing_res.tmalign_res.d0B = chain_align_res.d0B;
        cur_hing_res.tmalign_res.d0u = chain_align_res.d0u;
        cur_hing_res.tmalign_res.d0a = chain_align_res.d0a;
        cur_hing_res.tmalign_res.d0_out = chain_align_res.d0_out;
        cur_hing_res.tmalign_res.do_vec = chain_align_res.do_vec;
        cur_hing_res.tmalign_res.L_ali = chain_align_res.L_ali;
        cur_hing_res.tmalign_res.Liden = chain_align_res.Liden;
        cur_hing_res.tmalign_res.TM_ali = chain_align_res.TM_ali;
        cur_hing_res.tmalign_res.rmsd_ali = chain_align_res.rmsd_ali;
        se_main(xt, ya, seqx, seqy,
            cur_hing_res.tmalign_res,
            xlen, ylen, sequence, align_opts, 0, hinge+1);
        chain_align_res.d0A = cur_hing_res.tmalign_res.d0A;
        chain_align_res.d0B = cur_hing_res.tmalign_res.d0B;
        chain_align_res.d0u = cur_hing_res.tmalign_res.d0u;
        chain_align_res.d0a = cur_hing_res.tmalign_res.d0a;
        chain_align_res.do_vec = cur_hing_res.tmalign_res.do_vec;
        chain_align_res.Liden = cur_hing_res.tmalign_res.Liden;

        int new_ali=0;
        for (int r=0;r<cur_hing_res.tmalign_res.seqM.size();r++)
        {
           new_ali+=(cur_hing_res.tmalign_res.seqM[r]==hinge+'1');
        }
        if (cur_hing_res.tmalign_res.n_ali8 - chain_align_res.n_ali8<5)
        {
           new_ali=0;
        }

        if (new_ali>=5)
        {
            update_global_res(chain_align_res, cur_hing_res, tu_vec);
        }

        if (new_ali<5) break;
    }
}

void run_flexalign(
    const ParsedChain& chain1_data,
    const ParsedChain& chain2_data,
    const AlignCommonInput& common_inputs,
    const FlexalignParams& flex_params,
    FlexAlignResult& flexalign_res,
    std::ostream& os)
{
    CoordArray xa = chain1_data.chain_coords;
    CoordArray ya = chain2_data.chain_coords;
    std::string seqx = chain1_data.chain_seq;
    std::string seqy = chain2_data.chain_seq;
    std::string secx = chain1_data.chain_sec;
    std::string secy = chain2_data.chain_sec;
    int xa_len = chain1_data.chain_len;
    int ya_len = chain2_data.chain_len;
    int mol_type = chain1_data.cur_complex_mol_list + chain2_data.cur_complex_mol_list;

    FlexAlignResult flexalign_main_res;
    align_with_flexalign_main(
        xa, ya, seqx, seqy, secx, secy, xa_len, ya_len,
        common_inputs, mol_type, flex_params.hinge_opt, flexalign_main_res, os);
    double best_global_max_TM = (flexalign_main_res.TM1 > flexalign_main_res.TM2) ? flexalign_main_res.TM1 : flexalign_main_res.TM2; // Best max(TM1,TM2) threshold (shared for min_resid_num pruning)

    switch (flex_params.mode)
    {
        case FLEX_USBCAT:
            if (best_global_max_TM >= flex_params.TMpass)
            {
                flexalign_res = flexalign_main_res;
                flexalign_res.hingeNum = flexalign_res.tu_vec.size();
            }
            else
            {
                flexalign_usbcat_main(
                    xa, ya, seqx, seqy, secx, secy, xa_len, ya_len, common_inputs,
                    flexalign_main_res, mol_type, best_global_max_TM, flex_params.sparse_val,
                    flex_params, flexalign_res, os);
            }
            break;
        case FLEX_BEST:
            flexalign_res = flexalign_main_res;
            break;
        default:
            PrintErrorAndQuit("ERROR! Unknown FlexAlignMode.");
            break;
    }

    return;
}

std::vector<USBCAT_AFP> extract_initial_afps(
    const CoordArray& xa,
    const CoordArray& ya,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double min_resid_num,
    const double cur_rmsdCut,
    const double cur_badRmsd,
    const int sparse_val)
{

        std::vector<USBCAT_AFP> initial_afps;
        int step = sparse_val + 1;

        
        CoordArray r1(usb_cat_para.fragLen), r2(usb_cat_para.fragLen);

        for (int i = 0; i <= xlen - usb_cat_para.fragLen; i += step)
        {
            for (int j = 0; j <= ylen - usb_cat_para.fragLen; j += step)
            {
                int cur_resid_num = cur_diag_resid_num(i, j, xlen, ylen, usb_cat_para.fragLen);
                if (cur_resid_num < min_resid_num)
                    continue;

                double dist1 = tables.chain1_dist_table[i][usb_cat_para.fragLen - 1];
                double dist2 = tables.chain2_dist_table[j][usb_cat_para.fragLen - 1];

                if (std::fabs(dist1 - dist2) > 2.0 * cur_rmsdCut)
                    continue;

                for (int k = 0; k < usb_cat_para.fragLen; k++)
                {
                    r1[k][0] = xa[i + k][0];
                    r1[k][1] = xa[i + k][1];
                    r1[k][2] = xa[i + k][2];
                    r2[k][0] = ya[j + k][0];
                    r2[k][1] = ya[j + k][1];
                    r2[k][2] = ya[j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
            
                Kabsch(r1, r2, usb_cat_para.fragLen, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / usb_cat_para.fragLen);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    USBCAT_AFP afp;
                    afp.i = i;
                    afp.j = j;
                    afp.len = usb_cat_para.fragLen;
                    afp.score = afp_score(usb_cat_para.resScore, usb_cat_para.fragLen, rmsd_tmp, cur_badRmsd);
                    initial_afps.push_back(afp);
                }
            }
        }
    return initial_afps;
}

std::vector<USBCAT_AFP> merge_cur_afps_group(
    const std::vector<USBCAT_AFP>& cur_group,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_rmsdCut,
    const double cur_badRmsd)
{
    std::vector<USBCAT_AFP> merged_afps;
    int afp_num = cur_group.size();
    std::vector<bool> merged_status(afp_num, false);

    int max_merge_len = std::min(xlen, ylen);
    CoordArray r1_merge(max_merge_len), r2_merge(max_merge_len);

    for (int idx = 0; idx < afp_num; idx++)
    {
        if (merged_status[idx]) continue;

        USBCAT_AFP curr_afp = cur_group[idx];
        for (int nxt_idx = idx + 1; nxt_idx < afp_num; nxt_idx++)
        {
            USBCAT_AFP nxt_afp = cur_group[nxt_idx];
            if (nxt_afp.i > curr_afp.i + curr_afp.len)
                break;  // Sorted by i, subsequent nxt entries cannot overlap

            if (nxt_afp.i + nxt_afp.len > curr_afp.i + curr_afp.len)
            {
                int new_len = (nxt_afp.i + nxt_afp.len) - curr_afp.i;
                for (int k = 0; k < new_len; k++)
                {
                    r1_merge[k][0] = xa[curr_afp.i + k][0];
                    r1_merge[k][1] = xa[curr_afp.i + k][1];
                    r1_merge[k][2] = xa[curr_afp.i + k][2];
                    r2_merge[k][0] = ya[curr_afp.j + k][0];
                    r2_merge[k][1] = ya[curr_afp.j + k][1];
                    r2_merge[k][2] = ya[curr_afp.j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
                Kabsch(r1_merge, r2_merge, new_len, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / new_len);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    curr_afp.len = new_len;
                    curr_afp.score = afp_score(usb_cat_para.resScore, new_len, rmsd_tmp, cur_badRmsd);
                    merged_status[nxt_idx] = true;
                }
            }
        }
        merged_afps.push_back(curr_afp);
    }
    return merged_afps;
}

void find_candidate_prevs(const USBCAT_AFP& cur_afp,
                                 const std::vector<int>& afp_aft_index,
                                 const std::vector<int>& afp_bef_index,
                                 const int xlen, const int ylen,
                                 const USBCATParams& usb_cat_para,
                                 std::vector<int>& candidate_prevs)
{
    int x_right = cur_afp.i - usb_cat_para.fragLen;
    int x_left2 = std::max(0, x_right - usb_cat_para.misCut);
    int x_left1 = std::max(0, cur_afp.i - usb_cat_para.maxGapFrag);

    int y_right = cur_afp.j - usb_cat_para.fragLen;
    int y_left2 = std::max(0, y_right - usb_cat_para.misCut);
    int y_left1 = std::max(0, cur_afp.j - usb_cat_para.maxGapFrag);

    for (int st = 0; st < 2; st++)
    {
        int x_left, x_right_b, y_left, y_right_b;
        if (st == 0)
        {
            x_left = std::max(x_left1, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left2, 0);
            y_right_b = std::min(y_right, ylen - 1);
        }
        else
        {
            x_left = std::max(x_left2, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left1, 0);
            y_right_b = std::min(y_left2 - 1, ylen - 1);
        }

        //Check whether the window overlaps the valid range [0, ylen-1]
        if (y_left >= ylen || y_right_b < 0) continue;

        for (int prev_i = x_left; prev_i <= x_right_b; prev_i++)
        {
            int s1 = afp_aft_index[prev_i * ylen + y_left];
            int s2 = afp_bef_index[prev_i * ylen + y_right_b];
            if (s1 != -1 && s2 != -1 && s1 <= s2)
                for (int s = s1; s <= s2; s++)
                    candidate_prevs.push_back(s);
        }
    }
}

void split_candidate_blocks(
    std::vector<AFPBlock>& candidate_blocks,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const double cur_local_badRmsd)
{
    const size_t max_blocks = (size_t)(usb_cat_para.max_hinge_num + 1);
    bool can_split_more = candidate_blocks.size() < max_blocks; 
    while (can_split_more)
    {
        bool splitted = false; 
        double max_rmsd = 0.0;
        int target_b = -1;
    
        for (size_t block_idx = 0; block_idx < candidate_blocks.size(); block_idx++)
        {
            if (candidate_blocks[block_idx].afps.size() > 2)
            {
                double cur_rmsd = calc_block_rmsd(candidate_blocks[block_idx].afps, xa, ya);
                if (cur_rmsd > max_rmsd)
                {
                    max_rmsd = cur_rmsd;
                    target_b = block_idx;
                }
            }
        }
        
        if (max_rmsd >= cur_local_badRmsd && target_b != -1)
        {
            AFPBlock& target_blk = candidate_blocks[target_b]; 
            double max_drmsd = 0;
            int cut_afp_idx = 0;
            for (size_t idx = 1; idx < target_blk.afps.size(); idx++)
            {
                if (target_blk.drmsd[idx] > max_drmsd)
                {
                    max_drmsd = target_blk.drmsd[idx];
                    cut_afp_idx = idx;
                }
            }
            
            if (cut_afp_idx > 0)
            {
                AFPBlock right_blk;
                right_blk.afps.assign(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                right_blk.drmsd.assign(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                right_blk.drmsd[0] = 0.0;
                target_blk.afps.erase(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                target_blk.drmsd.erase(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                candidate_blocks.insert(candidate_blocks.begin() + target_b + 1, right_blk);
                splitted = true;
            }
        }
        can_split_more = splitted && (candidate_blocks.size() < max_blocks);
    }
}

void build_global_align_result(
    std::vector<RegionAlignResult>& region_align_res,
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const int region_num,
    const std::string& seqx,
    const std::string& seqy,
    GlobalAlignResult& global_align_res)
{
    for (int region_idx = 0; region_idx < region_num; region_idx++)
    {
        int L1_sub = chain1_bounds[region_idx + 1] - chain1_bounds[region_idx];
        int L2_sub = chain2_bounds[region_idx + 1] - chain2_bounds[region_idx];

        if (!region_align_res[region_idx].valid)
        {
            // Fill gaps if region was invalid or bypassed
            for (int i = 0; i < L1_sub; i++)
            {
                global_align_res.seqxA += seqx[chain1_bounds[region_idx] + i];
                global_align_res.seqyA += '-';
                global_align_res.seqM += ' ';
            }
            for (int i = 0; i < L2_sub; i++)
            {
                global_align_res.seqxA += '-';
                global_align_res.seqyA += seqy[chain2_bounds[region_idx] + i];
                global_align_res.seqM += ' ';
            }
            continue;
        }

        FlexAlignResult& flexalign_res = region_align_res[region_idx].align_re;
        if (flexalign_res.tu_vec.empty())
        {
            std::vector<double> tu_tmp(12);
            t_u2tu(flexalign_res.t0, flexalign_res.u0, tu_tmp);
            flexalign_res.tu_vec.push_back(tu_tmp);
        }

        int base_tu_idx = (int)global_align_res.tu_vec.size();
        for (size_t m = 0; m < flexalign_res.tu_vec.size(); m++)
            global_align_res.tu_vec.push_back(flexalign_res.tu_vec[m]);

        int rx = chain1_bounds[region_idx];
        int current_global_idx = base_tu_idx;

        for (size_t i = 0; i < flexalign_res.seqxA.length(); i++)
        {
            char c = flexalign_res.seqM[i];
            if (c != ' ' && c != '.' && c != ':')
            {
                int local_hinge_idx = -1;
                if (c >= '0' && c <= '9')
                    local_hinge_idx = c - '0';
                else if (c >= 'a' && c <= 'z')
                    local_hinge_idx = c - 'a' + 10;
                else if (c >= 'A' && c <= 'Z')
                    local_hinge_idx = c - 'A' + 36;

                if (local_hinge_idx >= 0 && local_hinge_idx < (int)flexalign_res.tu_vec.size())
                    current_global_idx = base_tu_idx + local_hinge_idx;
            }

            if (flexalign_res.seqxA[i] != '-')
            {
                global_align_res.res_tu[rx] = current_global_idx;
                rx++;
            }

            if (flexalign_res.seqxA[i] != '-' && flexalign_res.seqyA[i] != '-')
            {
                if (c != ' ' && c != '.' && c != ':')
                {
                    char global_c;
                    if (current_global_idx < 10)
                        global_c = '0' + current_global_idx;
                    else if (current_global_idx < 36)
                        global_c = 'a' + (current_global_idx - 10);
                    else if (current_global_idx < 62)
                        global_c = 'A' + (current_global_idx - 36);
                    else
                        global_c = '*';
                    flexalign_res.seqM[i] = global_c;
                }
                else
                    flexalign_res.seqM[i] = c;
            }
            else
                flexalign_res.seqM[i] = ' ';
        }

        global_align_res.seqM += flexalign_res.seqM;
        global_align_res.seqxA += flexalign_res.seqxA;
        global_align_res.seqyA += flexalign_res.seqyA;
    }
}

void recompute_global_metrics(
    const GlobalAlignResult& sd,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const int xlen,
    const int ylen,
    const AlignCommonInput& common_inputs,
    const int mol_type,
    const double d0_out,
    FlexAlignResult& flexalign_res)
{
    const UserOptions& opts = common_inputs.user_options;

    double cur_d0A, cur_d0B, cur_d0a, cur_d0u;
    calc_cur_d0(xlen, ylen, opts, mol_type, cur_d0A, cur_d0B, cur_d0a, cur_d0u);
    reset_flexalign_result(flexalign_res, sd.seqxA.length());

    int i_res = 0, j_res = 0;
    int cached_matrix_idx = -1;
    Vec3 t_k;
    RotMat u_k;
    for (size_t r = 0; r < sd.seqxA.length(); r++)
    {
        bool x_valid = (sd.seqxA[r] != '-');
        bool y_valid = (sd.seqyA[r] != '-');

        double d = -1.0;
        if (x_valid && y_valid)
        {
            int matrix_idx = sd.res_tu[i_res];
            if (matrix_idx >= 0 && matrix_idx < (int)sd.tu_vec.size())
            {
                if (matrix_idx != cached_matrix_idx)
                {
                    tu2t_u(sd.tu_vec[matrix_idx], t_k, u_k);
                    cached_matrix_idx = matrix_idx;
                }

                Vec3 x_rot;
                transform(t_k, u_k, xa[i_res], x_rot);
                double dist2 = dist(x_rot, ya[j_res]);
                d = std::sqrt(dist2);

                flexalign_res.TM2 += 1.0 / (1.0 + dist2 / (cur_d0B * cur_d0B));
                flexalign_res.TM1 += 1.0 / (1.0 + dist2 / (cur_d0A * cur_d0A));
                if (opts.a_opt)
                    flexalign_res.TM3 += 1.0 / (1.0 + dist2 / (cur_d0a * cur_d0a));
                if (opts.u_opt)
                    flexalign_res.TM4 += 1.0 / (1.0 + dist2 / (cur_d0u * cur_d0u));
                if (opts.d_opt)
                    flexalign_res.TM5 += 1.0 / (1.0 + dist2 / (opts.d0_scale * opts.d0_scale));

                flexalign_res.n_ali++;

                if (d <= d0_out)
                {
                    flexalign_res.rmsd0 += dist2;
                    flexalign_res.n_ali8++;
                    if (seqx[i_res] == seqy[j_res])
                        flexalign_res.Liden += 1.0;
                }
            }
        }
        flexalign_res.do_vec.push_back(d);

        if (x_valid)
            i_res++;
        if (y_valid)
            j_res++;
    }

    flexalign_res.TM2 /= xlen;
    flexalign_res.TM1 /= ylen;
    if (opts.a_opt)
        flexalign_res.TM3 /= (xlen + ylen) * 0.5;
    if (opts.u_opt)
        flexalign_res.TM4 /= opts.Lnorm_ass;
    if (opts.d_opt)
        flexalign_res.TM5 /= ylen;
    if (flexalign_res.n_ali8 > 0)
        flexalign_res.rmsd0 = std::sqrt(flexalign_res.rmsd0 / flexalign_res.n_ali8);
    else
        flexalign_res.rmsd0 = 0.0;

    flexalign_res.seqM = sd.seqM;
    flexalign_res.seqxA = sd.seqxA;
    flexalign_res.seqyA = sd.seqyA;
    flexalign_res.tu_vec = sd.tu_vec;
    flexalign_res.d0A = cur_d0A;
    flexalign_res.d0B = cur_d0B;
    flexalign_res.d0a = cur_d0a;
    flexalign_res.d0u = cur_d0u;
    if (!flexalign_res.tu_vec.empty())
        tu2t_u(flexalign_res.tu_vec[0], flexalign_res.t0, flexalign_res.u0);

    flexalign_res.TM_ali = flexalign_res.TM1;
    flexalign_res.rmsd_ali = flexalign_res.rmsd0;
    flexalign_res.L_ali = flexalign_res.n_ali;
}

