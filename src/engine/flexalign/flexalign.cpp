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

    extract_aligned_resi(rasmol_state, seqM, seqxA, seqyA, resi_vec1, resi_vec2, mm_opt);

    read_rasmol_structure1(rasmol_state, xname, tu_vec, mm_opt, ter_opt, mirror_opt);

    read_rasmol_structure2(rasmol_state, yname, mm_opt, ter_opt);

    write_rasmol_to_file(rasmol_state, fname_super, mm_opt);
}

