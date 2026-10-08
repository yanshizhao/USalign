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
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &flexalign_res){
    flexalign_res.hingeNum = flexalign_main(
        xa, ya, seqx, seqy, secx, secy, xlen, ylen, common_inputs, mol_type, hinge_opt, ss_opt, flexalign_res);

    bool refine_fallback = (hinge_opt) && (flexalign_res.hingeNum <= 1) && (flexalign_res.n_ali8 < 0.6 * getmin(xlen, ylen));
    if (refine_fallback)
    {
        FlexAlignResult refine_res;
        refine_res.tu_vec.push_back(flexalign_res.tu_vec[0]);
        tu2t_u(flexalign_res.tu_vec[0], refine_res.t0, refine_res.u0);

        refine_res.hingeNum = flexalign_main(
            xa, ya, seqx, seqy, secx, secy, xlen, ylen, common_inputs, mol_type, hinge_opt, ss_opt, refine_res);

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
    FlexAlignResult& flexalign_res){
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
        global_best_align, best_global_max_TM);

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

int Flexalign(AlignCommonInput& common_inputs, const FlexalignParams& flex_params, FlexAlignResult& flex_result)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    std::vector<std::vector<std::string> > PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> > PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    int read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;  // whether to read residue index

    const PdbParseOptions pdb_parse_opts1{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1};
    const PdbParseOptions pdb_parse_opts2{
        user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2};

    // loop over file names
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1, user_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            ParsedChain chain1_data;
            if (!parse_chain(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1,
                    chain_i, user_opts.mol_opt, user_opts.mirror_opt, read_resi,
                    user_opts.atom_opt, chain1_data))
                continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j = j_start; j < parsed_input.chain2_list.size(); j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size() == 0)
                {
                    user_opts.yname = parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2, user_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ParsedChain chain2_data;
                    if (!parse_chain(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2,
                            chain_j, user_opts.mol_opt, 0, read_resi,
                            user_opts.atom_opt, chain2_data))
                        continue;

                    if (user_opts.byresi_opt)
                        extract_aln_from_resi(parsed_input.sequence, chain1_data.chain_seq, chain2_data.chain_seq,
                            chain1_data.resi_vec, chain2_data.resi_vec, user_opts.byresi_opt);

                    flex_result = FlexAlignResult();
                    run_flexalign(chain1_data, chain2_data, common_inputs, flex_params, flex_result);

                    if (user_opts.outfmt_opt==0) print_version(USALIGN_VERSION);
                    output_flexalign_results(chain1_data, chain2_data, flex_result, user_opts, false);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j

            PDB_lines1[chain_i].clear();
        } // chain_i

        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i

    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
    }
    return 0;
}

