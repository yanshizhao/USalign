#include "flexalign.h"

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