#include "SOIalign.h"

void soi_egs(DoubleMatrix& score, const int xlen, const int ylen, std::vector<int>& invmap,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const int mm_opt){
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
    while(1)
    {
        max_score=0;
        maxi=maxj=-1;
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

// sequence order independent alignment
int SOIalign(AlignCommonInput& common_inputs, const SoiAlignParams& soi_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;

    // declare previously global variables
    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    xlen, ylen;         // chain length
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    std::string secx;                // for the secondary structure
    std::string secy;
    std::string seqx, seqy;         // for the protein sequence
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;
                               // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure
                               // --> superpose xa onto ya
    SoiAlignContext soi_context;
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

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
            xlen=parse_chain_for_soi(PDB_lines1[chain_i], user_opts.xname,
                mol_vec1[chain_i], user_opts, soi_params, user_opts.mirror_opt,
                xa, seqx, secx, resi_vec1, soi_context.k_nearest1, soi_context.chain1_bonds);
            if (xlen<3) continue;

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j=j_start;j<parsed_input.chain2_list.size();j++)
            {
                if (user_opts.dirpair_opt.size() && i!=j) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    user_opts.yname=parsed_input.chain2_list[j];
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
                    ylen=parse_chain_for_soi(PDB_lines2[chain_j], user_opts.yname,
                        mol_vec2[chain_j], user_opts, soi_params, false,
                        ya, seqy, secy, resi_vec2, soi_context.k_nearest2, soi_context.chain2_bonds);
                    if (ylen<3) continue;

                    ChainPairAlignResult chain_align_res = { 0};
                    chain_align_res.d0_out = 5.0;
                    ChainPairAlignOptions soi_align_opts = soi_params.align;
                    soi_align_opts.fast_opt = (getmin(xlen,ylen)>1500)?true:user_opts.fast_opt;
                    soi_align_opts.mol_type = mol_vec1[chain_i]+mol_vec2[chain_j];
                    soi_context.pair_distances.assign(ylen+1, 0.0);

                    // entry function for structure alignment
                    soi_align_pair(chain_align_res, xa, ya, soi_context.k_nearest1, soi_context.k_nearest2, soi_params.closeK_opt,
                        seqx, seqy, secx, secy, xlen, ylen, soi_align_opts,
                        parsed_input.sequence, soi_context.pair_distances, soi_context.chain1_bonds, soi_context.chain2_bonds,
                        soi_params.mm_opt, user_opts.outfmt_opt);

                    // print chain_align_res
                    output_pair_alignment(common_inputs,
                        xlen, ylen, chain_align_res,
                        PDB_lines1[chain_i], PDB_lines2[chain_j],
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        resi_vec1, resi_vec2, &soi_context.pair_distances);

                    resi_vec2.clear();
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
            resi_vec1.clear();
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

