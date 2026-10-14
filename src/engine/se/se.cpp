#include "se.h"

int se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt, const int hinge)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqyA = chain_align_res.seqyA;
    std::vector<int> &invmap = chain_align_res.invmap;
    const int mol_type = align_opts.mol_type;
    SeSearchState search;
    CharMatrix path;          // for dynamic programming (char: 1/0)
    DoubleMatrix val;          // for dynamic programming

    std::vector<int> m1;
    std::vector<int> m2;
    if (outfmt_opt<2)
    {
        m1.resize(xlen); //alignd index in x
        m2.resize(ylen); //alignd index in y
    }

    /***********************/
    // allocate memory
    /***********************/
    path.assign(xlen+1, std::vector<char>(ylen+1));
    val.assign(xlen+1, std::vector<double>(ylen+1));
    std::vector<int> invmap0(ylen+1);
    int j;
    if (hinge==0) for (j=0;j<=ylen;j++) invmap0[j]=-1;
    else for (j=0;j<ylen;j++) invmap0[j]=invmap[j];
    std::vector<char> seqM_char;
    if (hinge)
    {
        seqM_char.assign(ylen,hinge+'0');
        j=-1;
        for (int r=0;r<seqM.size();r++)
        {
            j+=seqyA[r]!='-';
            if (seqM[r]!=' ') seqM_char[j]=seqM[r];
        }
    }

    se_prepare_search_parameters(xlen, ylen, align_opts, mol_type, chain_align_res, search);
    se_build_initial_alignment(path, val, xa, ya, xlen, ylen, sequence, align_opts, search, chain_align_res, hinge);
    se_score_alignment(xa, ya, xlen, ylen, align_opts, search, invmap0, chain_align_res, m1, m2, outfmt_opt, hinge);
    if (outfmt_opt>=2)
    {
        return 0;
    }

    se_extract_alignment_strings(xa, ya, seqx, seqy, xlen, ylen, m1, m2, seqM_char, chain_align_res, hinge);
    return 0; // zero for no exception
}
