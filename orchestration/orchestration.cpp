#include "orchestration.h"

#include "se.h"

// ---- Single chain-pair structure alignment dispatcher (cp / se / TMalign_main) ----
void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (align_opts.cp_opt) CPalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
    else if (align_opts.se_opt)
    {
        chain_align_res.invmap.assign(ylen + 1, 0);
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);
        se_main(xa, ya, seqx, seqy,
            chain_align_res,
            xlen, ylen, sequence, align_opts, outfmt_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0; r2<ylen; r2++)
            {
                r1 = chain_align_res.invmap[r2];
                if (r1 < 0) continue;
                chain_align_res.L_ali += 1;
                chain_align_res.Liden += (seqx[r1] == seqy[r2]);
            }
        }
    }
    else TMalign_main(
        xa, ya, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
}
