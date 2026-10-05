#pragma once
#ifndef SOIalign_h
#define SOIalign_h 1

#include "TMalign.h"

inline void assign_sec_bond(IntPairArray& secx_bond, const std::string& secx, const int xlen)
{
    int i,j,starti=-1,endi=-1;
    char ss,prev_ss=0;
    for (i=0; i<xlen; i++) { ss=secx[i]; secx_bond[i][0]=secx_bond[i][1]=-1;
        if (ss!=prev_ss && !(ss=='C' && prev_ss=='T') && !(ss=='T' && prev_ss=='C')) {
            if (starti>=0) { endi=i; for (j=starti;j<endi;j++) { secx_bond[j][0]=starti; secx_bond[j][1]=endi; } }
            if (ss=='H' || ss=='E' || ss=='<' || ss=='>') starti=i; else starti=-1;
        } prev_ss=secx[i];
    }
    if (starti>=0) { endi=i; for (j=starti;j<endi;j++) { secx_bond[j][0]=starti; secx_bond[j][1]=endi; } }
    for (i=0;i<xlen;i++) if (secx_bond[i][1]-secx_bond[i][0]==1) secx_bond[i][0]=secx_bond[i][1]=-1;
}

inline void getCloseK(const CoordArray& xa, const int xlen, const int closeK_opt, CoordArray& xk)
{
    DoubleMatrix score;
    score.assign(xlen+1, std::vector<double>(xlen+1, 0));
    std::vector<std::pair<double,int> > close_idx_vec(xlen, std::make_pair(0,0));
    int i,j,k;
    for(i=0;i<xlen;i++) {
        score[i+1][i+1]=0;
        for(j=i+1;j<xlen;j++) score[j+1][i+1]=score[i+1][j+1]=dist(xa[i], xa[j]);
    }
    for(i=0;i<xlen;i++) {
        for(j=0;j<xlen;j++) { close_idx_vec[j].first=score[i+1][j+1]; close_idx_vec[j].second=j; }
        std::sort(close_idx_vec.begin(), close_idx_vec.end());
        for(k=0;k<closeK_opt;k++) {
            j=close_idx_vec[k % xlen].second;
            xk[i*closeK_opt+k][0]=xa[j][0]; xk[i*closeK_opt+k][1]=xa[j][1]; xk[i*closeK_opt+k][2]=xa[j][2];
        }
    }
}

// check if pairing i to j conform to sequantiality within the SSE



inline bool sec2sq(const int i, const int j,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const std::vector<int>& fwdmap, const std::vector<int>& invmap)
{
    if (i<0 || j<0) return true;
    int ii,jj;
    if (secx_bond[i][0]>=0) {
        for (ii=secx_bond[i][0];ii<secx_bond[i][1];ii++) {
            jj=fwdmap[ii]; if (jj>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    if (secy_bond[j][0]>=0) {
        for (jj=secy_bond[j][0];jj<secy_bond[j][1];jj++) {
            ii=invmap[jj]; if (ii>=0 && (i-ii)*(j-jj)<=0) return false;
        }
    }
    return true;
}

void soi_egs(DoubleMatrix& score, const int xlen, const int ylen, std::vector<int>& invmap,
    const IntPairArray& secx_bond, const IntPairArray& secy_bond, const int mm_opt)
{
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


/* entry function for se
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * */


inline int soi_se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt);



inline void SOI_super2score(const CoordArray& xt, const CoordArray& ya, const int xlen,
    const int ylen, DoubleMatrix& score, double d0, double score_d8)
{
    int i,j; double d02=d0*d0,score_d82=score_d8*score_d8,d2;
    for(i=0;i<xlen;i++) for(j=0;j<ylen;j++){
        d2=dist(xt[i],ya[j]);
        if(d2>score_d82)score[i+1][j+1]=0;else score[i+1][j+1]=1./(1+d2/d02);
    }
}

//heuristic run of dynamic programing iteratively to find the best alignment
//input: initial rotation matrix t, u
//       vectors x and y, d0
//output: best alignment that maximizes the TMscore, will be stored in invmap


inline double SOI_iter(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap0,
    int iteration_max, double local_d0_search,
    double Lnorm, double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt, const bool init_invmap=false)
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
        soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);

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
inline void get_SOI_initial_assign(CoordArray& xk, CoordArray& yk, const int closeK_opt,
    DoubleMatrix& score, CharMatrix& path, DoubleMatrix& val, const int xlen, const int ylen,
    Vec3& t, RotMat& u, std::vector<int>& invmap,
    double d0, double score_d8,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
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
    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);

}

inline void SOI_assign2super(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, CoordArray& xa, CoordArray& ya,
    const int xlen, const int ylen, Vec3& t, RotMat& u, std::vector<int>& invmap,
    double local_d0_search, double Lnorm, double d0, double score_d8)
{
    int i;
    int j;
    int k;
    double rmsd;

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
    TMscore8_search(r1, r2, xtm, ytm, xt, k, t, u,
        40, 8, rmsd, local_d0_search, Lnorm, score_d8, d0);
    do_rotation(xa, xt, xlen, t, u);
}

/* entry function for TM-align with circular permutation
 * i_opt, a_opt, u_opt, d_opt, TMcut are not implemented yet */





struct SoiAlignSearchState
{
    double D0_MIN;
    double Lnorm;
    double score_d8;
    double d0;
    double d0_search;
    double dcu0;
    double local_d0_search;
    int simplify_step;
    int score_sum_method;
    int iteration_max;
    double TM;
    double TMmax;
    Vec3 t;
    RotMat u;
    DoubleMatrix score;
    DoubleMatrix scoret;
    CharMatrix path;
    DoubleMatrix val;
    CoordArray xtm;
    CoordArray ytm;
    CoordArray xt;
    CoordArray yt;
    CoordArray r1;
    CoordArray r2;
    std::vector<int> fwdmap0;
    std::vector<int> invmap0;
};

inline void soi_align_initial_seq(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen, const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &local_d0_search = soi_search_state.local_d0_search;
    DoubleMatrix &score = soi_search_state.score;
    DoubleMatrix &scoret = soi_search_state.scoret;
    CharMatrix &path = soi_search_state.path;
    DoubleMatrix &val = soi_search_state.val;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &yt = soi_search_state.yt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double &TMmax = soi_search_state.TMmax;
    double &TM = soi_search_state.TM;
    int &iteration_max = soi_search_state.iteration_max;
    int i;
    int j;
    CPalign_main(xa_c, ya_c, seqx, seqy, secx, secy,
        chain_align_res,
        xlen, ylen, sequence, align_opts);
    chain_align_res.do_vec.clear();
    if (mm_opt==6)
    {
        i=0;
        j=0;
        for (int r=0;r<seqxA.size();r++)
        {
            if (seqxA[r]=='*') // circular permutation point
            {
                for (int jj=0;jj<j;jj++) if (invmap0[jj]>=0)
                    invmap0[jj]+=xlen - i;
                i=0;
                continue;
            }
            if (seqyA[r]!='-')
            {
                if (seqxA[r]!='-') invmap0[j]=i;
                j++;
            }
            if (seqxA[r]!='-') i++;
        }
        for (j=0;j<ylen;j++)
        {
            i=invmap0[j];
            if (i>=0) fwdmap0[i]=j;
        }
    }
    do_rotation(xa_c, xt, xlen, t0, u0);
    SOI_super2score(xt, ya_c, xlen, ylen, score, d0, score_d8);
    for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];
    TMmax=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa_c, ya_c,
        xlen, ylen, t0, u0, invmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt, true);
    TM   =SOI_iter(r2, r1, ytm, xtm, yt, scoret, path, val, ya_c, xa_c,
        ylen, xlen, t0, u0, fwdmap0, iteration_max,
        local_d0_search, Lnorm, d0, score_d8, secy_bond, secx_bond, mm_opt, true);
    if (TM>TMmax)
    {
        TMmax = TM;
        for (j=0; j<ylen; j++) invmap0[j]=-1;
        for (i=0; i<xlen; i++) 
        {
            j=fwdmap0[i];
            if (j>=0) invmap0[j]=i;
        }
    }
}

inline void soi_align_initial_closek(CoordArray& xa_c, CoordArray& ya_c,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    int xlen, int ylen, SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
{
    std::vector<int> &invmap = chain_align_res.invmap;
    Vec3 &t = soi_search_state.t;
    RotMat &u = soi_search_state.u;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &local_d0_search = soi_search_state.local_d0_search;
    DoubleMatrix &score = soi_search_state.score;
    DoubleMatrix &scoret = soi_search_state.scoret;
    CharMatrix &path = soi_search_state.path;
    DoubleMatrix &val = soi_search_state.val;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &yt = soi_search_state.yt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double &TMmax = soi_search_state.TMmax;
    double &TM = soi_search_state.TM;
    int &iteration_max = soi_search_state.iteration_max;
    int i;
    int j;
    /***************************************************************/
    // initial alignment with sequence order independent alignment
    /***************************************************************/
    if (closeK_opt>=3)
    {
        get_SOI_initial_assign(xk, yk, closeK_opt, score, path, val,
            xlen, ylen, t, u, invmap, d0, score_d8,
            secx_bond, secy_bond, mm_opt);
        for (i=0;i<xlen;i++) for (j=0;j<ylen;j++) scoret[j+1][i+1]=score[i+1][j+1];

        SOI_assign2super(r1, r2, xtm, ytm, xt, xa_c, ya_c,
            xlen, ylen, t, u, invmap, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r1, r2, xtm, ytm, xt, score, path, val, xa_c, ya_c,
            xlen, ylen, t, u, invmap, iteration_max,
            local_d0_search, Lnorm, d0, score_d8, secx_bond, secy_bond, mm_opt);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j = 0; j<ylen; j++) invmap0[j] = invmap[j];
        }

        for (i=0;i<xlen;i++) fwdmap0[i]=-1;
        if (mm_opt==6) NWDP_TM(scoret, path, val, ylen, xlen, -0.6, fwdmap0);
        soi_egs(scoret, ylen, xlen, fwdmap0, secy_bond, secx_bond, mm_opt);
        SOI_assign2super(r2, r1, ytm, xtm, yt, ya_c, xa_c,
            ylen, xlen, t, u, fwdmap0, local_d0_search, Lnorm, d0, score_d8);
        TM=SOI_iter(r2, r1, ytm, xtm, yt, scoret, path, val, ya_c, xa_c, ylen, xlen, t, u,
            fwdmap0, iteration_max, local_d0_search, Lnorm, d0, score_d8,secy_bond, secx_bond, mm_opt);
        if (TM>TMmax)
        {
            TMmax = TM;
            for (j=0; j<ylen; j++) invmap0[j]=-1;
            for (i=0; i<xlen; i++) 
            {
                j=fwdmap0[i];
                if (j>=0) invmap0[j]=i;
            }
        }
    }
}

inline void soi_align_final_tmscore_inv(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res);

inline void soi_align_final_tmscore(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM2 = chain_align_res.TM2;
    double &d0B = chain_align_res.d0B;
    double &rmsd0 = chain_align_res.rmsd0;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    const int mol_type = align_opts.mol_type;
    const bool fast_opt = align_opts.fast_opt;
    double &D0_MIN = soi_search_state.D0_MIN;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &d0_search = soi_search_state.d0_search;
    double &local_d0_search = soi_search_state.local_d0_search;
    int &simplify_step = soi_search_state.simplify_step;
    int &score_sum_method = soi_search_state.score_sum_method;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    Vec3 &t = soi_search_state.t;
    RotMat &u = soi_search_state.u;
    double &TM = soi_search_state.TM;
    int i;
    int j;

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alignment
    simplify_step=1;
    if (fast_opt) simplify_step=40;
    score_sum_method=8;
    TM = detailed_search_standard(r1, r2, xtm, ytm, xt, xa_c, ya_c, ylen,
        invmap0, t, u, simplify_step, score_sum_method, local_d0_search,
        false, Lnorm, score_d8, d0);
    
    double rmsd;
    simplify_step=1;
    score_sum_method=0;

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    double d;
    std::vector<int> m1(xlen); //alignd index in x
    std::vector<int> m2(ylen); //alignd index in y
    copy_t_u(t, u, t0, u0);
    
    //****************************************//
    //              Final TMscore 1           //
    //****************************************//

    do_rotation(xa_c, xt, xlen, t, u);
    k=0;
    n_ali=0;
    for (i=0; i<xlen; i++)
    {
        j=fwdmap0[i];
        if(j>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xt[i], ya_c[j]));
            if (d <= score_d8)
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa_c[i][0];
                xtm[k][1]=xa_c[i][1];
                xtm[k][2]=xa_c[i][2];

                ytm[k][0]=ya_c[j][0];
                ytm[k][1]=ya_c[j][1];
                ytm[k][2]=ya_c[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya_c[j][0];
                r2[k][1] = ya_c[j][1];
                r2[k][2] = ya_c[j][2];

                k++;
            }
            else fwdmap0[i]=-1;
        }
    }
    n_ali8=k;

    Kabsch(r1, r2, n_ali8, 0, rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    rmsd0 = sqrt(rmsd0 / n_ali8);

    //normalized by length of structure A
    parameter_set4final(xlen+0.0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0B=d0;
    local_d0_search = d0_search;
    TM2 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t, u, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    soi_align_final_tmscore_inv(xa_c, ya_c, xlen, ylen, align_opts, soi_search_state, chain_align_res);
}

inline void soi_align_final_tmscore_inv(CoordArray& xa_c, CoordArray& ya_c,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM1 = chain_align_res.TM1;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    double &d0_0 = chain_align_res.d0_0;
    double &TM_0 = chain_align_res.TM_0;
    double &d0A = chain_align_res.d0A;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    double &d0_out = chain_align_res.d0_out;
    int &n_ali8 = chain_align_res.n_ali8;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;
    const int a_opt = align_opts.a_opt;
    const bool u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    const int mol_type = align_opts.mol_type;
    double &D0_MIN = soi_search_state.D0_MIN;
    double &Lnorm = soi_search_state.Lnorm;
    double &score_d8 = soi_search_state.score_d8;
    double &d0 = soi_search_state.d0;
    double &d0_search = soi_search_state.d0_search;
    double &local_d0_search = soi_search_state.local_d0_search;
    int &simplify_step = soi_search_state.simplify_step;
    int &score_sum_method = soi_search_state.score_sum_method;
    CoordArray &xtm = soi_search_state.xtm;
    CoordArray &ytm = soi_search_state.ytm;
    CoordArray &xt = soi_search_state.xt;
    CoordArray &r1 = soi_search_state.r1;
    CoordArray &r2 = soi_search_state.r2;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    double rmsd;
    double Lnorm_0=ylen;
    int i;
    int j;
    int k;
    double d;
    std::vector<int> m1(xlen);
    std::vector<int> m2(ylen);
    //****************************************//
    //              Final TMscore 2           //
    //****************************************//

    do_rotation(xa_c, xt, xlen, t0, u0);
    k=0;
    for (j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
            d=sqrt(dist(xt[i], ya_c[j]));
            if (d <= score_d8)
            {
                m1[k]=i;
                m2[k]=j;

                xtm[k][0]=xa_c[i][0];
                xtm[k][1]=xa_c[i][1];
                xtm[k][2]=xa_c[i][2];

                ytm[k][0]=ya_c[j][0];
                ytm[k][1]=ya_c[j][1];
                ytm[k][2]=ya_c[j][2];

                r1[k][0] = xt[i][0];
                r1[k][1] = xt[i][1];
                r1[k][2] = xt[i][2];
                r2[k][0] = ya_c[j][0];
                r2[k][1] = ya_c[j][1];
                r2[k][2] = ya_c[j][2];

                k++;
            }
            else invmap[j]=invmap0[j]=-1;
        }
    }

    //normalized by length of structure B
    parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
    d0A=d0;
    d0_0=d0A;
    local_d0_search = d0_search;
    TM1 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0, simplify_step,
        score_sum_method, rmsd, local_d0_search, Lnorm, score_d8, d0);
    TM_0 = TM1;

    if (a_opt>0)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, D0_MIN, Lnorm, d0, d0_search, mol_type);
        d0a=d0;
        d0_0=d0a;
        local_d0_search = d0_search;

        TM3 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0, d0_search, mol_type);
        d0u=d0;
        d0_0=d0u;
        Lnorm_0=Lnorm_ass;
        local_d0_search = d0_search;
        TM4 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, Lnorm, d0, d0_search);
        d0_out=d0_scale;
        d0_0=d0_scale;
        local_d0_search = d0_search;
        TM5 = TMscore8_search(r1, r2, xtm, ytm, xt, n_ali8, t0, u0,
            simplify_step, score_sum_method, rmsd, local_d0_search, Lnorm,
            score_d8, d0);
        TM_0=TM5;
    }
}


inline void soi_derive_alignment_strings(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy, int xlen, int ylen,
    SoiAlignSearchState& soi_search_state, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Liden = chain_align_res.Liden;
    const double d0_out = chain_align_res.d0_out;
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    CoordArray &xt = soi_search_state.xt;
    std::vector<int> &fwdmap0 = soi_search_state.fwdmap0;
    std::vector<int> &invmap0 = soi_search_state.invmap0;
    int i;
    int j;
    int k;
    double d;
    // derive alignment from superposition
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap0[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    
    do_rotation(xa_c, xt, xlen, t0, u0);

    Liden=0;
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap0[j];
        dist_list[j]=-1;
        if (i<0) continue;
        d=sqrt(dist(xt[i], ya_c[j]));
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        dist_list[j]=d;
        seqxA[j]=seqx[i];
        Liden+=(seqx[i]==seqy[j]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap0[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }

    // clean up

}

inline int SOIalign_main(CoordArray& xa_c, CoordArray& ya_c,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    const std::string &seqx, const std::string &seqy, const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
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
        xlen, ylen, sequence, align_opts, soi_search_state, chain_align_res, secx_bond, secy_bond, mm_opt);
    soi_align_initial_closek(xa_c, ya_c, xk, yk, closeK_opt,
        xlen, ylen, soi_search_state, chain_align_res, secx_bond, secy_bond, mm_opt);

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
        std::cout << "There is no alignment between the two structures! "
             << "Program stop with no result!" << std::endl;
        TM1=TM2=TM3=TM4=TM5=0;
        init_rt_mat(t0, u0);
        return 1;
    }

    soi_align_final_tmscore(xa_c, ya_c, xlen, ylen, align_opts, soi_search_state, chain_align_res);
    soi_derive_alignment_strings(xa_c, ya_c, seqx, seqy, xlen, ylen, soi_search_state, chain_align_res, dist_list);

    return 0;
}




inline void soi_build_score_seed_invmap(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts,
    ChainPairAlignResult& chain_align_res, IntPairArray& secx_bond, IntPairArray& secy_bond,
    const int mm_opt, DoubleMatrix& score)
{
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const int mol_type = align_opts.mol_type;
    double D0_MIN;
    double Lnorm;
    double score_d8,d0,d0_search,dcu0;
    CharMatrix  path;
    DoubleMatrix val;
    int i;
    int j;

    /***********************/
    // allocate memory
    /***********************/
    score.assign(xlen+1, std::vector<double>(ylen+1));
    path.assign( xlen+1, std::vector<char>(ylen+1));
    val.assign(  xlen+1, std::vector<double>(ylen+1));


    // set d0
    parameter_set4search(xlen, ylen, D0_MIN, Lnorm,
        score_d8, d0, d0_search, dcu0); // set score_d8
    parameter_set4final(xlen, D0_MIN, Lnorm,
        d0B, d0_search, mol_type); // set d0B
    parameter_set4final(ylen, D0_MIN, Lnorm,
        d0A, d0_search, mol_type); // set d0A
    if (a_opt)
        parameter_set4final((xlen+ylen)*0.5, D0_MIN, Lnorm,
            d0a, d0_search, mol_type); // set d0a
    if (u_opt)
    {
        parameter_set4final(Lnorm_ass, D0_MIN, Lnorm,
            d0u, d0_search, mol_type); // set d0u
        if (u_opt==2)
        {
            parameter_set4search(Lnorm_ass, Lnorm_ass, D0_MIN, Lnorm,
                score_d8, d0, d0_search, dcu0); // set score_d8
        }
    }

    // perform alignment
    for(j=0; j<ylen; j++) invmap[j]=-1;
    double d02=d0*d0;
    double score_d82=score_d8*score_d8;
    double d2;
    for(i=0; i<xlen; i++)
    {
        for(j=0; j<ylen; j++)
        {
            d2=dist(xa[i], ya[j]);
            if (d2>score_d82) score[i+1][j+1]=0;
            else score[i+1][j+1]=1./(1+ d2/d02);
        }
    }
    if (mm_opt==6) NWDP_TM(score, path, val, xlen, ylen, -0.6, invmap);

    soi_egs(score, xlen, ylen, invmap, secx_bond, secy_bond, mm_opt);
}

inline bool soi_se_score_alignment(CoordArray& xa, CoordArray& ya,
    int xlen, int ylen, const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    const DoubleMatrix& score, ChainPairAlignResult& chain_align_res, std::vector<double>& dist_list)
{
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0u = chain_align_res.d0u;
    double &d0a = chain_align_res.d0a;
    double &rmsd0 = chain_align_res.rmsd0;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    std::vector<int> m1;
    std::vector<int> m2;
    int i;
    int j;
    double d;
    if (outfmt_opt<2)
    {
        m1.resize(xlen);
        m2.resize(ylen);
    }

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(j=0; j<ylen; j++)
    {
        i=invmap[j];
        dist_list[j]=-1;
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xa[i], ya[j]));
            dist_list[j]=d;
            if (score[i+1][j+1]>0)
            {
                if (outfmt_opt<2)
                {
                    m1[k]=i;
                    m2[k]=j;
                }
                k++;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);
    if (outfmt_opt>=2)
    {
        return true;
    }
    return false;
}

inline void soi_se_extract_strings(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    ChainPairAlignResult& chain_align_res)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    double &Liden = chain_align_res.Liden;
    const double d0_out = chain_align_res.d0_out;
    std::vector<int> &invmap = chain_align_res.invmap;
    int i;
    int j;
    int k = chain_align_res.n_ali8;
    double d;

    // extract aligned sequence
    int ali_len=xlen+ylen;
    for (j=0;j<ylen;j++) ali_len-=(invmap[j]>=0);
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');

    std::vector<int> fwdmap(xlen+1, -1);
    for (j=0;j<ylen;j++)
    {
        seqyA[j]=seqy[j];
        i=invmap[j];
        if (i<0) continue;
        d=sqrt(dist(xa[i], ya[j]));
        if (d<d0_out) seqM[j]=':';
        else seqM[j]='.';
        fwdmap[i]=j;
        seqxA[j]=seqx[i];
        Liden+=(seqxA[k]==seqyA[k]);
    }
    k=0;
    for (i=0;i<xlen;i++)
    {
        j=fwdmap[i];
        if (j>=0) continue;
        seqxA[ylen+k]=seqx[i];
        k++;
    }
}

inline int soi_se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const ChainPairAlignOptions& align_opts, const int outfmt_opt,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt)
{
    DoubleMatrix score;
    soi_build_score_seed_invmap(xa, ya, xlen, ylen, align_opts, chain_align_res,
        secx_bond, secy_bond, mm_opt, score);
    if (soi_se_score_alignment(xa, ya, xlen, ylen, align_opts, outfmt_opt,
        score, chain_align_res, dist_list)) return 0;
    soi_se_extract_strings(xa, ya, seqx, seqy, xlen, ylen, chain_align_res);
    return 0;
}


inline void soi_align_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    CoordArray& xk, CoordArray& yk, const int closeK_opt,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    std::vector<double>& dist_list,
    IntPairArray& secx_bond, IntPairArray& secy_bond, const int mm_opt,
    int outfmt_opt)
{
    chain_align_res.invmap.assign(ylen + 1, 0);
    if (align_opts.se_opt)
    {
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);
        soi_se_main(xa, ya, seqx, seqy, chain_align_res,
            xlen, ylen, align_opts, outfmt_opt, dist_list, secx_bond, secy_bond, mm_opt);
        if (outfmt_opt >= 2)
        {
            chain_align_res.Liden=chain_align_res.L_ali=0;
            int r1;
            int r2;
            for (r2=0;r2<ylen;r2++)
            {
                r1=chain_align_res.invmap[r2];
                if (r1<0) continue;
                chain_align_res.L_ali+=1;
                chain_align_res.Liden+=(seqx[r1]==seqy[r2]);
            }
        }
    }
    else SOIalign_main(xa, ya, xk, yk, closeK_opt,
        seqx, seqy, secx, secy, chain_align_res,
        xlen, ylen, sequence, align_opts, dist_list, secx_bond, secy_bond, mm_opt);
}

#endif
