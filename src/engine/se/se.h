#pragma once
#include "common.h"
#include "engine_common.h"
#include "TMalign.h"

/* entry function for se
 * outfmt_opt>=2 should not parse sequence alignment
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * if hinge>0, append to original invmap */
struct SeSearchState
{
    double D0_MIN;
    double Lnorm;
    double score_d8;
    double d0;
    double d0_search;
    double dcu0;
};


inline void se_prepare_search_parameters(int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, int mol_type, ChainPairAlignResult& chain_align_res, SeSearchState& se_search_state)
{
    double &D0_MIN = se_search_state.D0_MIN;
    double &Lnorm = se_search_state.Lnorm;
    double &score_d8 = se_search_state.score_d8;
    double &d0 = se_search_state.d0;
    double &d0_search = se_search_state.d0_search;
    double &dcu0 = se_search_state.dcu0;
    double &d0A = chain_align_res.d0A;
    double &d0B = chain_align_res.d0B;
    double &d0a = chain_align_res.d0a;
    double &d0u = chain_align_res.d0u;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const double Lnorm_ass = align_opts.Lnorm;
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

}

inline void se_build_initial_alignment(CharMatrix& path, DoubleMatrix& val,
    CoordArray& xa, CoordArray& ya, int xlen, int ylen,
    const std::vector<std::string>& sequence, const ChainPairAlignOptions& align_opts,
    const SeSearchState& se_search_state, ChainPairAlignResult& chain_align_res, const int hinge)
{
    const double d0 = se_search_state.d0;
    const bool i_opt = align_opts.i_opt;
    std::vector<int> &invmap = chain_align_res.invmap;
    int j;
    if (hinge==0) for(j=0; j<ylen; j++) invmap[j]=-1;
    if (!i_opt) NWDP_SE(path, val, xa, ya, xlen, ylen, d0*d0, 0, invmap, hinge);
    else
    {
        int i1 = -1;// in C version, index starts from zero, not from one
        int i2 = -1;
        int L1 = sequence[0].size();
        int L2 = sequence[1].size();
        int L = std::min(L1, L2);// Get positions for aligned residues
        for (int kk1 = 0; kk1 < L; kk1++)
        {
            if (sequence[0][kk1] != '-') i1++;
            if (sequence[1][kk1] != '-')
            {
                i2++;
                if (i2 >= ylen || i1 >= xlen) kk1 = L;
                else if (sequence[0][kk1] != '-') invmap[i2] = i1;
            }
        }
    }
    
}

inline void se_score_alignment(CoordArray& xa, CoordArray& ya, int xlen, int ylen,
    const ChainPairAlignOptions& align_opts, const SeSearchState& se_search_state, const std::vector<int>& invmap0,
    ChainPairAlignResult& chain_align_res, std::vector<int>& m1, std::vector<int>& m2,
    const int outfmt_opt, const int hinge)
{
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    double &rmsd0 = chain_align_res.rmsd0;
    int &n_ali = chain_align_res.n_ali;
    int &n_ali8 = chain_align_res.n_ali8;
    std::vector<int> &invmap = chain_align_res.invmap;
    const double Lnorm_ass = align_opts.Lnorm;
    const double d0_scale = align_opts.d0_scale;
    const double d0A = chain_align_res.d0A;
    const double d0B = chain_align_res.d0B;
    const double d0a = chain_align_res.d0a;
    const double d0u = chain_align_res.d0u;
    const double score_d8 = se_search_state.score_d8;
    const bool a_opt = align_opts.a_opt;
    const int u_opt = align_opts.u_opt;
    const bool d_opt = align_opts.d_opt;
    const bool i_opt = align_opts.i_opt;
    double d;
    if (hinge==0) rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    else
    {
        TM2*=xlen;
        TM1*=ylen;
        TM3*=(xlen+ylen)*0.5;
        TM4*=Lnorm_ass;
        TM5*=ylen;
        rmsd0=rmsd0*rmsd0*n_ali8;
    }
    int k=0;
    n_ali=0;
    n_ali8=0;
    for(int i=0,j=0; j<ylen; j++)
    {
        i=invmap[j];
        if(i>=0)//aligned
        {
            n_ali++;
            d=sqrt(dist(xa[i], ya[j]));
            if (d <= score_d8 || i_opt || invmap0[j]==i)
            {
                if (outfmt_opt<2)
                {
                    m1[k]=i;
                    m2[k]=j;
                }
                k++;
                if (invmap0[j]==i) continue;
                TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
                TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
                if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
                if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
                if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
                rmsd0+=d*d;
            }
            else if (hinge) invmap[j]=-1;
        }
    }
    n_ali8=k;
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);

}

inline void se_extract_alignment_strings(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, int xlen, int ylen,
    const std::vector<int>& m1, const std::vector<int>& m2, const std::vector<char>& seqM_char,
    ChainPairAlignResult& chain_align_res, const int hinge)
{
    std::string &seqM = chain_align_res.seqM;
    std::string &seqxA = chain_align_res.seqxA;
    std::string &seqyA = chain_align_res.seqyA;
    std::vector<double> &do_vec = chain_align_res.do_vec;
    double &Liden = chain_align_res.Liden;
    const double d0_out = chain_align_res.d0_out;
    const int n_ali8 = chain_align_res.n_ali8;
    int j;
    double d;
    int ali_len=xlen+ylen; //maximum length of alignment
    seqxA.assign(ali_len,'-');
    seqM.assign( ali_len,' ');
    seqyA.assign(ali_len,'-');
    do_vec.clear();
    do_vec.assign(ali_len,0);
    
    int kk = 0;
    int i_old = 0;
    int j_old = 0;
    d=0;
    Liden=0;
    for(int k=0; k<n_ali8; k++)
    {
        for(int i=i_old; i<m1[k]; i++)
        {
            //align x to gap
            seqxA[kk]=seqx[i];
            seqyA[kk]='-';
            seqM[kk]=' ';                    
            kk++;
        }

        for(int j=j_old; j<m2[k]; j++)
        {
            //align y to gap
            seqxA[kk]='-';
            seqyA[kk]=seqy[j];
            seqM[kk]=' ';
            kk++;
        }

        seqxA[kk]=seqx[m1[k]];
        seqyA[kk]=seqy[m2[k]];
        Liden+=(seqxA[kk]==seqyA[kk]);
        d=sqrt(dist(xa[m1[k]], ya[m2[k]]));
        if(d<d0_out) seqM[kk]=':';
        else         seqM[kk]='.';
        do_vec[kk]=d;
        kk++;  
        i_old=m1[k]+1;
        j_old=m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        seqxA[kk]=seqx[i];
        seqyA[kk]='-';
        seqM[kk]=' ';
        kk++;
    }    
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        seqxA[kk]='-';
        seqyA[kk]=seqy[j];
        seqM[kk]=' ';
        kk++;
    }
    seqxA=seqxA.substr(0,kk);
    seqyA=seqyA.substr(0,kk);
    seqM =seqM.substr(0,kk);
    if (hinge)
    {
        j=-1;
        for (int r=0;r<seqM.size();r++)
        {
            j+=seqyA[r]!='-';
            if (seqM[r]!=' ') seqM[r]=seqM_char[j];
        }
    }

}

