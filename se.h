#pragma once
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
    const ChainPairAlignOptions& opt, int mol_type, ChainPairAlignResult& res, SeSearchState& search)
{
    double &D0_MIN = search.D0_MIN;
    double &Lnorm = search.Lnorm;
    double &score_d8 = search.score_d8;
    double &d0 = search.d0;
    double &d0_search = search.d0_search;
    double &dcu0 = search.dcu0;
    double &d0A = res.d0A;
    double &d0B = res.d0B;
    double &d0a = res.d0a;
    double &d0u = res.d0u;
    const bool a_opt = opt.a_opt;
    const int u_opt = opt.u_opt;
    const double Lnorm_ass = opt.Lnorm;
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

}

inline void se_build_initial_alignment(CharMatrix& path, DoubleMatrix& val,
    CoordArray& xa, CoordArray& ya, int xlen, int ylen,
    const std::vector<std::string>& sequence, const ChainPairAlignOptions& opt,
    const SeSearchState& search, ChainPairAlignResult& res, const int hinge)
{
    const double d0 = search.d0;
    const bool i_opt = opt.i_opt;
    std::vector<int> &invmap = res.invmap;
    int j;
    // perform alignment
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
    const ChainPairAlignOptions& opt, const SeSearchState& search, const std::vector<int>& invmap0,
    ChainPairAlignResult& res, std::vector<int>& m1, std::vector<int>& m2,
    const int outfmt_opt, const int hinge)
{
    double &TM1 = res.TM1;
    double &TM2 = res.TM2;
    double &TM3 = res.TM3;
    double &TM4 = res.TM4;
    double &TM5 = res.TM5;
    double &rmsd0 = res.rmsd0;
    int &n_ali = res.n_ali;
    int &n_ali8 = res.n_ali8;
    std::vector<int> &invmap = res.invmap;
    const double Lnorm_ass = opt.Lnorm;
    const double d0_scale = opt.d0_scale;
    const double d0A = res.d0A;
    const double d0B = res.d0B;
    const double d0a = res.d0a;
    const double d0u = res.d0u;
    const double score_d8 = search.score_d8;
    const bool a_opt = opt.a_opt;
    const int u_opt = opt.u_opt;
    const bool d_opt = opt.d_opt;
    const bool i_opt = opt.i_opt;
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
    ChainPairAlignResult& res, const int hinge)
{
    std::string &seqM = res.seqM;
    std::string &seqxA = res.seqxA;
    std::string &seqyA = res.seqyA;
    std::vector<double> &do_vec = res.do_vec;
    double &Liden = res.Liden;
    const double d0_out = res.d0_out;
    const int n_ali8 = res.n_ali8;
    int j;
    double d;
    // extract aligned sequence
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

inline int se_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    ChainPairAlignResult& res,
    const int xlen, const int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& opt, const int outfmt_opt, const int hinge=0)
{
    std::string &seqM = res.seqM;
    std::string &seqyA = res.seqyA;
    std::vector<int> &invmap = res.invmap;
    const int mol_type = opt.mol_type;
    SeSearchState search;
    CharMatrix path;          // for dynamic programming (char: 1/0)
    DoubleMatrix val;          // for dynamic programming

    std::vector<int> m1;
    std::vector<int> m2;
    double d;
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
    int i;
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

    se_prepare_search_parameters(xlen, ylen, opt, mol_type, res, search);
    se_build_initial_alignment(path, val, xa, ya, xlen, ylen, sequence, opt, search, res, hinge);
    se_score_alignment(xa, ya, xlen, ylen, opt, search, invmap0, res, m1, m2, outfmt_opt, hinge);
    if (outfmt_opt>=2)
    {
        if (hinge) seqM_char.clear();    

        return 0;
    }

    se_extract_alignment_strings(xa, ya, seqx, seqy, xlen, ylen, m1, m2, seqM_char, res, hinge);
    return 0; // zero for no exception
    se_extract_alignment_strings(xa, ya, seqx, seqy, xlen, ylen, m1, m2, seqM_char, res, hinge);
    return 0; // zero for no exception
}


int se_main(
    CoordArray& xa, CoordArray& ya, const std::string &seqx, const std::string &seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    std::string &seqM, std::string &seqxA, std::string &seqyA, std::vector<double> &do_vec,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen, const std::vector<std::string> &sequence,
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool a_opt, const int u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, std::vector<int>& invmap, const int hinge=0)
{
    ChainPairAlignResult res = { 0};
    res.TM1 = TM1;
    res.TM2 = TM2;
    res.TM3 = TM3;
    res.TM4 = TM4;
    res.TM5 = TM5;
    res.d0_0 = d0_0;
    res.TM_0 = TM_0;
    res.d0A = d0A;
    res.d0B = d0B;
    res.d0u = d0u;
    res.d0a = d0a;
    res.d0_out = d0_out;
    res.seqM = seqM;
    res.seqxA = seqxA;
    res.seqyA = seqyA;
    res.do_vec = do_vec;
    res.rmsd0 = rmsd0;
    res.L_ali = L_ali;
    res.Liden = Liden;
    res.TM_ali = TM_ali;
    res.rmsd_ali = rmsd_ali;
    res.n_ali = n_ali;
    res.n_ali8 = n_ali8;
    res.invmap = invmap;
    ChainPairAlignOptions opt;
    opt.i_opt = i_opt;
    opt.a_opt = a_opt;
    opt.u_opt = u_opt;
    opt.d_opt = d_opt;
    opt.fast_opt = false;
    opt.se_opt = false;
    opt.cp_opt = false;
    opt.Lnorm = Lnorm_ass;
    opt.d0_scale = d0_scale;
    opt.TMcut = -1;
    opt.parallel_threads = 1;
    opt.ss_opt = 0;
    opt.mol_type = mol_type;
    int rc = se_main(xa, ya, seqx, seqy, res, xlen, ylen, sequence, opt, outfmt_opt, hinge);
    TM1 = res.TM1;
    TM2 = res.TM2;
    TM3 = res.TM3;
    TM4 = res.TM4;
    TM5 = res.TM5;
    d0_0 = res.d0_0;
    TM_0 = res.TM_0;
    d0A = res.d0A;
    d0B = res.d0B;
    d0u = res.d0u;
    d0a = res.d0a;
    d0_out = res.d0_out;
    seqM = res.seqM;
    seqxA = res.seqxA;
    seqyA = res.seqyA;
    do_vec = res.do_vec;
    rmsd0 = res.rmsd0;
    L_ali = res.L_ali;
    Liden = res.Liden;
    TM_ali = res.TM_ali;
    rmsd_ali = res.rmsd_ali;
    n_ali = res.n_ali;
    n_ali8 = res.n_ali8;
    invmap = res.invmap;
    return rc;
}

