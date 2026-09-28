#pragma once
#ifndef TMalign_HwRMSD_h
#define TMalign_HwRMSD_h 1
#include "NWalign.h"
#include "se.h"

double Kabsch_Superpose(CoordArray& r1, CoordArray& r2, CoordArray& xt,
    CoordArray& xa, CoordArray& ya, int xlen, int ylen, std::vector<int>& invmap,
    int& L_ali, Vec3& t, RotMat& u)
{
    L_ali = 0;
    int i;
    int j;
    for (j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            r1[L_ali][0]  = xa[i][0];
            r1[L_ali][1]  = xa[i][1];
            r1[L_ali][2]  = xa[i][2];

            r2[L_ali][0]  = ya[j][0];
            r2[L_ali][1]  = ya[j][1];
            r2[L_ali][2]  = ya[j][2];

            L_ali++;
        }
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");
    }

    double RMSD = 0;
    Kabsch(r1, r2, L_ali, 1, RMSD, t, u);
    RMSD = sqrt( RMSD/(1.0*L_ali) );

    for (i=0; i<xlen; i++)
    {
        xt[i][0] = xa[i][0];
        xt[i][1] = xa[i][1];
        xt[i][2] = xa[i][2];
    }
    do_rotation(xa, xt, xlen, t,u);
    return RMSD;
}

void parse_alignment_into_invmap(const std::string seqxA_tmp,
    const std::string seqyA_tmp, const int xlen, const int ylen, std::vector<int>& invmap_tmp)
{
    if (seqxA_tmp.size()==0) return;
    int i1=-1;
    int i2=-1;
    int j = 0;
    int L = std::min(seqxA_tmp.size(), seqyA_tmp.size());
    for (j = 0; j < ylen; j++) invmap_tmp[j] = -1;
    for (j = 0; j<L; j++)
    {
        if (seqxA_tmp[j] != '-') i1++;
        if (seqyA_tmp[j] != '-')
        {
            i2++;
            if (i2 >= ylen || i1 >= xlen) j = L;
            else if (seqxA_tmp[j] != '-') invmap_tmp[i2] = i1;
        }
    }
    return;
}


int HwRMSD_main(CoordArray& xa, CoordArray& ya, const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy, Vec3& t0, RotMat& u0,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &TM_0, double &d0A, double &d0B, double &d0u,
    double &d0a, double &d0_out, std::string &seqM, std::string &seqxA, std::string &seqyA,
    double &rmsd0, int &L_ali, double &Liden, double &TM_ali,
    double &rmsd_ali, int &n_ali, int &n_ali8, const int xlen, const int ylen,
    const std::vector<std::string>&sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt,
    const int a_opt, const bool u_opt, const bool d_opt, const int mol_type,
    std::vector<int>& invmap, const int glocal=0, const int iter_opt=10,
    const int seq_opt=3, const double early_opt=0.01)
{
    /***********************/
    // allocate memory    
    /***********************/
    Vec3 t; RotMat u;
    CoordArray xt;
    CoordArray r1, r2;
    int minlen = std::min(xlen, ylen);
    xt.resize(xlen);
    r1.resize(minlen);
    r2.resize(minlen);
    ChainPairAlignResult chain_align_res = {};
    chain_align_res.d0_out = d0_out;
    chain_align_res.invmap.assign(ylen+1, 0);
    std::vector<int>& invmap_tmp = chain_align_res.invmap;

    int i;
    int j;
    int i1;
    double& TM1_tmp = chain_align_res.TM1;
    double& TM2_tmp = chain_align_res.TM2;
    double& TM3_tmp = chain_align_res.TM3;
    double& TM4_tmp = chain_align_res.TM4;
    double& TM5_tmp = chain_align_res.TM5;
    double& TM_ali_tmp = chain_align_res.TM_ali;
    std::string& seqxA_tmp = chain_align_res.seqxA;
    std::string& seqyA_tmp = chain_align_res.seqyA;
    std::string& seqM_tmp = chain_align_res.seqM;
    double& rmsd0_tmp = chain_align_res.rmsd0;
    int& L_ali_tmp = chain_align_res.L_ali;
    int& n_ali_tmp = chain_align_res.n_ali;
    int& n_ali8_tmp = chain_align_res.n_ali8;
    double& Liden_tmp = chain_align_res.Liden;
    double& rmsd_ali_tmp = chain_align_res.rmsd_ali;
    double max_TM=0;
    double cur_TM=0;

    ChainPairAlignOptions align_opts;
    align_opts.i_opt = (i_opt==3);
    align_opts.a_opt = a_opt;
    align_opts.u_opt = u_opt;
    align_opts.d_opt = d_opt;
    align_opts.fast_opt = false;
    align_opts.se_opt = false;
    align_opts.cp_opt = false;
    align_opts.Lnorm = Lnorm_ass;
    align_opts.d0_scale = d0_scale;
    align_opts.TMcut = -1;
    align_opts.parallel_threads = 1;
    align_opts.ss_opt = 0;
    align_opts.mol_type = mol_type;

    TM1=TM2=TM1_tmp=TM2_tmp=L_ali=-1;

    if (i_opt)
    {
        seqxA_tmp=sequence[0];
        seqyA_tmp=sequence[1];
    }
    else if (seq_opt==2) NWalign_main(secx.c_str(), secy.c_str(), xlen, ylen,
            seqxA_tmp, seqyA_tmp, mol_type, invmap_tmp, 1, glocal);
    else NWalign_main(seqx.c_str(), seqy.c_str(), xlen, ylen,
            seqxA_tmp, seqyA_tmp, mol_type, invmap_tmp, 1, glocal);
    int total_iter=(i_opt==3 || iter_opt<1)?1:iter_opt;

    for (int iter=0;iter<total_iter;iter++)
    {
        n_ali_tmp=n_ali8_tmp=0;
        if (iter==1 && !i_opt && seq_opt==3) NWalign_main(secx.c_str(), secy.c_str(), xlen,
            ylen, seqxA_tmp, seqyA_tmp, mol_type, invmap_tmp, 1, glocal);

        parse_alignment_into_invmap(seqxA_tmp, seqyA_tmp, xlen, ylen, invmap_tmp);

        Kabsch_Superpose(r1, r2, xt, xa, ya, xlen, ylen, invmap_tmp,
            L_ali, t, u);

        se_main(xt, ya, seqx.c_str(), seqy.c_str(), chain_align_res,
            xlen, ylen, sequence, align_opts, 1);
        d0A = chain_align_res.d0A;
        d0B = chain_align_res.d0B;
        d0u = chain_align_res.d0u;
        d0a = chain_align_res.d0a;

        if (n_ali8_tmp==0)
        {
            if (xlen>=ylen) seqxA_tmp = seqx;
            if (xlen<=ylen) seqyA_tmp = seqy;
            if (xlen<ylen)
            {
                seqxA_tmp.clear();
                for (i1=0;i1<static_cast<int>((ylen-xlen)/2);i1++) seqxA_tmp+='-';
                seqxA_tmp += seqx;
                for (i1=seqxA_tmp.size();i1<ylen;i1++) seqxA_tmp+='-';
            }
            if (xlen>ylen)
            {
                seqyA_tmp.clear();
                for (i1=0;i1<static_cast<int>((xlen-ylen)/2);i1++) seqyA_tmp+='-';
                seqyA_tmp += seqy;
                for (i1=seqyA_tmp.size();i1<xlen;i1++) seqyA_tmp+='-';
            }
        
            parse_alignment_into_invmap(seqxA_tmp, seqyA_tmp, xlen, ylen, invmap_tmp);

            Kabsch_Superpose(r1, r2, xt, xa, ya, xlen, ylen, invmap_tmp,
                L_ali, t, u);

            se_main(xt, ya, seqx.c_str(), seqy.c_str(), chain_align_res,
                xlen, ylen, sequence, align_opts, 1);
        }

        if (TM1_tmp>TM1 && TM2_tmp>TM2)
        {
            for (i=0; i<3; i++)
            {
                t0[i]=t[i];
                for (j=0;j<3;j++) u0[i][j]=u[i][j];
            }
            TM1=TM1_tmp;
            TM2=TM2_tmp;
            TM3=TM3_tmp;
            TM4=TM4_tmp;
            TM5=TM5_tmp;

            TM_0=TM1;
            if (a_opt>0) TM_0=TM3;
            else if (u_opt) TM_0=TM4;
            else if (d_opt) TM_0=TM5;

            seqxA =seqxA_tmp;
            seqM  =seqM_tmp;
            seqyA =seqyA_tmp;
            for (j=0; j<ylen; j++) invmap[j]=invmap_tmp[j];

            rmsd0 =rmsd0_tmp;
            Liden =Liden_tmp;
            n_ali =n_ali_tmp;
            n_ali8=n_ali8_tmp;

            if (i_opt && L_ali==-1)
            {
                L_ali=L_ali_tmp;
                TM_ali=TM_ali_tmp;
                rmsd_ali=rmsd_ali_tmp;
            }
        }
        else
        {
            if (iter>=2) break;
            seqxA_tmp  = seqxA;
            seqyA_tmp  = seqyA;
            for (j=0; j<ylen; j++) invmap_tmp[j]=invmap[j];
            rmsd0_tmp  = 0;
            Liden_tmp  = 0;
            n_ali_tmp  = 0;
            n_ali8_tmp = 0;
        }

        if (iter>=2 && early_opt>0)
        {
            cur_TM=(TM1+TM2)/2;
            if (cur_TM-max_TM<early_opt) break;
            max_TM=cur_TM;
        }
    }

    return 0;
}

#endif
