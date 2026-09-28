#pragma once
#include "TMalign.h"

int score_fun8(const CoordArray& xa, const CoordArray& ya, int n_ali, double d, std::vector<int>& i_ali,
    double &score1, int score_sum_method, const double Lnorm,
    const double score_d8, const double d0,
    double GDT_list_tmp[5], double &maxsub_tmp)
{
    double score_sum=0;
    double di;
    double d_tmp=d*d;
    double d02=d0*d0;
    double score_d8_cut = score_d8*score_d8;

    int i;
    int n_cut;
    int inc=0;

    while(1)
    {
        for (i=0;i<5;i++) GDT_list_tmp[i]=0;
        maxsub_tmp=0;

        n_cut=0;
        score_sum=0;
        for(i=0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if(di<d_tmp)
            {
                i_ali[n_cut]=i;
                n_cut++;
            }
            if(score_sum_method==8)
            {
                if(di<=score_d8_cut) score_sum += 1/(1+di/d02);
            }
            else score_sum += 1/(1+di/d02);

            if (di<64)
            {
                GDT_list_tmp[4]+=1;
                if (di<16)
                {
                    GDT_list_tmp[3]+=1;
                    if (di<12.25)
                    {
                        maxsub_tmp+=1/(1+di/12.25);
                        if (di<4)
                        {
                            GDT_list_tmp[2]+=1;
                            if (di<1)
                            {
                                GDT_list_tmp[1]+=1;
                                if (di<0.25)
                                    GDT_list_tmp[0]+=1;
                            }
                        }
                    }
                }
            }
        }
        if(n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc=(d+inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    score1=score_sum/Lnorm;
    return n_cut;
}

int score_fun8_standard(const CoordArray& xa, const CoordArray& ya, int n_ali, double d,
    std::vector<int>& i_ali, double &score1, int score_sum_method,
    double score_d8, double d0, double GDT_list_tmp[5], double &maxsub_tmp)
{
    double score_sum = 0;
    double di;
    double d_tmp = d*d;
    double d02 = d0*d0;
    double score_d8_cut = score_d8*score_d8;

    int i;
    int n_cut;
    int inc = 0;
    while (1)
    {
        for (i=0;i<5;i++) GDT_list_tmp[i]=0;
        maxsub_tmp=0;
        n_cut = 0;
        score_sum = 0;
        for (i = 0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if (di<d_tmp)
            {
                i_ali[n_cut] = i;
                n_cut++;
            }
            if (score_sum_method == 8)
            {
                if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02);
            }
            else
            {
                score_sum += 1 / (1 + di / d02);
            }

            if (di<64)
            {
                GDT_list_tmp[4]+=1;
                if (di<16)
                {
                    GDT_list_tmp[3]+=1;
                    if (di<12.25)
                    {
                        maxsub_tmp+=1/(1+di/12.25);
                        if (di<4)
                        {
                            GDT_list_tmp[2]+=1;
                            if (di<1)
                            {
                                GDT_list_tmp[1]+=1;
                                if (di<0.25)
                                    GDT_list_tmp[0]+=1;
                            }
                        }
                    }
                }
            }
        }
        if (n_cut<3 && n_ali>3)
        {
            inc++;
            double dinc = (d + inc*0.5);
            d_tmp = dinc * dinc;
        }
        else break;
    }

    score1 = score_sum / n_ali;
    return n_cut;
}

double TMscore8_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, int Lali, Vec3& t0, RotMat& u0, int simplify_step,
    int score_sum_method, double &Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, double GDT_list[5], double &maxsub)
{
    double GDT_list_tmp[5]={0,0,0,0,0};
    double maxsub_tmp=0;
    int i;
    int m;
    double score_max;
    double score;
    double rmsd;
    const int kmax=Lali;
    std::vector<int> k_ali(kmax);
    int ka;
    int k;
    Vec3 t;
    RotMat u;
    double d;


    //iterative parameters
    int n_it=20;            //maximum number of iterations
    int n_init_max=6; //maximum number of different fragment length
    std::vector<int> L_ini(n_init_max);  //fragment lengths, Lali, Lali/2, Lali/4 ... 4
    int L_ini_min=4;
    if(Lali<L_ini_min) L_ini_min=Lali;

    int n_init=0;
    int i_init;
    for(i=0; i<n_init_max-1; i++)
    {
        n_init++;
        L_ini[i]=static_cast<int>(Lali/pow(2.0, static_cast<double>(i)));
        if(L_ini[i]<=L_ini_min)
        {
            L_ini[i]=L_ini_min;
            break;
        }
    }
    if(i==n_init_max-1)
    {
        n_init++;
        L_ini[i]=L_ini_min;
    }

    score_max=-1;
    //find the maximum score starting from local structures superposition
    std::vector<int> i_ali(kmax);
    int n_cut;
    int L_frag; //fragment length
    int iL_max; //maximum starting postion for the fragment

    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;

        i=0;
        while(1)
        {
            //extract the fragment starting from position i
            ka=0;
            for(k=0; k<L_frag; k++)
            {
                int kk=k+i;
                r1[k][0]=xtm[kk][0];
                r1[k][1]=xtm[kk][1];
                r1[k][2]=xtm[kk][2];

                r2[k][0]=ytm[kk][0];
                r2[k][1]=ytm[kk][1];
                r2[k][2]=ytm[kk][2];

                k_ali[ka]=kk;
                ka++;
            }

            //extract rotation matrix based on the fragment
            Kabsch(r1, r2, L_frag, 1, rmsd, t, u);
            if (simplify_step != 1)
                Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);

            //get subsegment of this fragment
            d = local_d0_search - 1;
            n_cut=score_fun8(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, Lnorm, score_d8, d0,
                GDT_list_tmp, maxsub_tmp);
            if(score>score_max)
            {
                score_max=score;

                //save the rotation matrix
                for(k=0; k<3; k++)
                {
                    t0[k]=t[k];
                    u0[k][0]=u[k][0];
                    u0[k][1]=u[k][1];
                    u0[k][2]=u[k][2];
                }
            }
            if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
            for (k=0;k<5;k++)
                if (GDT_list_tmp[k]>GDT_list[k])
                    GDT_list[k]=GDT_list_tmp[k];

            //try to extend the alignment iteratively
            d = local_d0_search + 1;
            for(int it=0; it<n_it; it++)
            {
                ka=0;
                for(k=0; k<n_cut; k++)
                {
                    m=i_ali[k];
                    r1[k][0]=xtm[m][0];
                    r1[k][1]=xtm[m][1];
                    r1[k][2]=xtm[m][2];

                    r2[k][0]=ytm[m][0];
                    r2[k][1]=ytm[m][1];
                    r2[k][2]=ytm[m][2];

                    k_ali[ka]=m;
                    ka++;
                }
                //extract rotation matrix based on the fragment
                Kabsch(r1, r2, n_cut, 1, rmsd, t, u);

                do_rotation(xtm, xt, Lali, t, u);
                n_cut=score_fun8(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, Lnorm, score_d8, d0);
                if(score>score_max)
                {
                    score_max=score;

                    //save the rotation matrix
                    for(k=0; k<3; k++)
                    {
                        t0[k]=t[k];
                        u0[k][0]=u[k][0];
                        u0[k][1]=u[k][1];
                        u0[k][2]=u[k][2];
                    }
                }
                if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
                for (k=0;k<5;k++)
                    if (GDT_list_tmp[k]>GDT_list[k])
                        GDT_list[k]=GDT_list_tmp[k];

                //check if it converges
                if(n_cut==ka)
                {
                    for(k=0; k<n_cut; k++)
                    {
                        if(i_ali[k]!=k_ali[k]) break;
                    }
                    if(k==n_cut) break;
                }
            } //for iteration

            if(i<iL_max)
            {
                i=i+simplify_step; //shift the fragment
                if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
            }
            else if(i>=iL_max) break;
        }//while(1)
        //end of one fragment
    }//for(i_init
    return score_max;
}


double TMscore8_search_standard(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, int Lali,
    Vec3& t0, RotMat& u0, int simplify_step, int score_sum_method,
    double &Rcomm, double local_d0_search, double score_d8, double d0,
    double GDT_list[5], double &maxsub)
{
    double GDT_list_tmp[5]={0,0,0,0,0};
    double maxsub_tmp=0;
    int i;
    int m;
    double score_max;
    double score;
    double rmsd;
    const int kmax = Lali;
    std::vector<int> k_ali(kmax);
    int ka;
    int k;
    Vec3 t;
    RotMat u;
    double d;

    //iterative parameters
    int n_it = 20;
    int n_init_max = 6;
    std::vector<int> L_ini(n_init_max);
    int L_ini_min = 4;
    if (Lali<L_ini_min) L_ini_min = Lali;

    int n_init = 0;
    int i_init;
    for (i = 0; i<n_init_max - 1; i++)
    {
        n_init++;
        L_ini[i] = static_cast<int>(Lali / pow(2.0, static_cast<double>(i)));
        if (L_ini[i] <= L_ini_min)
        {
            L_ini[i] = L_ini_min;
            break;
        }
    }
    if (i == n_init_max - 1)
    {
        n_init++;
        L_ini[i] = L_ini_min;
    }

    score_max = -1;
    std::vector<int> i_ali(kmax);
    int n_cut;
    int L_frag;
    int iL_max;

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

        i = 0;
        while (1)
        {
            ka = 0;
            for (k = 0; k<L_frag; k++)
            {
                int kk = k + i;
                r1[k][0] = xtm[kk][0];
                r1[k][1] = xtm[kk][1];
                r1[k][2] = xtm[kk][2];

                r2[k][0] = ytm[kk][0];
                r2[k][1] = ytm[kk][1];
                r2[k][2] = ytm[kk][2];

                k_ali[ka] = kk;
                ka++;
            }

            Kabsch(r1, r2, L_frag, 1, rmsd, t, u);

            if (simplify_step != 1)
                Rcomm = 0;
            do_rotation(xtm, xt, Lali, t, u);

            d = local_d0_search - 1;
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, score_d8, d0, GDT_list_tmp, maxsub_tmp);

            if (score>score_max)
            {
                score_max = score;
                for (k = 0; k<3; k++)
                {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }
            if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
            for (k=0;k<5;k++)
                if (GDT_list_tmp[k]>GDT_list[k])
                    GDT_list[k]=GDT_list_tmp[k];

            d = local_d0_search + 1;
            for (int it = 0; it<n_it; it++)
            {
                ka = 0;
                for (k = 0; k<n_cut; k++)
                {
                    m = i_ali[k];
                    r1[k][0] = xtm[m][0];
                    r1[k][1] = xtm[m][1];
                    r1[k][2] = xtm[m][2];

                    r2[k][0] = ytm[m][0];
                    r2[k][1] = ytm[m][1];
                    r2[k][2] = ytm[m][2];

                    k_ali[ka] = m;
                    ka++;
                }

                Kabsch(r1, r2, n_cut, 1, rmsd, t, u);

                do_rotation(xtm, xt, Lali, t, u);
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, score_d8, d0, GDT_list_tmp, maxsub_tmp);
                if (score>score_max)
                {
                    score_max = score;
                    for (k = 0; k<3; k++)
                    {
                        t0[k] = t[k];
                        u0[k][0] = u[k][0];
                        u0[k][1] = u[k][1];
                        u0[k][2] = u[k][2];
                    }
                }
                if (maxsub_tmp>maxsub) maxsub=maxsub_tmp;
                for (k=0;k<5;k++)
                    if (GDT_list_tmp[k]>GDT_list[k])
                        GDT_list[k]=GDT_list_tmp[k];

                if (n_cut == ka)
                {
                    for (k = 0; k<n_cut; k++)
                    {
                        if (i_ali[k] != k_ali[k]) break;
                    }
                    if (k == n_cut) break;
                }
            }

            if (i<iL_max)
            {
                i = i + simplify_step;
                if (i>iL_max) i = iL_max;
            }
            else if (i >= iL_max) break;
        }
    }
    return score_max;
}


double detailed_search_standard( CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, const CoordArray& x, const CoordArray& y,
    int ylen, std::vector<int>& invmap0, Vec3& t, RotMat& u,
    int simplify_step, int score_sum_method, double local_d0_search,
    const bool& bNormalize, double Lnorm, double score_d8, double d0,
    double GDT_list[5], double &maxsub)
{
    int i;
    int j;
    int k;
    double tmscore;
    double rmsd;

    k=0;
    for(i=0; i<ylen; i++)
    {
        j=invmap0[i];
        if(j>=0) //aligned
        {
            xtm[k][0]=x[j][0];
            xtm[k][1]=x[j][1];
            xtm[k][2]=x[j][2];

            ytm[k][0]=y[i][0];
            ytm[k][1]=y[i][1];
            ytm[k][2]=y[i][2];
            k++;
        }
    }

    //detailed search 40-->1
    tmscore = TMscore8_search_standard( r1, r2, xtm, ytm, xt, k, t, u,
        simplify_step, score_sum_method, rmsd, local_d0_search, score_d8, d0,
        GDT_list, maxsub);
    if (bNormalize)
        tmscore = tmscore * k / Lnorm;
    return tmscore;
}




void output_TMscore_results(
    const std::string xname, const std::string yname,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const Vec3& t, const RotMat& u,
    const double TM1, const double TM2,
    const double TM3, const double TM4, const double TM5,
    const double rmsd, const double d0_out,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int n_ali8,
    const double d0A, const double d0B,
    const double Lnorm_ass, const double d0_scale, 
    const double d0a, const double d0u, const std::string& fname_matrix,
    const int outfmt_opt, const int ter_opt, const std::string& fname_super,
    const int a_opt, const bool u_opt, const bool d_opt, const int mirror_opt,
    int L_lt_d, double GDT_list[5], double maxsub, const int split_opt,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2)
{
    if (outfmt_opt<=0)
    {
        fcout("\nStructure1: %s%s    Length=%5d\n",
            xname, chainID1, xlen);
        fcout("Structure2: %s%s    Length=%5d (by which all scores are normalized)\n",
            yname, chainID2, ylen);

        fcout("Number of residues in common=%5d\n", n_ali8);
        fcout("RMSD of  the common residues=%9.3f\n\n", rmsd);
        fcout("TM-score    = %6.4f  (d0= %.2f)\n", TM1, d0A);
        fcout("MaxSub-score= %6.4f  (d0= 3.50)\n", maxsub/ylen);

        double gdt_ts_score=0;
        double gdt_ha_score=0;
        int i;
        for (i=0;i<4;i++)
        {
            gdt_ts_score+=GDT_list[i+1];
            gdt_ha_score+=GDT_list[i];
        }
        gdt_ts_score/=(4*ylen);
        gdt_ha_score/=(4*ylen);
        fcout("GDT-TS-score= %6.4f %%(d<1)=%6.4f %%(d<2)=%6.4f %%(d<4)=%6.4f %%(d<8)=%6.4f\n",
            gdt_ts_score, GDT_list[1]/ylen, GDT_list[2]/ylen,
                          GDT_list[3]/ylen, GDT_list[4]/ylen);
        fcout("GDT-HA-score= %6.4f %%(d<0.5)=%6.4f %%(d<1)=%6.4f %%(d<2)=%6.4f %%(d<4)=%6.4f\n",
            gdt_ha_score, GDT_list[0]/ylen, GDT_list[1]/ylen,
                          GDT_list[2]/ylen, GDT_list[3]/ylen);

        if (a_opt==1)
            fcout("TM-score    = %5.4f  (if normalized by average length of two structures, i.e., LN= %.1f, d0= %.2f)\n", TM3, (xlen+ylen)*0.5, d0a);
        if (u_opt)
            fcout("TM-score    = %5.4f  (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
        if (d_opt)
            fcout("TM-score    = %5.5f  (if scaled by user-specified d0= %.2f, and LN= %d)\n", TM5, d0_scale, ylen);
    

        std::cout << "\n -------- rotation matrix to rotate Chain-1 to Chain-2 ------\n";
        std::cout << " i          t(i)         u(i,1)         u(i,2)         u(i,3)\n";
        fcout(" 1 %17.10f %14.10f %14.10f %14.10f\n",clean_fmt(t[0]),clean_fmt(u[0][0]),clean_fmt(u[0][1]),clean_fmt(u[0][2]));
        fcout(" 2 %17.10f %14.10f %14.10f %14.10f\n",clean_fmt(t[1]),clean_fmt(u[1][0]),clean_fmt(u[1][1]),clean_fmt(u[1][2]));
        fcout(" 3 %17.10f %14.10f %14.10f %14.10f\n",clean_fmt(t[2]),clean_fmt(u[2][0]),clean_fmt(u[2][1]),clean_fmt(u[2][2]));

        //output alignment
        std::string seq_scale=seqM;
        for (i=0;i<seqM.size();i++)
        {
            L_lt_d+=seqM[i]==':';
            seq_scale[i]=(i+1)%10+'0';
        }
        fcout("\nSuperposition in the TM-score: Length(d<%3.1f)= %d\n", d0_out, L_lt_d);
        fcout("(\":\" denotes the residue pairs of distance <%4.1f Angstrom)\n", d0_out);
        std::cout << seqxA << "\n";
        std::cout << seqM << "\n";
        std::cout << seqyA << "\n";
        std::cout << seq_scale.c_str() << "\n";
        seq_scale.clear();
    }
    else if (outfmt_opt==1)
    {
        fcout(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            xname, chainID1, xlen, d0B, Liden/xlen, TM2);
        std::cout << seqxA << "\n";
        fcout(">%s%s\tL=%d\td0=%.2f\tseqID=%.3f\tTM-score=%.5f\n",
            yname, chainID2, ylen, d0A, Liden/ylen, TM1);
        std::cout << seqyA << "\n";

        fcout("# Lali=%d\tRMSD=%.2f\tseqID_ali=%.3f\n",
            n_ali8, rmsd, (n_ali8>0)?Liden/n_ali8:0);

        if(a_opt)
            fcout("# TM-score=%.5f (normalized by average length of two structures: L=%.1f\td0=%.2f)\n", TM3, (xlen+ylen)*0.5, d0a);

        if(u_opt)
            fcout("# TM-score=%.5f (normalized by user-specified L=%.2f\td0=%.2f)\n", TM4, Lnorm_ass, d0u);

        if(d_opt)
            fcout("# TM-score=%.5f (scaled by user-specified d0=%.2f\tL=%d)\n", TM5, d0_scale, ylen);

        std::cout << "$$$$\n";
    }
    else if (outfmt_opt==2)
    {
        fcout("%s%s\t%s%s\t%.4f\t%.4f\t%.2f\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            TM2, TM1, rmsd, Liden/xlen, Liden/ylen, (n_ali8>0)?Liden/n_ali8:0,
            xlen, ylen, n_ali8);
    }
    std::cout << std::endl;

    if (!fname_matrix.empty()) 
        output_rotation_matrix(fname_matrix, t, u);
    if (!fname_super.empty())
        output_pymol(xname, yname, fname_super, t, u, ter_opt,
            0, split_opt, mirror_opt, seqM, seqxA, seqyA,
            resi_vec1, resi_vec2, chainID1, chainID2);
}

struct TMscoreWorkspace
{
    RotMat u;                     //Kabsch translation vector and rotation matrix
    Vec3 t;
    CoordArray xtm;               // for TMscore search engine
    CoordArray ytm;
    CoordArray xt;                //for saving the superposed version of r_1 or xtm
    CoordArray r1;                // for Kabsch rotation
    CoordArray r2;
    std::vector<int> invmap;
    std::vector<int> invmap0;
    std::vector<int> m1;          //alignd index in x
    std::vector<int> m2;          //alignd index in y
    double D0_MIN;                //for d0
    double Lnorm;                 //normalization length
    double score_d8;              //for TMscore search
    double d0;
    double d0_search;
    double dcu0;
    double local_d0_search;
    double TM;
    double TMmax;
    double rmsd;
    int simplify_step;
    int score_sum_method;
};

int tmscore_initial_alignment(TMscoreWorkspace& ws, const std::vector<std::string>& sequence,
    CoordArray& xa, CoordArray& ya, ChainPairAlignResult& chain_align_res,
    const int mol_type, const int a_opt, const double TMcut)
{
    const int xlen = xa.size();
    const int ylen = ya.size();

    //************************************************//
    //    Stick to the initial alignment              //
    //************************************************//
    for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
        ws.invmap[j] = -1;

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
            else if (sequence[0][kk1] != '-') ws.invmap[i2] = i1;
        }
    }

    //--------------- 2. Align proteins from original alignment
    double prevD0_MIN = ws.D0_MIN;// stored for later use
    int prevLnorm = ws.Lnorm;
    double prevd0 = ws.d0;
    chain_align_res.TM_ali = standard_TMscore(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, xa, ya, ylen,
        ws.invmap, chain_align_res.L_ali, chain_align_res.rmsd_ali, ws.D0_MIN, ws.Lnorm, ws.d0,
        ws.score_d8, ws.t, ws.u, mol_type);
    ws.D0_MIN = prevD0_MIN;
    ws.Lnorm = prevLnorm;
    ws.d0 = prevd0;
    ws.TM = detailed_search_standard(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, xa, ya, ylen,
        ws.invmap, ws.t, ws.u, 40, 8, ws.local_d0_search, true, ws.Lnorm, ws.score_d8, ws.d0);
    if (ws.TM > ws.TMmax)
    {
        ws.TMmax = ws.TM;
        for (int i = 0; i<ylen; i++) ws.invmap0[i] = ws.invmap[i];
    }

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for(int i=0; i<ylen; i++)
    {
        if(ws.invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        std::cout << "There is no alignment between the two structures! "
             << "Program stop with no result!" << std::endl;
        return 1;
    }

    // last TM-score pre-termination
    if (TMcut>0)
    {
        double TMtmp=approx_TM(xlen, ylen, a_opt,
            xa, ya, chain_align_res.t0, chain_align_res.u0, ws.invmap0, mol_type);

        if (TMtmp<0.6*TMcut)
        {
            chain_align_res.TM1=chain_align_res.TM2=chain_align_res.TM3=chain_align_res.TM4=chain_align_res.TM5=TMtmp;
            return 7;
        }
    }

    return 0;
}

void tmscore_detailed_search(TMscoreWorkspace& ws, CoordArray& xa, CoordArray& ya,
    ChainPairAlignResult& chain_align_res, const bool fast_opt,
    double GDT_list[5], double &maxsub)
{
    const int xlen = xa.size();
    const int ylen = ya.size();

    //********************************************************************//
    //    Detailed TMscore search engine --> prepare for final TMscore    //
    //********************************************************************//
    //run detailed TMscore search engine for the best alignment, and
    //extract the best rotation matrix (t, u) for the best alignment
    ws.simplify_step=1;
    if (fast_opt) ws.simplify_step=40;
    ws.score_sum_method=8;
    ws.TM = detailed_search_standard(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, xa, ya, ylen,
        ws.invmap0, ws.t, ws.u, ws.simplify_step, ws.score_sum_method, ws.local_d0_search,
        false, ws.Lnorm, ws.score_d8, ws.d0,
        GDT_list, maxsub);

    //select pairs with dis<d8 for final TMscore computation and output alignment
    int k=0;
    do_rotation(xa, ws.xt, xlen, ws.t, ws.u);
    k=0;
    for(int j=0; j<ylen; j++)
    {
        int i=ws.invmap0[j];
        if(i>=0)//aligned
        {
            chain_align_res.n_ali++;
            ws.m1[k]=i;
            ws.m2[k]=j;

            ws.xtm[k][0]=xa[i][0];
            ws.xtm[k][1]=xa[i][1];
            ws.xtm[k][2]=xa[i][2];

            ws.ytm[k][0]=ya[j][0];
            ws.ytm[k][1]=ya[j][1];
            ws.ytm[k][2]=ya[j][2];

            ws.r1[k][0] = ws.xt[i][0];
            ws.r1[k][1] = ws.xt[i][1];
            ws.r1[k][2] = ws.xt[i][2];
            ws.r2[k][0] = ya[j][0];
            ws.r2[k][1] = ya[j][1];
            ws.r2[k][2] = ya[j][2];

            k++;
        }
    }
    chain_align_res.n_ali8=k;

    Kabsch(ws.r1, ws.r2, chain_align_res.n_ali8, 0, chain_align_res.rmsd0, ws.t, ws.u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
    chain_align_res.rmsd0 = sqrt(chain_align_res.rmsd0 / chain_align_res.n_ali8);
}

void tmscore_final_scores(TMscoreWorkspace& ws, ChainPairAlignResult& chain_align_res,
    const double Lnorm_ass, const double d0_scale, const int xlen, const int ylen,
    const int a_opt, const bool u_opt, const bool d_opt, const int mol_type,
    double GDT_list[5], double &maxsub)
{
    //****************************************//
    //              Final TMscore             //
    //    Please set parameters for output    //
    //****************************************//
    double rmsd;
    ws.simplify_step=1;
    ws.score_sum_method=0;
    double Lnorm_0=ylen;


    //normalized by length of structure A
    parameter_set4final(Lnorm_0, ws.D0_MIN, ws.Lnorm, ws.d0, ws.d0_search, mol_type);
    chain_align_res.d0A=ws.d0;
    chain_align_res.d0_0=chain_align_res.d0A;
    ws.local_d0_search = ws.d0_search;
    chain_align_res.TM1 = TMscore8_search(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, chain_align_res.n_ali8, chain_align_res.t0, chain_align_res.u0, ws.simplify_step,
        ws.score_sum_method, rmsd, ws.local_d0_search, ws.Lnorm, ws.score_d8, ws.d0,
        GDT_list, maxsub);
    chain_align_res.TM_0 = chain_align_res.TM1;

    if (a_opt>0)
    {
        //normalized by average length of structures A, B
        Lnorm_0=(xlen+ylen)*0.5;
        parameter_set4final(Lnorm_0, ws.D0_MIN, ws.Lnorm, ws.d0, ws.d0_search, mol_type);
        chain_align_res.d0a=ws.d0;
        chain_align_res.d0_0=chain_align_res.d0a;
        ws.local_d0_search = ws.d0_search;

        chain_align_res.TM3 = TMscore8_search(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, chain_align_res.n_ali8, chain_align_res.t0, chain_align_res.u0,
            ws.simplify_step, ws.score_sum_method, rmsd, ws.local_d0_search, ws.Lnorm,
            ws.score_d8, ws.d0);
        chain_align_res.TM_0=chain_align_res.TM3;
    }
    if (u_opt)
    {
        //normalized by user assigned length
        parameter_set4final(Lnorm_ass, ws.D0_MIN, ws.Lnorm,
            ws.d0, ws.d0_search, mol_type);
        chain_align_res.d0u=ws.d0;
        chain_align_res.d0_0=chain_align_res.d0u;
        Lnorm_0=Lnorm_ass;
        ws.local_d0_search = ws.d0_search;
        chain_align_res.TM4 = TMscore8_search(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, chain_align_res.n_ali8, chain_align_res.t0, chain_align_res.u0,
            ws.simplify_step, ws.score_sum_method, rmsd, ws.local_d0_search, ws.Lnorm,
            ws.score_d8, ws.d0);
        chain_align_res.TM_0=chain_align_res.TM4;
    }
    if (d_opt)
    {
        //scaled by user assigned d0
        parameter_set4scale(ylen, d0_scale, ws.Lnorm, ws.d0, ws.d0_search);
        chain_align_res.d0_out=d0_scale;
        chain_align_res.d0_0=d0_scale;
        //Lnorm_0=ylen;
        ws.local_d0_search = ws.d0_search;
        chain_align_res.TM5 = TMscore8_search(ws.r1, ws.r2, ws.xtm, ws.ytm, ws.xt, chain_align_res.n_ali8, chain_align_res.t0, chain_align_res.u0,
            ws.simplify_step, ws.score_sum_method, rmsd, ws.local_d0_search, ws.Lnorm,
            ws.score_d8, ws.d0);
        chain_align_res.TM_0=chain_align_res.TM5;
    }
}

void tmscore_derive_alignment(TMscoreWorkspace& ws, CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy, ChainPairAlignResult& chain_align_res)
{
    const int xlen = xa.size();
    const int ylen = ya.size();

    // derive alignment from superposition
    int ali_len=xlen+ylen; //maximum length of alignment
    chain_align_res.seqxA.assign(ali_len,'-');
    chain_align_res.seqM.assign( ali_len,' ');
    chain_align_res.seqyA.assign(ali_len,'-');

    do_rotation(xa, ws.xt, xlen, chain_align_res.t0, chain_align_res.u0);

    int kk=0;
    int i_old=0;
    int j_old=0;
    double d=0;
    for(int k=0; k<chain_align_res.n_ali8; k++)
    {
        for(int i=i_old; i<ws.m1[k]; i++)
        {
            //align x to gap
            chain_align_res.seqxA[kk]=seqx[i];
            chain_align_res.seqyA[kk]='-';
            chain_align_res.seqM[kk]=' ';                    
            kk++;
        }

        for(int j=j_old; j<ws.m2[k]; j++)
        {
            //align y to gap
            chain_align_res.seqxA[kk]='-';
            chain_align_res.seqyA[kk]=seqy[j];
            chain_align_res.seqM[kk]=' ';
            kk++;
        }

        chain_align_res.seqxA[kk]=seqx[ws.m1[k]];
        chain_align_res.seqyA[kk]=seqy[ws.m2[k]];
        chain_align_res.Liden+=(chain_align_res.seqxA[kk]==chain_align_res.seqyA[kk]);
        d=sqrt(dist(ws.xt[ws.m1[k]], ya[ws.m2[k]]));
        //if(d<d0_out) seqM[kk]=':';
        //else         seqM[kk]='.';
        if(d<5) chain_align_res.seqM[kk]=':';
        kk++;
        i_old=ws.m1[k]+1;
        j_old=ws.m2[k]+1;
    }

    //tail
    for(int i=i_old; i<xlen; i++)
    {
        //align x to gap
        chain_align_res.seqxA[kk]=seqx[i];
        chain_align_res.seqyA[kk]='-';
        chain_align_res.seqM[kk]=' ';
        kk++;
    }
    for(int j=j_old; j<ylen; j++)
    {
        //align y to gap
        chain_align_res.seqxA[kk]='-';
        chain_align_res.seqyA[kk]=seqy[j];
        chain_align_res.seqM[kk]=' ';
        kk++;
    }
    chain_align_res.seqxA=chain_align_res.seqxA.substr(0,kk);
    chain_align_res.seqyA=chain_align_res.seqyA.substr(0,kk);
    chain_align_res.seqM =chain_align_res.seqM.substr(0,kk);
}

int TMscore_main(CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    ChainPairAlignResult& chain_align_res,
    const std::vector<std::string> sequence, const double Lnorm_ass,
    const double d0_scale, const int a_opt,
    const bool u_opt, const bool d_opt, const bool fast_opt,
    const int mol_type, double GDT_list[5], double &maxsub,
    const double TMcut)
{
    const int xlen = xa.size();
    const int ylen = ya.size();

    /***********************/
    // allocate memory
    /***********************/
    TMscoreWorkspace ws;
    int minlen = std::min(xlen, ylen);
    ws.xtm.resize(minlen);
    ws.ytm.resize(minlen);
    ws.xt.resize(xlen);
    ws.r1.resize(minlen);
    ws.r2.resize(minlen);
    ws.invmap.resize(ylen+1);
    ws.invmap0.assign(ylen+1, -1);
    ws.m1.resize(xlen);
    ws.m2.resize(ylen);

    /***********************/
    //    parameter set   
    /***********************/
    parameter_set4search(xlen, ylen, ws.D0_MIN, ws.Lnorm,
        ws.score_d8, ws.d0, ws.d0_search, ws.dcu0);
    ws.simplify_step    = 40; //for simplified search engine
    ws.score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8
    ws.TMmax = -1;
    ws.local_d0_search = ws.d0_search;

    int rc = tmscore_initial_alignment(ws, sequence, xa, ya, chain_align_res,
        mol_type, a_opt, TMcut);
    if (rc) return rc;

    tmscore_detailed_search(ws, xa, ya, chain_align_res, fast_opt, GDT_list, maxsub);

    tmscore_final_scores(ws, chain_align_res, Lnorm_ass, d0_scale, xlen, ylen,
        a_opt, u_opt, d_opt, mol_type, GDT_list, maxsub);

    tmscore_derive_alignment(ws, xa, ya, seqx, seqy, chain_align_res);

    return 0; // zero for no exception
}
