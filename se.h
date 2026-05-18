#include "TMalign.h"

/* entry function for se
 * outfmt_opt>=2 should not parse sequence alignment 
 * u_opt corresponds to option -L
 *       if u_opt==2, use d0 from Lnorm_ass for alignment
 * if hinge>0, append to original invmap */
// C++ string overload (forward bridge)
// NOTE: This overload must NOT be called from cross-scope call sites where seqx
// and seqy are declared in different stack frames (e.g., seqy in outer for(iter),
// seqx in inner for(tm_idx)).  The const std::string& vs const char* parameter
// difference changes se_main's stack frame layout, causing a crash in deeply
// nested contexts like mTMalign hinge recovery.  Callers in such contexts must
// pass seqx.c_str() / seqy.c_str() explicitly to reach the const char* version.
int se_main(
    double **xa, double **ya, const std::string &seqx, const std::string &seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA, vector<double> &do_vec,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen, const vector<string> &sequence,
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool a_opt, const int u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, int *invmap, const int hinge=0)
{
    return se_main(xa, ya, seqx.c_str(), seqy.c_str(),
        TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA, do_vec,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence,
        Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt, mol_type,
        outfmt_opt, invmap, hinge);
}

int se_main(
    double **xa, double **ya, const char *seqx, const char *seqy,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0,
    double &d0A, double &d0B, double &d0u, double &d0a, double &d0_out,
    string &seqM, string &seqxA, string &seqyA, vector<double> &do_vec,
    double &rmsd0, int &L_ali, double &Liden,
    double &TM_ali, double &rmsd_ali, int &n_ali, int &n_ali8,
    const int xlen, const int ylen, const vector<string> &sequence,
    const double Lnorm_ass, const double d0_scale, const bool i_opt,
    const bool a_opt, const int u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, int *invmap, const int hinge=0)
{
    double D0_MIN;        //for d0
    double Lnorm;         //normalization length
    double score_d8;
    double d0;
    double d0_search;
    double dcu0; //for TMscore search
    double **score;       // Input score table for dynamic programming
    bool   **path;        // for dynamic programming  
    double **val;         // for dynamic programming  

    int *m1=nullptr;
    int *m2=nullptr;
    double d;
    if (outfmt_opt<2)
    {
        m1=new int[xlen]; //alignd index in x
        m2=new int[ylen]; //alignd index in y
    }

    /***********************/
    // allocate memory    
    /***********************/
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path, xlen+1, ylen+1);
    NewArray(&val, xlen+1, ylen+1);
    int *invmap0          = new int[ylen+1];
    int i;
    int j;
    if (hinge==0) for (j=0;j<=ylen;j++) invmap0[j]=-1;
    else for (j=0;j<ylen;j++) invmap0[j]=invmap[j];
    vector<char> seqM_char;
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
    if (hinge==0) for(j=0; j<ylen; j++) invmap[j]=-1;
    if (!i_opt) NWDP_SE(path, val, xa, ya, xlen, ylen, d0*d0, 0, invmap, hinge);
    else
    {
        int i1 = -1;// in C version, index starts from zero, not from one
        int i2 = -1;
        int L1 = sequence[0].size();
        int L2 = sequence[1].size();
        int L = min(L1, L2);// Get positions for aligned residues
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
            d=sqrt(dist(&xa[i][0], &ya[j][0]));
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

    if (outfmt_opt>=2)
    {
        if (hinge) seqM_char.clear();    
        delete []invmap0;
        DeleteArray(&score, xlen+1);
        DeleteArray(&path, xlen+1);
        DeleteArray(&val, xlen+1);
        return 0;
    }

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
        d=sqrt(dist(&xa[m1[k]][0], &ya[m2[k]][0]));
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

    // free memory
    delete [] invmap0;
    delete [] m1;
    delete [] m2;
    DeleteArray(&score, xlen+1);
    DeleteArray(&path, xlen+1);
    DeleteArray(&val, xlen+1);
    return 0; // zero for no exception
}
