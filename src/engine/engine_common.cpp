#include "engine_common.h"
#include "se.h"

void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    //parameter initialization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    D0_MIN=0.5;
    dcu0=4.25;                       //update 3.85-->4.25

    Lnorm=getmin(xlen, ylen);        //normalize TMscore by this in searching
    if (Lnorm<=19)                    //update 15-->19
        d0=0.168;                   //update 0.5-->0.168
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search)
{
    D0_MIN=0.3;

    Lnorm=len;            //normalize TMscore by this in searching
    if(Lnorm<=11) d0=0.3;
    else if(Lnorm>11&&Lnorm<=15) d0=0.4;
    else if(Lnorm>15&&Lnorm<=19) d0=0.5;
    else if(Lnorm>19&&Lnorm<=23) d0=0.6;
    else if(Lnorm>23&&Lnorm<30)  d0=0.7;
    else d0=(0.6*pow((Lnorm*1.0-0.5), 1.0/2)-2.5);

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type)
{
    if (mol_type>0) // RNA
    {
        parameter_set4final_C3prime(len, D0_MIN, Lnorm,
            d0, d0_search);
        return;
    }
    D0_MIN=0.5;

    Lnorm=len;            //normalize TMscore by this in searching
    if (Lnorm<=21) d0=0.5;
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search)
{
    d0=d_s;
    Lnorm=len;            //normalize TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void init_gotoh_mat(GotohMatrix& m, const int xlen, const int ylen, const int gapopen,
    const int gapext, const int glocal, const int alt_init)
{
    IntMatrix& S = m.S;
    IntMatrix& JumpH = m.JumpH;
    IntMatrix& JumpV = m.JumpV;
    IntMatrix& P = m.P;
    IntMatrix& H = m.H;
    IntMatrix& V = m.V;
    // fill first row/colum of JumpH,jumpV and path matrix P
    int i;
    int j;
    for (i=0;i<xlen+1;i++)
        for (j=0;j<ylen+1;j++)
            H[i][j]=V[i][j]=P[i][j]=JumpH[i][j]=JumpV[i][j]=0;
    for (i=0;i<xlen+1;i++)
    {
        if (glocal<2) P[i][0]=4; // -
        JumpV[i][0]=i;
    }
    for (j=0;j<ylen+1;j++)
    {
        if (glocal<1) P[0][j]=2; // |
        JumpH[0][j]=j;
    }
    if (glocal<2) for (i=1;i<xlen+1;i++) S[i][0]=gapopen+gapext*(i-1);
    if (glocal<1) for (j=1;j<ylen+1;j++) S[0][j]=gapopen+gapext*(j-1);
    if (alt_init==0)
    {
        for (i=1;i<xlen+1;i++) H[i][0]=gapopen+gapext*(i-1);
        for (j=1;j<ylen+1;j++) V[0][j]=gapopen+gapext*(j-1);
    }
    else
    {
        if (glocal<2) for (i=1;i<xlen+1;i++) V[i][0]=gapopen+gapext*(i-1);
        if (glocal<1) for (j=1;j<ylen+1;j++) H[0][j]=gapopen+gapext*(j-1);
        for (i=0;i<xlen+1;i++) H[i][0]=-99999; // INT_MIN cause bug on ubuntu
        for (j=0;j<ylen+1;j++) V[0][j]=-99999; // INT_MIN;
    }
}

void find_highest_align_score( GotohMatrix& m,
    int &aln_score, const int xlen,const int ylen)
{
    const IntMatrix& S = m.S;
    IntMatrix& P = m.P;
    // locate the cell with highest alignment score
    int max_aln_i=xlen;
    int max_aln_j=ylen;
    int i;
    int j;
    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (S[i][j]>=aln_score)
            {
                max_aln_i=i;
                max_aln_j=j;
                aln_score=S[i][j];
            }
        }
    }

    // reset all path after [max_aln_i][max_aln_j]
    for (i=max_aln_i+1;i<xlen+1;i++) for (j=0;j<ylen+1;j++) P[i][j]=0;
    for (i=0;i<xlen+1;i++) for (j=max_aln_j+1;j<ylen+1;j++) P[i][j]=0;
}

int calculate_score_gotoh(const int xlen,const int ylen, GotohMatrix& m, const int gapopen,const int gapext,
    const int glocal, const int alt_init)
{
    IntMatrix& S = m.S;
    IntMatrix& JumpH = m.JumpH;
    IntMatrix& JumpV = m.JumpV;
    IntMatrix& P = m.P;
    IntMatrix& H = m.H;
    IntMatrix& V = m.V;

    // fill first row/colum of JumpH,jumpV and path matrix P
    int i;
    int j;
    init_gotoh_mat(m, xlen, ylen,
        gapopen, gapext, glocal, alt_init);

    // fill S and P
    int diag_score;
    int left_score;
    int up_score;
    for (i=1;i<xlen+1;i++)
    {
        for (j=1;j<ylen+1;j++)
        {
            // penalty of consective deletion
            if (glocal<1 || i<xlen || glocal>=3)
            {
                H[i][j]=std::max(S[i][j-1]+gapopen,H[i][j-1]+gapext);
                JumpH[i][j]=(H[i][j]==H[i][j-1]+gapext)?(JumpH[i][j-1]+1):1;
            }
            else
            {
                H[i][j]=std::max(S[i][j-1],H[i][j-1]);
                JumpH[i][j]=(H[i][j]==H[i][j-1])?(JumpH[i][j-1]+1):1;
            }
            // penalty of consective insertion
            if (glocal<2 || j<ylen || glocal>=3)
            {
                V[i][j]=std::max(S[i-1][j]+gapopen,V[i-1][j]+gapext);
                JumpV[i][j]=(V[i][j]==V[i-1][j]+gapext)?(JumpV[i-1][j]+1):1;
            }
            else
            {
                V[i][j]=std::max(S[i-1][j],V[i-1][j]);
                JumpV[i][j]=(V[i][j]==V[i-1][j])?(JumpV[i-1][j]+1):1;
            }

            diag_score=S[i-1][j-1]+S[i][j]; // match-mismatch '\'
            left_score=H[i][j];     // deletion       '-'
            up_score  =V[i][j];     // insertion      '|'

            if (diag_score>=left_score && diag_score>=up_score)
            {
                S[i][j]=diag_score;
                P[i][j]+=1;
            }
            if (up_score>=diag_score && up_score>=left_score)
            {
                S[i][j]=up_score;
                P[i][j]+=2;
            }
            if (left_score>=diag_score && left_score>=up_score)
            {
                S[i][j]=left_score;
                P[i][j]+=4;
            }
            if (glocal>=3 && S[i][j]<0)
            {
                S[i][j]=0;
                P[i][j]=0;
                H[i][j]=0;
                V[i][j]=0;
                JumpH[i][j]=0;
                JumpV[i][j]=0;
            }
        }
    }
    int aln_score=S[xlen][ylen];

    // re-fill first row/column of path matrix P for back-tracing
    for (i=1;i<xlen+1;i++) if (glocal<3 || P[i][0]>0) P[i][0]=2; // |
    for (j=1;j<ylen+1;j++) if (glocal<3 || P[0][j]>0) P[0][j]=4; // -

    // calculate alignment score and alignment path for swalign
    if (glocal>=3)
        find_highest_align_score(m,aln_score,xlen,ylen);

    // release memory
    return aln_score; // final alignment score
}

int NWalign_main(const std::string &seqx, const std::string &seqy, const int xlen,
    const int ylen, std::string & seqxA, std::string & seqyA, const int mol_type,
    std::vector<int>& invmap, const int invmap_only, const int glocal)
{
    GotohMatrix m;
    m.init(xlen, ylen);
    IntMatrix& S = m.S;

    int aln_score;
    int gapopen=gapopen_blosum62;
    int gapext =gapext_blosum62;
    int i;
    int j;
    if (mol_type>0) // RNA or DNA
    {
        gapopen=gapopen_blastn;
        gapext =gapext_blastn;
        if (glocal==3)
        {
            gapopen=-5;
            gapext =-2;
        }
    }

    for (i=0;i<xlen+1;i++)
    {
        for (j=0;j<ylen+1;j++)
        {
            if (i*j==0) S[i][j]=0;
            else S[i][j]=BLOSUM[seqx[i-1]][seqy[j-1]];
        }
    }

    aln_score=calculate_score_gotoh(xlen, ylen, m,
        gapopen, gapext, glocal);

    seqxA.clear();
    seqyA.clear();

    if (glocal<3) trace_back_gotoh(seqx, seqy, m,
            seqxA, seqyA, xlen, ylen, invmap, invmap_only);
    else trace_back_sw(seqx, seqy, m, seqxA, seqyA,
            xlen, ylen, invmap, invmap_only);

    return aln_score; // aligment score
}

void get_seqID(const std::vector<int>& invmap, const std::string& seqx, const std::string& seqy,
    const int ylen, double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    int i;
    int j;
    for (j=0;j<ylen;j++)
    {
        i=invmap[j];
        if (i<0) continue;
        L_ali+=1;
        Liden+=(seqx[i]==seqy[j]);
    }
}

void get_seqID(const std::string& seqxA, const std::string& seqyA,
    std::string &seqM,double &Liden,int &L_ali)
{
    Liden=0;
    L_ali=0;
    for (int i=0;i<seqxA.size();i++)
    {
        if (seqxA[i]==seqyA[i] && seqxA[i]!='-')
        {
            Liden++;
            seqM+=':';
        }
        else seqM+=' ';
        L_ali+=(seqxA[i]!='-' && seqyA[i]!='-');
    }
}

void output_NWalign_results(
    const std::string xname, const std::string yname,
    const std::string& chainID1, const std::string& chainID2,
    const int xlen, const int ylen, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int L_ali, const int aln_score, const int outfmt_opt)
{
    if (outfmt_opt<=0)
    {
        fcout("\nName of Chain_1: %s%s\n", xname, chainID1);
        fcout("Name of Chain_2: %s%s\n", yname, chainID2);
        fcout("Length of Chain_1: %d residues\n", xlen);
        fcout("Length of Chain_2: %d residues\n\n", ylen);

        fcout("Aligned length= %d, Alignment score= %d, Seq_ID=n_identical/n_aligned= %4.3f\n",
            L_ali, aln_score, Liden/L_ali);
        fcout("Seq_ID= %6.5f (if normalized by length of Chain_1\n", Liden/xlen);
        fcout("Seq_ID= %6.5f (if normalized by length of Chain_2\n", Liden/ylen);
        std::cout << "(You should use Seq_ID normalized by length of the reference structure)\n";

        //output alignment
        fcout("\n(\":\" denotes pairs with identical residue type)\n");
        std::cout << seqxA << "\n";
        std::cout << seqM << "\n";
        std::cout << seqyA << "\n";
    }
    else if (outfmt_opt==1)
    {
        fcout(">%s%s\tL=%d\tseqID=%.3f\n",
            xname, chainID1, xlen, Liden/xlen);
        std::cout << seqxA << "\n";
        fcout(">%s%s\tL=%d\tseqID=%.3f\n",
            yname, chainID2, ylen, Liden/ylen);
        std::cout << seqyA << "\n";
        fcout("#score=%d\tLali=%d\tseqID_ali=%.3f\n", aln_score, L_ali, Liden/L_ali);
        std::cout << "$$$$\n";
    }
    else if (outfmt_opt==2)
    {
        fcout("%s%s\t%s%s\t%4.3f\t%4.3f\t%4.3f\t%d\t%d\t%d",
            xname, chainID1, yname, chainID2,
            Liden/xlen, Liden/ylen, Liden/L_ali,
            xlen, ylen, L_ali);
    }
    std::cout << std::endl;
}

int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=resi_vec1.size();
    int ylen=resi_vec2.size();
    if (byresi_opt==4 || byresi_opt==5 || byresi_opt==7) // global or glocal sequence alignment
    {
        std::vector<int> invmap(ylen+1);
        int glocal=0;
        if (byresi_opt==5 || byresi_opt==7) glocal=2;
        int mol_type=0;
        for (i1=0;i1<xlen;i1++)
            if ('a'<seqx[i1] && seqx[i1]<'z') mol_type++;
            else mol_type--;
        for (i2=0;i2<ylen;i2++)
            if ('a'<seqy[i2] && seqy[i2]<'z') mol_type++;
            else mol_type--;
        NWalign_main(seqx, seqy, xlen, ylen, sequence[0],sequence[1],
            mol_type, invmap, 0, glocal);
    }

    std::map<std::string,std::string> chainID_map1;
    std::map<std::string,std::string> chainID_map2;
    if (byresi_opt==3)
    {
        std::vector<std::string> chainID_vec;
        std::string chainID;
        std::ostringstream ss;
        int i;
        for (i=0;i<xlen;i++)
        {
            chainID=resi_vec1[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map1[chainID]=ss.str();
                ss.str("");
            }
        }
        chainID_vec.clear();
        for (i=0;i<ylen;i++)
        {
            chainID=resi_vec2[i].substr(5);
            if (!chainID_vec.size()|| chainID_vec.back()!=chainID)
            {
                chainID_vec.push_back(chainID);
                ss<<chainID_vec.size();
                chainID_map2[chainID]=ss.str();
                ss.str("");
            }
        }
    }
    std::string chainID1="";
    std::string chainID2="";
    std::string chainID1_prev="";
    std::string chainID2_prev="";
    while(i1<xlen && i2<ylen)
    {
        if (byresi_opt==2)
        {
            chainID1=resi_vec1[i1].substr(5);
            chainID2=resi_vec2[i2].substr(5);
        }
        else if (byresi_opt==3)
        {
            chainID1=chainID_map1[resi_vec1[i1].substr(5)];
            chainID2=chainID_map2[resi_vec2[i2].substr(5)];
        }

        if (chainID1==chainID2)
        {
            if (safe_stoi(resi_vec1[i1].substr(0,4).c_str())<
                safe_stoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
            }
            else if (safe_stoi(resi_vec1[i1].substr(0,4).c_str())>
                     safe_stoi(resi_vec2[i2].substr(0,4).c_str()))
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
            }
            chainID1_prev=chainID1;
            chainID2_prev=chainID2;
        }
        else
        {
            if (chainID1_prev==chainID1 && chainID2_prev!=chainID2)
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+='-';
                chainID1_prev=chainID1;
            }
            else if (chainID1_prev!=chainID1 && chainID2_prev==chainID2)
            {
                sequence[0]+='-';
                sequence[1]+=seqy[i2++];
                chainID2_prev=chainID2;
            }
            else
            {
                sequence[0]+=seqx[i1++];
                sequence[1]+=seqy[i2++];
                chainID1_prev=chainID1;
                chainID2_prev=chainID2;
            }
        }

    }
    return sequence[0].size();
}

int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const std::vector<int> xlen_vec, const std::vector<int> ylen_vec,
    const int chain_i, const int chain_j, const int byresi_opt)
{
    sequence.clear();
    sequence.push_back("");
    sequence.push_back("");

    int i1=0; // positions in resi_vec1
    int i2=0; // positions in resi_vec2
    int xlen=xlen_vec[chain_i];
    int ylen=ylen_vec[chain_j];
    int i;
    int j;
    for (i=0;i<chain_i;i++) i1+=xlen_vec[i];
    for (j=0;j<chain_j;j++) i2+=ylen_vec[j];

    if (byresi_opt==7)
    {
        std::vector<int> invmap(ylen+1);
        int glocal=2;
        int mol_type=0;

        for (i=0;i<xlen;i++)
            if ('a'<seqx[i] && seqx[i]<'z') mol_type++;
            else mol_type--;
        for (i=0;i<ylen;i++)
            if ('a'<seqy[i] && seqy[i]<'z') mol_type++;
            else mol_type--;
        NWalign_main(seqx, seqy, xlen, ylen, sequence[0],sequence[1],
            mol_type, invmap, 0, glocal);
        return sequence[0].size();
    }

    i=j=0;
    while(i<xlen && j<ylen)
    {
        if (safe_stoi(resi_vec1[i+i1].substr(0,4).c_str())<
            safe_stoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+=seqx[i++];
            sequence[1]+='-';
        }
        else if (safe_stoi(resi_vec1[i+i1].substr(0,4).c_str())>
                 safe_stoi(resi_vec2[j+i2].substr(0,4).c_str()))
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j++];
        }
        else
        {
            sequence[0]+=seqx[i++];
            sequence[1]+=seqy[j++];
        }
    }
    if (i<xlen && j==ylen)
    {
        for (;i<xlen;i++)
        {
            sequence[0]+=seqx[i];
            sequence[1]+='-';
        }
    }
    else if (i==xlen && j<ylen)
    {
        for (;j<ylen;j++)
        {
            sequence[0]+='-';
            sequence[1]+=seqy[j];
        }
    }
    return sequence[0].size();
}


void align_chain_pair(ChainPairAlignResult& chain_align_res,
    CoordArray& xa, CoordArray& ya,
    const std::string& seqx, const std::string& seqy,
    const std::string& secx, const std::string& secy,
    int xlen, int ylen,
    const ChainPairAlignOptions& align_opts,
    const std::vector<std::string>& sequence,
    int outfmt_opt)
{
    if (align_opts.cp_opt)
    {
        CPalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, align_opts);
    }
    else if (align_opts.se_opt)
    {
        chain_align_res.invmap.assign(ylen + 1, 0);
        init_rt_mat(chain_align_res.t0, chain_align_res.u0);

        se_main(xa, ya, seqx, seqy, chain_align_res, xlen, ylen, sequence, align_opts, outfmt_opt);
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
    else
    {
        TMalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, align_opts);
    }
}
