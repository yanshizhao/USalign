#include "TMalign.h"

int score_fun8(const CoordArray& xa, const CoordArray& ya, int n_ali, double d, std::vector<int>& i_ali,
    double &score1, int score_sum_method, const double Lnorm,
    const double score_d8, const double d0){
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
    double score_d8, double d0){
    double score_sum = 0;
    double di, d_tmp = d*d, d02 = d0*d0, score_d8_cut = score_d8*score_d8;
    int i, n_cut, inc = 0;
    while (1)
    {
        n_cut = 0; score_sum = 0;
        for (i = 0; i<n_ali; i++)
        {
            di = dist(xa[i], ya[i]);
            if (di<d_tmp) { i_ali[n_cut] = i; n_cut++; }
            if (score_sum_method == 8)
            { if (di <= score_d8_cut) score_sum += 1 / (1 + di / d02); }
            else score_sum += 1 / (1 + di / d02);
        }
        if (n_cut<3 && n_ali>3) { inc++; double dinc = (d + inc*0.5); d_tmp = dinc * dinc; }
        else break;
    }
    score1 = score_sum / n_ali;
    return n_cut;
}
double get_score_fast( CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    const CoordArray& x, const CoordArray& y, int ylen, std::vector<int>& invmap,
    double d0, double d0_search, Vec3& t, RotMat& u){
    double rms,tmscore,tmscore1,tmscore2;
    int i,j,k;
    k=0;
    for(j=0; j<ylen; j++) {
        i=invmap[j];
        if(i>=0) {
            r1[k][0]=x[i][0]; r1[k][1]=x[i][1]; r1[k][2]=x[i][2];
            r2[k][0]=y[j][0]; r2[k][1]=y[j][1]; r2[k][2]=y[j][2];
            xtm[k][0]=x[i][0]; xtm[k][1]=x[i][1]; xtm[k][2]=x[i][2];
            ytm[k][0]=y[j][0]; ytm[k][1]=y[j][1]; ytm[k][2]=y[j][2];
            k++;
        }
        else if(i!=-1) PrintErrorAndQuit("Wrong map!\n");
    }

    Kabsch(r1, r2, k, 1, rms, t, u);
    
    double di; const int len=k;
    std::vector<double> dis(len);
    double d00=d0_search,d002=d00*d00,d02=d0*d0;
    int n_ali=k;
    Vec3 xrot;
    tmscore=0;
    for(k=0; k<n_ali; k++) {
        transform(t, u, xtm[k], xrot);
        di=dist(xrot, ytm[k]); dis[k]=di;
        tmscore += 1/(1+di/d02);
    }
    double d002t=d002;
    std::vector<double> dis_vec(dis.begin(), dis.begin()+n_ali);
    std::sort(dis_vec.begin(), dis_vec.end());
    if (d002t<dis_vec[2]) d002t=dis_vec[2];
    dis_vec.clear();
    while(1) {
        j=0;
        for(k=0; k<n_ali; k++) {
            if(dis[k]<=d002t) {
                r1[j][0]=xtm[k][0]; r1[j][1]=xtm[k][1]; r1[j][2]=xtm[k][2];
                r2[j][0]=ytm[k][0]; r2[j][1]=ytm[k][1]; r2[j][2]=ytm[k][2];
                j++;
            }
        }
        if(j<3 && n_ali>3) d002t += 0.5; else break;
    }
    if(n_ali!=j) {
    
        Kabsch(r1, r2, j, 1, rms, t, u);
    
        tmscore1=0;
        for(k=0; k<n_ali; k++) {
            transform(t, u, xtm[k], xrot);
            di=dist(xrot, ytm[k]); dis[k]=di;
            tmscore1 += 1/(1+di/d02);
        }
        d002t=d002+1;
        std::vector<double> dis_vec2(dis.begin(), dis.begin()+n_ali);
        std::sort(dis_vec2.begin(), dis_vec2.end());
        if (d002t<dis_vec2[2]) d002t=dis_vec2[2];
        dis_vec2.clear();
        while(1) {
            j=0;
            for(k=0; k<n_ali; k++) {
                if(dis[k]<=d002t) {
                    r1[j][0]=xtm[k][0]; r1[j][1]=xtm[k][1]; r1[j][2]=xtm[k][2];
                    r2[j][0]=ytm[k][0]; r2[j][1]=ytm[k][1]; r2[j][2]=ytm[k][2];
                    j++;
                }
            }
            if(j<3 && n_ali>3) d002t += 0.5; else break;
        }
        
        Kabsch(r1, r2, j, 1, rms, t, u);
        
        tmscore2=0;
        for(k=0; k<n_ali; k++) {
            transform(t, u, xtm[k], xrot);
            di=dist(xrot, ytm[k]);
            tmscore2 += 1/(1+di/d02);
        }
    } else { tmscore1=tmscore; tmscore2=tmscore; }
    if(tmscore1>=tmscore) tmscore=tmscore1;
    if(tmscore2>=tmscore) tmscore=tmscore2;
    return tmscore;
}
char sec_str(double dis13, double dis14, double dis15,
            double dis24, double dis25, double dis35){
    char s='C';
    
    double delta=2.1;
    if (fabs(dis15-6.37)<delta && fabs(dis14-5.18)<delta && 
        fabs(dis25-5.18)<delta && fabs(dis13-5.45)<delta &&
        fabs(dis24-5.45)<delta && fabs(dis35-5.45)<delta)
    {
        s='H'; //helix                        
        return s;
    }

    delta=1.42;
    if (fabs(dis15-13  )<delta && fabs(dis14-10.4)<delta &&
        fabs(dis25-10.4)<delta && fabs(dis13-6.1 )<delta &&
        fabs(dis24-6.1 )<delta && fabs(dis35-6.1 )<delta)
    {
        s='E'; //strand
        return s;
    }

    if (dis15 < 8) s='T'; //turn
    return s;
}
void make_sec(const CoordArray& x, int len, std::string& sec){
    int j1;
    int j2;
    int j3;
    int j4;
    int j5;
    double d13;
    double d14;
    double d15;
    double d24;
    double d25;
    double d35;
    sec.assign(len, 'C');
    for(int i=0; i<len; i++)
    {
        j1=i-2;
        j2=i-1;
        j3=i;
        j4=i+1;
        j5=i+2;

        if(j1>=0 && j5<len)
        {
            d13=sqrt(dist(x[j1], x[j3]));
            d14=sqrt(dist(x[j1], x[j4]));
            d15=sqrt(dist(x[j1], x[j5]));
            d24=sqrt(dist(x[j2], x[j4]));
            d25=sqrt(dist(x[j2], x[j5]));
            d35=sqrt(dist(x[j3], x[j5]));
            sec[i]=sec_str(d13, d14, d15, d24, d25, d35);
        }
    }
}
void sec_str(int len, const std::vector<std::vector<bool> >&bp,
    int a, int b,int &c, int &d){
    int i;
    
    for (i=0;i<len;i++)
    {
        if (a+i<len-3 && b-i>0)
        {
            if (a+i<b-i && bp[a+i][b-i]) continue;
            break;
        }
    }
    c=a+i-1;d=b-i+1;
}
void make_sec(const std::string& seq, const CoordArray& x, int len, std::string& sec, const std::string atom_opt){
    int ii;
    int jj;
    int i;
    int j;

    float lb=12.5;
    float ub=15.0;
    if     (atom_opt==" C4'") {lb=14.0;ub=16.0;}
    else if(atom_opt==" C5'") {lb=16.0;ub=18.0;}
    else if(atom_opt==" O3'") {lb=13.5;ub=16.5;}
    else if(atom_opt==" O5'") {lb=15.5;ub=18.5;}
    else if(atom_opt==" P  ") {lb=16.5;ub=21.0;}

    float dis;
    std::vector<bool> bp_tmp(len,false);
    std::vector<std::vector<bool> > bp(len,bp_tmp);
    bp_tmp.clear();
    sec.assign(len, '.');
    for (i=0; i<len; i++)
    {
        for (j=i+1; j<len; j++)
        {
            if (((seq[i]=='u'||seq[i]=='t')&&(seq[j]=='a'             ))||
                ((seq[i]=='a'             )&&(seq[j]=='u'||seq[j]=='t'))||
                ((seq[i]=='g'             )&&(seq[j]=='c'||seq[j]=='u'))||
                ((seq[i]=='c'||seq[i]=='u')&&(seq[j]=='g'             )))
            {
                dis=sqrt(dist(x[i], x[j]));
                bp[j][i]=bp[i][j]=(dis>lb && dis<ub);
            }
        }
    }

    std::vector<int> A0_var,B0_var,C0_var,D0_var;
    for (i=0; i<len-2; i++)
    {
        for (j=i+3; j<len; j++)
        {
            if (!bp[i][j]) continue;
            if (i>0 && j+1<len && bp[i-1][j+1]) continue;
            if (!bp[i+1][j-1]) continue;
            sec_str(len, bp, i,j,ii,jj);
            if (jj<i || j<ii)
            {
                ii=i;
                jj=j;
            }
            A0_var.push_back(i);
            B0_var.push_back(j);
            C0_var.push_back(ii);
            D0_var.push_back(jj);
        }
    }

    for (i=0;i<A0_var.size();i++)
    {
        for (j=0;;j++)
        {
            if(A0_var[i]+j>C0_var[i]) break;
            sec[A0_var[i]+j]='<';
            sec[D0_var[i]+j]='>';
        }
    }
    // clean up
    A0_var.clear();
    B0_var.clear();
    C0_var.clear();
    D0_var.clear();
    bp.clear();
}
bool get_initial5( CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CharMatrix& path, DoubleMatrix& val,
    const CoordArray& x, const CoordArray& y, int xlen, int ylen, std::vector<int>& y2x,
    double d0, double d0_search, const bool fast_opt, const double D0_MIN){
    double GL;
    double rmsd;
    Vec3 t;
    RotMat u;

    double d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    double d02 = d01 * d01;

    double GLmax = 0;
    int aL = getmin(xlen, ylen);
    std::vector<int> invmap(ylen + 1);

    // jump on sequence1
    int n_jump1 = 0;
    if (xlen > 250)
        n_jump1 = 45;
    else if (xlen > 200)
        n_jump1 = 35;
    else if (xlen > 150)
        n_jump1 = 25;
    else
        n_jump1 = 15;
    if (n_jump1 > (xlen / 3))
        n_jump1 = xlen / 3;

    // jump on sequence2
    int n_jump2 = 0;
    if (ylen > 250)
        n_jump2 = 45;
    else if (ylen > 200)
        n_jump2 = 35;
    else if (ylen > 150)
        n_jump2 = 25;
    else
        n_jump2 = 15;
    if (n_jump2 > (ylen / 3))
        n_jump2 = ylen / 3;

    // fragment sizes for superposition
    std::array<int, 2> n_frag = { 20, 100 };
    if (n_frag[0] > (aL / 3))
        n_frag[0] = aL / 3;
    if (n_frag[1] > (aL / 2))
        n_frag[1] = aL / 2;

    // fast mode: skip more positions
    if (fast_opt)
    {
        n_jump1 *= 5;
        n_jump2 *= 5;
    }

    std::vector<double*> xv(xlen), yv(ylen);
    for (int _i = 0; _i < xlen; _i++) xv[_i] = (double*)x[_i].data();
    for (int _j = 0; _j < ylen; _j++) yv[_j] = (double*)y[_j].data();

    bool flag = false;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i < m1; i += n_jump1)
        {
            for (int j = 0; j < m2; j += n_jump2)
            {
                // extract fragment
                for (int k = 0; k < n_frag[i_frag]; k++)
                {
                    r1[k][0] = x[k + i][0];
                    r1[k][1] = x[k + i][1];
                    r1[k][2] = x[k + i][2];

                    r2[k][0] = y[k + j][0];
                    r2[k][1] = y[k + j][1];
                    r2[k][2] = y[k + j][2];
                }

                // superpose the two structures and rotate it
                Kabsch(r1, r2, n_frag[i_frag], 1, rmsd, t, u);

                double gap_open = 0.0;
                NWDP_TM(path, val, x, y, xlen, ylen,
                    t, u, d02, gap_open, invmap);
                GL = get_score_fast(r1, r2, xtm, ytm, x, y, ylen,
                    invmap, d0, d0_search, t, u);
                if (GL > GLmax)
                {
                    GLmax = GL;
                    for (int ii = 0; ii < ylen; ii++) y2x[ii] = invmap[ii];
                    flag = true;
                }
            }
        }
    }

    return flag;
}
void find_max_frag(const CoordArray& x, int len, int &start_max,
    int &end_max, double dcu0, const bool fast_opt){
    int r_min, fra_min=4;
    if (fast_opt) fra_min=8;
    int start, Lfr_max=0;
    r_min= static_cast<int>(len*1.0/3.0);
    if(r_min > fra_min) r_min=fra_min;
    int inc=0;
    double dcu0_cut=dcu0*dcu0;
    double dcu_cut=dcu0_cut;
    while(Lfr_max < r_min) {
        Lfr_max=0;
        int j=1; start=0;
        for(int i=1; i<len; i++) {
            if(dist(x[i-1], x[i]) < dcu_cut) {
                j++;
                if(i==(len-1)) {
                    if(j > Lfr_max) { Lfr_max=j; start_max=start; end_max=i; }
                    j=1;
                }
            } else {
                if(j>Lfr_max) { Lfr_max=j; start_max=start; end_max=i-1; }
                j=1; start=i;
            }
        }
        if(Lfr_max<r_min) {
            inc++;
            double dinc=pow(1.1, static_cast<double>(inc)) * dcu0;
            dcu_cut= dinc*dinc;
        }
    }
}