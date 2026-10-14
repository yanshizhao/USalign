#include "TMalign.h"
#include "common.h"

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
    if (dis_vec.size()>2)
    {
        if (d002t < dis_vec[2])
            d002t = dis_vec[2];
    }
    else if (d002t < dis_vec.back())
        d002t = dis_vec.back();
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
        if (dis_vec2.size()>2)
        {
            if (d002t < dis_vec2[2])
                d002t = dis_vec2[2];
        }
        else if (d002t < dis_vec2.back())
            d002t = dis_vec2.back();
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
    A0_var.clear();
    B0_var.clear();
    C0_var.clear();
    D0_var.clear();
    bp.clear();
}
bool get_initial5(const CoordArray& x, const CoordArray& y, int xlen, int ylen,
    std::vector<int>& y2x, double d0, double d0_search, const bool fast_opt,
    const double D0_MIN, int parallel_threads){
    double d01 = d0 + 1.5;
    if (d01 < D0_MIN) d01 = D0_MIN;
    double d02 = d01 * d01;

    int aL = getmin(xlen, ylen);

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

    std::vector<std::array<int, 3> > frag_tasks;
    for (int i_frag = 0; i_frag < 2; i_frag++)
    {
        int m1 = xlen - n_frag[i_frag] + 1;
        int m2 = ylen - n_frag[i_frag] + 1;

        for (int i = 0; i < m1; i += n_jump1)
        {
            for (int j = 0; j < m2; j += n_jump2)
            {
                std::array<int, 3> task = { {i_frag, i, j} };
                frag_tasks.push_back(task);
            }
        }
    }
    const int n_tasks = static_cast<int>(frag_tasks.size());

    double GLmax = 0;
    bool flag = false;
    int best_idx = INT_MAX;
    std::vector<int> best_inv(ylen + 1);

#ifdef _OPENMP
    #pragma omp parallel num_threads(parallel_threads) if(parallel_threads>1)
#endif
    {
        CoordArray r1(aL), r2(aL), xtm(aL), ytm(aL);
        std::vector<int> invmap(ylen + 1);
        CharMatrix path(xlen + 1, std::vector<char>(ylen + 1, 0));
        DoubleMatrix val(xlen + 1, std::vector<double>(ylen + 1, 0));
        Vec3 t;
        RotMat u;
        double rmsd;
        double GLmax_l = 0;
        int best_idx_l = INT_MAX;
        std::vector<int> best_inv_l(ylen + 1);

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (int task = 0; task < n_tasks; task++)
        {
            int i_frag = frag_tasks[task][0];
            int i = frag_tasks[task][1];
            int j = frag_tasks[task][2];

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
            double GL = get_score_fast(r1, r2, xtm, ytm, x, y, ylen,
                invmap, d0, d0_search, t, u);
            if (GL > GLmax_l ||
                (GL == GLmax_l && best_idx_l != INT_MAX && task < best_idx_l))
            {
                GLmax_l = GL;
                best_idx_l = task;
                for (int ii = 0; ii < ylen; ii++) best_inv_l[ii] = invmap[ii];
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        if (GLmax_l > GLmax ||
            (GLmax_l == GLmax && best_idx != INT_MAX && best_idx_l != INT_MAX && best_idx_l < best_idx))
        {
            GLmax = GLmax_l;
            best_idx = best_idx_l;
            for (int ii = 0; ii < ylen; ii++) best_inv[ii] = best_inv_l[ii];
            flag = true;
        }
    }

    if (flag)
        for (int ii = 0; ii < ylen; ii++) y2x[ii] = best_inv[ii];
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
void parse_file_into_cache(const std::string& fname,
    const ChainParseOptions& parse_opts,
    std::vector<ParsedChain>& all_chains,
    std::map<std::string, std::vector<int>>& cache,
    const std::string& key)
{
    if (cache.count(key)) return;
    std::vector<std::vector<std::string>> PDB_lines;
    std::vector<int> mol_vec;
    std::vector<std::string> chainID_list;
    int nchain = get_PDB_lines(fname, PDB_lines, chainID_list, mol_vec, parse_opts.infmt_opt,
        parse_opts.pdb_parse_opts);
    if (nchain == 0)
    {
        std::cerr << "Warning! Cannot parse file: " << fname
            << ". Chain number 0." << std::endl;
        return;
    }
    std::vector<int> indices;
    for (int c = 0; c < nchain; c++) {
        int len = (int)PDB_lines[c].size();
        if (parse_opts.mol_opt == "RNA") mol_vec[c] = 1;
        else if (parse_opts.mol_opt == "protein") mol_vec[c] = -1;
        if (!len)
        {
            std::cerr << "Warning! Cannot parse file: " << fname
                << ". Chain length 0." << std::endl;
            indices.push_back(-1);
            continue;
        }
        else if (len < 3)
        {
            std::cerr << "Sequence is too short <3!: " << fname << std::endl;
            indices.push_back(-1);
            continue;
        }
        int idx = (int)all_chains.size();
        all_chains.emplace_back();
        auto& chain = all_chains.back();
        chain.filename = fname; chain.chain_len = len;
        chain.chain_id = chainID_list[c]; chain.cur_complex_mol_list = mol_vec[c];
        chain.chain_coords.reserve(len);
        std::string seq;
        chain.chain_len = read_PDB(PDB_lines[c], chain.chain_coords, seq,
            chain.resi_vec, parse_opts.read_resi);
        if (parse_opts.mirror_opt)
            for (int r = 0; r < chain.chain_len; r++)
                chain.chain_coords[r][2] = -chain.chain_coords[r][2];
        chain.chain_seq = seq;
        if (mol_vec[c] > 0)
            make_sec(seq, chain.chain_coords, chain.chain_len, chain.chain_sec, parse_opts.pdb_parse_opts.atom_opt);
        else
            make_sec(chain.chain_coords, chain.chain_len, chain.chain_sec);
        if (parse_opts.retain_pdb_lines) chain.pdb_lines = std::move(PDB_lines[c]);
        else PDB_lines[c].clear();
        indices.push_back(idx);
    }
    PDB_lines.clear();
    cache[key] = indices;
}

static inline std::string make_cache_key(const std::string& fname, bool struct1)
{
    return (struct1 ? std::string("1:") : std::string("2:")) + fname;
}

// TMalign, RNAalign, CPalign, TMscore
#ifdef _OPENMP
int run_batch_parallel(
    const std::vector<std::string>& chain1_list, const std::vector<std::string>& chain2_list,
    const std::vector<std::string>& sequence, const std::string& dir_opt, const std::string& dir1_opt,
    const std::string& dir2_opt, const std::string& dirpair_opt, const ChainParseOptions& parse_opts1, const ChainParseOptions& parse_opts2,
    const ChainPairAlignOptions& align_opts, const OutputOptions& out_opts, bool do_opt, int parallel_threads)
{
    int i, j, chain_i;

    // ---- Phase 1: pre-parse all unique files, build task list ----
    std::vector<ParsedChain> all_chains;
    std::map<std::string, std::vector<int>> file_to_idx;
    std::vector<PairTask> tasks;

    for (i = 0; i < (int)chain1_list.size(); i++) {
        const std::string key1 = make_cache_key(chain1_list[i], true);
        parse_file_into_cache(chain1_list[i], parse_opts1, all_chains, file_to_idx, key1);

        auto& c1_indices = file_to_idx[key1];
        for (chain_i = 0; chain_i < (int)c1_indices.size(); chain_i++) {
            int c1_idx = c1_indices[chain_i];
            if (c1_idx < 0) continue;
            int j_start = (dir_opt.size() > 0) * (i + 1);

            for (j = j_start; j < (int)chain2_list.size(); j++) {
                if (dirpair_opt.size() && j != i) continue;
                const std::string key2 = make_cache_key(chain2_list[j], false);

                parse_file_into_cache(chain2_list[j], parse_opts2, all_chains, file_to_idx, key2);
                auto& c2_indices = file_to_idx[key2];

                for (int c2_i = 0; c2_i < (int)c2_indices.size(); c2_i++) {
                    int c2_idx = c2_indices[c2_i];
                    if (c2_idx < 0) continue;
                    tasks.push_back({c1_idx, c2_idx, (int)tasks.size()});

                }
            }
        }
    }

    // ---- Phase 2: parallel pair processing ----
    std::vector<std::string> out_lines(tasks.size());
    #pragma omp parallel for schedule(dynamic, 8) num_threads(parallel_threads)
    for (int t = 0; t < (int)tasks.size(); t++) {
        auto& task = tasks[t];
        auto& c1 = all_chains[task.chain1_idx];
        auto& c2 = all_chains[task.chain2_idx];
        CoordArray xa_c = c1.chain_coords;
        CoordArray ya_c = c2.chain_coords;

        bool force_fast = (std::min(c1.chain_len, c2.chain_len) > 1500) ? true : align_opts.fast_opt;
        ChainPairAlignResult chain_align_res = { 0};
        chain_align_res.d0_out = 5.0;

        ChainPairAlignOptions pair_align_opts = align_opts;
        pair_align_opts.fast_opt = force_fast;
        pair_align_opts.parallel_threads = 1;
        pair_align_opts.ss_opt = 0;
        pair_align_opts.mol_type = c1.cur_complex_mol_list + c2.cur_complex_mol_list;

        align_chain_pair(chain_align_res, xa_c, ya_c, c1.chain_seq, c2.chain_seq,
            c1.chain_sec, c2.chain_sec, c1.chain_len, c2.chain_len,
            pair_align_opts, sequence, out_opts.outfmt_opt);

        std::stringstream ss;
        std::string xname_out = c1.filename.substr( dir1_opt.size() + dir_opt.size() + dirpair_opt.size());
        std::string yname_out = c2.filename.substr(dir2_opt.size() + dir_opt.size() + dirpair_opt.size());

        if (out_opts.outfmt_opt == 0) print_version(USALIGN_VERSION, ss);

        int right_num=0;
        if (align_opts.cp_opt) output_cp(xname_out, yname_out, chain_align_res.seqxA, out_opts.outfmt_opt, right_num, ss);

        const ChainPairDesc pair_desc{xname_out, yname_out, c1.chain_id, c2.chain_id, c1.chain_len, c2.chain_len};
        output_results(
            pair_desc, chain_align_res, out_opts, c1.resi_vec, c2.resi_vec, ss);

        if (do_opt || (align_opts.cp_opt && out_opts.outfmt_opt <= 0))
            output_do_distance_table(ss, chain_align_res.seqxA, chain_align_res.seqyA,
                c1.pdb_lines, c2.pdb_lines, chain_align_res.do_vec, right_num);

        out_lines[task.order] = ss.str();
    }

    // ---- Phase 3: serial output in original order ----
    for (int t = 0; t < (int)tasks.size(); t++)
        std::cout << out_lines[t];

    return 0;
}
#endif

void fill_tmalign_params(TMalignParams& params, const AlignCommonInput& common_inputs)
{
    const UserOptions& user_opts = common_inputs.user_options;

    params.read_resi = ((user_opts.byresi_opt == 0) && user_opts.o_opt) ? 2 : user_opts.byresi_opt;
    fill_chain_pair_align_options(params.align, user_opts, common_inputs.control_options);
}

int TMalign(AlignCommonInput& common_inputs, const TMalignParams& tm_params)
{
    UserOptions& user_opts = common_inputs.user_options;
    ParsedInput& parsed_input = common_inputs.parsed_input;
    ControlOptions& ctrl_opts = common_inputs.control_options;

    std::vector<std::vector<std::string> >PDB_lines1; // text of chain1
    std::vector<std::vector<std::string> >PDB_lines2; // text of chain2
    std::vector<int> mol_vec1;              // molecule type of chain1, RNA if >0
    std::vector<int> mol_vec2;              // molecule type of chain2, RNA if >0
    std::vector<std::string> chainID_list1;      // list of chainID1
    std::vector<std::string> chainID_list2;      // list of chainID2
    int    i,j;                // file index
    int    chain_i,chain_j;    // chain index
    int    r;                  // residue index
    int    xlen, ylen;         // chain length
    int    xchainnum=0,ychainnum=0;// number of chains in a PDB file
    std::string secx;                // for the secondary structure
    std::string secy;
    CoordArray xa;                  // for input vectors xa[0...xlen-1][0..2] and
    CoordArray ya;                  // ya[0...ylen-1][0..2], in general,
                               // ya is regarded as native structure
                               // --> superpose xa onto ya
    std::vector<std::string> resi_vec1;  // residue index for chain1
    std::vector<std::string> resi_vec2;  // residue index for chain2

#ifdef _OPENMP
    // === Parallel batch mode ===
    if (ctrl_opts.parallel_threads > 1 && (parsed_input.chain1_list.size() > 1 || parsed_input.chain2_list.size() > 1)) {
        ChainParseOptions parse_opts1;
        parse_opts1.pdb_parse_opts = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
            user_opts.split_opt, user_opts.het_opt,
            user_opts.chain2parse1, user_opts.model2parse1);
        parse_opts1.infmt_opt = user_opts.infmt1_opt;
        parse_opts1.mol_opt = user_opts.mol_opt;
        parse_opts1.read_resi = tm_params.read_resi;
        parse_opts1.mirror_opt = (user_opts.mirror_opt != 0);
        parse_opts1.retain_pdb_lines = ctrl_opts.cp_opt;

        ChainParseOptions parse_opts2;
        parse_opts2.pdb_parse_opts = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
            user_opts.split_opt, user_opts.het_opt,
            user_opts.chain2parse2, user_opts.model2parse2);
        parse_opts2.infmt_opt = user_opts.infmt2_opt;
        parse_opts2.mol_opt = user_opts.mol_opt;
        parse_opts2.read_resi = tm_params.read_resi;
        parse_opts2.mirror_opt = false;
        parse_opts2.retain_pdb_lines = ctrl_opts.cp_opt;

        ChainPairAlignOptions align_opts = {0};
        fill_chain_pair_align_options(align_opts, user_opts, ctrl_opts);

        OutputOptions out_opts;
        fill_out_opts_from_user_opts(out_opts, user_opts);

        return run_batch_parallel(
            parsed_input.chain1_list, parsed_input.chain2_list,
            parsed_input.sequence,
            user_opts.dir_opt, user_opts.dir1_opt, user_opts.dir2_opt, user_opts.dirpair_opt,
            parse_opts1, parse_opts2, align_opts, out_opts,
            false, ctrl_opts.parallel_threads);
    }
#endif  // _OPENMP

    const PdbParseOptions pdb_parse_opts1 = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse1, user_opts.model2parse1);
    const PdbParseOptions pdb_parse_opts2 = make_pdb_parse_opts(user_opts.ter_opt, user_opts.atom_opt, parsed_input.normalize_atom_name,
        user_opts.split_opt, user_opts.het_opt,
        user_opts.chain2parse2, user_opts.model2parse2);

    // loop over file names (original serial code)
    for (i=0;i<parsed_input.chain1_list.size();i++)
    {
        // parse chain 1
        user_opts.xname=parsed_input.chain1_list[i];
        xchainnum=get_PDB_lines(user_opts.xname, PDB_lines1, chainID_list1, mol_vec1, user_opts.infmt1_opt,
            pdb_parse_opts1);
        if (!xchainnum)
        {
            std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                <<". Chain number 0."<<std::endl;
            continue;
        }
        for (chain_i=0;chain_i<xchainnum;chain_i++)
        {
            xlen=PDB_lines1[chain_i].size();
            if (user_opts.mol_opt=="RNA") mol_vec1[chain_i]=1;
            else if (user_opts.mol_opt=="protein") mol_vec1[chain_i]=-1;
            if (!xlen)
            {
                std::cerr<<"Warning! Cannot parse file: "<<user_opts.xname
                    <<". Chain length 0."<<std::endl;
                continue;
            }
            else if (xlen<3)
            {
                std::cerr<<"Sequence is too short <3!: "<<user_opts.xname<<std::endl;
                continue;
            }
            xa.clear();
            xa.reserve(xlen);
            std::string seqx;
            secx.resize(xlen + 1);
            xlen = read_PDB(PDB_lines1[chain_i], xa, seqx,
                resi_vec1, tm_params.read_resi);
            if (user_opts.mirror_opt) for (r=0;r<xlen;r++) xa[r][2]=-xa[r][2];
            if (mol_vec1[chain_i]>0) make_sec(seqx, xa, xlen, secx, user_opts.atom_opt);
            else make_sec(xa, xlen, secx); // secondary structure assignment

            int j_start = (user_opts.dir_opt.size() > 0) * (i + 1);
            for (j=j_start;j<parsed_input.chain2_list.size();j++)
            {
                if (user_opts.dirpair_opt.size() && j!=i) continue;
                // parse chain 2
                if (PDB_lines2.size()==0)
                {
                    user_opts.yname=parsed_input.chain2_list[j];
                    ychainnum=get_PDB_lines(user_opts.yname, PDB_lines2, chainID_list2, mol_vec2, user_opts.infmt2_opt,
                        pdb_parse_opts2);
                    if (!ychainnum)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname<<". Chain number 0."<<std::endl;
                        continue;
                    }
                }
                for (chain_j=0;chain_j<ychainnum;chain_j++)
                {
                    ylen=PDB_lines2[chain_j].size();
                    if (user_opts.mol_opt=="RNA") mol_vec2[chain_j]=1;
                    else if (user_opts.mol_opt=="protein") mol_vec2[chain_j]=-1;
                    if (!ylen)
                    {
                        std::cerr<<"Warning! Cannot parse file: "<<user_opts.yname
                            <<". Chain length 0."<<std::endl;
                        continue;
                    }
                    else if (ylen<3)
                    {
                        std::cerr<<"Sequence is too short <3!: "<<user_opts.yname<<std::endl;
                        continue;
                    }
                    ya.clear();
                    ya.reserve(ylen);
                    std::string seqy;
                    secy.resize(ylen + 1);
                    ylen = read_PDB(PDB_lines2[chain_j], ya, seqy,
                        resi_vec2, tm_params.read_resi);
                    if (mol_vec2[chain_j]>0)
                         make_sec(seqy, ya, ylen, secy, user_opts.atom_opt);
                    else make_sec(ya, ylen, secy);

                    if (user_opts.byresi_opt) extract_aln_from_resi(parsed_input.sequence,
                        seqx, seqy, resi_vec1, resi_vec2, user_opts.byresi_opt);

                    ChainPairAlignResult chain_align_res = { 0};
                    chain_align_res.d0_out = 5.0;
                    ChainPairAlignOptions align_opts = tm_params.align;
                    align_opts.fast_opt = (getmin(xlen, ylen) > 1500) ? true : user_opts.fast_opt;
                    align_opts.mol_type = mol_vec1[chain_i] + mol_vec2[chain_j];
                    align_chain_pair(chain_align_res, xa, ya, seqx, seqy, secx, secy,
                        xlen, ylen, align_opts, parsed_input.sequence,
                        user_opts.outfmt_opt);

                    output_pair_alignment(common_inputs,
                        xlen, ylen, chain_align_res,
                        PDB_lines1[chain_i], PDB_lines2[chain_j],
                        chainID_list1[chain_i], chainID_list2[chain_j],
                        resi_vec1, resi_vec2);
                } // chain_j
                if (parsed_input.chain2_list.size()>1)
                {
                    user_opts.yname.clear();
                    PDB_lines2.clear();
                    chainID_list2.clear();
                    mol_vec2.clear();
                }
            } // j
            PDB_lines1[chain_i].clear();
            resi_vec1.clear();
        } // chain_i
        user_opts.xname.clear();
        PDB_lines1.clear();
        chainID_list1.clear();
        mol_vec1.clear();
    } // i
    if (parsed_input.chain2_list.size()==1)
    {
        user_opts.yname.clear();
    }
    return 0;
}

// Output the check warning for the forced molecule type from -mol
// Params: file_path - structure file path
//         chain_num - total number of chains
//         prot_count, na_count - protein / RNA chain counts
//         filtered_count - number of filtered chains
//         mol_opt - value of the -mol option (protein / RNA)

int TMalign_main(CoordArray& xa_c, CoordArray& ya_c,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string> sequence,
    const ChainPairAlignOptions& align_opts,
    std::ostream& os)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double &TM1 = chain_align_res.TM1;
    double &TM2 = chain_align_res.TM2;
    double &TM3 = chain_align_res.TM3;
    double &TM4 = chain_align_res.TM4;
    double &TM5 = chain_align_res.TM5;
    TMscoreAlignContext align_context;

    /***********************/
    // allocate memory
    /***********************/
    int minlen = std::min(xlen, ylen);
    align_context.score.assign(xlen+1, std::vector<double>(ylen+1));
    align_context.path.assign( xlen+1, std::vector<char>(ylen+1));
    align_context.val.assign(  xlen+1, std::vector<double>(ylen+1));
    align_context.coords.xtm.resize(minlen);
    align_context.coords.ytm.resize(minlen);
    align_context.coords.xt.resize(xlen);
    align_context.coords.r1.resize(minlen);
    align_context.coords.r2.resize(minlen);

    /***********************/
    //    parameter set
    /***********************/
    parameter_set4search(xlen, ylen, align_context.search_params.D0_MIN, align_context.search_params.Lnorm,
        align_context.search_params.score_d8, align_context.search_params.d0, align_context.search_params.d0_search, align_context.search_params.dcu0);
    align_context.search_params.simplify_step    = 40; //for simplified search engine
    align_context.search_params.score_sum_method = 8;  //for scoring method, whether only sum over pairs with dis<score_d8

    align_context.res_map.invmap0.assign(ylen+1, -1);
    align_context.res_map.invmap.assign(ylen+1, 0);
    align_context.transform.TMmax = -1;

    align_context.ddcc = 0.4;
    if (align_context.search_params.Lnorm <= 40) align_context.ddcc = 0.1;   //Lnorm was setted in parameter_set4search
    align_context.search_params.local_d0_search = align_context.search_params.d0_search;

    int rc = tmalign_build_initial(xa_c, ya_c, secx, secy,
        xlen, ylen, sequence, align_opts, align_context, chain_align_res);
    if (rc) return rc;

    //*******************************************************************//
    //    The alignment will not be changed any more in the following    //
    //*******************************************************************//
    //check if the initial alignment is generated appropriately
    bool flag=false;
    for(int i=0; i<ylen; i++)
    {
        if(align_context.res_map.invmap0[i]>=0)
        {
            flag=true;
            break;
        }
    }
    if(!flag)
    {
        os << "There is no alignment between the two structures! "
             << "Program stop with no result!" << std::endl;
        TM1=TM2=TM3=TM4=TM5=0;
        init_rt_mat(t0, u0);
        return 1;
    }

    // last TM-score pre-termination
    if (align_opts.TMcut>0)
    {
        double TMtmp=approx_TM(xlen, ylen, align_opts.a_opt,
            xa_c, ya_c, t0, u0, align_context.res_map.invmap0, align_opts.mol_type);

        if (TMtmp<0.6*align_opts.TMcut)
        {
            TM1=TM2=TM3=TM4=TM5=TMtmp;


            return 7;
        }
    }

    tmalign_select_aligned_pairs(xa_c, ya_c, xlen, ylen, align_opts, align_context, chain_align_res);
    tmalign_final_tmscore(xlen, ylen, align_opts, align_context, chain_align_res);
    tmalign_extract_strings(xa_c, ya_c, seqx, seqy, xlen, ylen,
        chain_align_res.d0_out, chain_align_res.n_ali8, align_context, chain_align_res);

    return 0; // zero for no exception
}

int CPalign_main(CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    ChainPairAlignResult& chain_align_res,
    const int xlen, const int ylen,
    const std::vector<std::string>& sequence,
    const ChainPairAlignOptions& align_opts,
    std::ostream& os)
{
    Vec3& t0 = chain_align_res.t0;
    RotMat& u0 = chain_align_res.u0;
    double& TM3 = chain_align_res.TM3;
    double& TM4 = chain_align_res.TM4;
    double& TM5 = chain_align_res.TM5;
    double& d0_0 = chain_align_res.d0_0;
    double& TM_0 = chain_align_res.TM_0;
    double& d0A = chain_align_res.d0A;
    double& d0B = chain_align_res.d0B;
    double& d0u = chain_align_res.d0u;
    double& d0a = chain_align_res.d0a;
    double& d0_out = chain_align_res.d0_out;
    std::string& seqM = chain_align_res.seqM;
    std::string& seqxA = chain_align_res.seqxA;
    std::string& seqyA = chain_align_res.seqyA;
    std::vector<double>& do_vec = chain_align_res.do_vec;
    double& rmsd0 = chain_align_res.rmsd0;
    int& L_ali = chain_align_res.L_ali;
    double& Liden = chain_align_res.Liden;
    double& TM_ali = chain_align_res.TM_ali;
    double& rmsd_ali = chain_align_res.rmsd_ali;
    int& n_ali = chain_align_res.n_ali;
    int& n_ali8 = chain_align_res.n_ali8;
    const double d0_scale = align_opts.d0_scale;
    const double Lnorm_cp=getmin(xlen,ylen);
    double TM4_cp;
    ChainPairAlignOptions cp_trial_opts;
    cp_trial_opts.i_opt = 0;
    cp_trial_opts.a_opt = 0;
    cp_trial_opts.u_opt = true;
    cp_trial_opts.d_opt = false;
    cp_trial_opts.fast_opt = true;
    cp_trial_opts.se_opt = false;
    cp_trial_opts.cp_opt = false;
    cp_trial_opts.Lnorm = Lnorm_cp;
    cp_trial_opts.d0_scale = d0_scale;
    cp_trial_opts.TMcut = -1;
    cp_trial_opts.parallel_threads = 1;
    cp_trial_opts.ss_opt = 0;
    cp_trial_opts.mol_type = align_opts.mol_type;

    std::string seqx_cp;
    std::string secx_cp;
    CoordArray xa_cp;
    std::string seqxA_cp,seqyA_cp;
    int cp_point=0;
    int cp_aln_best=0;

    build_cp_doubled_template(xa, seqx, secx, xlen, xa_cp, seqx_cp, secx_cp);

    ChainPairAlignResult cp_trial_res = chain_align_res;
    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy, cp_trial_res, xlen*2, ylen, sequence, cp_trial_opts, os);

    t0 = cp_trial_res.t0;
    u0 = cp_trial_res.u0;
    TM3 = cp_trial_res.TM3;
    TM4_cp = cp_trial_res.TM4;
    TM5 = cp_trial_res.TM5;
    d0_0 = cp_trial_res.d0_0;
    TM_0 = cp_trial_res.TM_0;
    d0A = cp_trial_res.d0A;
    d0B = cp_trial_res.d0B;
    d0u = cp_trial_res.d0u;
    d0a = cp_trial_res.d0a;
    d0_out = cp_trial_res.d0_out;
    seqM = cp_trial_res.seqM;
    seqxA_cp = cp_trial_res.seqxA;
    seqyA_cp = cp_trial_res.seqyA;
    do_vec = cp_trial_res.do_vec;
    rmsd0 = cp_trial_res.rmsd0;
    L_ali = cp_trial_res.L_ali;
    Liden = cp_trial_res.Liden;
    TM_ali = cp_trial_res.TM_ali;
    rmsd_ali = cp_trial_res.rmsd_ali;
    n_ali = cp_trial_res.n_ali;
    n_ali8 = cp_trial_res.n_ali8;

    strip_gaps_from_pair_alignment(seqxA_cp, seqyA_cp, seqxA, seqyA);

    find_cp_point(seqyA, xlen, cp_aln_best, cp_point);
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    seqxA_cp.clear();
    seqyA_cp.clear();
    rmsd0=Liden=n_ali=n_ali8=0;

    TMalign_main(xa, ya, seqx, seqy, secx, secy, chain_align_res, xlen, ylen, sequence, cp_trial_opts, os);
    seqxA_cp = chain_align_res.seqxA;
    seqyA_cp = chain_align_res.seqyA;

    if (n_ali8>=cp_aln_best || TM4>=TM4_cp) cp_point=0;

    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    rmsd0=Liden=n_ali=n_ali8=0;
    shift_cp_template_by_point(xa_cp, seqx_cp, secx_cp, xlen, cp_point);
    seqx_cp[xlen]=0;
    secx_cp[xlen]=0;

    if (cp_point!=0)
    {
        ChainPairAlignResult cp_best_res = chain_align_res;
        cp_best_res.n_ali8 = cp_aln_best;
        TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy,
            cp_best_res, xlen, ylen, sequence, cp_trial_opts, os);
        t0 = cp_best_res.t0;
        u0 = cp_best_res.u0;
        TM3 = cp_best_res.TM3;
        TM4_cp = cp_best_res.TM4;
        TM5 = cp_best_res.TM5;
        d0_0 = cp_best_res.d0_0;
        TM_0 = cp_best_res.TM_0;
        d0A = cp_best_res.d0A;
        d0B = cp_best_res.d0B;
        d0u = cp_best_res.d0u;
        d0a = cp_best_res.d0a;
        d0_out = cp_best_res.d0_out;
        seqM = cp_best_res.seqM;
        seqxA_cp = cp_best_res.seqxA;
        seqyA_cp = cp_best_res.seqyA;
        do_vec = cp_best_res.do_vec;
        rmsd0 = cp_best_res.rmsd0;
        L_ali = cp_best_res.L_ali;
        Liden = cp_best_res.Liden;
        TM_ali = cp_best_res.TM_ali;
        rmsd_ali = cp_best_res.rmsd_ali;
        n_ali = cp_best_res.n_ali;
        cp_aln_best = cp_best_res.n_ali8;
        if (n_ali8>=cp_aln_best || TM4>=TM4_cp)
        {
            cp_point=0;
            restore_cp_template(xa_cp, seqx_cp, secx_cp, xa, seqx, secx, xlen);
        }
    }

    TMalign_main(xa_cp, ya, seqx_cp, seqy, secx_cp, secy, chain_align_res, xlen, ylen, sequence, align_opts, os);
    seqxA_cp = chain_align_res.seqxA;
    seqyA_cp = chain_align_res.seqyA;

    if (cp_point>0)
    {
        apply_cp_marker_to_alignment(seqxA_cp, seqyA_cp, xlen, cp_point, seqM, seqxA, seqyA);
    }
    else
    {
        seqxA=seqxA_cp;
        seqyA=seqyA_cp;
    }

    return cp_point;
}

void tmscore_search_pos(CoordArray& xtm, const CoordArray& ytm,
    int Lali, int L_frag, int ii,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0)
{
    CoordArray r1_l(Lali), r2_l(Lali), xt_l(Lali);
    std::vector<int> i_ali_l(Lali), k_ali_l(Lali);
    Vec3 t_l, t_l_best; RotMat u_l, u_l_best;
    double score_l, score_max_l = -1;
    int n_cut_l, ka_l; double rmsd_l, d_l;

    // ① Extract fragment of L_frag residues starting at ii
    ka_l = 0;
    for (int kk = 0; kk < L_frag; kk++) {
        int idx = kk + ii;
        r1_l[kk][0] = xtm[idx][0];
        r1_l[kk][1] = xtm[idx][1];
        r1_l[kk][2] = xtm[idx][2];
        r2_l[kk][0] = ytm[idx][0];
        r2_l[kk][1] = ytm[idx][1];
        r2_l[kk][2] = ytm[idx][2];
        k_ali_l[ka_l++] = idx;
    }

    // ② Kabsch rotation from the fragment
    Kabsch(r1_l, r2_l, L_frag, 1, rmsd_l, t_l, u_l);
    do_rotation(xtm, xt_l, Lali, t_l, u_l);

    // ③ Initial score
    d_l = local_d0_search - 1;
    if (use_standard)
        n_cut_l = score_fun8_standard(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
            score_sum_method, score_d8, d0);
    else
        n_cut_l = score_fun8(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
            score_sum_method, Lnorm, score_d8, d0);
    if (score_l > score_max_l) {
        score_max_l = score_l;
        t_l_best = t_l; u_l_best = u_l;
    }

    // ④ Iterative refinement (≤ n_it iterations)
    d_l = local_d0_search + 1;
    for (int it = 0; it < n_it; it++) {
        ka_l = 0;
        for (int kk = 0; kk < n_cut_l; kk++) {
            int m = i_ali_l[kk];
            r1_l[kk][0] = xtm[m][0];
            r1_l[kk][1] = xtm[m][1];
            r1_l[kk][2] = xtm[m][2];
            r2_l[kk][0] = ytm[m][0];
            r2_l[kk][1] = ytm[m][1];
            r2_l[kk][2] = ytm[m][2];
            k_ali_l[ka_l++] = m;
        }
        Kabsch(r1_l, r2_l, n_cut_l, 1, rmsd_l, t_l, u_l);
        do_rotation(xtm, xt_l, Lali, t_l, u_l);
        if (use_standard)
            n_cut_l = score_fun8_standard(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
                score_sum_method, score_d8, d0);
        else
            n_cut_l = score_fun8(xt_l, ytm, Lali, d_l, i_ali_l, score_l,
                score_sum_method, Lnorm, score_d8, d0);
        if (score_l > score_max_l) {
            score_max_l = score_l;
            t_l_best = t_l; u_l_best = u_l;
        }
        if (n_cut_l == ka_l) {
            int kk;
            for (kk = 0; kk < n_cut_l; kk++)
                if (i_ali_l[kk] != k_ali_l[kk]) break;
            if (kk == n_cut_l) break;
        }
    }

    // ⑤ Update global best (protected by critical section for thread safety)
    #pragma omp critical(tmscore_update)
    if (score_max_l > score_max) {
        score_max = score_max_l;
        for (int kk = 0; kk < 3; kk++) {
            t0[kk] = t_l_best[kk];
            for (int j = 0; j < 3; j++)
                u0[kk][j] = u_l_best[kk][j];
        }
    }
}

void tmscore_search_serial(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt,
    int Lali, int L_frag, int iL_max, int simplify_step,
    double local_d0_search, int score_sum_method, double Lnorm,
    double score_d8, double d0, int n_it,
    double& Rcomm,
    bool use_standard,
    double& score_max, Vec3& t0, RotMat& u0)
{
    int i = 0, k, m, ka, n_cut;
    double rmsd, d, score;
    std::vector<int> k_ali(Lali), i_ali(Lali);
    Vec3 t; RotMat u;

    while (1) {
        ka = 0;
        for (k = 0; k < L_frag; k++) {
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
        if (simplify_step != 1) Rcomm = 0;
        do_rotation(xtm, xt, Lali, t, u);

        d = local_d0_search - 1;
        if (use_standard)
            n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, score_d8, d0);
        else
            n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
                score_sum_method, Lnorm, score_d8, d0);
        if (score > score_max) {
            score_max = score;
            for (k = 0; k < 3; k++) {
                t0[k] = t[k];
                u0[k][0] = u[k][0];
                u0[k][1] = u[k][1];
                u0[k][2] = u[k][2];
            }
        }

        d = local_d0_search + 1;
        for (int it = 0; it < n_it; it++) {
            ka = 0;
            for (k = 0; k < n_cut; k++) {
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
            if (use_standard)
                n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, score_d8, d0);
            else
                n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
                    score_sum_method, Lnorm, score_d8, d0);
            if (score > score_max) {
                score_max = score;
                for (k = 0; k < 3; k++) {
                    t0[k] = t[k];
                    u0[k][0] = u[k][0];
                    u0[k][1] = u[k][1];
                    u0[k][2] = u[k][2];
                }
            }
            if (n_cut == ka) {
                for (k = 0; k < n_cut; k++)
                    if (i_ali[k] != k_ali[k]) break;
                if (k == n_cut) break;
            }
        }

        if (i < iL_max) {
            i = i + simplify_step;
            if (i > iL_max) i = iL_max;
        } else if (i >= iL_max) break;
    }
}

double TMscore8_search(CoordArray& r1, CoordArray& r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& xt, int Lali, Vec3& t0, RotMat& u0, int simplify_step,
    int score_sum_method, double &Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0, int parallel_threads)
{
    int i;
    double score_max;

    int n_it=20;
    int n_init_max=6;
    std::array<int, 6> L_ini;
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
    int L_frag;
    int iL_max;

    for(i_init=0; i_init<n_init; i_init++)
    {
        L_frag=L_ini[i_init];
        iL_max=Lali-L_frag;

#ifdef _OPENMP
        if (parallel_threads > 1) {
            if (simplify_step != 1) Rcomm = 0;
            int n_pos = (iL_max + simplify_step - 1) / simplify_step + 1;
            #pragma omp parallel for num_threads(parallel_threads)
            for (int pos = 0; pos < n_pos; pos++) {
                int ii = pos * simplify_step;
                if (ii > iL_max) ii = iL_max;
                tmscore_search_pos(xtm, ytm, Lali, L_frag, ii,
                    local_d0_search, score_sum_method, Lnorm, score_d8, d0, n_it,
                    false, score_max, t0, u0);
            }
        } else
#endif
        {
            tmscore_search_serial(r1, r2, xtm, ytm, xt,
                Lali, L_frag, iL_max, simplify_step,
                local_d0_search, score_sum_method, Lnorm,
                score_d8, d0, n_it,
                Rcomm, false, score_max, t0, u0);
        }
    }
    return score_max;
}

double TMscore8_search_standard(CoordArray& r1, CoordArray& r2,
    CoordArray& xtm, CoordArray& ytm, CoordArray& xt, int Lali,
    Vec3& t0, RotMat& u0, int simplify_step, int score_sum_method,
    double &Rcomm, double local_d0_search, double score_d8, double d0, int parallel_threads)
{
    int i;
    double score_max;

    int n_it = 20;
    int n_init_max = 6;
    std::array<int, 6> L_ini;
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
    int L_frag;
    int iL_max;

    for (i_init = 0; i_init<n_init; i_init++)
    {
        L_frag = L_ini[i_init];
        iL_max = Lali - L_frag;

#ifdef _OPENMP
        if (parallel_threads > 1) {
            if (simplify_step != 1) Rcomm = 0;
            int n_pos = (iL_max + simplify_step - 1) / simplify_step + 1;
            #pragma omp parallel for num_threads(parallel_threads)
            for (int pos = 0; pos < n_pos; pos++) {
                int ii = pos * simplify_step;
                if (ii > iL_max) ii = iL_max;
                tmscore_search_pos(xtm, ytm, Lali, L_frag, ii,
                    local_d0_search, score_sum_method, 0.0, score_d8, d0, n_it,
                    true, score_max, t0, u0);
            }
        } else
#endif
        {
            tmscore_search_serial(r1, r2, xtm, ytm, xt,
                Lali, L_frag, iL_max, simplify_step,
                local_d0_search, score_sum_method, 0.0,
                score_d8, d0, n_it,
                Rcomm, true, score_max, t0, u0);
        }
    }
    return score_max;
}