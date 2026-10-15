#pragma once
/* Functions for the core TMalign algorithm, including the entry function
 * flexalign_main */

#include "TMalign.h"
#include "se.h"
#include "common.h"

#define MAX_SEC_STRUCT_OPT 1

inline void t_u2tu(const Vec3& t0, const RotMat& u0, std::vector<double> &tu_tmp)
{
    for (int i=0;i<3;i++) tu_tmp[i]=t0[i];
    int k=3;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    {
        tu_tmp[k]=u0[i][j];
        k++;
    }
}

inline void tu2t_u(std::vector<double> tu_tmp, Vec3& t0, RotMat& u0)
{
    for (int i=0;i<3;i++) t0[i]=tu_tmp[i];
    int k=3;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
    {
        u0[i][j]=tu_tmp[k];
        k++;
    }
}

// extract rotation matrix based on TMscore8
void output_flexalign_matrix(const std::string& fname_matrix,
    const DoubleMatrix&tu_vec, std::ostream& os = std::cout);

struct FlexalignRasmolState
{
    std::stringstream buf;
    std::stringstream buf_all;
    std::stringstream buf_atm;
    std::stringstream buf_all_atm;
    std::stringstream buf_all_atm_lig;
    std::stringstream buf_tm;
    std::string line;
    Vec3 x;
    Vec3 x1;
    bool after_ter;
    std::string asym_id;
    std::map<std::string,int> resi2hinge_dict;
    std::string resi;
    int read_resi;
    std::map<std::string,int> atom_site;
    int atom_site_pos;
    std::vector<std::string> line_vec;
    std::string atom;
    std::string AA;
    std::string inscode;
    std::string model_index;
    bool is_mmcif;
    int ca_idx1;
    int lig_idx1;
    std::vector<int> idx_vec;
    int ca_idx2;
    int lig_idx2;
    std::vector<std::string> resi_aln1;
    std::vector<std::string> resi_aln2;
    int i1;
    int i2;
    std::string chain1_sele;
    std::string chain2_sele;
};

inline void build_resi_hinge_map(FlexalignRasmolState& rasmol_state,
    const std::string& seqM, const std::string& seqxA,
    const std::vector<std::string>& resi_vec1, const int xlen)
{
    int r;
    int i;
    int j=-1;
    char hinge_char=0;
    int ali_len=seqM.size();
    for (r=0;r<seqxA.size();r++)
    {
        if (seqxA[r]=='-') continue;
        j++;
        hinge_char=seqM[r];
        if (hinge_char==' ')
        {
            for (i=1;i<ali_len;i++)
            {
                if (r-i>=0 && seqM[r-i]!=' ')
                    hinge_char=seqM[r-i];
                else if (r+i<xlen && seqM[r+i]!=' ')
                    hinge_char=seqM[r+i];
                if (hinge_char!=' ') break;
            }
        }
        rasmol_state.resi2hinge_dict[resi_vec1[j]]=hinge_char-'0';
    }
}

inline void build_rasmol_headers(FlexalignRasmolState& rasmol_state,
    const std::string& xname, const std::string& yname,
    const std::string& chainID1, const std::string& chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden, const int mm_opt)
{
    rasmol_state.buf_tm<<"REMARK US-align"
        <<"\nREMARK Structure 1:"<<std::setw(11)<<std::left<<xname+chainID1<<" Size= "<<xlen
        <<"\nREMARK Structure 2:"<<std::setw(11)<<yname+chainID2<<std::right<<" Size= "<<ylen
        <<" (TM-score is normalized by "<<std::setw(4)<<ylen<<", d0="
        <<std::setiosflags(std::ios::fixed)<<std::setprecision(2)<<std::setw(6)<<d0A<<")"
        <<"\nREMARK Aligned length="<<std::setw(4)<<n_ali8<<", RMSD="
        <<std::setw(6)<<setiosflags(std::ios::fixed)<<std::setprecision(2)<<rmsd
        <<", TM-score="<<std::setw(7)<<setiosflags(std::ios::fixed)<<std::setprecision(5)<<TM1
        <<", ID="<<std::setw(5)<<setiosflags(std::ios::fixed)<<std::setprecision(3)
        <<((n_ali8>0)?Liden/n_ali8:0)<<std::endl;
    std::string rasmol_CA_header="load inline\nselect *A\nwireframe .45\nselect *B\nwireframe .20\nselect all\ncolor white\n";
    std::string rasmol_cartoon_header="load inline\nselect all\ncartoon\nselect *A\ncolor blue\nselect *B\ncolor red\nselect ligand\nwireframe 0.25\nselect solvent\nspacefill 0.25\nselect all\nexit\n"+rasmol_state.buf_tm.str();
    if (!mm_opt) rasmol_state.buf<<rasmol_CA_header;
    rasmol_state.buf_all<<rasmol_CA_header;
    if (!mm_opt) rasmol_state.buf_atm<<rasmol_cartoon_header;
    rasmol_state.buf_all_atm<<rasmol_cartoon_header;
    rasmol_state.buf_all_atm_lig<<rasmol_cartoon_header;
}

inline void extract_aligned_resi(FlexalignRasmolState& rasmol_state,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2,
    const int mm_opt)
{
    int i;
    rasmol_state.i1=-1;
    rasmol_state.i2=-1;
    if (!mm_opt)
    {
        for (i=0;i<seqM.size();i++)
        {
            rasmol_state.i1+=(seqxA[i]!='-');
            rasmol_state.i2+=(seqyA[i]!='-');
            if (seqM[i]==' ') continue;
            rasmol_state.resi_aln1.push_back(resi_vec1[rasmol_state.i1].substr(0,4));
            rasmol_state.resi_aln2.push_back(resi_vec2[rasmol_state.i2].substr(0,4));
            if (seqM[i]!=':') continue;
            rasmol_state.buf    <<"select "<<rasmol_state.resi_aln1.back()<<":A,"
                   <<rasmol_state.resi_aln2.back()<<":B\ncolor red\n";
            rasmol_state.buf_all<<"select "<<rasmol_state.resi_aln1.back()<<":A,"
                   <<rasmol_state.resi_aln2.back()<<":B\ncolor red\n";
        }
        rasmol_state.buf<<"select all\nexit\n"<<rasmol_state.buf_tm.str();
    }
    rasmol_state.buf_all<<"select all\nexit\n"<<rasmol_state.buf_tm.str();
}

void read_rasmol_structure1(FlexalignRasmolState& rasmol_state,
    const std::string& xname, const DoubleMatrix& tu_vec,
    const int mm_opt, const int ter_opt, const int mirror_opt);

void read_rasmol_structure2(FlexalignRasmolState& rasmol_state,
    const std::string& yname, const int mm_opt, const int ter_opt);

inline void write_rasmol_to_file(const FlexalignRasmolState& rasmol_state,
    const std::string& fname_super, const int mm_opt)
{
    std::ofstream fp;
    if (!mm_opt)
    {
        fp.open(fname_super);
        fp << rasmol_state.buf.str();
        fp.close();
    }
    fp.open(fname_super + "_all");
    fp << rasmol_state.buf_all.str();
    fp.close();
    if (!mm_opt)
    {
        fp.open(fname_super + "_atm");
        fp << rasmol_state.buf_atm.str();
        fp.close();
    }
    fp.open(fname_super + "_all_atm");
    fp << rasmol_state.buf_all_atm.str();
    fp.close();
    fp.open(fname_super + "_all_atm_lig");
    fp << rasmol_state.buf_all_atm_lig.str();
    fp.close();
}

void output_flexalign_rasmol(const std::string xname, const std::string yname,
    const std::string fname_super,const DoubleMatrix&tu_vec,
    const int ter_opt,
    const int mm_opt, const int split_opt, const int mirror_opt,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>&resi_vec1, const std::vector<std::string>&resi_vec2,
    const std::string chainID1, const std::string chainID2,
    const int xlen, const int ylen, const double d0A, const int n_ali8,
    const double rmsd, const double TM1, const double Liden);

struct FlexalignPymolState
{
    std::stringstream buf;
    std::map<std::string, int> resi2hinge_dict;
    std::string chain1_sele;
    std::string chain2_sele;
    std::string resi1_sele;
    std::string resi2_sele;
    std::string resi1_bond;
    std::string resi2_bond;
    int infmt;
};

struct FlexalignPymolReadState
{
    std::stringstream& buf;
    std::map<std::string, int>& resi2hinge_dict;
    std::ifstream fin;
    CompressedStream fin_gz;
    std::string line;
    std::string resi;
    std::string asym_id = "."; // this is similar to chainID, except that
                          // chainID is char while asym_id is a string
                          // with possibly multiple char
    Vec3 t;
    RotMat u;
    int compress_type;
    int infmt;
    int read_resi;

    FlexalignPymolReadState(FlexalignPymolState& pymol_state)
        : buf(pymol_state.buf), resi2hinge_dict(pymol_state.resi2hinge_dict) {}
};

inline bool read_pymol_line(FlexalignPymolReadState& read_state)
{
    if (read_state.compress_type)
    {
        if (!read_state.fin_gz.good())
            return false;
        getline(read_state.fin_gz, read_state.line);
    }
    else
    {
        if (!read_state.fin.good())
            return false;
        getline(read_state.fin, read_state.line);
    }
    return true;
}

inline void read_pymol_line_or_quit(FlexalignPymolReadState& read_state,
    const std::string& xname)
{
    if (!read_pymol_line(read_state))
        PrintErrorAndQuit("ERROR! Unexpected end of " + xname);
}

inline int hinge_char_to_index(char hinge_char)
{
    if (hinge_char >= '0' && hinge_char <= '9')
        return hinge_char - '0';
    if (hinge_char >= 'a' && hinge_char <= 'z')
        return hinge_char - 'a' + 10;
    if (hinge_char >= 'A' && hinge_char <= 'Z')
        return hinge_char - 'A' + 36;
    return 0;
}

inline char find_nearest_hinge_char(const std::string& seqM,
    const int aln_pos, const int xlen)
{
    int ali_len = seqM.size();
    for (int i = 1; i < ali_len; i++)
    {
        if (aln_pos - i >= 0 && seqM[aln_pos - i] != ' ')
            return seqM[aln_pos - i];
        if (aln_pos + i < xlen && seqM[aln_pos + i] != ' ')
            return seqM[aln_pos + i];
    }
    return ' ';
}

inline void build_pymol_hinge_map(FlexalignPymolState& pymol_state,
    const std::string& seqM, const std::string& seqxA,
    const std::vector<std::string>& resi_vec1)
{
    int resi_idx = -1;
    const int xlen = resi_vec1.size();
    for (int aln_pos = 0; aln_pos < seqxA.size(); aln_pos++)
    {
        if (seqxA[aln_pos] == '-')
            continue;
        resi_idx++;
        char hinge_char = seqM[aln_pos];
        if (hinge_char == ' ')
            hinge_char = find_nearest_hinge_char(seqM, aln_pos, xlen);
        pymol_state.resi2hinge_dict[resi_vec1[resi_idx]] = hinge_char_to_index(hinge_char);
    }
}

inline void read_pymol_pdb_atom_line(FlexalignPymolReadState& read_state,
    const DoubleMatrix& tu_vec, const int mirror_opt)
{
    std::stringstream& buf = read_state.buf;
    read_state.infmt = 0;
    Vec3 x;  // before transform
    x[0] = atof(read_state.line.substr(30, 8).c_str());
    x[1] = atof(read_state.line.substr(38, 8).c_str());
    x[2] = atof(read_state.line.substr(46, 8).c_str());
    if (mirror_opt)
        x[2] = -x[2];
    if (read_state.read_resi == 1)
        read_state.resi = read_state.line.substr(22, 5);
    else
        read_state.resi = read_state.line.substr(22, 5) + read_state.line[21];
    int hinge = 0;
    if (read_state.resi2hinge_dict.count(read_state.resi))
        hinge = read_state.resi2hinge_dict[read_state.resi];
    tu2t_u(tu_vec[hinge], read_state.t, read_state.u);
    Vec3 x1; // after transform
    transform(read_state.t, read_state.u, x, x1);
    buf << read_state.line.substr(0, 30) << std::setiosflags(std::ios::fixed)
        << std::setprecision(3)
        << std::setw(8) << x1[0] << std::setw(8) << x1[1] << std::setw(8) << x1[2]
        << read_state.line.substr(54) << '\n';
}

void read_pymol_mmcif_loop(FlexalignPymolReadState& read_state,
    const std::string& xname,
    const DoubleMatrix& tu_vec, const int mirror_opt);

void read_pymol_structure(FlexalignPymolState& pymol_state,
    const std::string& xname, const std::vector<std::string>& resi_vec1,
    const DoubleMatrix& tu_vec, const int mirror_opt, const int ter_opt);

inline void build_pymol_chain_suffix(FlexalignPymolState& pymol_state,
    const std::string& chainID1, const std::string& chainID2,
    const int mm_opt, const int split_opt, const int ter_opt)
{
    int i;
    if (!mm_opt)
    {
        if (split_opt == 2 && ter_opt >= 1) // align one chain from model 1
        {
            pymol_state.chain1_sele = " and c. " + chainID1.substr(1);
            pymol_state.chain2_sele = " and c. " + chainID2.substr(1);
        }
        else if (split_opt == 2 && ter_opt == 0) // align one chain from each model
        {
            for (i = 1; i < chainID1.size(); i++)
                if (chainID1[i] == ',')
                    break;
            pymol_state.chain1_sele = " and c. " + chainID1.substr(i + 1);
            for (i = 1; i < chainID2.size(); i++)
                if (chainID2[i] == ',')
                    break;
            pymol_state.chain2_sele = " and c. " + chainID2.substr(i + 1);
        }
    }
}

void build_pymol_resi_suffix(FlexalignPymolState& pymol_state,
    const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
    const std::vector<std::string>& resi_vec1, const std::vector<std::string>& resi_vec2,
    const int mm_opt);

void write_pymol_to_file(const FlexalignPymolState& pymol_state,
    const std::string& fname_super,
    const std::string& yname, const int mm_opt, const int ter_opt);

void output_flexalign_pymol(const std::string xname, const std::string yname,
                            const std::string fname_super, const std::vector<std::vector<double> > &tu_vec,
                            const int ter_opt,
                            const int mm_opt, const int split_opt, const int mirror_opt,
                            const std::string& seqM, const std::string& seqxA, const std::string& seqyA,
                            const std::vector<std::string> &resi_vec1, const std::vector<std::string> &resi_vec2,
                            const std::string chainID1, const std::string chainID2);

struct FlexAlignResult : ChainPairAlignResult
{
    std::vector<std::vector<double> > tu_vec;
    int hingeNum;

    FlexAlignResult()
    {
        TM1 = TM2 = TM3 = TM4 = TM5 = -1.0;
        d0_0 = TM_0 = d0A = d0B = d0u = d0a = 0.0;
        d0_out = 5.0;
        rmsd0 = Liden = TM_ali = rmsd_ali = 0.0;
        L_ali = n_ali = n_ali8 = 0;
        hingeNum = 0;
        init_rt_mat(t0, u0);
    }
};

enum FlexAlignMode
{
    FLEX_BEST = 0,
    FLEX_USBCAT = 1
};

struct FlexalignParams
{
    FlexAlignMode mode;
    int hinge_opt;
    bool hinge_set;
    double TMpass;
    int sparse_val;
};

void output_flexalign_results(
    const ParsedChain& chain1_data, const ParsedChain& chain2_data,
    const FlexAlignResult& flexalign_res,
    const UserOptions& opts,
    const int mm_opt, std::ostream& os = std::cout);


inline void reassign_segments_by_distance(CoordArray& xa, CoordArray& ya,
    CoordArray& xt, std::vector<int>& invmap, DoubleMatrix& tu_vec,
    Vec3& t0, RotMat& u0, const int xlen, const int ylen,
    std::string &seqM, const std::string &seqyA,
    std::vector<char>& seqM_char, std::vector<double>& di_vec)
{
    for (int hinge = (int)tu_vec.size() - 1; hinge >= 0; hinge--)
    {
        tu2t_u(tu_vec[hinge], t0, u0);
        do_rotation(xa, xt, xlen, t0, u0);
        for (int j = 0; j < ylen; j++)
        {
            int i = invmap[j];
            if (i < 0) continue;
            double d = sqrt(dist(xt[i], ya[j]));
            if (di_vec[j] < 0 || d <= di_vec[j])
            {
                di_vec[j] = d;
                seqM_char[j] = hinge + '0';
            }
        }
    }

    int y_idx = -1;
    for (size_t r = 0; r < seqM.size(); r++)
    {
        if (seqyA[r] == '-') continue;
        y_idx++;
        seqM[r] = seqM_char[y_idx];
    }
}

inline void smooth_singleton_insert(std::string &seqM, std::vector<char>& seqM_char,
    const std::string &seqyA, const int n_segments)
{
    for (int hinge=n_segments-1; hinge>=0; hinge--)
    {
        char cur_seg = hinge + '0';
        char left_seg, right_seg;
        int  y_idx = -1;

        for (size_t r=0; r<seqM.size(); r++)
        {
            if (seqyA[r]=='-') continue;
            y_idx++;

            if (seqM_char[y_idx] != cur_seg) continue;

            left_seg  = (r > 0) ? seqM[r-1] : ' ';
            right_seg = (r < seqM.size()-1) ? seqM[r+1] : ' ';

            if (right_seg == cur_seg || right_seg == ' ') continue;
            if (left_seg  == cur_seg || left_seg  == ' ') continue;
            if (left_seg != right_seg) continue;

            char neighbor_seg = (left_seg != ' ') ? left_seg : right_seg;
            seqM[r] = seqM_char[y_idx] = neighbor_seg;
        }
    }
}

inline char find_nearest_seg_char(const std::string &seqM, int r, int step)
{
    for (int i = r + step; i >= 0 && i < (int)seqM.size(); i += step) {
        if (seqM[i] != ' ') return seqM[i];
    }
    return ' ';
}

inline char find_nearest_seg_and_dist(const std::string &seqM, int r, int step, int &dist)
{
    for (int i = r + step; i >= 0 && i < (int)seqM.size(); i += step) {
        if (seqM[i] != ' ') {
            dist = (step > 0) ? (i - r) : (r - i);
            return seqM[i];
        }
    }
    dist = 0;
    return ' ';
}

inline void smooth_singleton_at_end(std::string &seqM, std::vector<char>& seqM_char,
    const std::string &seqyA, const int n_segments)
{
    for (int hinge=n_segments-1; hinge>=0; hinge--)
    {
        char cur_seg = hinge + '0';
        int  y_idx   = -1;

        for (size_t r=0; r<seqM.size(); r++)
        {
            if (seqyA[r]=='-') continue;
            y_idx++;

            if (seqM[r] != cur_seg) continue;

            if (r > 0 && seqM[r-1]==' ' && r+1 < seqM.size() && seqM[r+1]==' ')
                continue;

            char left_hinge  = find_nearest_seg_char(seqM, r, -1);
            if (left_hinge == cur_seg) continue;

            char right_hinge = find_nearest_seg_char(seqM, r, +1);
            if (right_hinge == cur_seg) continue;

            if (left_hinge!=' ' && right_hinge!=' ' && left_hinge!=right_hinge)
                continue;

            char neighbor_seg = (right_hinge != ' ') ? right_hinge : left_hinge;
            seqM[r] = seqM_char[y_idx] = neighbor_seg;
        }
    }
}

inline void smooth_dimer_insert(std::string &seqM, std::vector<char>& seqM_char,
    const std::string &seqyA, const int n_segments)
{
    for (int hinge=n_segments-1; hinge>=0; hinge--)
    {
        char cur_seg = hinge + '0';
        int  y_idx   = -1;

        for (size_t r=0; r+1<seqM.size(); r++)
        {
            if (seqyA[r]=='-') continue;
            y_idx++;

            if (seqM[r] != cur_seg || seqM[r+1] != cur_seg) continue;

            if (r+2 < seqM.size() && (seqM[r+2]==' ' || seqM[r+2]==cur_seg)) continue;
            if (r > 0 && (seqM[r-1]==' ' || seqM[r-1]==cur_seg)) continue;
            if (r+2 < seqM.size() && r > 0 && seqM[r-1] != seqM[r+2]) continue;

            char neighbor_seg = (r > 0) ? seqM[r-1] : seqM[r+2];
            seqM[r]   = seqM_char[y_idx]       = neighbor_seg;
            seqM[r+1] = seqM_char[y_idx + 1]   = neighbor_seg;
        }
    }
}

inline void smooth_disconnected_singleton(std::string &seqM, std::vector<char>& seqM_char,
    const std::string &seqyA, const int n_segments)
{
    for (int hinge=n_segments-1; hinge>=0; hinge--)
    {
        char cur_seg = hinge + '0';
        int  y_idx   = -1;

        for (size_t r=0; r<seqM.size(); r++)
        {
            if (seqyA[r]=='-') continue;
            y_idx++;
            if (seqM[r] != cur_seg) continue;

            int left_dist, right_dist;
            char left_hinge  = find_nearest_seg_and_dist(seqM, r, -1, left_dist);
            if (left_hinge == cur_seg) continue;

            char right_hinge = find_nearest_seg_and_dist(seqM, r, +1, right_dist);
            if (right_hinge == cur_seg) continue;

            if (right_hinge == ' ')
                seqM[r] = seqM_char[y_idx] = left_hinge;
            else if (left_hinge == ' ')
                seqM[r] = seqM_char[y_idx] = right_hinge;
            else
                seqM[r] = seqM_char[y_idx] = (left_dist < right_dist) ? left_hinge : right_hinge;
        }
    }
}

inline void smooth_segment_assignment(std::string &seqM, std::vector<char>& seqM_char,
    const std::string &seqyA, const int n_segments)
{
    smooth_singleton_insert(seqM, seqM_char, seqyA, n_segments);
    smooth_singleton_at_end(seqM, seqM_char, seqyA, n_segments);
    smooth_dimer_insert(seqM, seqM_char, seqyA, n_segments);
    smooth_disconnected_singleton(seqM, seqM_char, seqyA, n_segments);
}

inline void recompute_flexalign_scores(const std::string &seqM, const std::string &seqxA,
    const std::string &seqyA, std::vector<char>& seqM_char, std::vector<double>& di_vec,
    const std::vector<int>& invmap, DoubleMatrix& tu_vec, CoordArray& xa,
    CoordArray& ya, CoordArray& xt, Vec3& t0, RotMat& u0,
    const int xlen, const int ylen,
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &rmsd0, double &Liden,
    const double d0A, const double d0B, const double d0a, const double d0u,
    const double d0_scale, const double Lnorm_ass,
    const int a_opt, const bool u_opt, const bool d_opt, const int n_ali8)
{
    for (int hinge=(int)tu_vec.size()-1; hinge>=0; hinge--)
    {
        tu2t_u(tu_vec[hinge],t0,u0);
        do_rotation(xa, xt, xlen, t0, u0);
        char cur_seg = hinge + '0';
        for (int j=0; j<ylen; j++)
        {
            int i = invmap[j];
            if (i < 0) continue;
            if (seqM_char[j] != cur_seg) continue;
            double d = sqrt(dist(xt[i], ya[j]));
            if (di_vec[j] < 0 || d <= di_vec[j])
                di_vec[j] = d;
        }
    }

    rmsd0=TM1=TM2=TM3=TM4=TM5=0;
    Liden=0;
    for (int r=0; r<(int)seqM.size(); r++)
    {
        if (seqM[r]!=' ')
            Liden+=(seqxA[r]==seqyA[r]);
    }

    for(int j=0; j<ylen; j++)
    {
        int i=invmap[j];
        if(i<0) continue;
        {
            double d=di_vec[j];
            TM2+=1/(1+(d/d0B)*(d/d0B)); // chain_1
            TM1+=1/(1+(d/d0A)*(d/d0A)); // chain_2
            if (a_opt) TM3+=1/(1+(d/d0a)*(d/d0a)); // -a
            if (u_opt) TM4+=1/(1+(d/d0u)*(d/d0u)); // -u
            if (d_opt) TM5+=1/(1+(d/d0_scale)*(d/d0_scale)); // -d
            rmsd0+=d*d;
        }
    }
    TM2/=xlen;
    TM1/=ylen;
    TM3/=(xlen+ylen)*0.5;
    TM4/=Lnorm_ass;
    TM5/=ylen;
    if (n_ali8) rmsd0=sqrt(rmsd0/n_ali8);
}

inline void remove_unused_segments(const std::string &seqM, DoubleMatrix& tu_vec)
{
    for (int hinge = (int)tu_vec.size() - 1; hinge > 0; hinge--)
    {
        int afp_len = 0;
        for (size_t r = 0; r < seqM.size(); r++) afp_len += (seqM[r] == hinge + '0');
        if (afp_len) break;
        tu_vec.pop_back(); // remove unnecessary afp
    }
}


inline void score_with_cur_transform(
    CoordArray& xa, CoordArray& ya, CoordArray& xt,
    const std::string &seqx, const std::string &seqy,
    const int xlen, const int ylen,
    const std::vector<std::string> &sequence,
    const ChainPairAlignOptions& align_opts,
    const int hinge,
    ChainPairAlignResult& chain_align_res)
{
    do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
    ChainPairAlignOptions score_opts = align_opts;
    score_opts.fast_opt = false;
    score_opts.ss_opt = 0;
    se_main(xt, ya, seqx, seqy, chain_align_res, xlen, ylen, sequence, score_opts, 0, hinge);
}

#define NUM_SPLIT_VARIANTS 2

struct AlignSubstruct {
    std::string seqx_hinge_trimmed;
    std::string seqy_hinge_trimmed;
    std::string secx_hinge_trimmed;
    std::string secy_hinge_trimmed;
    CoordArray xa_hinge_trimmed;
    CoordArray ya_hinge_trimmed;
    int xlen_hinge_trimmed;
    int ylen_hinge_trimmed;
};

struct SubstructAlignRes {
    double TM1;
    double TM2;
    double TM3;
    double TM4;
    double TM5;
    double rmsd0;
    double d0_0;
    double TM_0;
    double d0A;
    double d0B;
    double d0u;
    double d0a;
    double d0_out;
    std::string seqM;
    std::string seqxA;
    std::string seqyA;
    int n_ali;
    int n_ali8;
    int L_ali;
    double Liden;
    double TM_ali;
    double rmsd_ali;
    std::vector<int> invmap;
    Vec3 t_u_result;
    RotMat u_u_result;
};

inline void extract_split_substruct(
        const CoordArray& xa, const CoordArray& ya,
        const std::string &seqx, const std::string &seqy,
        const std::string &secx, const std::string &secy,
        const std::string &seqxA, const std::string &seqyA,
        const int xlen, const int ylen, const int n_ali8,
        AlignSubstruct sub[NUM_SPLIT_VARIANTS])
{
    AlignSubstruct &sub0 = sub[0];
    AlignSubstruct &sub1 = sub[1];
    sub0.seqx_hinge_trimmed.clear();
    sub0.seqy_hinge_trimmed.clear();
    sub0.secx_hinge_trimmed.resize(xlen + 1);
    sub0.secy_hinge_trimmed.resize(ylen + 1);
    sub0.xa_hinge_trimmed.resize(xlen);
    sub0.ya_hinge_trimmed.resize(ylen);
    sub0.xlen_hinge_trimmed = n_ali8;
    sub0.ylen_hinge_trimmed = ylen - n_ali8;
    sub1.seqx_hinge_trimmed.clear();
    sub1.seqy_hinge_trimmed.clear();
    sub1.secx_hinge_trimmed.resize(xlen + 1);
    sub1.secy_hinge_trimmed.resize(ylen + 1);
    sub1.xa_hinge_trimmed.resize(xlen);
    sub1.ya_hinge_trimmed.resize(ylen);
    sub1.xlen_hinge_trimmed = xlen - n_ali8;
    sub1.ylen_hinge_trimmed = n_ali8;
    int i = -1;
    int j = -1;
    int a0 = 0;
    int b0 = 0;
    int a1 = 0;
    int b1 = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        i += (seqxA[r] != '-');
        j += (seqyA[r] != '-');
        if (seqxA[r] != '-' && seqyA[r] != '-')
        {
            sub0.seqx_hinge_trimmed += seqx[i];
            sub0.secx_hinge_trimmed[a0] = secx[i];
            sub0.xa_hinge_trimmed[a0][0] = xa[i][0];
            sub0.xa_hinge_trimmed[a0][1] = xa[i][1];
            sub0.xa_hinge_trimmed[a0][2] = xa[i][2];
            a0++;
            sub1.seqy_hinge_trimmed += seqy[j];
            sub1.secy_hinge_trimmed[b1] = secy[j];
            sub1.ya_hinge_trimmed[b1][0] = ya[j][0];
            sub1.ya_hinge_trimmed[b1][1] = ya[j][1];
            sub1.ya_hinge_trimmed[b1][2] = ya[j][2];
            b1++;
        }
        else if (seqxA[r] != '-')
        {
            sub1.seqx_hinge_trimmed += seqx[i];
            sub1.secx_hinge_trimmed[a1] = secx[i];
            sub1.xa_hinge_trimmed[a1][0] = xa[i][0];
            sub1.xa_hinge_trimmed[a1][1] = xa[i][1];
            sub1.xa_hinge_trimmed[a1][2] = xa[i][2];
            a1++;
        }
        else if (seqyA[r] != '-')
        {
            sub0.seqy_hinge_trimmed += seqy[j];
            sub0.secy_hinge_trimmed[b0] = secy[j];
            sub0.ya_hinge_trimmed[b0][0] = ya[j][0];
            sub0.ya_hinge_trimmed[b0][1] = ya[j][1];
            sub0.ya_hinge_trimmed[b0][2] = ya[j][2];
            b0++;
        }
    }
}

inline void score_split_substruct(
        AlignSubstruct& sub,
        CoordArray& xa, CoordArray& ya,
        const std::string &seqx, const std::string &seqy,
        std::vector<double>& do_vec,
        const int xlen, const int ylen,
        const std::vector<std::string> &sequence,
        const ChainPairAlignOptions& align_opts,
        Vec3& t0, RotMat& u0,
        SubstructAlignRes& sub_align_res,
        std::ostream& os = std::cout)
{
    ChainPairAlignResult chain_align_res = { 0};
    chain_align_res.d0_out = 5.0;
    chain_align_res.t0 = t0;
    chain_align_res.u0 = u0;
    chain_align_res.invmap.assign(ylen + 1, -1);
    chain_align_res.do_vec = do_vec;
    TMalign_main(sub.xa_hinge_trimmed, sub.ya_hinge_trimmed, sub.seqx_hinge_trimmed, sub.seqy_hinge_trimmed, sub.secx_hinge_trimmed, sub.secy_hinge_trimmed,
        chain_align_res, sub.xlen_hinge_trimmed, sub.ylen_hinge_trimmed, sequence, align_opts, os);
    t0 = chain_align_res.t0;
    u0 = chain_align_res.u0;
    do_vec = chain_align_res.do_vec;
    sub_align_res.t_u_result = chain_align_res.t0;
    sub_align_res.u_u_result = chain_align_res.u0;
    chain_align_res.TM1 = chain_align_res.TM2 = chain_align_res.TM3 = chain_align_res.TM4 = chain_align_res.TM5 = chain_align_res.rmsd0 = 0;
    chain_align_res.seqM = "";
    chain_align_res.seqxA = "";
    chain_align_res.seqyA = "";
    chain_align_res.n_ali = chain_align_res.n_ali8 = 0;
    CoordArray xt;
    xt.resize(xlen);
    score_with_cur_transform(
        xa, ya, xt, seqx, seqy, xlen, ylen,
        sequence, align_opts, 1, chain_align_res);
    do_vec = chain_align_res.do_vec;
    sub_align_res.TM1 = chain_align_res.TM1;
    sub_align_res.TM2 = chain_align_res.TM2;
    sub_align_res.TM3 = chain_align_res.TM3;
    sub_align_res.TM4 = chain_align_res.TM4;
    sub_align_res.TM5 = chain_align_res.TM5;
    sub_align_res.seqM = chain_align_res.seqM;
    sub_align_res.seqxA = chain_align_res.seqxA;
    sub_align_res.seqyA = chain_align_res.seqyA;
    sub_align_res.n_ali = chain_align_res.n_ali;
    sub_align_res.n_ali8 = chain_align_res.n_ali8;
    sub_align_res.rmsd0 = chain_align_res.rmsd0;
    sub_align_res.d0_0 = chain_align_res.d0_0;
    sub_align_res.TM_0 = chain_align_res.TM_0;
    sub_align_res.d0A = chain_align_res.d0A;
    sub_align_res.d0B = chain_align_res.d0B;
    sub_align_res.d0u = chain_align_res.d0u;
    sub_align_res.d0a = chain_align_res.d0a;
    sub_align_res.d0_out = chain_align_res.d0_out;
    sub_align_res.L_ali = chain_align_res.L_ali;
    sub_align_res.Liden = chain_align_res.Liden;
    sub_align_res.TM_ali = chain_align_res.TM_ali;
    sub_align_res.rmsd_ali = chain_align_res.rmsd_ali;
    sub_align_res.invmap = chain_align_res.invmap;
}

inline void fill_sub_struct_align_result(
        const SubstructAlignRes sub_align_res[NUM_SPLIT_VARIANTS],
        DoubleMatrix& tu_vec,
        ChainPairAlignResult& chain_align_res)
{
    double TM_h=(sub_align_res[0].TM1>sub_align_res[0].TM2)?sub_align_res[0].TM1:sub_align_res[0].TM2;
    double TM  =(sub_align_res[1].TM1>sub_align_res[1].TM2)?sub_align_res[1].TM1:sub_align_res[1].TM2;
    const SubstructAlignRes* best_align;
    if (TM_h>TM)
        best_align=&sub_align_res[0];
    else
        best_align=&sub_align_res[1];
    chain_align_res.t0 = best_align->t_u_result;
    chain_align_res.u0 = best_align->u_u_result;
    t_u2tu(chain_align_res.t0, chain_align_res.u0, tu_vec[0]);
    chain_align_res.TM1=best_align->TM1;
    chain_align_res.TM2=best_align->TM2;
    chain_align_res.TM3=best_align->TM3;
    chain_align_res.TM4=best_align->TM4;
    chain_align_res.TM5=best_align->TM5;
    chain_align_res.seqM=best_align->seqM;
    chain_align_res.seqxA=best_align->seqxA;
    chain_align_res.seqyA=best_align->seqyA;
    chain_align_res.rmsd0=best_align->rmsd0;
    chain_align_res.n_ali=best_align->n_ali;
    chain_align_res.n_ali8=best_align->n_ali8;
    chain_align_res.L_ali=best_align->L_ali;
    chain_align_res.Liden=best_align->Liden;
    chain_align_res.TM_ali=best_align->TM_ali;
    chain_align_res.rmsd_ali=best_align->rmsd_ali;
    chain_align_res.d0_0=best_align->d0_0;
    chain_align_res.TM_0=best_align->TM_0;
    chain_align_res.d0A=best_align->d0A;
    chain_align_res.d0B=best_align->d0B;
    chain_align_res.d0u=best_align->d0u;
    chain_align_res.d0a=best_align->d0a;
    chain_align_res.d0_out=best_align->d0_out;
    chain_align_res.invmap = best_align->invmap;
}

inline void refine_via_split_alignment(
        CoordArray& xa,
        CoordArray& ya,
        DoubleMatrix& tu_vec,
        const std::string &seqx,
        const std::string &seqy,
        const std::string &secx,
        const std::string &secy,
        const int xlen,
        const int ylen,
        const std::vector<std::string> &sequence,
        const ChainPairAlignOptions& align_opts,
        ChainPairAlignResult& chain_align_res,
        std::ostream& os = std::cout)
{
    AlignSubstruct sub_struct[NUM_SPLIT_VARIANTS];
    extract_split_substruct(xa, ya, seqx, seqy, secx, secy,
        chain_align_res.seqxA, chain_align_res.seqyA, xlen, ylen, chain_align_res.n_ali8, sub_struct);
    SubstructAlignRes sub_align_res[NUM_SPLIT_VARIANTS];
    score_split_substruct(sub_struct[0], xa, ya, seqx, seqy, chain_align_res.do_vec,
        xlen, ylen, sequence, align_opts,
        chain_align_res.t0, chain_align_res.u0, sub_align_res[0], os);
    score_split_substruct(sub_struct[1], xa, ya, seqx, seqy, chain_align_res.do_vec,
        xlen, ylen, sequence, align_opts,
        chain_align_res.t0, chain_align_res.u0, sub_align_res[1], os);
    fill_sub_struct_align_result(sub_align_res, tu_vec, chain_align_res);
}

static inline void extract_unalign_subcoords(
        const std::string& seqxA, const std::string& seqyA,
        const std::string& seqx, const std::string& seqy,
        const std::string& secx, const std::string& secy,
        const CoordArray& xa, const CoordArray& ya,
        const int xlen_hinge_trimmed, const int ylen_hinge_trimmed,
        CoordArray& xa_hinge_trimmed, CoordArray& ya_hinge_trimmed,
        std::string& seqx_hinge_trimmed, std::string& seqy_hinge_trimmed,
        std::string& secx_hinge_trimmed, std::string& secy_hinge_trimmed,
        std::vector<int>& r1toi, std::vector<int>& r2toj)
{
    xa_hinge_trimmed.resize(xlen_hinge_trimmed);
    ya_hinge_trimmed.resize(ylen_hinge_trimmed);
    secx_hinge_trimmed.resize(xlen_hinge_trimmed + 1);
    secy_hinge_trimmed.resize(ylen_hinge_trimmed + 1);
    r1toi.assign(xlen_hinge_trimmed, 0);
    r2toj.assign(ylen_hinge_trimmed, 0);
    int i = -1, j = -1, r1 = 0, r2 = 0;
    for (size_t r = 0; r < seqxA.size(); r++)
    {
        i += (seqxA[r] != '-');
        j += (seqyA[r] != '-');
        if (seqyA[r] == '-')
        {
            seqx_hinge_trimmed.push_back(seqx[i]);
            secx_hinge_trimmed[r1] = secx[i];
            xa_hinge_trimmed[r1][0] = xa[i][0];
            xa_hinge_trimmed[r1][1] = xa[i][1];
            xa_hinge_trimmed[r1][2] = xa[i][2];
            r1toi[r1] = i;
            r1++;
        }
        if (seqxA[r] == '-')
        {
            seqy_hinge_trimmed.push_back(seqy[j]);
            secy_hinge_trimmed[r2] = secy[j];
            ya_hinge_trimmed[r2][0] = ya[j][0];
            ya_hinge_trimmed[r2][1] = ya[j][1];
            ya_hinge_trimmed[r2][2] = ya[j][2];
            r2toj[r2] = j;
            r2++;
        }
    }
}

struct CurHingeResult {
    ChainPairAlignResult tmalign_res;
};

// Copy main alignment output state into CurHingeResult
static inline void fill_hing_res(
        CurHingeResult& cur_hing_res,
        const double& TM1, const double& TM2, const double& TM3,
        const double& TM4, const double& TM5, const double& rmsd0,
        const std::string& seqM, const std::string& seqxA,
        const std::string& seqyA, const int& n_ali, const int& n_ali8,
        const std::vector<int>& invmap, const int ylen)
{
    cur_hing_res.tmalign_res.TM1 = TM1;  cur_hing_res.tmalign_res.TM2 = TM2;  cur_hing_res.tmalign_res.TM3 = TM3;
    cur_hing_res.tmalign_res.TM4 = TM4;  cur_hing_res.tmalign_res.TM5 = TM5;  cur_hing_res.tmalign_res.rmsd0 = rmsd0;
    cur_hing_res.tmalign_res.seqM = seqM;
    cur_hing_res.tmalign_res.seqxA = seqxA;
    cur_hing_res.tmalign_res.seqyA = seqyA;
    cur_hing_res.tmalign_res.n_ali = n_ali;
    cur_hing_res.tmalign_res.n_ali8 = n_ali8;
    cur_hing_res.tmalign_res.invmap.assign(ylen + 1, -1);
    for (int j = 0; j < ylen + 1; j++)
        cur_hing_res.tmalign_res.invmap[j] = invmap[j];
}

static inline void update_global_res(ChainPairAlignResult& chain_align_res,
        const CurHingeResult& cur_hing_res,
        DoubleMatrix& tu_vec)
{
    const ChainPairAlignResult& cur_res = cur_hing_res.tmalign_res;
    chain_align_res.TM1 = cur_res.TM1;
    chain_align_res.TM2 = cur_res.TM2;
    chain_align_res.TM3 = cur_res.TM3;
    chain_align_res.TM4 = cur_res.TM4;
    chain_align_res.TM5 = cur_res.TM5;
    chain_align_res.rmsd0 = cur_res.rmsd0;
    chain_align_res.seqM  = cur_res.seqM;
    chain_align_res.seqxA = cur_res.seqxA;
    chain_align_res.seqyA = cur_res.seqyA;
    chain_align_res.n_ali  = cur_res.n_ali;
    chain_align_res.n_ali8 = cur_res.n_ali8;
    std::vector<double> tu_tmp(12, 0);
    t_u2tu(chain_align_res.t0, chain_align_res.u0, tu_tmp);
    tu_vec.push_back(tu_tmp);
    chain_align_res.invmap = cur_res.invmap;
}

inline void search_hinge_regions(
        CoordArray& xa, CoordArray& ya, CoordArray& xt, DoubleMatrix& tu_vec,
        const std::string &seqx, const std::string &seqy, const std::string &secx,
        const std::string &secy,
        const int xlen, const int ylen, const int hinge_opt, const std::vector<std::string> &sequence,
        const ChainPairAlignOptions& align_opts,
        ChainPairAlignResult& chain_align_res,
        std::ostream& os = std::cout)
{
    int r;
    for (int r2=0;r2<chain_align_res.seqM.size();r2++)
    {
        if (chain_align_res.seqM[r2]=='1')
            chain_align_res.seqM[r2]='0';
    }

    int minlen = std::min(xlen, ylen);
    for (int hinge=0;hinge<hinge_opt;hinge++)
    {
        if (minlen-chain_align_res.n_ali8<5) break;

        int xlen_hinge_trimmed=xlen - chain_align_res.n_ali8;
        int ylen_hinge_trimmed=ylen - chain_align_res.n_ali8;
        CoordArray xa_hinge_trimmed;
        CoordArray ya_hinge_trimmed;
        std::string seqx_hinge_trimmed, seqy_hinge_trimmed, secx_hinge_trimmed, secy_hinge_trimmed;
        std::vector<int> r1toi, r2toj;

        extract_unalign_subcoords(chain_align_res.seqxA, chain_align_res.seqyA, seqx, seqy, secx, secy, xa, ya,
            xlen_hinge_trimmed, ylen_hinge_trimmed, xa_hinge_trimmed, ya_hinge_trimmed, seqx_hinge_trimmed, seqy_hinge_trimmed, secx_hinge_trimmed, secy_hinge_trimmed, r1toi, r2toj);

        CurHingeResult cur_hing_res;
        cur_hing_res.tmalign_res.t0 = chain_align_res.t0;
        cur_hing_res.tmalign_res.u0 = chain_align_res.u0;
        cur_hing_res.tmalign_res.do_vec = chain_align_res.do_vec;
        TMalign_main(xa_hinge_trimmed, ya_hinge_trimmed, seqx_hinge_trimmed, seqy_hinge_trimmed,
            secx_hinge_trimmed, secy_hinge_trimmed, cur_hing_res.tmalign_res,
            xlen_hinge_trimmed, ylen_hinge_trimmed, sequence, align_opts, os);
        chain_align_res.t0 = cur_hing_res.tmalign_res.t0;
        chain_align_res.u0 = cur_hing_res.tmalign_res.u0;
        chain_align_res.do_vec = cur_hing_res.tmalign_res.do_vec;

        do_rotation(xa, xt, xlen, chain_align_res.t0, chain_align_res.u0);
        fill_hing_res(cur_hing_res, chain_align_res.TM1, chain_align_res.TM2, chain_align_res.TM3, chain_align_res.TM4, chain_align_res.TM5, chain_align_res.rmsd0,
            chain_align_res.seqM, chain_align_res.seqxA, chain_align_res.seqyA, chain_align_res.n_ali, chain_align_res.n_ali8, chain_align_res.invmap, ylen);

        cur_hing_res.tmalign_res.d0_0 = chain_align_res.d0_0;
        cur_hing_res.tmalign_res.TM_0 = chain_align_res.TM_0;
        cur_hing_res.tmalign_res.d0A = chain_align_res.d0A;
        cur_hing_res.tmalign_res.d0B = chain_align_res.d0B;
        cur_hing_res.tmalign_res.d0u = chain_align_res.d0u;
        cur_hing_res.tmalign_res.d0a = chain_align_res.d0a;
        cur_hing_res.tmalign_res.d0_out = chain_align_res.d0_out;
        cur_hing_res.tmalign_res.do_vec = chain_align_res.do_vec;
        cur_hing_res.tmalign_res.L_ali = chain_align_res.L_ali;
        cur_hing_res.tmalign_res.Liden = chain_align_res.Liden;
        cur_hing_res.tmalign_res.TM_ali = chain_align_res.TM_ali;
        cur_hing_res.tmalign_res.rmsd_ali = chain_align_res.rmsd_ali;
        se_main(xt, ya, seqx, seqy,
            cur_hing_res.tmalign_res,
            xlen, ylen, sequence, align_opts, 0, hinge+1);
        chain_align_res.d0A = cur_hing_res.tmalign_res.d0A;
        chain_align_res.d0B = cur_hing_res.tmalign_res.d0B;
        chain_align_res.d0u = cur_hing_res.tmalign_res.d0u;
        chain_align_res.d0a = cur_hing_res.tmalign_res.d0a;
        chain_align_res.do_vec = cur_hing_res.tmalign_res.do_vec;
        chain_align_res.Liden = cur_hing_res.tmalign_res.Liden;

        int new_ali=0;
        for (r=0;r<cur_hing_res.tmalign_res.seqM.size();r++)
        {
           new_ali+=(cur_hing_res.tmalign_res.seqM[r]==hinge+'1');
        }
        if (cur_hing_res.tmalign_res.n_ali8 - chain_align_res.n_ali8<5)
        {
           new_ali=0;
        }

        if (new_ali>=5)
        {
            update_global_res(chain_align_res, cur_hing_res, tu_vec);
        }

        if (new_ali<5) break;
    }
}

inline void denoise_segments(
        CoordArray& xa, CoordArray& ya, CoordArray& xt,
        DoubleMatrix& tu_vec, const int xlen, const int ylen,
        const ChainPairAlignOptions& align_opts,
        ChainPairAlignResult& chain_align_res)
{
    std::vector<char> seqM_char(ylen,' ');
    std::vector<double> di_vec(ylen,-1);
    reassign_segments_by_distance(xa, ya, xt, chain_align_res.invmap, tu_vec, chain_align_res.t0, chain_align_res.u0,
        xlen, ylen, chain_align_res.seqM, chain_align_res.seqyA, seqM_char, di_vec);
    smooth_segment_assignment(chain_align_res.seqM, seqM_char, chain_align_res.seqyA, tu_vec.size());
    recompute_flexalign_scores(chain_align_res.seqM, chain_align_res.seqxA, chain_align_res.seqyA, seqM_char, di_vec, chain_align_res.invmap,
        tu_vec, xa, ya, xt, chain_align_res.t0, chain_align_res.u0, xlen, ylen,
        chain_align_res.TM1, chain_align_res.TM2, chain_align_res.TM3, chain_align_res.TM4, chain_align_res.TM5, chain_align_res.rmsd0, chain_align_res.Liden,
        chain_align_res.d0A, chain_align_res.d0B, chain_align_res.d0a, chain_align_res.d0u, align_opts.d0_scale, align_opts.Lnorm,
        align_opts.a_opt, align_opts.u_opt, align_opts.d_opt, chain_align_res.n_ali8);
    remove_unused_segments(chain_align_res.seqM, tu_vec);
}

struct USBCAT_AFP
{
    int i; 
    int j; 
    int len;
    double score;
};





inline bool parse_chain(
    const std::string& filename,
    const std::vector<std::vector<std::string>>& PDB_lines,
    const std::vector<std::string>& chainID_list,
    const std::vector<int>& mol_vec,
    int chain_idx,
    const std::string& mol_opt,
    int mirror_opt,
    int read_resi,
    const std::string& atom_opt,
    ParsedChain& chain)
{
    chain.filename = filename;
    int len = PDB_lines[chain_idx].size();
    if (!len)     { chain.chain_len = 0; std::cerr<<"Warning! Cannot parse file: "<<filename<<". Chain length 0."<<std::endl; return false; }
    if (len < 3)  { chain.chain_len = len; std::cerr<<"Sequence is too short <3!: "<<filename<<std::endl; return false; }

    chain.cur_complex_mol_list = mol_vec[chain_idx];
    if (mol_opt == "RNA")     chain.cur_complex_mol_list = 1;
    else if (mol_opt == "protein") chain.cur_complex_mol_list = -1;

    chain.chain_id = chainID_list[chain_idx];
    chain.chain_coords.clear();
    chain.chain_coords.reserve(len);
    chain.chain_seq.clear();
    chain.chain_sec.resize(len + 1, '\0');

    chain.chain_len = read_PDB(PDB_lines[chain_idx], chain.chain_coords, chain.chain_seq,
                       chain.resi_vec, read_resi);

    if (mirror_opt) {
        for (int r = 0; r < chain.chain_len; r++)
            chain.chain_coords[r][2] = -chain.chain_coords[r][2];
    }

    if (chain.cur_complex_mol_list > 0)
        make_sec(chain.chain_seq, chain.chain_coords, chain.chain_len, chain.chain_sec, atom_opt);
    else
        make_sec(chain.chain_coords, chain.chain_len, chain.chain_sec);

    return true;
}



void run_flexalign_main(
    CoordArray& xa, CoordArray& ya, const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen,
    const AlignCommonInput& common_inputs,
    const int mol_type, const int hinge_opt, const int ss_opt, FlexAlignResult &flexalign_res,
    std::ostream& os = std::cout);


inline void align_with_flexalign_main(
    CoordArray& xa, CoordArray& ya,
    const std::string &seqx, const std::string &seqy,
    const std::string &secx, const std::string &secy,
    int xlen, int ylen,
    const AlignCommonInput& common_inputs,
    int mol_type, int hinge_opt,
    FlexAlignResult& flexalign_res,
    std::ostream& os = std::cout)
{
    double TM_best_max = -1.0;
    for (int cur_ss_opt = 0; cur_ss_opt <= MAX_SEC_STRUCT_OPT; cur_ss_opt++)
    {
        FlexAlignResult cur_flexalign_res;
        run_flexalign_main(xa, ya, seqx, seqy, secx, secy,
            xlen, ylen, common_inputs, mol_type, hinge_opt, cur_ss_opt, cur_flexalign_res, os);
        double cur_max_TM = (cur_flexalign_res.TM1 > cur_flexalign_res.TM2) ? cur_flexalign_res.TM1 : cur_flexalign_res.TM2;
        if (cur_max_TM > TM_best_max)
        {
            TM_best_max = cur_max_TM;
            flexalign_res = cur_flexalign_res;
        }
    }
}


void flexalign_usbcat_main(
    CoordArray& xa,
    CoordArray& ya,
    const std::string &seqx,
    const std::string &seqy,
    const std::string &secx,
    const std::string &secy,
    int xlen,
    int ylen,
    const AlignCommonInput& common_inputs,
    const FlexAlignResult& flexalign_main_res,
    int mol_type, double best_global_max_TM, // Best max(TM1,TM2) threshold (shared for min_resid_num pruning)
    int sparse_val,
    const FlexalignParams& flex_params,
    FlexAlignResult& flexalign_res,
    std::ostream& os = std::cout);





inline void run_flexalign(
    const ParsedChain& chain1_data,
    const ParsedChain& chain2_data,
    const AlignCommonInput& common_inputs,
    const FlexalignParams& flex_params,
    FlexAlignResult& flexalign_res,
    std::ostream& os = std::cout)
{
    CoordArray xa = chain1_data.chain_coords;
    CoordArray ya = chain2_data.chain_coords;
    std::string seqx = chain1_data.chain_seq;
    std::string seqy = chain2_data.chain_seq;
    std::string secx = chain1_data.chain_sec;
    std::string secy = chain2_data.chain_sec;
    int xa_len = chain1_data.chain_len;
    int ya_len = chain2_data.chain_len;
    int mol_type = chain1_data.cur_complex_mol_list + chain2_data.cur_complex_mol_list;

    FlexAlignResult flexalign_main_res;
    align_with_flexalign_main(
        xa, ya, seqx, seqy, secx, secy, xa_len, ya_len,
        common_inputs, mol_type, flex_params.hinge_opt, flexalign_main_res, os);
    double best_global_max_TM = (flexalign_main_res.TM1 > flexalign_main_res.TM2) ? flexalign_main_res.TM1 : flexalign_main_res.TM2; // Best max(TM1,TM2) threshold (shared for min_resid_num pruning)

    switch (flex_params.mode)
    {
        case FLEX_USBCAT:
            if (best_global_max_TM >= flex_params.TMpass)
            {
                flexalign_res = flexalign_main_res;
                flexalign_res.hingeNum = flexalign_res.tu_vec.size();
            }
            else
            {
                flexalign_usbcat_main(
                    xa, ya, seqx, seqy, secx, secy, xa_len, ya_len, common_inputs,
                    flexalign_main_res, mol_type, best_global_max_TM, flex_params.sparse_val,
                    flex_params, flexalign_res, os);
            }
            break;
        case FLEX_BEST:
            flexalign_res = flexalign_main_res;
            break;
        default:
            PrintErrorAndQuit("ERROR! Unknown FlexAlignMode.");
            break;
    }

    return;
}

struct USBCATParams
{
    int    fragLen;
    int    max_gap;
    int    misCut;
    int    maxGapFrag;
    int    max_hinge_num;
    double resScore;
    double gap_penalty_fac;
    double disCut;
    double disSmooth;
    double twist_pen;
    double max_penalty;
    double afp_dis_cut;
};

inline void fill_usbcat_params(USBCATParams& usb_cat_para, int hinge_opt)
{
    usb_cat_para.fragLen      = 12;
    usb_cat_para.resScore     = 3.0;
    usb_cat_para.gap_penalty_fac      = -0.5;
    usb_cat_para.disCut       = 5.0;
    usb_cat_para.disSmooth    = 4.0;
    usb_cat_para.twist_pen    = -usb_cat_para.fragLen * usb_cat_para.resScore - 1;
    usb_cat_para.max_gap      = 40;
    usb_cat_para.max_penalty  = -5.0 * usb_cat_para.fragLen / 8.0;
    usb_cat_para.misCut       = 2 * usb_cat_para.fragLen;
    usb_cat_para.maxGapFrag   = usb_cat_para.fragLen + usb_cat_para.max_gap;
    usb_cat_para.afp_dis_cut  = usb_cat_para.fragLen * usb_cat_para.fragLen * (usb_cat_para.disCut * usb_cat_para.disCut);
    usb_cat_para.max_hinge_num   = hinge_opt;
}

// Local distance table: one per chain; tbl[a][d] = Euclidean distance between residue a and a+d (d <= max_residue_gap)
struct LocalDistTables
{
    int    max_residue_gap;      // Max residue offset (local window, = max_gap+2*fragLen+1 = 65)
    std::vector<std::vector<double> > chain1_dist_table;   // Local distance table for chain 1 (xa)
    std::vector<std::vector<double> > chain2_dist_table;   // Local distance table for chain 2 (ya)
};

inline void build_dist_table(const CoordArray& coords, const int max_residue_gap, DoubleMatrix& dist_table)
{
    int resid_num = (int)coords.size();
    dist_table.assign(resid_num, std::vector<double>(max_residue_gap, 0.0));
    for (int i = 0; i < resid_num; i++)
    {
        for (int j = i; j < std::min(resid_num, i + max_residue_gap); j++)
            dist_table[i][j - i] = std::sqrt(dist(coords[i], coords[j]));
    }
}

inline void build_local_dist_tables(const CoordArray& xa, const CoordArray& ya, const USBCATParams& usb_cat_para, LocalDistTables& dist_tables)
{
    dist_tables.max_residue_gap = usb_cat_para.max_gap + 2 * usb_cat_para.fragLen + 1;
    build_dist_table(xa, dist_tables.max_residue_gap, dist_tables.chain1_dist_table);
    build_dist_table(ya, dist_tables.max_residue_gap, dist_tables.chain2_dist_table);
}


inline double dtable_dist(const LocalDistTables& t, bool is_chain1, int a, int b)
{
    const std::vector<std::vector<double> >& tbl = is_chain1 ? t.chain1_dist_table : t.chain2_dist_table;
    if (a >= b) return tbl[b][a - b];
    return tbl[a][b - a];
}


inline double cal_apf_drmsd_sq(const USBCAT_AFP& prv, const USBCAT_AFP& cur,
                             const LocalDistTables& tables, const USBCATParams& usb_cat_para)
{
    double rms_sq = 0;
    for (int i_idx = 0; i_idx < usb_cat_para.fragLen; i_idx++)
    {
        for (int j_idx = 0; j_idx < usb_cat_para.fragLen; j_idx++)
        {
            double dist1 = dtable_dist(tables, true,  cur.i + i_idx, prv.i + j_idx);
            double dist2 = dtable_dist(tables, false, cur.j + i_idx, prv.j + j_idx);
            rms_sq += (dist1 - dist2) * (dist1 - dist2);
        }
    }
    return rms_sq;
}


inline double cal_twist_penalty(double drmsd, const USBCATParams& usb_cat_para)
{
    return usb_cat_para.twist_pen * std::sqrt((drmsd - usb_cat_para.disCut + usb_cat_para.disSmooth) / usb_cat_para.disSmooth);
}



static const double INF_DRMSD = 1e9;


struct AfpEdgeInfo
{
    double drmsd;      // Normalized dRMSD (INF_DRMSD for hard twist edges)
    int    hinge_inc;  // Whether this edge counts as 1 hinge (only for hard twist edges)
    double penalty;    // DP penalty: 0 / soft twist penalty / twist_pen
};


inline AfpEdgeInfo eval_cur_afp_edge(const USBCAT_AFP& prv, const USBCAT_AFP& cur,
                              const LocalDistTables& tables,
                              const USBCATParams& usb_cat_para)
{
    double drmsd_sq = cal_apf_drmsd_sq(prv, cur, tables, usb_cat_para);
    AfpEdgeInfo cur_edge;
    if (drmsd_sq >= usb_cat_para.afp_dis_cut)
    {
        cur_edge.drmsd     = INF_DRMSD;
        cur_edge.hinge_inc = 1;
        cur_edge.penalty   = usb_cat_para.twist_pen;
    }
    else
    {
        cur_edge.drmsd     = std::sqrt(drmsd_sq / (usb_cat_para.fragLen * usb_cat_para.fragLen));
        cur_edge.hinge_inc = 0;
        cur_edge.penalty   = (cur_edge.drmsd > usb_cat_para.disCut - usb_cat_para.disSmooth)
                             ? cal_twist_penalty(cur_edge.drmsd, usb_cat_para)
                             : 0.0;
    }
    return cur_edge;
}



inline double get_gap_max_penalty(const USBCAT_AFP& prev_afp, const USBCAT_AFP& cur_afp,
                                 const USBCATParams& usb_cat_para)
{
    int gap_x = cur_afp.i - (prev_afp.i + prev_afp.len);
    int gap_y = cur_afp.j - (prev_afp.j + prev_afp.len);
    int max_gap = std::max(gap_x, gap_y);

    int max_overlap_num = 0;
    if (gap_x < 0 || gap_y < 0)
        max_overlap_num = (gap_x < gap_y) ? -gap_x : -gap_y;

    double gp = usb_cat_para.gap_penalty_fac * max_overlap_num;
    if (max_gap > 0)
        gp += usb_cat_para.gap_penalty_fac * max_gap;
    if (gp < usb_cat_para.max_penalty)
        gp = usb_cat_para.max_penalty;
    return gp;
}


inline double calc_block_rmsd(const std::vector<USBCAT_AFP>& afp_list,
                              const CoordArray& xa, const CoordArray& ya)
{
    std::vector<int> r1, r2;
    for (size_t a = 0; a < afp_list.size(); a++)
    {
        for (int l = 0; l < afp_list[a].len; l++)
        {
            r1.push_back(afp_list[a].i + l);
            r2.push_back(afp_list[a].j + l);
        }
    }
    int n = (int)r1.size();
    if (n < 3)
        return 0.0;
    CoordArray p1(n), p2(n);
    for (int i = 0; i < n; i++)
    {
        p1[i][0] = xa[r1[i]][0];
        p1[i][1] = xa[r1[i]][1];
        p1[i][2] = xa[r1[i]][2];
        p2[i][0] = ya[r2[i]][0];
        p2[i][1] = ya[r2[i]][1];
        p2[i][2] = ya[r2[i]][2];
    }
    double rms_sq_sum;
    Vec3 t_tmp;
    RotMat u_tmp;
    Kabsch(p1, p2, n, 0, rms_sq_sum, t_tmp, u_tmp);
    return std::sqrt(rms_sq_sum / n);
}

struct RegionMeta
{
    int original_region_idx;
    double drmsd;
    int region_x_len; 
    int region_y_len;
};

inline bool region_meta_drmsd_desc(const RegionMeta& a, const RegionMeta& b)
{ return a.drmsd > b.drmsd; }

struct RegionBoundsAllChain
{
    std::vector<int> chain1_bounds, chain2_bounds;
    std::vector<RegionMeta> region_meta;
};

struct RegionBoundPool
{
    std::vector<RegionBoundsAllChain> region_bounds;
};


inline bool afp_less_by_i(const USBCAT_AFP& a, const USBCAT_AFP& b)
{ 
    return a.i < b.i; 
}

inline bool afp_less_by_ij(const USBCAT_AFP& a, const USBCAT_AFP& b)
{ 
    if (a.i == b.i) 
        return a.j < b.j; 
    return a.i < b.i;
}

struct AFPBlock
{
    std::vector<USBCAT_AFP> afps;
    std::vector<double> drmsd;
};

struct RegionBounds
{
    int region_x_start; 
    int region_x_end;
    int region_y_start; 
    int region_y_end;
};


struct RegionAlignResult
{
    bool valid = false;
    FlexAlignResult align_re;  // Full result
};


inline double afp_score(double resScore, int len, double rmsd, double badRmsd)
{
    double t = rmsd / badRmsd;
    return resScore * len * (1.0 - t * t);
}


inline int cur_diag_resid_num(int i, int j, int xlen, int ylen, int fragLen)
{
    return std::min(i, j)
         + std::min(xlen - (i + fragLen - 1), ylen - (j + fragLen))
         + fragLen;
}


inline std::vector<USBCAT_AFP> extract_initial_afps(
    const CoordArray& xa,
    const CoordArray& ya,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double min_resid_num,
    const double cur_rmsdCut,
    const double cur_badRmsd,
    const int sparse_val)
{

        std::vector<USBCAT_AFP> initial_afps;
        int step = sparse_val + 1;

        
        CoordArray r1(usb_cat_para.fragLen), r2(usb_cat_para.fragLen);

        for (int i = 0; i <= xlen - usb_cat_para.fragLen; i += step)
        {
            for (int j = 0; j <= ylen - usb_cat_para.fragLen; j += step)
            {
                int cur_resid_num = cur_diag_resid_num(i, j, xlen, ylen, usb_cat_para.fragLen);
                if (cur_resid_num < min_resid_num)
                    continue;

                double dist1 = tables.chain1_dist_table[i][usb_cat_para.fragLen - 1];
                double dist2 = tables.chain2_dist_table[j][usb_cat_para.fragLen - 1];

                if (std::fabs(dist1 - dist2) > 2.0 * cur_rmsdCut)
                    continue;

                for (int k = 0; k < usb_cat_para.fragLen; k++)
                {
                    r1[k][0] = xa[i + k][0];
                    r1[k][1] = xa[i + k][1];
                    r1[k][2] = xa[i + k][2];
                    r2[k][0] = ya[j + k][0];
                    r2[k][1] = ya[j + k][1];
                    r2[k][2] = ya[j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
            
                Kabsch(r1, r2, usb_cat_para.fragLen, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / usb_cat_para.fragLen);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    USBCAT_AFP afp;
                    afp.i = i;
                    afp.j = j;
                    afp.len = usb_cat_para.fragLen;
                    afp.score = afp_score(usb_cat_para.resScore, usb_cat_para.fragLen, rmsd_tmp, cur_badRmsd);
                    initial_afps.push_back(afp);
                }
            }
        }
    return initial_afps;
}

// Greedy merge of all AFPs along one diagonal
// Input: cur_group (AFPs on the same diagonal, sorted ascending by i)
// Output: merged AFP list
inline std::vector<USBCAT_AFP> merge_cur_afps_group(
    const std::vector<USBCAT_AFP>& cur_group,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_rmsdCut,
    const double cur_badRmsd)
{
    std::vector<USBCAT_AFP> merged_afps;
    int afp_num = cur_group.size();
    std::vector<bool> merged_status(afp_num, false);

    int max_merge_len = std::min(xlen, ylen);
    CoordArray r1_merge(max_merge_len), r2_merge(max_merge_len);

    for (int idx = 0; idx < afp_num; idx++)
    {
        if (merged_status[idx]) continue;

        USBCAT_AFP curr_afp = cur_group[idx];
        for (int nxt_idx = idx + 1; nxt_idx < afp_num; nxt_idx++)
        {
            USBCAT_AFP nxt_afp = cur_group[nxt_idx];
            if (nxt_afp.i > curr_afp.i + curr_afp.len)
                break;  // Sorted by i, subsequent nxt entries cannot overlap

            if (nxt_afp.i + nxt_afp.len > curr_afp.i + curr_afp.len)
            {
                int new_len = (nxt_afp.i + nxt_afp.len) - curr_afp.i;
                for (int k = 0; k < new_len; k++)
                {
                    r1_merge[k][0] = xa[curr_afp.i + k][0];
                    r1_merge[k][1] = xa[curr_afp.i + k][1];
                    r1_merge[k][2] = xa[curr_afp.i + k][2];
                    r2_merge[k][0] = ya[curr_afp.j + k][0];
                    r2_merge[k][1] = ya[curr_afp.j + k][1];
                    r2_merge[k][2] = ya[curr_afp.j + k][2];
                }

                double rms_sum_sq;
                Vec3 t_tmp;
                RotMat u_tmp;
                Kabsch(r1_merge, r2_merge, new_len, 0, rms_sum_sq, t_tmp, u_tmp);
                double rmsd_tmp = std::sqrt(rms_sum_sq / new_len);

                if (rmsd_tmp < cur_rmsdCut)
                {
                    curr_afp.len = new_len;
                    curr_afp.score = afp_score(usb_cat_para.resScore, new_len, rmsd_tmp, cur_badRmsd);
                    merged_status[nxt_idx] = true;
                }
            }
        }
        merged_afps.push_back(curr_afp);
    }
    return merged_afps;
}

inline std::vector<USBCAT_AFP> merge_afps_group(
    const std::vector<USBCAT_AFP>& initial_afps,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_rmsdCut,
    const double cur_badRmsd)
{
    int max_diagonal_idx = xlen + ylen + 1;
    std::vector<std::vector<USBCAT_AFP>> afps_groups(max_diagonal_idx);
    for (size_t k = 0; k < initial_afps.size(); k++)
    {
        int diagonal_idx = initial_afps[k].i - initial_afps[k].j + ylen;
        afps_groups[diagonal_idx].push_back(initial_afps[k]);
    }

    std::vector<USBCAT_AFP> merged_afps;
    for (int diagonal_idx = 0; diagonal_idx < max_diagonal_idx; diagonal_idx++)
    {
        if (afps_groups[diagonal_idx].empty()) continue;

        std::vector<USBCAT_AFP> cur_group = afps_groups[diagonal_idx];
        std::sort(cur_group.begin(), cur_group.end(), afp_less_by_i);

        std::vector<USBCAT_AFP> cur_merged_result = merge_cur_afps_group(
            cur_group, xa, ya, usb_cat_para, xlen, ylen, cur_rmsdCut, cur_badRmsd);
        merged_afps.insert(merged_afps.end(), cur_merged_result.begin(), cur_merged_result.end());
    }

    std::sort(merged_afps.begin(), merged_afps.end(), afp_less_by_ij);

    if (merged_afps.empty())
        return std::vector<USBCAT_AFP>();

    return merged_afps;
}

inline void propagate_index_table(std::vector<int>& table, int xlen, int ylen, int direction)
{
    for (int idx_x = 0; idx_x < xlen; idx_x++)
    {
        int curr_val = -1;
        if (direction > 0)
        {
            for (int idx_y = 0; idx_y < ylen; idx_y++)
            {
                int pos = idx_x * ylen + idx_y;
                if (table[pos] != -1) curr_val = table[pos];
                else                  table[pos] = curr_val;
            }
        }
        else
        {
            for (int idx_y = ylen - 1; idx_y >= 0; idx_y--)
            {
                int pos = idx_x * ylen + idx_y;
                if (table[pos] != -1) curr_val = table[pos];
                else                  table[pos] = curr_val;
            }
        }
    }
}

inline void build_index_tables(
    const std::vector<std::vector<std::pair<int, int>>>& xchain_idx_map_afp,
    std::vector<int>& afp_aft_index,
    std::vector<int>& afp_bef_index,
    int xlen, int ylen)
{
    for (int idx_x = 0; idx_x < xlen; idx_x++)
    {
        if (xchain_idx_map_afp[idx_x].empty()) continue;
        int cur_idx_afp_num = (int)xchain_idx_map_afp[idx_x].size();
        for (int entry_idx = 0; entry_idx < cur_idx_afp_num; entry_idx++)
        {
            int y_start = xchain_idx_map_afp[idx_x][entry_idx].first;
            afp_aft_index[idx_x * ylen + y_start] = xchain_idx_map_afp[idx_x][entry_idx].second;
            afp_bef_index[idx_x * ylen + y_start] = xchain_idx_map_afp[idx_x][entry_idx].second;
        }
    }
    propagate_index_table(afp_bef_index, xlen, ylen, +1);
    propagate_index_table(afp_aft_index, xlen, ylen, -1);
}

inline void find_candidate_prevs(const USBCAT_AFP& cur_afp,
                                 const std::vector<int>& afp_aft_index,
                                 const std::vector<int>& afp_bef_index,
                                 const int xlen, const int ylen,
                                 const USBCATParams& usb_cat_para,
                                 std::vector<int>& candidate_prevs)
{
    int x_right = cur_afp.i - usb_cat_para.fragLen;
    int x_left2 = std::max(0, x_right - usb_cat_para.misCut);
    int x_left1 = std::max(0, cur_afp.i - usb_cat_para.maxGapFrag);

    int y_right = cur_afp.j - usb_cat_para.fragLen;
    int y_left2 = std::max(0, y_right - usb_cat_para.misCut);
    int y_left1 = std::max(0, cur_afp.j - usb_cat_para.maxGapFrag);

    for (int st = 0; st < 2; st++)
    {
        int x_left, x_right_b, y_left, y_right_b;
        if (st == 0)
        {
            x_left = std::max(x_left1, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left2, 0);
            y_right_b = std::min(y_right, ylen - 1);
        }
        else
        {
            x_left = std::max(x_left2, 0);
            x_right_b = std::min(x_right, xlen - 1);
            y_left = std::max(y_left1, 0);
            y_right_b = std::min(y_left2 - 1, ylen - 1);
        }

        //Check whether the window overlaps the valid range [0, ylen-1]
        if (y_left >= ylen || y_right_b < 0) continue;

        for (int prev_i = x_left; prev_i <= x_right_b; prev_i++)
        {
            int s1 = afp_aft_index[prev_i * ylen + y_left];
            int s2 = afp_bef_index[prev_i * ylen + y_right_b];
            if (s1 != -1 && s2 != -1 && s1 <= s2)
                for (int s = s1; s <= s2; s++)
                    candidate_prevs.push_back(s);
        }
    }
}

struct DpState
{
    double score;      // best chain score ending at this AFP (inclusive)
    int    prev;       // best predecessor AFP index, -1 = chain start
    int    hinge_num;  // hinges consumed along the best chain
};

inline void update_dp_state(int afp_idx,
                            const std::vector<int>& candidate_prevs,
                            const std::vector<USBCAT_AFP>& merged_afps,
                            const LocalDistTables& tables,
                            const USBCATParams& usb_cat_para,
                            std::vector<DpState>& dp)
{
    double curr_afp_score = merged_afps[afp_idx].score;
    for (size_t cand_idx = 0; cand_idx < candidate_prevs.size(); cand_idx++)
    {
        int prev_afp = candidate_prevs[cand_idx];
        int hinge_num = dp[prev_afp].hinge_num;
        if (hinge_num > usb_cat_para.max_hinge_num)
            continue;

        double gp = get_gap_max_penalty(merged_afps[prev_afp], merged_afps[afp_idx], usb_cat_para);
        AfpEdgeInfo cur_edge = eval_cur_afp_edge(merged_afps[prev_afp], merged_afps[afp_idx], tables, usb_cat_para);
        hinge_num += cur_edge.hinge_inc;

        if (hinge_num > usb_cat_para.max_hinge_num) continue;
        double cur_score = dp[prev_afp].score + curr_afp_score + cur_edge.penalty + gp;
        if (cur_score > dp[afp_idx].score)
        {
            dp[afp_idx].score     = cur_score;
            dp[afp_idx].prev      = prev_afp;
            dp[afp_idx].hinge_num = hinge_num;
        }
    }
}

inline std::vector<int> find_best_path(
    const std::vector<USBCAT_AFP>& merged_afps,
    const std::vector<int>& afp_aft_index,
    const std::vector<int>& afp_bef_index,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen, const int ylen)
{
    int n_afps = (int)merged_afps.size();
    if (n_afps == 0)
        return std::vector<int>();

    std::vector<DpState> dp(n_afps);

    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
        dp[afp_idx] = {merged_afps[afp_idx].score, -1, 0};
    std::vector<int> candidate_prevs;
    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
    {
        candidate_prevs.clear();
        find_candidate_prevs(merged_afps[afp_idx], afp_aft_index, afp_bef_index,
                             xlen, ylen, usb_cat_para, candidate_prevs);
        update_dp_state(afp_idx, candidate_prevs, merged_afps, tables,
                        usb_cat_para, dp);
    }

    int best_end_idx = 0;
    for (int afp_idx = 1; afp_idx < n_afps; afp_idx++)
        if (dp[afp_idx].score > dp[best_end_idx].score)
            best_end_idx = afp_idx;

    std::vector<int> best_afp_path;
    for (int cur_idx = best_end_idx; cur_idx != -1; cur_idx = dp[cur_idx].prev)
        best_afp_path.push_back(cur_idx);
    std::reverse(best_afp_path.begin(), best_afp_path.end());
    return best_afp_path;
}

inline std::vector<AFPBlock> build_candidate_blocks_list(
    const std::vector<int>& path,
    const std::vector<USBCAT_AFP>& merged_afps,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para)
{
    std::vector<AFPBlock> candidate_blocks;
    if (path.empty())
        return candidate_blocks;

    AFPBlock curr_block;
    curr_block.afps.push_back(merged_afps[path[0]]);
    curr_block.drmsd.push_back(0.0);

    for (size_t k = 1; k < path.size(); k++)
    {
        USBCAT_AFP curr_afp = merged_afps[path[k]];
        USBCAT_AFP prv = merged_afps[path[k - 1]];
        double drmsd = eval_cur_afp_edge(prv, curr_afp, tables, usb_cat_para).drmsd;

        if (drmsd >= usb_cat_para.disCut)
        {  //Open a new block
            candidate_blocks.push_back(curr_block);
            curr_block.afps.clear();
            curr_block.drmsd.clear();
            curr_block.afps.push_back(curr_afp);
            curr_block.drmsd.push_back(0.0);
        }
        else
        {
            curr_block.afps.push_back(curr_afp);
            curr_block.drmsd.push_back(drmsd);
        }
    }
    if (!curr_block.afps.empty()) candidate_blocks.push_back(curr_block);

    return candidate_blocks;
}

inline void build_prev_lookup_tables(
    const std::vector<USBCAT_AFP>& merged_afps,
    std::vector<int>& afp_aft_index,
    std::vector<int>& afp_bef_index,
    const int xlen, const int ylen)
{
    int n_afps = (int)merged_afps.size();
    std::vector<std::vector<std::pair<int, int> > > xchain_idx_map_afp(xlen);
    for (int afp_idx = 0; afp_idx < n_afps; afp_idx++)
        xchain_idx_map_afp[merged_afps[afp_idx].i].push_back(std::make_pair(merged_afps[afp_idx].j, afp_idx));
    afp_aft_index.assign(xlen * ylen, -1);
    afp_bef_index.assign(xlen * ylen, -1);
    build_index_tables(xchain_idx_map_afp, afp_aft_index, afp_bef_index, xlen, ylen);
}

inline void split_candidate_blocks(
    std::vector<AFPBlock>& candidate_blocks,
    const CoordArray& xa,
    const CoordArray& ya,
    const USBCATParams& usb_cat_para,
    const double cur_local_badRmsd)
{
    const size_t max_blocks = (size_t)(usb_cat_para.max_hinge_num + 1);
    bool can_split_more = candidate_blocks.size() < max_blocks; 
    while (can_split_more)
    {
        bool splitted = false; 
        double max_rmsd = 0.0;
        int target_b = -1;
    
        for (size_t block_idx = 0; block_idx < candidate_blocks.size(); block_idx++)
        {
            if (candidate_blocks[block_idx].afps.size() > 2)
            {
                double cur_rmsd = calc_block_rmsd(candidate_blocks[block_idx].afps, xa, ya);
                if (cur_rmsd > max_rmsd)
                {
                    max_rmsd = cur_rmsd;
                    target_b = block_idx;
                }
            }
        }
        
        if (max_rmsd >= cur_local_badRmsd && target_b != -1)
        {
            AFPBlock& target_blk = candidate_blocks[target_b]; 
            double max_drmsd = 0;
            int cut_afp_idx = 0;
            for (size_t idx = 1; idx < target_blk.afps.size(); idx++)
            {
                if (target_blk.drmsd[idx] > max_drmsd)
                {
                    max_drmsd = target_blk.drmsd[idx];
                    cut_afp_idx = idx;
                }
            }
            
            if (cut_afp_idx > 0)
            {
                AFPBlock right_blk;
                right_blk.afps.assign(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                right_blk.drmsd.assign(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                right_blk.drmsd[0] = 0.0;
                target_blk.afps.erase(target_blk.afps.begin() + cut_afp_idx, target_blk.afps.end());
                target_blk.drmsd.erase(target_blk.drmsd.begin() + cut_afp_idx, target_blk.drmsd.end());
                candidate_blocks.insert(candidate_blocks.begin() + target_b + 1, right_blk);
                splitted = true;
            }
        }
        can_split_more = splitted && (candidate_blocks.size() < max_blocks);
    }
}

inline void remove_single_elem_block(
    std::vector<AFPBlock>& candidate_blocks,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen)
{
    for (int block_idx = 0; block_idx < (int)candidate_blocks.size(); block_idx++)
    {
        if (candidate_blocks[block_idx].afps.size() <= 1)
        {
            int x_right = (block_idx < (int)candidate_blocks.size() - 1) ? candidate_blocks[block_idx + 1].afps.front().i : xlen;
            int x_left = (block_idx > 0) ? candidate_blocks[block_idx - 1].afps.back().i + candidate_blocks[block_idx - 1].afps.back().len : 0;


            int y_right = (block_idx < (int)candidate_blocks.size() - 1) ? candidate_blocks[block_idx + 1].afps.front().j : ylen;
            int y_left = (block_idx > 0) ? candidate_blocks[block_idx - 1].afps.back().j + candidate_blocks[block_idx - 1].afps.back().len : 0;

            int span = std::min(x_right - x_left, y_right - y_left);
            if (span < 2 * usb_cat_para.fragLen)
            {
                candidate_blocks.erase(candidate_blocks.begin() + block_idx);
                block_idx--; //Element shifts left after erase; step back one to re-check the same index
            }
        }
    }
}

inline void merge_adjacent_blocks(
    std::vector<AFPBlock>& candidate_blocks,
    const CoordArray& xa,
    const CoordArray& ya,
    const double cur_local_badRmsd)
{
    bool merged = true;
    while (merged && candidate_blocks.size() > 1)
    {
        merged = false;
        double min_rmsd = 1e9;
        int min_b = -1;
        for (size_t block_idx = 0; block_idx < candidate_blocks.size() - 1; block_idx++)
        {
            std::vector<USBCAT_AFP> temp_merged = candidate_blocks[block_idx].afps;
            temp_merged.insert(temp_merged.end(),
                              candidate_blocks[block_idx + 1].afps.begin(),
                              candidate_blocks[block_idx + 1].afps.end());
            double cur_rmsd = calc_block_rmsd(temp_merged, xa, ya);
            if (cur_rmsd < min_rmsd)
            {
                min_rmsd = cur_rmsd;
                min_b = (int)block_idx;
            }
        }
        if (min_rmsd < cur_local_badRmsd && min_b != -1)
        {
            candidate_blocks[min_b].afps.insert(candidate_blocks[min_b].afps.end(),
                                               candidate_blocks[min_b + 1].afps.begin(),
                                               candidate_blocks[min_b + 1].afps.end());
            candidate_blocks.erase(candidate_blocks.begin() + min_b + 1);
            merged = true;
        }
    }
}

inline void calc_region_from_block(
    const AFPBlock& block,
    int& resi_idx_x,
    int& resi_idx_y,
    RegionBounds& out_region)
{
    int region_x_start = -1, region_x_end = -1, region_y_start = -1, region_y_end = -1;
    for (size_t afp_idx = 0; afp_idx < block.afps.size(); afp_idx++)
    {
        USBCAT_AFP afp = block.afps[afp_idx];
        int overlap_len_x = resi_idx_x - afp.i;
        int overlap_len_y = resi_idx_y - afp.j;
        int gap_len = std::max(std::max(overlap_len_x, overlap_len_y), 0);
        if (gap_len >= afp.len)
            continue;

        int eff_i = afp.i + gap_len;
        int eff_j = afp.j + gap_len;
        int eff_L = afp.len - gap_len;
        if (region_x_start == -1)
        {
            region_x_start = eff_i;
            region_y_start = eff_j;
        }
        region_x_end = eff_i + eff_L;
        region_y_end = eff_j + eff_L;
        resi_idx_x = region_x_end;
        resi_idx_y = region_y_end;
    }
    out_region = {region_x_start, region_x_end, region_y_start, region_y_end};
}

inline void build_usbcat_regions(
    const std::vector<AFPBlock>& candidate_blocks,
    std::vector<RegionBounds>& usbcat_regions_list)
{
    int resi_idx_x = 0, resi_idx_y = 0;
    RegionBounds cur_region_bounds;
    for (size_t block_idx = 0; block_idx < candidate_blocks.size(); block_idx++)
    {
        calc_region_from_block(candidate_blocks[block_idx], resi_idx_x, resi_idx_y, cur_region_bounds);
        bool valid_region = (cur_region_bounds.region_x_start != -1)
            && (cur_region_bounds.region_x_end - cur_region_bounds.region_x_start >= 4)
            && (cur_region_bounds.region_y_end - cur_region_bounds.region_y_start >= 4);
        if (valid_region)
        {
            usbcat_regions_list.push_back(cur_region_bounds);
        }
    }
}

inline void fill_region_meta(
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int region_num,
    std::vector<RegionMeta>& region_meta);

inline void build_domains_bounds(
    const std::vector<AFPBlock>& candidate_blocks_list,
    const CoordArray& xa,
    const CoordArray& ya,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int xlen,
    const int ylen,
    const double cur_local_badRmsd,
    RegionBoundsAllChain& bounds)
{
    std::vector<AFPBlock> candidate_blocks = candidate_blocks_list;

    split_candidate_blocks(candidate_blocks, xa, ya, usb_cat_para, cur_local_badRmsd);

    remove_single_elem_block(candidate_blocks, usb_cat_para, xlen, ylen);

    merge_adjacent_blocks(candidate_blocks, xa, ya, cur_local_badRmsd);

    std::vector<RegionBounds> usbcat_regions_list;
    build_usbcat_regions(candidate_blocks, usbcat_regions_list);
    if (usbcat_regions_list.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    std::vector<int> chain1_bounds, chain2_bounds;
    chain1_bounds.push_back(0);
    chain2_bounds.push_back(0);
    for (size_t region_idx = 0; region_idx < usbcat_regions_list.size() - 1; region_idx++)
    {
        chain1_bounds.push_back((usbcat_regions_list[region_idx].region_x_end + usbcat_regions_list[region_idx + 1].region_x_start) / 2);
        chain2_bounds.push_back((usbcat_regions_list[region_idx].region_y_end + usbcat_regions_list[region_idx + 1].region_y_start) / 2);
    }
    chain1_bounds.push_back(xlen);
    chain2_bounds.push_back(ylen);

    bounds.chain1_bounds = chain1_bounds;
    bounds.chain2_bounds = chain2_bounds;

    int region_num = (int)chain1_bounds.size() - 1;
    fill_region_meta(chain1_bounds, chain2_bounds, tables, usb_cat_para, region_num, bounds.region_meta);
}

inline void fill_region_bound_pools(
    const RegionBoundsAllChain& bounds_default,
    const RegionBoundsAllChain& bounds_strict,
    RegionBoundPool& region_bound_pool)
{
    region_bound_pool.region_bounds.push_back(bounds_default);
    if (bounds_strict.chain1_bounds != bounds_default.chain1_bounds || bounds_strict.chain2_bounds != bounds_default.chain2_bounds)
    {
        region_bound_pool.region_bounds.push_back(bounds_strict);
    }
}

inline void generate_bounds(const CoordArray& xa, const CoordArray& ya,
                            const LocalDistTables& tables, const USBCATParams& usb_cat_para,
                            const int xlen, const int ylen, const double min_resid_num,
                            const double cur_rmsdCut, const double cur_badRmsd,
                            const double cur_local_badRmsd, const int sparse_val,
                            RegionBoundsAllChain& bounds)
{
    std::vector<USBCAT_AFP> initial_afps =
        extract_initial_afps(xa, ya, tables, usb_cat_para, xlen, ylen, min_resid_num,
                             cur_rmsdCut, cur_badRmsd, sparse_val);

    std::vector<USBCAT_AFP> merged_afps =
        merge_afps_group(initial_afps, xa, ya, usb_cat_para, xlen, ylen,
                            cur_rmsdCut, cur_badRmsd);
    if (merged_afps.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    std::vector<int> afp_aft_index;
    std::vector<int> afp_bef_index;
    build_prev_lookup_tables(merged_afps, afp_aft_index, afp_bef_index, xlen, ylen);

    std::vector<int> best_afp_path = find_best_path(merged_afps, afp_aft_index, afp_bef_index,
                                            tables, usb_cat_para, xlen, ylen);

    std::vector<AFPBlock> candidate_blocks =
        build_candidate_blocks_list(best_afp_path, merged_afps, tables, usb_cat_para);
    if (candidate_blocks.empty())
    {
        bounds = RegionBoundsAllChain{};
        return;
    }

    build_domains_bounds(candidate_blocks, xa, ya, tables, usb_cat_para, xlen, ylen,
                         cur_local_badRmsd, bounds);
}

struct RegionPdbData
{
    CoordArray xa;
    CoordArray ya;
    std::string seqx, secx, seqy, secy;
};

struct GlobalAlignResult
{
    std::string seqM, seqxA, seqyA;
    std::vector<std::vector<double> > tu_vec;
    std::vector<int> res_tu;
};

inline void get_cur_region_pdb_data(
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& secx,
    const std::string& seqy,
    const std::string& secy,
    const int region_x_start,
    const int region_x_len,
    const int region_y_start,
    const int region_y_len,
    RegionPdbData& cur_reg_data)
{
    cur_reg_data.xa.resize(region_x_len);
    cur_reg_data.ya.resize(region_y_len);
    cur_reg_data.seqx.resize(region_x_len);
    cur_reg_data.secx.resize(region_x_len);
    cur_reg_data.seqy.resize(region_y_len);
    cur_reg_data.secy.resize(region_y_len);
    for (int i = 0; i < region_x_len; i++)
    {
        cur_reg_data.xa[i][0] = xa[region_x_start + i][0];
        cur_reg_data.xa[i][1] = xa[region_x_start + i][1];
        cur_reg_data.xa[i][2] = xa[region_x_start + i][2];
        cur_reg_data.seqx[i] = seqx[region_x_start + i];
        cur_reg_data.secx[i] = secx[region_x_start + i];
    }
    for (int i = 0; i < region_y_len; i++)
    {
        cur_reg_data.ya[i][0] = ya[region_y_start + i][0];
        cur_reg_data.ya[i][1] = ya[region_y_start + i][1];
        cur_reg_data.ya[i][2] = ya[region_y_start + i][2];
        cur_reg_data.seqy[i] = seqy[region_y_start + i];
        cur_reg_data.secy[i] = secy[region_y_start + i];
    }
}

inline void calc_region_rms_sq(
    const LocalDistTables& tables,
    const int region_x_start,
    const int region_y_start,
    const int region_len,
    double& rms_sq,
    int& count)
{
    for (int i = 0; i < region_len; i++)
    {
        int j_end = std::min((int)region_len, i + tables.max_residue_gap);

        for (int j = i + 2; j < j_end; j++)
        {
            double d1 = dtable_dist(tables, true, region_x_start + i, region_x_start + j);
            double d2 = dtable_dist(tables, false, region_y_start + i, region_y_start + j);
            rms_sq += (d1 - d2) * (d1 - d2);
            count++;
        }
    }
}

inline void fill_region_meta(
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const LocalDistTables& tables,
    const USBCATParams& usb_cat_para,
    const int region_num,
    std::vector<RegionMeta>& region_meta)
{
    for (int region_idx = 0; region_idx < region_num; region_idx++)
    {
        int region_x_start = chain1_bounds[region_idx], region_x_end = chain1_bounds[region_idx + 1];
        int region_y_start = chain2_bounds[region_idx], region_y_end = chain2_bounds[region_idx + 1];
        int region_x_len = region_x_end - region_x_start;
        int region_y_len = region_y_end - region_y_start;
        int region_len = std::min(region_x_len, region_y_len);

        double region_drmsd = 0.0;
        if (region_len < 2 * usb_cat_para.fragLen)
        {
            region_meta.push_back({region_idx, region_drmsd, region_x_len, region_y_len});
            continue;
        }

        double rms_sq = 0.0;
        int count = 0;
        calc_region_rms_sq(tables, region_x_start, region_y_start, region_len, rms_sq, count);
        if (count > 0)
            region_drmsd = std::sqrt(rms_sq / count);
        region_meta.push_back({region_idx, region_drmsd, region_x_len, region_y_len});
    }
}


inline int calc_local_hinge_opt(
    const int region_x_len,
    const int region_y_len,
    const USBCATParams& usb_cat_para,
    const bool hinge_set,
    const int remaining_hinges)
{
    bool long_enough = (std::min(region_x_len, region_y_len) >= 2 * usb_cat_para.fragLen);
    int local_hinge_opt = 0;
    if (hinge_set)
    {
        if (remaining_hinges > 0 && long_enough) local_hinge_opt = remaining_hinges;
    }
    else
    {
        if (long_enough) local_hinge_opt = 2;
    }
    return local_hinge_opt;
}

inline void align_cur_region(
    const RegionBoundsAllChain& bounds,
    const int region_indx,
    const bool hinge_set,
    const int remaining_hinges,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const AlignCommonInput& common_inputs,
    const int mol_type,
    bool& valid,
    FlexAlignResult& cur_region_align_res,
    std::ostream& os = std::cout)
{
    const RegionMeta& region_meta = bounds.region_meta[region_indx];
    int orig_region_idx = region_meta.original_region_idx;

    int region_x_start = bounds.chain1_bounds[orig_region_idx];
    int region_y_start = bounds.chain2_bounds[orig_region_idx];
    int region_x_len = region_meta.region_x_len;
    int region_y_len = region_meta.region_y_len;

    valid = true;
    if (region_x_len < 3 || region_y_len < 3) { valid = false; return; }

    int local_hinge_opt = calc_local_hinge_opt(region_x_len, region_y_len, usb_cat_para, hinge_set, remaining_hinges);

    RegionPdbData cur_reg_data;
    get_cur_region_pdb_data(xa, ya, seqx, secx, seqy, secy, region_x_start, region_x_len, region_y_start, region_y_len, cur_reg_data);

    align_with_flexalign_main(
        cur_reg_data.xa, cur_reg_data.ya, cur_reg_data.seqx, cur_reg_data.seqy, cur_reg_data.secx, cur_reg_data.secy,
        region_x_len, region_y_len, common_inputs, mol_type, local_hinge_opt, cur_region_align_res, os);
}

inline void run_region_align(
    std::vector<RegionAlignResult>& region_align_res,
    int& remaining_hinges,
    const RegionBoundsAllChain& cur_bound_pool,
    const int region_num,
    const bool hinge_set,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const AlignCommonInput& common_inputs,
    const int mol_type,
    std::ostream& os = std::cout)
{
    for (int region_indx = 0; region_indx < region_num; region_indx++)
    {
        FlexAlignResult cur_align_res;
        bool region_valid = false;
        align_cur_region(
            cur_bound_pool, region_indx, hinge_set, remaining_hinges,
            usb_cat_para, xa, ya, seqx, seqy, secx, secy,
            common_inputs, mol_type, region_valid, cur_align_res, os);

        int orig_region_idx = cur_bound_pool.region_meta[region_indx].original_region_idx;
        if (!region_valid) continue;
        region_align_res[orig_region_idx].align_re = cur_align_res;
        region_align_res[orig_region_idx].valid = true;

        if (hinge_set && !cur_align_res.tu_vec.empty())
        {
            int consumed_hinges = (int)cur_align_res.tu_vec.size() - 1;
            if (consumed_hinges > 0)
            {
                remaining_hinges -= consumed_hinges;
                if (remaining_hinges < 0) remaining_hinges = 0;
            }
        }
    }
}

inline void build_global_align_result(
    std::vector<RegionAlignResult>& region_align_res,
    const std::vector<int>& chain1_bounds,
    const std::vector<int>& chain2_bounds,
    const int region_num,
    const std::string& seqx,
    const std::string& seqy,
    GlobalAlignResult& global_align_res)
{
    for (int region_idx = 0; region_idx < region_num; region_idx++)
    {
        int L1_sub = chain1_bounds[region_idx + 1] - chain1_bounds[region_idx];
        int L2_sub = chain2_bounds[region_idx + 1] - chain2_bounds[region_idx];

        if (!region_align_res[region_idx].valid)
        {
            // Fill gaps if region was invalid or bypassed
            for (int i = 0; i < L1_sub; i++)
            {
                global_align_res.seqxA += seqx[chain1_bounds[region_idx] + i];
                global_align_res.seqyA += '-';
                global_align_res.seqM += ' ';
            }
            for (int i = 0; i < L2_sub; i++)
            {
                global_align_res.seqxA += '-';
                global_align_res.seqyA += seqy[chain2_bounds[region_idx] + i];
                global_align_res.seqM += ' ';
            }
            continue;
        }

        FlexAlignResult& flexalign_res = region_align_res[region_idx].align_re;
        if (flexalign_res.tu_vec.empty())
        {
            std::vector<double> tu_tmp(12);
            t_u2tu(flexalign_res.t0, flexalign_res.u0, tu_tmp);
            flexalign_res.tu_vec.push_back(tu_tmp);
        }

        int base_tu_idx = (int)global_align_res.tu_vec.size();
        for (size_t m = 0; m < flexalign_res.tu_vec.size(); m++)
            global_align_res.tu_vec.push_back(flexalign_res.tu_vec[m]);

        int rx = chain1_bounds[region_idx];
        int current_global_idx = base_tu_idx;

        for (size_t i = 0; i < flexalign_res.seqxA.length(); i++)
        {
            char c = flexalign_res.seqM[i];
            if (c != ' ' && c != '.' && c != ':')
            {
                int local_hinge_idx = -1;
                if (c >= '0' && c <= '9')
                    local_hinge_idx = c - '0';
                else if (c >= 'a' && c <= 'z')
                    local_hinge_idx = c - 'a' + 10;
                else if (c >= 'A' && c <= 'Z')
                    local_hinge_idx = c - 'A' + 36;

                if (local_hinge_idx >= 0 && local_hinge_idx < (int)flexalign_res.tu_vec.size())
                    current_global_idx = base_tu_idx + local_hinge_idx;
            }

            if (flexalign_res.seqxA[i] != '-')
            {
                global_align_res.res_tu[rx] = current_global_idx;
                rx++;
            }

            if (flexalign_res.seqxA[i] != '-' && flexalign_res.seqyA[i] != '-')
            {
                if (c != ' ' && c != '.' && c != ':')
                {
                    char global_c;
                    if (current_global_idx < 10)
                        global_c = '0' + current_global_idx;
                    else if (current_global_idx < 36)
                        global_c = 'a' + (current_global_idx - 10);
                    else if (current_global_idx < 62)
                        global_c = 'A' + (current_global_idx - 36);
                    else
                        global_c = '*';
                    flexalign_res.seqM[i] = global_c;
                }
                else
                    flexalign_res.seqM[i] = c;
            }
            else
                flexalign_res.seqM[i] = ' ';
        }

        global_align_res.seqM += flexalign_res.seqM;
        global_align_res.seqxA += flexalign_res.seqxA;
        global_align_res.seqyA += flexalign_res.seqyA;
    }
}



inline void calc_cur_d0(
    const int xlen, const int ylen,
    const UserOptions& opts,
    const int mol_type,
    double& d0A, double& d0B, double& d0a, double& d0u)
{
    double dummy_D0_MIN, dummy_Lnorm, dummy_d0_search;

    parameter_set4final(ylen, dummy_D0_MIN, dummy_Lnorm, d0A, dummy_d0_search, mol_type);
    parameter_set4final(xlen, dummy_D0_MIN, dummy_Lnorm, d0B, dummy_d0_search, mol_type);
    parameter_set4final((xlen + ylen) * 0.5, dummy_D0_MIN, dummy_Lnorm, d0a, dummy_d0_search, mol_type);
    d0u = 0.0;
    if (opts.u_opt)
        parameter_set4final(opts.Lnorm_ass, dummy_D0_MIN, dummy_Lnorm, d0u, dummy_d0_search, mol_type);
}

inline void reset_flexalign_result(FlexAlignResult& flexalign_res, const size_t aln_len)
{
    flexalign_res.do_vec.clear();
    flexalign_res.do_vec.reserve(aln_len);
    flexalign_res.TM1 = 0.0;
    flexalign_res.TM2 = 0.0;
    flexalign_res.TM3 = 0.0;
    flexalign_res.TM4 = 0.0;
    flexalign_res.TM5 = 0.0;
    flexalign_res.rmsd0 = 0.0;
    flexalign_res.Liden = 0.0;
    flexalign_res.n_ali = 0;
    flexalign_res.n_ali8 = 0;
}

inline void recompute_global_metrics(
    const GlobalAlignResult& sd,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const int xlen,
    const int ylen,
    const AlignCommonInput& common_inputs,
    const int mol_type,
    const double d0_out,
    FlexAlignResult& flexalign_res)
{
    const UserOptions& opts = common_inputs.user_options;

    double cur_d0A, cur_d0B, cur_d0a, cur_d0u;
    calc_cur_d0(xlen, ylen, opts, mol_type, cur_d0A, cur_d0B, cur_d0a, cur_d0u);
    reset_flexalign_result(flexalign_res, sd.seqxA.length());

    int i_res = 0, j_res = 0;
    int cached_matrix_idx = -1;
    Vec3 t_k;
    RotMat u_k;
    for (size_t r = 0; r < sd.seqxA.length(); r++)
    {
        bool x_valid = (sd.seqxA[r] != '-');
        bool y_valid = (sd.seqyA[r] != '-');

        double d = -1.0;
        if (x_valid && y_valid)
        {
            int matrix_idx = sd.res_tu[i_res];
            if (matrix_idx >= 0 && matrix_idx < (int)sd.tu_vec.size())
            {
                if (matrix_idx != cached_matrix_idx)
                {
                    tu2t_u(sd.tu_vec[matrix_idx], t_k, u_k);
                    cached_matrix_idx = matrix_idx;
                }

                Vec3 x_rot;
                transform(t_k, u_k, xa[i_res], x_rot);
                double dist2 = dist(x_rot, ya[j_res]);
                d = std::sqrt(dist2);

                flexalign_res.TM2 += 1.0 / (1.0 + dist2 / (cur_d0B * cur_d0B));
                flexalign_res.TM1 += 1.0 / (1.0 + dist2 / (cur_d0A * cur_d0A));
                if (opts.a_opt)
                    flexalign_res.TM3 += 1.0 / (1.0 + dist2 / (cur_d0a * cur_d0a));
                if (opts.u_opt)
                    flexalign_res.TM4 += 1.0 / (1.0 + dist2 / (cur_d0u * cur_d0u));
                if (opts.d_opt)
                    flexalign_res.TM5 += 1.0 / (1.0 + dist2 / (opts.d0_scale * opts.d0_scale));

                flexalign_res.n_ali++;

                if (d <= d0_out)
                {
                    flexalign_res.rmsd0 += dist2;
                    flexalign_res.n_ali8++;
                    if (seqx[i_res] == seqy[j_res])
                        flexalign_res.Liden += 1.0;
                }
            }
        }
        flexalign_res.do_vec.push_back(d);

        if (x_valid)
            i_res++;
        if (y_valid)
            j_res++;
    }

    flexalign_res.TM2 /= xlen;
    flexalign_res.TM1 /= ylen;
    if (opts.a_opt)
        flexalign_res.TM3 /= (xlen + ylen) * 0.5;
    if (opts.u_opt)
        flexalign_res.TM4 /= opts.Lnorm_ass;
    if (opts.d_opt)
        flexalign_res.TM5 /= ylen;
    if (flexalign_res.n_ali8 > 0)
        flexalign_res.rmsd0 = std::sqrt(flexalign_res.rmsd0 / flexalign_res.n_ali8);
    else
        flexalign_res.rmsd0 = 0.0;

    flexalign_res.seqM = sd.seqM;
    flexalign_res.seqxA = sd.seqxA;
    flexalign_res.seqyA = sd.seqyA;
    flexalign_res.tu_vec = sd.tu_vec;
    flexalign_res.d0A = cur_d0A;
    flexalign_res.d0B = cur_d0B;
    flexalign_res.d0a = cur_d0a;
    flexalign_res.d0u = cur_d0u;
    if (!flexalign_res.tu_vec.empty())
        tu2t_u(flexalign_res.tu_vec[0], flexalign_res.t0, flexalign_res.u0);

    flexalign_res.TM_ali = flexalign_res.TM1;
    flexalign_res.rmsd_ali = flexalign_res.rmsd0;
    flexalign_res.L_ali = flexalign_res.n_ali;
}

inline void update_global_best_align(
    RegionBoundPool& region_bound_pool,
    const bool hinge_set,
    const int hinge_opt,
    const USBCATParams& usb_cat_para,
    const CoordArray& xa,
    const CoordArray& ya,
    const std::string& seqx,
    const std::string& seqy,
    const std::string& secx,
    const std::string& secy,
    const AlignCommonInput& common_inputs,
    const int mol_type,
    const int xlen,
    const int ylen,
    const double d0_out,
    FlexAlignResult& best_res,
    double& best_global_max_TM,
    std::ostream& os = std::cout)
{
    for (size_t pool_idx = 0; pool_idx < region_bound_pool.region_bounds.size(); pool_idx++)
    {
        RegionBoundsAllChain& cur_bound_pool = region_bound_pool.region_bounds[pool_idx];
        const std::vector<int>& chain1_bounds = cur_bound_pool.chain1_bounds;
        const std::vector<int>& chain2_bounds = cur_bound_pool.chain2_bounds;

        if (chain1_bounds.size() <= 2) continue;
        int region_num = (int)chain1_bounds.size() - 1;

        if (hinge_set)
            std::sort(cur_bound_pool.region_meta.begin(), cur_bound_pool.region_meta.end(), region_meta_drmsd_desc);

        int remaining_hinges = hinge_set ? std::max(0, hinge_opt + 1 - region_num) : 0;

        std::vector<RegionAlignResult> region_align_res(region_num);
        run_region_align(
            region_align_res, remaining_hinges,
            cur_bound_pool, region_num, hinge_set,
            usb_cat_para, xa, ya, seqx, seqy, secx, secy,
            common_inputs, mol_type, os);

        GlobalAlignResult global_align_res;
        global_align_res.res_tu.assign(xlen, -1);
        
        build_global_align_result(region_align_res, chain1_bounds, chain2_bounds, region_num, seqx, seqy, global_align_res);
        
        FlexAlignResult cur_res;
        recompute_global_metrics( global_align_res, xa, ya, seqx, seqy, xlen, ylen, common_inputs, mol_type, d0_out, cur_res);

        double cur_global_max_TM = (cur_res.TM1 > cur_res.TM2) ? cur_res.TM1 : cur_res.TM2;
        if (cur_global_max_TM > best_global_max_TM)
        {
            best_global_max_TM = cur_global_max_TM;
            best_res = cur_res;
        }
    }
}

void flexalign_usbcat_main(
    CoordArray& xa,
    CoordArray& ya,
    const std::string &seqx,
    const std::string &seqy,
    const std::string &secx,
    const std::string &secy,
    int xlen,
    int ylen,
    const AlignCommonInput& common_inputs,
    const FlexAlignResult& flexalign_main_res,
    int mol_type, double best_global_max_TM,
    int sparse_val,
    const FlexalignParams& flex_params,
    FlexAlignResult& flexalign_res,
    std::ostream& os);

void fill_flexalign_params(
    FlexalignParams& flex_params,
    const ControlOptions& control);
int Flexalign(AlignCommonInput& common_inputs, const FlexalignParams& flex_params, FlexAlignResult& flex_result);
