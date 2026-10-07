#pragma once

#include "io.h"
#include "BLOSUM.h"

/* Partial implementation of Needleman-Wunsch (NW) dynamic programming for
 * global alignment. The three NWDP_TM functions below are not complete
 * implementation of NW algorithm because gap jumping in the standard Gotoh
 * algorithm is not considered. Since the gap opening and gap extension is
 * the same, this is not a problem. This code was exploited in TM-align
 * because it is about 1.5 times faster than a complete NW implementation.
 * Nevertheless, if gap opening != gap extension shall be implemented in
 * the future, the Gotoh algorithm must be implemented. In rare scenarios,
 * it is also possible to have asymmetric alignment (i.e.
 * TMalign A.pdb B.pdb and TMalign B.pdb A.pdb have different TM_A and TM_B
 * values) caused by the NWPD_TM implement.
 */

// 3-state path encoding for fast backtracking
#define PATH_DIAG  1   // diagonal: residue i aligns with residue j
#define PATH_UP    2   // up:       residue i aligns with a gap
#define PATH_LEFT  3   // left:     gap aligns with residue j

// Forward declarations
inline void NWDP_TM(const DoubleMatrix& score, CharMatrix& path,
    DoubleMatrix& val, int len1, int len2, double gap_open, std::vector<int>& j2i);

inline void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const CoordArray& x, const CoordArray& y,
    int len1, int len2, const Vec3& t, const RotMat& u,
    double d02, double gap_open, std::vector<int>& j2i);

inline void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const std::string& secx, const std::string& secy,
    const int len1, const int len2, const double gap_open, std::vector<int>& j2i);

inline void NWDP_SE(CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y,
    int len1, int len2, double d02, double gap_open,
    std::vector<int>& j2i, const int hinge);


inline void NWDP_TM(const DoubleMatrix& score, CharMatrix& path,
    DoubleMatrix& val, int len1, int len2, double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            d = val[i - 1][j - 1] + score[i][j];
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack using stored path codes
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }  // path[i][j] == PATH_LEFT(left)
    }
}


inline void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const CoordArray& x, const CoordArray& y,
    int len1, int len2, const Vec3& t, const RotMat& u,
    double d02, double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;
    Vec3   xx;
    double dij;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        transform(t, u, x[i - 1], xx);

        for (j = 1; j <= len2; j++) {
            dij = dist(xx, y[j - 1]);
            d   = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
            h   = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack – no need to re-read val[] just to decide direction
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}


inline void NWDP_TM(CharMatrix& path, DoubleMatrix& val,
    const std::string& secx, const std::string& secy,
    const int len1, const int len2, const double gap_open, std::vector<int>& j2i)
{
    int    i, j;
    double h, v, d;

    for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
    for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            d = val[i - 1][j - 1] + 1.0 * (secx[i - 1] == secy[j - 1]);
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}


inline void NWDP_SE(CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y,
    int len1, int len2, double d02, double gap_open,
    std::vector<int>& j2i, const int hinge)
{
    int    i, j;
    double h, v, d;
    double dij;

    if (hinge == 0) {
        for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
        for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

        for (i = 1; i <= len1; i++) {
            for (j = 1; j <= len2; j++) {
                dij = dist(x[i - 1], y[j - 1]);
                d   = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
                h   = val[i - 1][j];
                if (path[i - 1][j] == PATH_DIAG) h += gap_open;
                v = val[i][j - 1];
                if (path[i][j - 1] == PATH_DIAG) v += gap_open;

                if (d >= h && d >= v) {
                    path[i][j] = PATH_DIAG;// diagonal
                    val[i][j] = d;
                } else if (v >= h) {
                    path[i][j] = PATH_LEFT;// left
                    val[i][j] = v;
                } else {
                    path[i][j] = PATH_UP;// up
                    val[i][j] = h;
                }
            }
        }

        i = len1; j = len2;
        while (i > 0 && j > 0) {
            if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
            else if (path[i][j] == PATH_UP) { i--; }
            else { j--; }
        }
        return;
    }

    for (i = 0; i <= len1; i++)
        for (j = 0; j <= len2; j++)
            val[i][j] = 0, path[i][j] = 0;

    for (j = 0; j < len2; j++) {
        i = j2i[j];
        if (i < 0) continue;
        path[i + 1][j + 1] = PATH_DIAG;// diagonal
        val[i + 1][j + 1] = 0;
    }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            dij = 0;
            if (path[i][j] == 0)
                dij = dist(x[i - 1], y[j - 1]);

            d = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
            h = val[i - 1][j];
            if (path[i - 1][j] == PATH_DIAG) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == PATH_DIAG) v += gap_open;

            if (d >= h && d >= v && val[i][j] == 0) {
                path[i][j] = PATH_DIAG;// diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = PATH_LEFT;// left
                val[i][j] = v;
            } else {
                path[i][j] = PATH_UP;// up
                val[i][j] = h;
            }
        }
    }

    for (j = 0; j <= len2; j++) j2i[j] = -1;

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == PATH_DIAG) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == PATH_UP) { i--; }
        else { j--; }
    }
}

const int gapopen_blosum62=-11;
const int gapext_blosum62=-1;

const int gapopen_blastn=-15; //-5;
const int gapext_blastn =-4;  //-2;

struct GotohMatrix
{
    IntMatrix S;
    IntMatrix JumpH;
    IntMatrix JumpV;
    IntMatrix P;
    IntMatrix H;
    IntMatrix V;

    void init(const int xlen, const int ylen)
    {
        S.assign(xlen+1, std::vector<int>(ylen+1));
        JumpH.assign(xlen+1, std::vector<int>(ylen+1));
        JumpV.assign(xlen+1, std::vector<int>(ylen+1));
        P.assign(xlen+1, std::vector<int>(ylen+1));
        H.assign(xlen+1, std::vector<int>(ylen+1));
        V.assign(xlen+1, std::vector<int>(ylen+1));
    }
};

// initialize matrix in gotoh algorithm
void init_gotoh_mat(GotohMatrix& m, const int xlen, const int ylen, const int gapopen,
    const int gapext, const int glocal=0, const int alt_init=1);

/* locate the cell with highest alignment score. reset path after
 * the cell to zero */
void find_highest_align_score( GotohMatrix& m,
    int &aln_score, const int xlen,const int ylen);

/* calculate dynamic programming matrix using gotoh algorithm
 * S     - cumulative scorefor each cell
 * P     - std::string representation for path
 *         0 :   uninitialized, for gaps at N- & C- termini when glocal>0
 *         1 : \ match-mismatch
 *         2 : | vertical gap (insertion)
 *         4 : - horizontal gap (deletion)
 * JumpH - horizontal long gap number.
 * JumpV - vertical long gap number.
 * all matrices are in the size of [len(seqx)+1]*[len(seqy)+1]
 *
 * glocal - global or local alignment
 *         0 : global alignment (Needleman-Wunsch dynamic programming)
 *         1 : glocal-query alignment
 *         2 : glocal-both alignment
 *         3 : local alignment (Smith-Waterman dynamic programming)
 *
 * alt_init - whether to adopt alternative matrix initialization
 *         1 : use wei zheng's matrix initialization
 *         0 : use yang zhang's matrix initialization, does NOT work
 *             for glocal alignment
 */
int calculate_score_gotoh(const int xlen,const int ylen, GotohMatrix& m, const int gapopen,const int gapext,
    const int glocal=0, const int alt_init=1);

inline void trace_back_gotoh(const std::string& seqx, const std::string& seqy,
    const GotohMatrix& m, std::string& seqxA, std::string& seqyA,
    const int xlen, const int ylen, std::vector<int>& invmap, const int invmap_only=1)
{
    const IntMatrix& JumpH = m.JumpH;
    const IntMatrix& JumpV = m.JumpV;
    const IntMatrix& P = m.P;
    int i;
    int j;
    int gaplen;
    std::string buf;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf.resize(std::max(xlen,ylen)+1);

    i=xlen;
    j=ylen;
    while(i+j)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqy[j],gaplen);
            seqyA=buf+seqyA;

            buf.assign(gaplen,'-');
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqx[i],gaplen);
            seqxA=buf+seqxA;

            buf.assign(gaplen,'-');
            seqyA=buf+seqyA;
        }
        else
        {
            if (i==0 && j!=0)
            {
                buf.assign(seqy.data(),j);
                seqyA=buf+seqyA;
                buf.assign(j,'-');
                seqxA=buf+seqxA;
                break;
            }
            if (i!=0 && j==0)
            {
                buf.assign(seqx.data(),i);
                seqxA=buf+seqxA;
                buf.assign(i,'-');
                seqyA=buf+seqyA;
                break;
            }
            i--;
            j--;
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
}


inline void trace_back_sw(const std::string& seqx, const std::string& seqy,
    const GotohMatrix& m, std::string& seqxA, std::string& seqyA,
    const int xlen, const int ylen, std::vector<int>& invmap, const int invmap_only=1)
{
    const IntMatrix& JumpH = m.JumpH;
    const IntMatrix& JumpV = m.JumpV;
    const IntMatrix& P = m.P;
    int i;
    int j;
    int gaplen;
    bool found_start_cell=false;
    std::string buf;

    if (invmap_only) for (j = 0; j < ylen; j++) invmap[j] = -1;
    if (invmap_only!=1) buf.resize(std::max(xlen,ylen)+1);

    i=xlen;
    j=ylen;
    for (i=xlen;i>=0;i--)
    {
        for (j=ylen;j>=0;j--)
        {
            if (P[i][j]!=0)
            {
                found_start_cell=true;
                break;
            }
        }
        if (found_start_cell) break;
    }

    // copy C terminal sequence
    if (invmap_only!=1)
    {
        buf.assign(ylen-j,'-');
        seqxA=buf;
        buf.assign(&seqx[i],xlen-i);
        seqxA+=buf;

        buf.assign(&seqy[j],ylen-j);
        seqyA+=buf;
        buf.assign(xlen-i,'-');
        seqyA+=buf;
    }

    if (i<0||j<0) return;

    // traceback aligned sequences
    while(P[i][j]!=0)
    {
        gaplen=0;
        if (P[i][j]>=4)
        {
            gaplen=JumpH[i][j];
            j-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqy[j],gaplen);
            seqyA=buf+seqyA;

            buf.assign(gaplen,'-');
            seqxA=buf+seqxA;
        }
        else if (P[i][j] % 4 >= 2)
        {
            gaplen=JumpV[i][j];
            i-=gaplen;
            if (invmap_only==1) continue;
            buf.assign(&seqx[i],gaplen);
            seqxA=buf+seqxA;

            buf.assign(gaplen,'-');
            seqyA=buf+seqyA;
        }
        else
        {
            i--;
            j--;
            if (invmap_only) invmap[j]=i;
            if (invmap_only!=1)
            {
                seqxA=seqx[i]+seqxA;
                seqyA=seqy[j]+seqyA;
            }
        }
    }
    // copy N terminal sequence
    if (invmap_only!=1)
    {
        buf.assign(j,'-');
        buf.append(seqx.data(),i);
        seqxA=buf+seqxA;

        buf.assign(seqy.data(),j);
        buf.append(i,'-');
        seqyA=buf+seqyA;
    }
}


/* entry function for NWalign
 * invmap_only - whether to return seqxA and seqyA or to return invmap
 *               0: only return seqxA and seqyA
 *               1: only return invmap
 *               2: return seqxA, seqyA and invmap */
int NWalign_main(const std::string &seqx, const std::string &seqy, const int xlen,
    const int ylen, std::string & seqxA, std::string & seqyA, const int mol_type,
    std::vector<int>& invmap, const int invmap_only=0, const int glocal=0);

void get_seqID(const std::vector<int>& invmap, const std::string& seqx, const std::string& seqy,
    const int ylen, double &Liden,int &L_ali);

void get_seqID(const std::string& seqxA, const std::string& seqyA,
    std::string &seqM,double &Liden,int &L_ali);

void output_NWalign_results(
    const std::string xname, const std::string yname,
    const std::string& chainID1, const std::string& chainID2,
    const int xlen, const int ylen, const std::string& seqM,
    const std::string& seqxA, const std::string& seqyA, const double Liden,
    const int L_ali, const int aln_score, const int outfmt_opt);

/* extract pairwise sequence alignment from residue index vectors,
 * assuming that "sequence" contains two empty strings.
 * return length of alignment, including gap. */
int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const int byresi_opt);

/* extract pairwise sequence alignment from residue index vectors,
 * return length of alignment, including gap. */
int extract_aln_from_resi(std::vector<std::string> &sequence, const std::string& seqx, const std::string& seqy,
    const std::vector<std::string> resi_vec1, const std::vector<std::string> resi_vec2,
    const std::vector<int> xlen_vec, const std::vector<int> ylen_vec,
    const int chain_i, const int chain_j, const int byresi_opt);