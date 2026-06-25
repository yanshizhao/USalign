#pragma once

#include <vector>
#include <string>

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

// ---------------------------------------------------------------------------
// NWDP_TM (score-matrix version)
// ---------------------------------------------------------------------------
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
            if (path[i - 1][j] == 1) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == 1) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = 1;  // diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = 3;  // left
                val[i][j] = v;
            } else {
                path[i][j] = 2;  // up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack using stored path codes
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == 1) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == 2) { i--; }
        else { j--; }  // path[i][j] == 3 (left)
    }
}

// ---------------------------------------------------------------------------
// NWDP_TM (coordinate version with rotation)  –  hot path, used in DP_iter
// ---------------------------------------------------------------------------
/* Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
 * Output: j2i[1:len2] ∈ {1:len1} ∪ {-1}
 * path[, ] ∈ {1=diagonal, 2=up,   3=left} */
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
            if (path[i - 1][j] == 1) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == 1) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = 1;  // diagonal
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = 3;  // left
                val[i][j] = v;
            } else {
                path[i][j] = 2;  // up
                val[i][j] = h;
            }
        }
    }

    // Fast backtrack – no need to re-read val[] just to decide direction
    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == 1) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == 2) { i--; }
        else { j--; }
    }
}

// ---------------------------------------------------------------------------
// NWDP_TM (secondary-structure version)
// ---------------------------------------------------------------------------
/* Input: secondary structure strings secx, secy, and gap_open
 * Output: j2i[1:len2] ∈ {1:len1} ∪ {-1} */
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
            if (path[i - 1][j] == 1) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == 1) v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = 1;
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = 3;
                val[i][j] = v;
            } else {
                path[i][j] = 2;
                val[i][j] = h;
            }
        }
    }

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == 1) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == 2) { i--; }
        else { j--; }
    }
}

// ---------------------------------------------------------------------------
// NWDP_SE  –  superposition-engine DP with optional hinge (flexible alignment)
// ---------------------------------------------------------------------------
inline void NWDP_SE(CharMatrix& path, DoubleMatrix& val,
    CoordArray& x, CoordArray& y,
    int len1, int len2, double d02, double gap_open,
    std::vector<int>& j2i, const int hinge)
{
    int    i, j;
    double h, v, d;
    double dij;

    if (hinge == 0) {
        // ---- ordinary DP (same as coordinate version, but no rotation) ----
        for (i = 0; i <= len1; i++) { val[i][0] = 0; path[i][0] = 0; }
        for (j = 0; j <= len2; j++) { val[0][j] = 0; path[0][j] = 0; j2i[j] = -1; }

        for (i = 1; i <= len1; i++) {
            for (j = 1; j <= len2; j++) {
                dij = dist(x[i - 1], y[j - 1]);
                d   = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
                h   = val[i - 1][j];
                if (path[i - 1][j] == 1) h += gap_open;
                v = val[i][j - 1];
                if (path[i][j - 1] == 1) v += gap_open;

                if (d >= h && d >= v) {
                    path[i][j] = 1;
                    val[i][j] = d;
                } else if (v >= h) {
                    path[i][j] = 3;
                    val[i][j] = v;
                } else {
                    path[i][j] = 2;
                    val[i][j] = h;
                }
            }
        }

        i = len1; j = len2;
        while (i > 0 && j > 0) {
            if (path[i][j] == 1) { j2i[j - 1] = i - 1;  i--;  j--; }
            else if (path[i][j] == 2) { i--; }
            else { j--; }
        }
        return;
    }

    // ---- hinge > 0 : masked DP with pre-set alignment anchors ----
    int L = (len2 > len1) ? len2 : len1;
    int int_min = L * (gap_open - 1);

    for (i = 0; i <= len1; i++)
        for (j = 0; j <= len2; j++)
            val[i][j] = 0, path[i][j] = 0;

    for (j = 0; j < len2; j++) {
        i = j2i[j];
        if (i < 0) continue;
        path[i + 1][j + 1] = 1;
        val[i + 1][j + 1] = 0;
    }

    for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
            dij = 0;
            if (path[i][j] == 0)
                dij = dist(x[i - 1], y[j - 1]);

            d = val[i - 1][j - 1] + 1.0 / (1.0 + dij / d02);
            h = val[i - 1][j];
            if (path[i - 1][j] == 1) h += gap_open;
            v = val[i][j - 1];
            if (path[i][j - 1] == 1) v += gap_open;

            if (d >= h && d >= v && val[i][j] == 0) {
                path[i][j] = 1;
                val[i][j] = d;
            } else if (v >= h) {
                path[i][j] = 3;
                val[i][j] = v;
            } else {
                path[i][j] = 2;
                val[i][j] = h;
            }
        }
    }

    for (j = 0; j <= len2; j++) j2i[j] = -1;

    i = len1; j = len2;
    while (i > 0 && j > 0) {
        if (path[i][j] == 1) { j2i[j - 1] = i - 1;  i--;  j--; }
        else if (path[i][j] == 2) { i--; }
        else { j--; }
    }
}
