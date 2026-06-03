with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Remove the NWDP bridge (static version)
old = "static void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,\n    int len1, int len2, CharMatrix& mask,\n    double t[3], double u[3][3], double d02, double gap_open, int j2i[])\n{\n    std::vector<double*> _xv(x.size()), _yv(y.size());\n    for (size_t i=0; i<x.size(); i++) { _xv[i]=(double*)x[i].data(); _yv[i]=(double*)y[i].data(); }\n    NWDP_TM_dimer(path, val, _xv.data(), _yv.data(), len1, len2, mask, t, u, d02, gap_open, j2i);\n}\n"

if old in content:
    content = content.replace(old, "")
    print("1: Removed NWDP bridge")
else:
    print("1 FAILED")
    idx = content.find("static void NWDP_TM_dimer")
    if idx != -1:
        end = content.find("}\n", idx)
        print(repr(content[idx:end+2]))
        content = content[:idx] + content[end+2:]
        print("1: Removed NWDP bridge (manual)")

# Revert NWDP_TM_dimer call to temp views
old2 = "                NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,\n                    t, u, d02, gap_open, invmap);"
new2 = "                {\n                    vector<double*> _xv2(x.size()), _yv2(y.size());\n                    for (size_t _i=0; _i<x.size(); _i++) _xv2[_i]=(double*)x[_i].data();\n                    for (size_t _i=0; _i<y.size(); _i++) _yv2[_i]=(double*)y[_i].data();\n                    NWDP_TM_dimer(path, val, _xv2.data(), _yv2.data(), xlen, ylen, mask,\n                        t, u, d02, gap_open, invmap);\n                }"

if old2 in content:
    content = content.replace(old2, new2)
    print("2: Reverted NWDP_TM_dimer to temp views")
else:
    print("2 FAILED - checking...")
    idx = content.find("NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,")
    if idx != -1:
        print(f"Found at {idx}: {repr(content[idx:idx+120])}")

# Also restore get_score_fast to use ALL bridge (remove temp views)
# First check current state
if "_xv(x.size()), _yv(y.size())" in content:
    # Temp views are present for get_score_fast
    
    # Remove get_score_fast temp views
    old3 = "                {\n                    vector<double*> _xv(x.size()), _yv(y.size());\n                    for (size_t _i=0; _i<x.size(); _i++) _xv[_i]=(double*)x[_i].data();\n                    for (size_t _i=0; _i<y.size(); _i++) _yv[_i]=(double*)y[_i].data();\n                    GL = get_score_fast(r1, r2, xtm, ytm, _xv.data(), _yv.data(), xlen, ylen,\n                        invmap, d0, d0_search, t, u);\n                }"
    new3 = "                GL = get_score_fast(r1, r2, xtm, ytm, x, y, xlen, ylen,\n                    invmap, d0, d0_search, t, u);"
    
    if old3 in content:
        content = content.replace(old3, new3)
        print("3: get_score_fast uses ALL bridge")
    else:
        print("3 FAILED")
else:
    print("3: get_score_fast already uses ALL bridge")

# Also add back get_score_fast ALL-CoordArray bridge in TMalign.h
with open("TMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content2 = f.read()

old4 = "\ndouble get_score_fast( double **r1, double **r2, CoordArray& xtm, CoordArray& ytm,\n    double **x, double **y,"
bridge4 = """
double get_score_fast( double **r1, double **r2, CoordArray& xtm, CoordArray& ytm,
    CoordArray& x, CoordArray& y, int xlen, int ylen, int invmap[],
    double d0, double d0_search, double t[3], double u[3][3])
{
    std::vector<double*> _xtmv(xtm.size()), _ytmv(ytm.size());
    std::vector<double*> _xv(x.size()), _yv(y.size());
    for (size_t i=0; i<xtm.size(); i++) _xtmv[i]=(double*)xtm[i].data();
    for (size_t i=0; i<ytm.size(); i++) _ytmv[i]=(double*)ytm[i].data();
    for (size_t i=0; i<x.size(); i++) _xv[i]=(double*)x[i].data();
    for (size_t i=0; i<y.size(); i++) _yv[i]=(double*)y[i].data();
    return get_score_fast(r1, r2, _xtmv.data(), _ytmv.data(), _xv.data(), _yv.data(),
        xlen, ylen, invmap, d0, d0_search, t, u);
}

double get_score_fast( double **r1, double **r2, CoordArray& xtm, CoordArray& ytm,
    double **x, double **y,"

if old4 in content2:
    content2 = content2.replace(old4, bridge4)
    print("4: Added get_score_fast ALL-CoordArray bridge")
else:
    print("4 FAILED")

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)
with open("TMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content2)

import subprocess
r = subprocess.run(["g++", "-O3", "-ffast-math", "-lm", "-static", "-o", "USalign_test.exe", "USalign.cpp"],
                   capture_output=True, text=True, timeout=120)
if r.returncode == 0:
    print("Compilation OK")
else:
    print("FAILED:")
    print(r.stderr[-500:])
