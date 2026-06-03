with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Change bridge to delegate to CharMatrix& version (line 2157) instead of char** version
old = "inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,\n    int len1, int len2, CharMatrix& mask,\n    double t[3], double u[3][3], double d02, double gap_open, int j2i[])\n{\n    std::vector<char*> _pv(len1+1);\n    std::vector<double*> _vv(len1+1);\n    for (int i=0; i<=len1; i++) { _pv[i]=path[i].data(); _vv[i]=val[i].data(); }\n    std::vector<char*> _mv(len1+1);\n    for (int i=0; i<=len1; i++) _mv[i]=mask[i].data();\n    std::vector<double*> _xv(x.size()), _yv(y.size());\n    for (size_t i=0; i<x.size(); i++) { _xv[i]=(double*)x[i].data(); _yv[i]=(double*)y[i].data(); }\n    NWDP_TM_dimer(_pv.data(), _vv.data(), _xv.data(), _yv.data(), len1, len2, _mv.data(), t, u, d02, gap_open, j2i);\n}"

new = "inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,\n    int len1, int len2, CharMatrix& mask,\n    double t[3], double u[3][3], double d02, double gap_open, int j2i[])\n{\n    std::vector<double*> _xv(x.size()), _yv(y.size());\n    for (size_t i=0; i<x.size(); i++) { _xv[i]=(double*)x[i].data(); _yv[i]=(double*)y[i].data(); }\n    NWDP_TM_dimer(path, val, _xv.data(), _yv.data(), len1, len2, mask, t, u, d02, gap_open, j2i);\n}"

if old in content:
    content = content.replace(old, new)
    print("Changed NWDP bridge to delegate to CharMatrix& version")
else:
    print("FAILED - trying to find bridge...")
    idx = content.find("NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y")
    if idx != -1:
        print(repr(content[idx:idx+500]))

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)

import subprocess
r = subprocess.run(["g++", "-O3", "-ffast-math", "-lm", "-static", "-o", "USalign_test.exe", "USalign.cpp"],
                   capture_output=True, text=True, timeout=120)
if r.returncode == 0:
    print("Compilation OK")
else:
    print("COMPILATION FAILED:")
    print(r.stderr[-500:])
