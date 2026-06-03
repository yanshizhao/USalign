with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Remove NWDP_TM_dimer temp views from get_initial5_dimer
old = "                {\n                    vector<double*> _xv2(x.size()), _yv2(y.size());\n                    for (size_t _i=0; _i<x.size(); _i++) _xv2[_i]=(double*)x[_i].data();\n                    for (size_t _i=0; _i<y.size(); _i++) _yv2[_i]=(double*)y[_i].data();\n                    NWDP_TM_dimer(path, val, _xv2.data(), _yv2.data(), xlen, ylen, mask,\n                        t, u, d02, gap_open, invmap);\n                }"
new = "                NWDP_TM_dimer(path, val, x, y, xlen, ylen, mask,\n                    t, u, d02, gap_open, invmap);"

if old in content:
    content = content.replace(old, new)
    print("1: Removed NWDP_TM_dimer temp views from caller")
else:
    print("1 FAILED")
    idx = content.find("_xv2")
    if idx != -1:
        print(repr(content[idx-20:idx+200]))

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)
print("Done")
