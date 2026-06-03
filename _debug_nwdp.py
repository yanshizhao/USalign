with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Add cerr inside NWDP_TM_dimer bridge to see if it's reached
old = "inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,\n    int len1, int len2, CharMatrix& mask,\n    double t[3], double u[3][3], double d02, double gap_open, int j2i[])\n{"
new = "inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,\n    int len1, int len2, CharMatrix& mask,\n    double t[3], double u[3][3], double d02, double gap_open, int j2i[])\n{\n    fprintf(stderr, \"DEBUG NWDP bridge called\n\");"

if old in content:
    content = content.replace(old, new)
    print("Added debug fprintf to NWDP bridge")
else:
    print("FAILED")
    idx = content.find("NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x")
    if idx != -1:
        print(repr(content[idx:idx+200]))
    else:
        print("Bridge not found at all!")
        # Check if the char** version bridge exists
        idx2 = content.find("NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x")
        if idx2 >= 0:
            print(f"Found at {idx2}")
        else:
            # Maybe the old bridge was replaced - check current state
            idx3 = content.find("CoordArray& x, CoordArray& y")
            if idx3 >= 0:
                print(f"'CoordArray& x, CoordArray& y' found at {idx3}")
                print(repr(content[idx3-100:idx3+200]))

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)
print("Done")
