with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Change inline to static for the NWDP bridge
old = "inline void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,"
new = "static void NWDP_TM_dimer(CharMatrix& path, DoubleMatrix& val, CoordArray& x, CoordArray& y,"

if old in content:
    content = content.replace(old, new)
    print("Changed inline → static for NWDP bridge")

    # Compile
    import subprocess
    r = subprocess.run(["g++", "-O3", "-ffast-math", "-lm", "-static", "-o", "USalign_test.exe", "USalign.cpp"],
                       capture_output=True, text=True, timeout=120)
    if r.returncode == 0:
        print("Compilation OK")
    else:
        print("COMPILATION FAILED:")
        print(r.stderr[-500:])
else:
    print("FAILED")

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)
