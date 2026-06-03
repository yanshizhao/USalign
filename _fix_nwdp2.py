with open("MMalign.h", "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

# Remove the broken fprintf line
old = "    fprintf(stderr, \"DEBUG NWDP bridge called\n\");\n"
if old in content:
    content = content.replace(old, "")
    print("Removed broken fprintf")
else:
    # The string might have different escaping
    print("Trying different pattern...")
    old2 = '    fprintf(stderr, "DEBUG NWDP bridge called'
    idx = content.find(old2)
    if idx != -1:
        # Find the end of this statement
        end = content.find(");", idx)
        if end != -1:
            content = content[:idx] + content[end+2:]
            print(f"Removed fprintf from byte {idx} to {end+2}")

with open("MMalign.h", "w", encoding="utf-8", errors="replace") as f:
    f.write(content)

# Verify compilation
import subprocess
r = subprocess.run(["g++", "-O3", "-ffast-math", "-lm", "-static", "-o", "USalign_test.exe", "USalign.cpp"], 
                   capture_output=True, text=True, timeout=120)
if r.returncode == 0:
    print("Compilation OK")
else:
    print("COMPILATION FAILED:")
    print(r.stderr[-500:])
