==============================================================================
   US-align: universal structure alignment of monomeric and complex proteins
   and nucleic acids

   References:
   * Chengxin Zhang, Lydia Freddolino, Yang Zhang
     (2026) Nat Protoc. 21, 517-541.
   * Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang
     (2022) Nat Methods. 19(9), 1109-1115.
   * Chengxin Zhang, Anna Marie Pyle (2022) iScience. 25(10), 105218.

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.

   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by J Yang
   2016/05/21: Several updates of this program were made by J Wu, including
              (1) fixed several compiling bugs
              (2) made I/O of C/C++ version consistent with the Fortran version
              (3) added outputs including full-atom and ligand structures
              (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
   2018/06/04: Several updates were made by C Zhang, including
              (1) Fixed bug in reading PDB files with negative residue index,
                  at the expense of the '-o' option now only being able to
                  output superposed structure instead of full rasmol script.
              (2) Implemented the fTM-align algorithm (by the '-fast' option)
                  as described in R Dong, S Pan, Z Peng, Y Zhang, J Yang
                  (2018) Nucleic acids research. gky430.
              (3) Included option to perform TM-align against a whole 
                  folder of PDB files. A full list of options not available
                  in the Fortran version can be explored by TMalign -h
   2018/07/27: Added the -byresi option for TM-score superposition without
               re-alignment as in TMscore and TMscore -c
   2018/08/07: Added the -dir option for all-against-all alignment within a
               folder
   2018/08/14: Added the -split option to split input structure by chain or
               by model
   2018/08/16: Added the -infmt1, -infmt2 options.
               TMalign can now read .gz and .bz2 compressed files.
   2018/10/20: C Zhang and S Gong updated the RNA alignment part of
               the program. Changes include:
               (1) new d0 calculation for RNA.
               (2) secondary structure assignment for RNA.
               (3) automatic detection of molecule type (protein vs RNA).
   2019/01/07: C Zhang added support for PDBx/mmCIF format.
   2019/02/09: Fixed asymmetric alignment bug where the output result is
               occasionally dependent on the order of input files.
   2019/03/17: Added the -cp option for circular permutation
   2019/03/27: Added the -mirror option for mirror structure alignment
   2019/04/25: The RNA-align algorithm was published by Bioinformatics
   2019/07/24: Fixed bug in displaying matching residues.
               Added GDT and MaxSub to TMscore program.
   2019/08/18: Prevent excessive circular permutation alignment by -cp.
   2020/05/19: Add back rasmol output by -rasmol
   2020/12/12: Fixed bug in double precision coordinate mmcif alignment
   2021/01/07: Fixed bug in TMscore -c for TM-score superposition between
               two oligomers
   2021/05/29: Remove unnecessary depedency on malloc.h, which prevent
               compilation on Mac OS
   2021/08/17: Complete implementation of MMalign for oligomer alignment
   2021/10/03: Support Windows
   2022/02/27: Add -seq (-byresi 4 & 5) for TM-score superimposition guided by
               sequence alignment.
   2022/04/12: Support AlphaFold CIF
   2022/05/11: Update -mm 4 output format
   2022/05/24: Limited support for sequence order independent alignment
   2022/05/30: Correct atom pair output for -mm 5
   2022/06/07: Sequence order semi-independent alignment
   2022/06/20: Ensure sequentiality within SSE in sequence order
               semi-independent alignment
   2022/06/22: Fix infinite loop for mal-formatted PDB
   2022/06/23: Fix -m for Windows. Add pymol plugin.
   2022/06/26: Add -full option for -mm 2 and 4 to show chain level alignment
   2022/09/24: Support -TMscore for complex when the chain order is different
   2023/06/09: Correct atom name justification in PDB file for standard amino
               acids and nucleotides
   2023/12/13: Refine chain assignment for -TMscore >=6 for large complex
   2023/12/22: Forbid chain assignment refinment for -chainmap
   2024/03/04: Add -chain1 and -chain2 options to specify chains to align
   2024/03/19: Chain mapping for trimer -TMscore
   2024/05/10: Better handling of alternative location indicator, e.g., 5e1n
   2024/06/02: Allow using -mm 4 with -se
   2024/07/05: Add -do to output distances for aligned residue pairs
   2024/07/30: Allow using -se with -byresi 6 7
   2024/10/30: Set default for -ter and -split based on -mm output
   2024/11/08: Add -chimerax for UCSF ChimeraX format output
   2026/03/28: Fix -mm 1 asymmetric alignment bug where oligomer alignment
	       output depends occasionally on input file order
   2026/05/27: Allow both -mm 1 and -dir/-dir1/-dir2
   2026/07/06: Add -threads N to set number of threads for OpenMP parallel
	       execution. For example, '-threads 4' uses 4 threads.
   2026/10/07: Output files now default to a "results" folder created next to
               the executable; change it with '-outdir <dir>' ('-outdir .' for
               the current directory). Source files reorganized into src/common/,
               src/engine/, orchestration/, src/output_funcs/ and app/; "make" now
               builds only USalign by default (list other tools explicitly);
               example inputs moved to the "example" folder.
   2026/10/08: src/app/usalign/usalign.cpp reduced to CLI routing only (parsing,
               validation, dispatch); each alignment mode's orchestration
               moved back into its engine under src/engine/<mode>/. New standalone
               "USFlexalign" executable for flexible alignment (-mm 7). The
               output module merged (tmalign_output.{h,cpp} -> output.{h,cpp}).
               Shared CLI helpers (build_file_lists / check_* families)
               deduplicated into src/common/common.h templates.
===============================================================================

================================
 Source code organization
================================
All C++ sources live under src/, organized into four layers, from the bottom up:

    src/common/       foundation shared by all other layers (common.h holds the
                      shared data structures, the CLI helper templates and the
                      single chain-pair dispatcher align_chain_pair, implemented
                      in src/common/orchestration.cpp)
    src/engine/       alignment engines (one sub-directory per engine):
                      tmalign/ mmalign/ mmdock/ mtmalign/ soialign/ flexalign/
                      tmscore/ hwrmsd/ se/ upgma/
    src/output_funcs/ formatting of alignment results into text / PDB / scripts
    src/app/          the USalign entry point (src/app/usalign/usalign.cpp +
                      src/app/usalign/usalign.h, routing only) and the other
                      standalone programs under src/app/tools/ (implementations
                      in tool_main/, declarations and the shared tools.h umbrella
                      header in tool_include/)

Each alignment mode (-mm 0/1/2/4/5/6/7) is driven by its engine: the chain
loop, per-pair alignment and result output live in the engine's own
.h/.cpp, while src/app/usalign/usalign.cpp only parses arguments, validates them and
dispatches to the engine entry point.

=========================
 How to install US-align
=========================
1. Linux / macOS
   1.1 Compile with make (Recommended):
           make clean
           make
       "make" builds the USalign program by default (see the section
       "Building the individual programs under app/" below for the other
       programs).
       On macOS the Makefile uses Homebrew GNU g++ automatically when it is
       installed (/opt/homebrew/bin/g++-* or /usr/local/bin/g++-*). Otherwise
       it falls back to Apple Clang and disables -fopenmp, because Apple Clang
       has no OpenMP support. USFlexalign builds with either compiler; USalign
       needs GNU g++ for OpenMP acceleration (brew install gcc).
   1.2 Manual compilation in bash (do not add -static; macOS cannot link
       statically):
           g++ -O3 -ffast-math -std=gnu++11 -fopenmp -Isrc/common -Isrc/engine -Isrc/engine/tmalign -Isrc/engine/flexalign -Isrc/engine/mmalign -Isrc/engine/mmdock -Isrc/engine/mtmalign -Isrc/engine/soialign -Isrc/engine/se -Isrc/engine/tmscore -Isrc/engine/hwrmsd -Isrc/engine/upgma -Isrc/output_funcs -Isrc/app/tools -o USalign src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp src/engine/tmalign/TMalign.cpp src/engine/mmalign/MMalign.cpp src/engine/mmdock/MMdock.cpp src/engine/mtmalign/mTMalign.cpp src/engine/flexalign/flexalign.cpp src/engine/tmscore/TMscore.cpp src/engine/hwrmsd/HwRMSD.cpp src/engine/soialign/SOIalign.cpp src/engine/se/se.cpp src/common/common.cpp src/engine/engine_common.cpp src/output_funcs/output.cpp -lm

2. Windows (MSYS2 / MinGW)
   2.1 Compile with make:
           mingw32-make clean
           mingw32-make
       Same targets as in 1.1: only USalign is built by default.
   2.2 Manual compilation in PowerShell (-static is required on Windows):
           g++ -static -O3 -ffast-math -std=gnu++11 -fopenmp -Isrc/common -Isrc/engine -Isrc/engine/tmalign -Isrc/engine/flexalign -Isrc/engine/mmalign -Isrc/engine/mmdock -Isrc/engine/mtmalign -Isrc/engine/soialign -Isrc/engine/se -Isrc/engine/tmscore -Isrc/engine/hwrmsd -Isrc/engine/upgma -Isrc/output_funcs -Isrc/app/tools -o USalign src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp src/engine/tmalign/TMalign.cpp src/engine/mmalign/MMalign.cpp src/engine/mmdock/MMdock.cpp src/engine/mtmalign/mTMalign.cpp src/engine/flexalign/flexalign.cpp src/engine/tmscore/TMscore.cpp src/engine/hwrmsd/HwRMSD.cpp src/engine/soialign/SOIalign.cpp src/engine/se/se.cpp src/common/common.cpp src/engine/engine_common.cpp src/output_funcs/output.cpp -lm

3. Note: the "make" command automatically includes the -fopenmp flag to enable OpenMP parallel acceleration for the USalign program. If you do not need OpenMP parallel acceleration, simply remove the -fopenmp flag during compilation.

Additional note:
1. The -static flag enables static linking. A statically linked executable has no external runtime dependencies and can be copied to and run directly on other compatible machines.
2. However, on Linux the "make" command does NOT pass -static: the resulting binary is dynamically linked and depends on the system libraries (e.g. libstdc++, libgomp and glibc), so copying it to another machine may fail if those libraries are missing or too old.
3. On Windows, the Makefile adds -static automatically when the OS environment
   variable is Windows_NT. This depends on the OS variable, not on whether you
   invoke make or mingw32-make. In a non-standard shell or toolchain, check that
   the actual compilation command contains -static.
4. If you need a statically linked Linux binary, add -static to the command in 1.2.
5. macOS does not support static linking.
6. Remove -fopenmp from all manual compilation commands if parallel acceleration is not required.

Compressed file support:
1. On Linux, macOS and WSL2, US-align reads gz and bz2 compressed files in
   addition to uncompressed files, provided that the "gunzip" and "bzcat"
   commands are respectively available.
2. Due to the lack of POSIX support, US-align natively compiled on Windows
   without WSL2 cannot read compressed files.

Compiler requirements:
1. US-align requires a compiler with C++11 support. It is known to be compilable by g++ version 6 or later (C++14 is the default standard since g++ 6), clang++ version 12.0.5 or later and mingw-w64 version 9.3 or later.
2. The Makefile sets -std=gnu++11 explicitly, so the build always compiles in C++11 mode regardless of the compiler's default standard. g++ 4.8.x/5.x defaults to C++98, but it will still work as long as it supports the gnu++11 mode. If you compile manually, remember to include -std=gnu++11 (see the manual compilation commands above).

Install via package managers:
1. If you are using Homebrew, run the following command. The executable is supposed to work on MacOS (x86_64, arm64) and Linux (x86_64).

```shell
brew install brewsci/bio/usalign
```

2. If you are using conda, run the following command. The executable is supposed to work on Linux (x86_64, aarch64) and MacOS (x86_64, arm64).

```shell
conda install -c bioconda usalign
```

3. The Homebrew and conda packages may not have the most recent version of US-align.

=====================
 How to use US-align
=====================
You can run the program without arguments to obtain a brief instruction

    ./USalign example/PDB1.pdb example/PDB2.pdb

(Example input files PDB1.pdb, PDB2.pdb and align.txt are in the "example"
folder.)

A full list of available options can be explored by:

    ./USalign -h

Flexible (hinge) alignment can also be run as a standalone program:

    ./USFlexalign example/PDB1.pdb example/PDB2.pdb

which is equivalent to "./USalign ... -mm 7" plus the flexible-alignment
options (-hinge / -afp / -TMpass).

Output files: by default, the files written by the program (-pymol, -m,
-rasmol, -chimerax and the upgma_tree files of -mm 4) are placed in a "results"
folder created next to the executable. Use

    -outdir <dir>    write output files to the given directory
    -outdir .        write output files to the current directory

to change this. Results that are not saved to a file are printed to stdout.

=========================
 How to view the results
=========================
The PyMOL (.pml) and ChimeraX (.cxc) scripts reference the structure files
by file name only, so enter the output directory before launching the viewer.
(The RasMol script embeds the coordinates inline and can be run from any
directory.)

    cd results
    pymol sup.pml
    chimerax sup_all_atm.cxc
    rasmol -script sup

Within an already-running PyMOL session, enter the output directory first and
then load the script:

    cd results
    @sup.pml

The scripts differ in the region and representation they display:

    sup              C-alpha trace of the aligned region
    sup_atm          full-atom representation of the aligned region
    sup_all          C-alpha trace of the whole chain
    sup_all_atm      full-atom representation of the whole chain
    sup_all_atm_lig  full-atom representation including all molecules and ligands

Note:
1. The option that outputs the superposed structure and PyMOL scripts, formerly
   "-o", is now named "-pymol". This is a backward-incompatible change.
2. With "-mm 2" on a multi-chain structure, "-pymol" writes the superposed
   structure of each chain to a separate file (e.g. "sup.0.pdb", "sup.1.pdb")
   and generates no scene script.
3. Structure 2 is written to a copy named "<name>_structure2.<ext>" in the
   output directory (decompressed first if it is a ".gz" or ".bz2" file), and
   the scripts reference this copy so that the superposition can be loaded
   directly.

========================================================
 Building the individual programs under app/
========================================================
"make" builds only USalign by default. Every other program has its own main()
under app/ and is built by listing its target name explicitly:

Target        Build command           Source
----------    --------------------    ------------------------------------------
USalign       make                    src/app/usalign/usalign.cpp (default target)
USFlexalign   make USFlexalign        src/app/tools/tool_main/usflexalign_main.cpp
TMalign       make TMalign            src/app/tools/tool_main/tmalign_main.cpp
TMscore       make TMscore            src/app/tools/tool_main/tmscore_main.cpp
MMalign       make MMalign            src/app/tools/tool_main/mmalign_main.cpp
se            make se                 src/app/tools/tool_main/se_main.cpp
HwRMSD        make HwRMSD             src/app/tools/tool_main/hwrmsd_main.cpp
qTMclust      make qTMclust           src/app/tools/tool_main/qtmclust_main.cpp
NWalign       make NWalign            src/app/tools/tool_main/nwalign_main.cpp
pdb2xyz       make pdb2xyz            src/app/tools/tool_main/pdb2xyz_main.cpp
pdb2fasta     make pdb2fasta          src/app/tools/tool_main/pdb2fasta_main.cpp
pdb2ss        make pdb2ss             src/app/tools/tool_main/pdb2ss_main.cpp
cif2pdb       make cif2pdb            src/app/tools/tool_main/cif2pdb_main.cpp
biounitasym   make biounitasym        src/app/tools/tool_main/biounitasym_main.cpp
pdbAtomName   make pdbAtomName        src/app/tools/tool_main/pdbatomname_main.cpp
addChainID    make addChainID         src/app/tools/tool_main/addchainid_main.cpp
xyz_sfetch    make xyz_sfetch         app/tools/xyz_sfetch_main.cpp

On Windows use "mingw32-make <target>". To compile a program manually, use the
same command as in the "Install" section but replace "src/app/usalign/usalign.cpp" with the
program's source file above, and include only the src/engine/common/output sources
that the program needs (the Makefile rule for each target lists them).

Run any program without arguments (or with "-h" where supported) to print its
usage.
