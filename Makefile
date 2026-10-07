CC=g++
MINGW=x86_64-w64-mingw32-g++
CFLAGS=-O3 -ffast-math -std=gnu++11
# Keep compiler temporaries inside the tree so the build does not depend on a
# writable system temp dir (e.g. C:\Windows\ when TEMP is unset under MSYS make).
BUILD_TMP:=$(CURDIR)/.build_tmp
$(shell mkdir -p $(BUILD_TMP))
export TMPDIR:=$(BUILD_TMP)
export TMP:=$(BUILD_TMP)
export TEMP:=$(BUILD_TMP)
INC=-Icommon -Iengine/tmalign -Iengine/flexalign -Iengine/mmalign -Iengine/soialign -Iengine/se -Iengine/tmscore -Iengine/hwrmsd -Iengine/upgma -Iorchestration -Ioutput_funcs -Iapp/tools
OMPFLAGS=-fopenmp
ifeq ($(OS),Windows_NT)
    # Windows (MSYS2, MinGW, PowerShell, cmd)
    STATIC=-static
    STATIC_MINGW=
else
    # Linux / macOS (native: no -static, cross-compile: STATIC_MINGW adds it)
    STATIC=
    STATIC_MINGW=-static
endif
LDFLAGS=$(STATIC) -lm
PROGRAM=qTMclust USalign TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta biounitasym pdb2ss NWalign HwRMSD cif2pdb pdbAtomName addChainID

ENGINE_CPP=engine/tmalign/TMalign.cpp engine/flexalign/flexalign.cpp engine/tmscore/TMscore.cpp engine/hwrmsd/HwRMSD.cpp engine/soialign/SOIalign.cpp
COMMON_CPP=common/common.cpp
ORCH_CPP=orchestration/orchestration.cpp
OUTPUT_CPP=output_funcs/output.cpp output_funcs/tmalign_output.cpp

all: USalign

qTMclust: app/tools/qtmclust_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

USalign: app/usalign.cpp engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ORCH_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $@ ${LDFLAGS}

USalign.exe: app/usalign.cpp engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ORCH_CPP) $(OUTPUT_CPP)
		${MINGW} ${STATIC_MINGW} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $@ ${LDFLAGS}

TMalign: app/tools/tmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

TMscore: app/tools/tmscore_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

MMalign: app/tools/mmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

se: app/tools/se_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2ss: app/tools/pdb2ss_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2xyz: app/tools/pdb2xyz_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

xyz_sfetch: app/tools/xyz_sfetch_main.cpp
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2fasta: app/tools/pdb2fasta_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

biounitasym: app/tools/biounitasym_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

NWalign: app/tools/nwalign_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

HwRMSD: app/tools/hwrmsd_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

cif2pdb: app/tools/cif2pdb_main.cpp
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdbAtomName: app/tools/pdbatomname_main.cpp
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

addChainID: app/tools/addchainid_main.cpp
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

clean:
		rm -f ${PROGRAM}
		rm -rf $(BUILD_TMP)
