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
INC=-Icommon -Iengine/tmalign -Iengine/flexalign -Iengine/mmalign -Iengine/mmdock -Iengine/mtmalign -Iengine/soialign -Iengine/se -Iengine/tmscore -Iengine/hwrmsd -Iengine/upgma -Ioutput_funcs -Iapp/tools -Iapp/tools/tool_include
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
PROGRAM=qTMclust USalign USFlexalign TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta biounitasym pdb2ss NWalign HwRMSD cif2pdb pdbAtomName addChainID

ENGINE_CPP=engine/tmalign/TMalign.cpp engine/flexalign/flexalign.cpp engine/tmscore/TMscore.cpp engine/hwrmsd/HwRMSD.cpp engine/soialign/SOIalign.cpp common/orchestration.cpp
USALIGN_ENGINE_CPP=engine/mmalign/MMalign.cpp engine/mmdock/MMdock.cpp engine/mtmalign/mTMalign.cpp
COMMON_CPP=common/common.cpp
OUTPUT_CPP=output_funcs/output.cpp

all: USalign

qTMclust: app/tools/tool_main/qtmclust_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

USalign: app/usalign/usalign.cpp engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $@ ${LDFLAGS}

USalign.exe: app/usalign/usalign.cpp engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${MINGW} ${STATIC_MINGW} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $@ ${LDFLAGS}

USFlexalign: app/tools/tool_main/usflexalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

TMalign: app/tools/tool_main/tmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

TMscore: app/tools/tool_main/tmscore_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

MMalign: app/tools/tool_main/mmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

se: app/tools/tool_main/se_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2ss: app/tools/tool_main/pdb2ss_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2xyz: app/tools/tool_main/pdb2xyz_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

xyz_sfetch: app/tools/tool_main/xyz_sfetch_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdb2fasta: app/tools/tool_main/pdb2fasta_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

biounitasym: app/tools/tool_main/biounitasym_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

NWalign: app/tools/tool_main/nwalign_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

HwRMSD: app/tools/tool_main/hwrmsd_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

cif2pdb: app/tools/tool_main/cif2pdb_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

pdbAtomName: app/tools/tool_main/pdbatomname_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

addChainID: app/tools/tool_main/addchainid_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $@ ${LDFLAGS}

clean:
		rm -f ${PROGRAM}
		rm -rf $(BUILD_TMP)
