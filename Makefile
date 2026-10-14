CC=g++
MINGW=x86_64-w64-mingw32-g++
CFLAGS=-O3 -ffast-math -std=gnu++11
OUT_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
# Keep compiler temporaries inside the tree so the build does not depend on a
# writable system temp dir (e.g. C:\Windows\ when TEMP is unset under MSYS make).
BUILD_TMP:=$(CURDIR)/.build_tmp
$(shell mkdir -p $(BUILD_TMP))
export TMPDIR:=$(BUILD_TMP)
export TMP:=$(BUILD_TMP)
export TEMP:=$(BUILD_TMP)
INC=-Isrc/common -Isrc/engine -Isrc/engine/tmalign -Isrc/engine/flexalign -Isrc/engine/mmalign -Isrc/engine/mmdock -Isrc/engine/mtmalign -Isrc/engine/soialign -Isrc/engine/se -Isrc/engine/tmscore -Isrc/engine/hwrmsd -Isrc/engine/upgma -Isrc/output_funcs -Isrc/app/tools -Isrc/app/tools/tool_include
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

ENGINE_CPP=src/engine/tmalign/TMalign.cpp src/engine/flexalign/flexalign.cpp src/engine/tmscore/TMscore.cpp src/engine/hwrmsd/HwRMSD.cpp src/engine/soialign/SOIalign.cpp src/engine/se/se.cpp
USALIGN_ENGINE_CPP=src/engine/mmalign/MMalign.cpp src/engine/mmdock/MMdock.cpp src/engine/mtmalign/mTMalign.cpp
COMMON_CPP=src/common/common.cpp
OUTPUT_CPP=src/output_funcs/output.cpp
ENGINE_COMMON_CPP=src/engine/engine_common.cpp

all: USalign

qTMclust: src/app/tools/tool_main/qtmclust_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

USalign: src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

USalign.exe: src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${MINGW} ${STATIC_MINGW} ${CFLAGS} ${INC} ${OMPFLAGS} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

USFlexalign: src/app/tools/tool_main/usflexalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

TMalign: src/app/tools/tool_main/tmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

TMscore: src/app/tools/tool_main/tmscore_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

MMalign: src/app/tools/tool_main/mmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

se: src/app/tools/tool_main/se_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

pdb2ss: src/app/tools/tool_main/pdb2ss_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

pdb2xyz: src/app/tools/tool_main/pdb2xyz_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

xyz_sfetch: src/app/tools/tool_main/xyz_sfetch_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

pdb2fasta: src/app/tools/tool_main/pdb2fasta_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

biounitasym: src/app/tools/tool_main/biounitasym_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

NWalign: src/app/tools/tool_main/nwalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

HwRMSD: src/app/tools/tool_main/hwrmsd_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

cif2pdb: src/app/tools/tool_main/cif2pdb_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

pdbAtomName: src/app/tools/tool_main/pdbatomname_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

addChainID: src/app/tools/tool_main/addchainid_main.cpp $(COMMON_CPP)
		${CC} ${CFLAGS} ${INC} $^ -o $(OUT_DIR)$@ ${LDFLAGS}

clean:
		rm -f $(addprefix $(OUT_DIR),${PROGRAM})
		rm -rf $(BUILD_TMP)
