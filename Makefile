ifeq ($(OS),Windows_NT)
    HOST_OS:=Windows_NT
else
    HOST_OS:=$(shell uname -s 2>/dev/null)
endif

# Use CXX/CXXFLAGS for C++ builds.  Keep command-line CC/CFLAGS overrides
# compatible with older build instructions.
ifeq ($(HOST_OS),Darwin)
    HOMEBREW_GCC_PREFIX:=$(shell command -v brew >/dev/null 2>&1 && brew --prefix gcc 2>/dev/null)
    HOMEBREW_GXX:=$(firstword $(wildcard $(HOMEBREW_GCC_PREFIX)/bin/g++-*))
endif

ifneq ($(filter default undefined,$(origin CXX)),)
    ifeq ($(origin CC),command line)
        CXX:=$(CC)
    else ifneq ($(HOMEBREW_GXX),)
        CXX:=$(HOMEBREW_GXX)
    else
        CXX:=g++
    endif
endif

MINGW?=x86_64-w64-mingw32-g++
FAST_MATH?=1
REPRODUCIBLE?=0
ifeq ($(REPRODUCIBLE),1)
    FAST_MATH:=0
endif
BASE_CXXFLAGS=-O3 -std=gnu++11
ifeq ($(FAST_MATH),1)
    BASE_CXXFLAGS+=-ffast-math
else ifneq ($(FAST_MATH),0)
    $(error FAST_MATH must be 0 or 1)
endif
CFLAGS?=$(BASE_CXXFLAGS)
CXXFLAGS?=$(CFLAGS)
OUT_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
INC=-Isrc/common -Isrc/engine -Isrc/engine/tmalign -Isrc/engine/flexalign -Isrc/engine/mmalign -Isrc/engine/mmdock -Isrc/engine/mtmalign -Isrc/engine/soialign -Isrc/engine/se -Isrc/engine/tmscore -Isrc/engine/hwrmsd -Isrc/engine/upgma -Isrc/output_funcs -Isrc/app/tools -Isrc/app/tools/tool_include

# OPENMP=auto enables OpenMP only when the compiler used by that target can
# compile and link a minimal program. OPENMP=0 is the portable serial build.
OPENMP?=auto
ifneq ($(origin OMPFLAGS),undefined)
    OMP_CXXFLAGS:=$(OMPFLAGS)
    OMP_LDFLAGS:=$(OMPFLAGS)
    MINGW_OMP_CXXFLAGS:=$(OMPFLAGS)
    MINGW_OMP_LDFLAGS:=$(OMPFLAGS)
else
    ifeq ($(OPENMP),0)
        OMP_CXXFLAGS:=
        OMP_LDFLAGS:=
        MINGW_OMP_CXXFLAGS:=
        MINGW_OMP_LDFLAGS:=
    else ifeq ($(OPENMP),1)
        OMP_CXXFLAGS:=-fopenmp
        OMP_LDFLAGS:=-fopenmp
        MINGW_OMP_CXXFLAGS:=-fopenmp
        MINGW_OMP_LDFLAGS:=-fopenmp
    else ifeq ($(OPENMP),auto)
        ifeq ($(strip $(MAKECMDGOALS)),clean)
            OMP_CXXFLAGS:=
            OMP_LDFLAGS:=
        else
            OPENMP_SUPPORTED:=$(shell printf 'int main(){return 0;}\n' | $(CXX) -x c++ -fopenmp -o /dev/null - >/dev/null 2>&1 && printf yes)
            ifeq ($(OPENMP_SUPPORTED),yes)
                OMP_CXXFLAGS:=-fopenmp
                OMP_LDFLAGS:=-fopenmp
            else
                OMP_CXXFLAGS:=
                OMP_LDFLAGS:=
                $(warning $(CXX) does not support -fopenmp; building USalign without OpenMP. Select an OpenMP-capable compiler or use OPENMP=1 to force it.)
            endif
        endif
        ifneq ($(filter %.exe,$(MAKECMDGOALS)),)
            MINGW_OPENMP_SUPPORTED:=$(shell printf 'int main(){return 0;}\n' | $(MINGW) -x c++ -fopenmp -o /dev/null - >/dev/null 2>&1 && printf yes)
            ifeq ($(MINGW_OPENMP_SUPPORTED),yes)
                MINGW_OMP_CXXFLAGS:=-fopenmp
                MINGW_OMP_LDFLAGS:=-fopenmp
            else
                MINGW_OMP_CXXFLAGS:=
                MINGW_OMP_LDFLAGS:=
                $(warning $(MINGW) does not support -fopenmp; building the Windows target without OpenMP.)
            endif
        endif
    else
        $(error OPENMP must be auto, 0, or 1)
    endif
endif
ifeq ($(OS),Windows_NT)
    EXEEXT=.exe
    STATIC=-static
    STATIC_MINGW=
else
    EXEEXT=
    STATIC=
    STATIC_MINGW=-static
endif
LDFLAGS?=$(STATIC)
LDLIBS?=-lm
PROGRAM=qTMclust USalign USFlexalign TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta biounitasym pdb2ss NWalign HwRMSD cif2pdb pdbAtomName addChainID
PROGRAM_TARGETS=$(addsuffix $(EXEEXT),$(PROGRAM))

ENGINE_CPP=src/engine/tmalign/TMalign.cpp src/engine/flexalign/flexalign.cpp src/engine/tmscore/TMscore.cpp src/engine/hwrmsd/HwRMSD.cpp src/engine/soialign/SOIalign.cpp src/engine/se/se.cpp
USALIGN_ENGINE_CPP=src/engine/mmalign/MMalign.cpp src/engine/mmdock/MMdock.cpp src/engine/mtmalign/mTMalign.cpp
COMMON_CPP=src/common/common.cpp
OUTPUT_CPP=src/output_funcs/output.cpp
ENGINE_COMMON_CPP=src/engine/engine_common.cpp
HEADERS=$(wildcard src/*.h src/*/*.h src/*/*/*.h src/*/*/*/*.h src/*.hpp src/*/*.hpp src/*/*/*.hpp src/*/*/*/*.hpp)

.PHONY: all clean

all: USalign$(EXEEXT)

qTMclust$(EXEEXT): src/app/tools/tool_main/qtmclust_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

USalign$(EXEEXT): src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(OMP_CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(OMP_LDFLAGS) $(LDLIBS)

ifneq ($(OS),Windows_NT)
USalign.exe: src/app/usalign/usalign.cpp src/engine/upgma/UPGMA.cpp $(ENGINE_CPP) $(USALIGN_ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(MINGW) $(CXXFLAGS) $(MINGW_OMP_CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(STATIC_MINGW) $(MINGW_OMP_LDFLAGS) $(LDLIBS)
endif

USFlexalign$(EXEEXT): src/app/tools/tool_main/usflexalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

TMalign$(EXEEXT): src/app/tools/tool_main/tmalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

TMscore$(EXEEXT): src/app/tools/tool_main/tmscore_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

MMalign$(EXEEXT): src/app/tools/tool_main/mmalign_main.cpp src/engine/mmalign/MMalign.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

se$(EXEEXT): src/app/tools/tool_main/se_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

pdb2ss$(EXEEXT): src/app/tools/tool_main/pdb2ss_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

pdb2xyz$(EXEEXT): src/app/tools/tool_main/pdb2xyz_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

xyz_sfetch$(EXEEXT): src/app/tools/tool_main/xyz_sfetch_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

pdb2fasta$(EXEEXT): src/app/tools/tool_main/pdb2fasta_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

biounitasym$(EXEEXT): src/app/tools/tool_main/biounitasym_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

NWalign$(EXEEXT): src/app/tools/tool_main/nwalign_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

HwRMSD$(EXEEXT): src/app/tools/tool_main/hwrmsd_main.cpp $(ENGINE_CPP) $(COMMON_CPP) $(ENGINE_COMMON_CPP) $(OUTPUT_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

cif2pdb$(EXEEXT): src/app/tools/tool_main/cif2pdb_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

pdbAtomName$(EXEEXT): src/app/tools/tool_main/pdbatomname_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

addChainID$(EXEEXT): src/app/tools/tool_main/addchainid_main.cpp $(COMMON_CPP) $(HEADERS)
		$(CXX) $(CXXFLAGS) $(INC) $(filter %.cpp,$^) -o "$(OUT_DIR)$@" $(LDFLAGS) $(LDLIBS)

ifeq ($(OS),Windows_NT)
.PHONY: $(PROGRAM)
$(PROGRAM): %: %.exe
endif

clean:
		$(RM) $(addprefix $(OUT_DIR),$(PROGRAM) $(addsuffix .exe,$(PROGRAM)))
