CXX = g++ -std=c++11
C99 = gcc -std=c99
LINK = g++
AR = ar
#DEBUG_FLAG=-g
CXXFLAGS = -O3 -Wall -fPIC $(DEBUG_FLAG)
CFLAGS = $(CXXFLAGS)
ARFLAGS = -rv

SRC_DIR = ./src/
INCLUDE_DIR = ./include/
OUT_DIR = ./build

# OBJS = $(OUT_DIR)/objs/world/cheaptrick.o $(OUT_DIR)/objs/world/world_common.o $(OUT_DIR)/objs/world/d4c.o $(OUT_DIR)/objs/world/dio.o $(OUT_DIR)/objs/world/world_fft.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/world_matlabfunctions.o $(OUT_DIR)/objs/world/stonemask.o $(OUT_DIR)/objs/world/synthesis.o $(OUT_DIR)/objs/world/synthesisrealtime.o
OBJS =  $(OUT_DIR)/objs/world/world_common.o $(OUT_DIR)/objs/world/world_fft.o $(OUT_DIR)/objs/world/harvest.o $(OUT_DIR)/objs/world/world_matlabfunctions.o
LIBS =

all: default test

###############################################################################################################
### Tests
###############################################################################################################
test: $(OUT_DIR)/test

test_OBJS=$(OUT_DIR)/objs/tools/audioio.o $(OUT_DIR)/objs/test/test.o
$(OUT_DIR)/test: $(OUT_DIR)/libworld.a $(test_OBJS)
	$(LINK) $(CXXFLAGS) -o $(OUT_DIR)/test $(test_OBJS) $(OUT_DIR)/libworld.a -lm


# $(OUT_DIR)/objs/test/test.o : tools/audioio.h $(INCLUDE_DIR)d4c.h $(INCLUDE_DIR)dio.h $(INCLUDE_DIR)harvest.hpp $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)cheaptrick.h $(INCLUDE_DIR)stonemask.h $(INCLUDE_DIR)synthesis.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/test/ctest.o : tools/audioio.h $(INCLUDE_DIR)d4c.h $(INCLUDE_DIR)dio.h $(INCLUDE_DIR)harvest.hpp $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)cheaptrick.h $(INCLUDE_DIR)stonemask.h $(INCLUDE_DIR)synthesis.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)macrodefinitions.h

$(OUT_DIR)/objs/test/test.o : tools/audioio.h $(INCLUDE_DIR)harvest.hpp $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)macrodefinitions.h

###############################################################################################################
### Library
###############################################################################################################
default: $(OUT_DIR)/libworld.a

$(OUT_DIR)/libworld.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libworld.a $(OBJS) $(LIBS)
	@echo Done.

# $(OUT_DIR)/objs/world/cheaptrick.o : $(INCLUDE_DIR)cheaptrick.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
$(OUT_DIR)/objs/world/world_common.o : $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/world/d4c.o : $(INCLUDE_DIR)d4c.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/world/dio.o : $(INCLUDE_DIR)dio.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
$(OUT_DIR)/objs/world/world_fft.o : $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)macrodefinitions.h
$(OUT_DIR)/objs/world/harvest.o : $(INCLUDE_DIR)harvest.hpp $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
$(OUT_DIR)/objs/world/world_matlabfunctions.o : $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/world/stonemask.o : $(INCLUDE_DIR)stonemask.h $(INCLUDE_DIR)world_fft.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/world/synthesis.o : $(INCLUDE_DIR)synthesis.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h
# $(OUT_DIR)/objs/world/synthesisrealtime.o : $(INCLUDE_DIR)synthesisrealtime.h $(INCLUDE_DIR)world_common.h $(INCLUDE_DIR)world_constantnumbers.h $(INCLUDE_DIR)world_matlabfunctions.h $(INCLUDE_DIR)macrodefinitions.h


###############################################################################################################
### Global rules
###############################################################################################################
$(OUT_DIR)/objs/test/%.o : test/%.c
	mkdir -p $(OUT_DIR)/objs/test
	$(C99) $(CFLAGS) -I$(INCLUDE_DIR) -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/test/%.o : test/%.cpp
	mkdir -p $(OUT_DIR)/objs/test
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -Itools -o "$@" -c "$<"

$(OUT_DIR)/objs/tools/%.o : tools/%.cpp
	mkdir -p $(OUT_DIR)/objs/tools
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o "$@" -c "$<"

$(OUT_DIR)/objs/world/%.o : src/%.cpp
	mkdir -p $(OUT_DIR)/objs/world
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -o "$@" -c "$<"

clean:
	@echo 'Removing all temporary binaries... '
	@$(RM) $(OUT_DIR)/libworld.a $(OBJS)
	@$(RM) $(test_OBJS) $(ctest_OBJS) $(OUT_DIR)/test $(OUT_DIR)/ctest
	@echo Done.

clear: clean

.PHONY: clean clear test default
.DELETE_ON_ERRORS:
