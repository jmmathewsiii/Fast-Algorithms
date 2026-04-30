# Compilers and flags
CXX       = g++
NVCC      = nvcc
CXXFLAGS  = -Wall -Wpedantic -O3 -Iinclude
NVCCFLAGS = -O3 -Iinclude
# LDFLAGS = -fsanitize=address

# Directories
SRC_DIR    = src
KERNEL_DIR = src/kernels
INC_DIR    = include
OBJ_DIR    = obj
KOBJ_DIR   = obj/kernels
BIN_DIR    = bin
OUT_DIR    = output
SCRIPT_DIR = scripts

# Shared objects (used by all binaries)
COMMON_OBJS = $(OBJ_DIR)/main.o \
              $(OBJ_DIR)/random.o \
              $(OBJ_DIR)/initialize.o \
              $(OBJ_DIR)/plotter.o \
              $(OBJ_DIR)/energy.o

# Targets
TARGET_SERIAL = $(BIN_DIR)/serial
TARGET_CUDA   = $(BIN_DIR)/cuda

# ===== Top-level build targets =====
all:    serial cuda
serial: $(TARGET_SERIAL)
cuda:   $(TARGET_CUDA)

# ===== Link rules =====
# The serial binary links both Direct and BH kernels and runs both back-to-back.
$(TARGET_SERIAL): $(COMMON_OBJS) \
                  $(OBJ_DIR)/integrator_cpu.o \
                  $(OBJ_DIR)/quadtree_cpu.o \
                  $(KOBJ_DIR)/direct_cpu.o \
                  $(KOBJ_DIR)/bh_cpu.o | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TARGET_CUDA): $(COMMON_OBJS) \
                $(OBJ_DIR)/integrator_gpu.o \
                $(OBJ_DIR)/quadtree_gpu.o \
                $(KOBJ_DIR)/direct_gpu.o \
                $(KOBJ_DIR)/bh_gpu.o | $(BIN_DIR)
	$(NVCC) $^ -o $@

# ===== Compile rules =====
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(KOBJ_DIR)/%.o: $(KERNEL_DIR)/%.cpp | $(KOBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu | $(OBJ_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(KOBJ_DIR)/%.o: $(KERNEL_DIR)/%.cu | $(KOBJ_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ===== Directories =====
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(KOBJ_DIR):
	mkdir -p $(KOBJ_DIR)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(SCRIPT_DIR):
	mkdir -p $(SCRIPT_DIR)

# ===== Clean =====
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

cleanall: clean
	rm -rf $(OUT_DIR)/* $(SCRIPT_DIR)/*

# ===== Run targets =====
run_serial: serial | $(OUT_DIR) $(SCRIPT_DIR)
	./$(TARGET_SERIAL)

run_cuda: cuda | $(OUT_DIR) $(SCRIPT_DIR)
	./$(TARGET_CUDA)

run_all: run_serial run_cuda

# ===== Phony targets =====
.PHONY: all serial cuda \
        clean cleanall \
        run_serial run_cuda run_all

# ===== Header dependencies =====
$(OBJ_DIR)/main.o:            $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h $(INC_DIR)/plotter.h $(INC_DIR)/integrator.h
$(OBJ_DIR)/random.o:          $(INC_DIR)/random.h
$(OBJ_DIR)/initialize.o:      $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h
$(OBJ_DIR)/plotter.o:         $(INC_DIR)/star.h $(INC_DIR)/plotter.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h
$(OBJ_DIR)/energy.o:          $(INC_DIR)/integrator.h $(INC_DIR)/star.h
$(OBJ_DIR)/integrator_cpu.o:  $(INC_DIR)/integrator.h $(INC_DIR)/force_cpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h
$(OBJ_DIR)/integrator_gpu.o:  $(INC_DIR)/integrator.h $(INC_DIR)/force_gpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_cpu.o:    $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_gpu.o:    $(INC_DIR)/quadtree_gpu.h
$(KOBJ_DIR)/direct_cpu.o:     $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_cpu.o:         $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h
$(KOBJ_DIR)/direct_gpu.o:     $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_gpu.o:         $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_gpu.h
