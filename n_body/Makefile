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

# Shared objects (used by all four binaries)
COMMON_OBJS = $(OBJ_DIR)/main.o \
              $(OBJ_DIR)/random.o \
              $(OBJ_DIR)/initialize.o \
              $(OBJ_DIR)/plotter.o

# Targets
TARGET_SERIAL_DIRECT = $(BIN_DIR)/serial_direct
TARGET_SERIAL_BH     = $(BIN_DIR)/serial_bh
TARGET_CUDA_DIRECT   = $(BIN_DIR)/cuda_direct
TARGET_CUDA_BH       = $(BIN_DIR)/cuda_bh

# ===== Top-level build targets =====
all:    serial cuda
serial: serial_direct serial_bh
cuda:   cuda_direct cuda_bh

serial_direct: $(TARGET_SERIAL_DIRECT)
serial_bh:     $(TARGET_SERIAL_BH)
cuda_direct:   $(TARGET_CUDA_DIRECT)
cuda_bh:       $(TARGET_CUDA_BH)

# ===== Link rules =====
$(TARGET_SERIAL_DIRECT): $(COMMON_OBJS) $(OBJ_DIR)/integrator_cpu.o $(KOBJ_DIR)/direct_cpu.o | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TARGET_SERIAL_BH): $(COMMON_OBJS) $(OBJ_DIR)/integrator_cpu.o $(KOBJ_DIR)/bh_cpu.o $(OBJ_DIR)/quadtree_cpu.o | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TARGET_CUDA_DIRECT): $(COMMON_OBJS) $(OBJ_DIR)/integrator_gpu.o $(KOBJ_DIR)/direct_gpu.o | $(BIN_DIR)
	$(NVCC) $^ -o $@

$(TARGET_CUDA_BH): $(COMMON_OBJS) $(OBJ_DIR)/integrator_gpu.o $(KOBJ_DIR)/bh_gpu.o $(OBJ_DIR)/quadtree_gpu.o | $(BIN_DIR)
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

# ===== Clean =====
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

cleanall: clean
	rm -rf $(OUT_DIR)/*

# ===== Run targets =====
run_serial_direct: serial_direct | $(OUT_DIR)
	./$(TARGET_SERIAL_DIRECT)

run_serial_bh: serial_bh | $(OUT_DIR)
	./$(TARGET_SERIAL_BH)

run_cuda_direct: cuda_direct | $(OUT_DIR)
	./$(TARGET_CUDA_DIRECT)

run_cuda_bh: cuda_bh | $(OUT_DIR)
	./$(TARGET_CUDA_BH)

run_all: run_serial_direct run_serial_bh run_cuda_direct run_cuda_bh

# ===== Phony targets =====
.PHONY: all serial cuda \
        serial_direct serial_bh cuda_direct cuda_bh \
        clean cleanall \
        run_serial_direct run_serial_bh run_cuda_direct run_cuda_bh run_all

# ===== Header dependencies =====
$(OBJ_DIR)/main.o:            $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h $(INC_DIR)/plotter.h $(INC_DIR)/integrator.h
$(OBJ_DIR)/random.o:          $(INC_DIR)/random.h
$(OBJ_DIR)/initialize.o:      $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h
$(OBJ_DIR)/plotter.o:         $(INC_DIR)/star.h $(INC_DIR)/plotter.h
$(OBJ_DIR)/integrator_cpu.o:  $(INC_DIR)/integrator.h $(INC_DIR)/force_cpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h
$(OBJ_DIR)/integrator_gpu.o:  $(INC_DIR)/integrator.h $(INC_DIR)/force_gpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_cpu.o:    $(INC_DIR)/quadtree.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_gpu.o:    $(INC_DIR)/quadtree_gpu.h
$(KOBJ_DIR)/direct_cpu.o:     $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_cpu.o:         $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree.h
$(KOBJ_DIR)/direct_gpu.o:     $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_gpu.o:         $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_gpu.h
