# Compilers and flags
CXX       = g++
NVCC      = nvcc
CXXFLAGS  = -Wall -Wpedantic -O3 -Iinclude
NVCCFLAGS = -O3 -Iinclude
# LDFLAGS = -fsanitize=address

# Directories
SRC_DIR     = src
KERNEL_DIR  = src/kernels
INC_DIR     = include
OBJ_DIR     = obj
KOBJ_DIR    = obj/kernels
BIN_DIR     = bin
OUT_DIR     = output
LYAP_DIR    = output/lyapunov
SCRIPT_DIR  = scripts
LYAP_SCRIPT = scripts/lyapunov

# Shared objects (used by all binaries)
COMMON_OBJS = $(OBJ_DIR)/random.o \
              $(OBJ_DIR)/initialize.o \
              $(OBJ_DIR)/plotter.o \
              $(OBJ_DIR)/energy.o \
              $(OBJ_DIR)/validate.o

# Targets
TARGET_SERIAL_DIRECT = $(BIN_DIR)/serial_direct
TARGET_SERIAL_BH     = $(BIN_DIR)/serial_bh
TARGET_CUDA_DIRECT   = $(BIN_DIR)/cuda_direct

# ===== Top-level build targets =====
all:           serial_direct serial_bh cuda_direct
serial_direct: $(TARGET_SERIAL_DIRECT)
serial_bh:     $(TARGET_SERIAL_BH)
cuda_direct:   $(TARGET_CUDA_DIRECT)

# ===== Link rules =====
$(TARGET_SERIAL_DIRECT): $(OBJ_DIR)/main.serial_direct.o $(COMMON_OBJS) \
                         $(OBJ_DIR)/integrator_cpu.o \
                         $(OBJ_DIR)/quadtree_cpu.o \
                         $(KOBJ_DIR)/direct_cpu.o \
                         $(KOBJ_DIR)/bh_cpu.o | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TARGET_SERIAL_BH): $(OBJ_DIR)/main.serial_bh.o $(COMMON_OBJS) \
                     $(OBJ_DIR)/integrator_cpu.o \
                     $(OBJ_DIR)/quadtree_cpu.o \
                     $(KOBJ_DIR)/direct_cpu.o \
                     $(KOBJ_DIR)/bh_cpu.o | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TARGET_CUDA_DIRECT): $(OBJ_DIR)/main.cuda_direct.o $(COMMON_OBJS) \
                       $(OBJ_DIR)/integrator_gpu.o \
                       $(OBJ_DIR)/quadtree_cpu.o \
                       $(OBJ_DIR)/quadtree_gpu.o \
                       $(KOBJ_DIR)/direct_gpu.o \
                       $(KOBJ_DIR)/bh_gpu.o | $(BIN_DIR)
	$(NVCC) $^ -o $@

# ===== Per-binary main.o (selects algorithm via preprocessor) =====
$(OBJ_DIR)/main.serial_direct.o: $(SRC_DIR)/main.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -DBUILD_SERIAL_DIRECT -c $< -o $@

$(OBJ_DIR)/main.serial_bh.o: $(SRC_DIR)/main.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -DBUILD_SERIAL_BH -c $< -o $@

$(OBJ_DIR)/main.cuda_direct.o: $(SRC_DIR)/main.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -DBUILD_CUDA_DIRECT -c $< -o $@

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

$(LYAP_DIR):
	mkdir -p $(LYAP_DIR)

$(SCRIPT_DIR):
	mkdir -p $(SCRIPT_DIR)

$(LYAP_SCRIPT):
	mkdir -p $(LYAP_SCRIPT)

# ===== Clean =====
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

cleanall: clean
	rm -rf $(OUT_DIR)/* $(SCRIPT_DIR)/*

# ===== Run targets =====
run-naive: serial_direct | $(OUT_DIR) $(LYAP_DIR) $(SCRIPT_DIR) $(LYAP_SCRIPT)
	./$(TARGET_SERIAL_DIRECT)

run-serial-BH: serial_bh | $(OUT_DIR) $(LYAP_DIR) $(SCRIPT_DIR) $(LYAP_SCRIPT)
	./$(TARGET_SERIAL_BH)

run-CUDA-direct: cuda_direct | $(OUT_DIR) $(LYAP_DIR) $(SCRIPT_DIR) $(LYAP_SCRIPT)
	./$(TARGET_CUDA_DIRECT)

# ===== Phony targets =====
.PHONY: all serial_direct serial_bh cuda_direct \
        clean cleanall \
        run-naive run-serial-BH run-CUDA-direct

# ===== Header dependencies =====
$(OBJ_DIR)/main.serial_direct.o: $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h $(INC_DIR)/plotter.h $(INC_DIR)/integrator.h $(INC_DIR)/validate.h
$(OBJ_DIR)/main.serial_bh.o:     $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h $(INC_DIR)/plotter.h $(INC_DIR)/integrator.h $(INC_DIR)/validate.h
$(OBJ_DIR)/main.cuda_direct.o:   $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h $(INC_DIR)/plotter.h $(INC_DIR)/integrator.h $(INC_DIR)/validate.h
$(OBJ_DIR)/validate.o:           $(INC_DIR)/validate.h $(INC_DIR)/star.h
$(OBJ_DIR)/random.o:             $(INC_DIR)/random.h
$(OBJ_DIR)/initialize.o:         $(INC_DIR)/star.h $(INC_DIR)/random.h $(INC_DIR)/initialize.h
$(OBJ_DIR)/plotter.o:            $(INC_DIR)/star.h $(INC_DIR)/plotter.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h
$(OBJ_DIR)/energy.o:             $(INC_DIR)/integrator.h $(INC_DIR)/star.h
$(OBJ_DIR)/integrator_cpu.o:     $(INC_DIR)/integrator.h $(INC_DIR)/force_cpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h $(INC_DIR)/validate.h
$(OBJ_DIR)/integrator_gpu.o:     $(INC_DIR)/integrator.h $(INC_DIR)/force_gpu.h $(INC_DIR)/plotter.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_cpu.o:       $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h $(INC_DIR)/star.h
$(OBJ_DIR)/quadtree_gpu.o:       $(INC_DIR)/quadtree_gpu.h
$(KOBJ_DIR)/direct_cpu.o:        $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_cpu.o:            $(INC_DIR)/force_cpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_cpu.h $(INC_DIR)/node.h
$(KOBJ_DIR)/direct_gpu.o:        $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h
$(KOBJ_DIR)/bh_gpu.o:            $(INC_DIR)/force_gpu.h $(INC_DIR)/star.h $(INC_DIR)/quadtree_gpu.h
