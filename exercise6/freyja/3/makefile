#SRC := main.c ns_differential.c integrator.c rk_classic.c lagrange_interpolation.c
CFLAGS = -O3 -Wall -g
COMPILER = clang 
LINK := -g -lm
INCLUDE_OPTIONS := 
TARGETS := poisson

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = inc






INCLUDES := $(shell find $(INC_DIR) -type d | sed s/^/-I/) $(INCLUDE_OPTIONS) 
SRC := $(shell find $(SRC_DIR) -name '*.c')
.PHONY:
all: depend directories $(TARGETS)
OBJ := $(subst $(SRC_DIR),$(OBJ_DIR),$(patsubst %.c,%.o,$(SRC)))
OBJDIRS := $(sort $(dir $(OBJ)))
SRCDIRS := $(sort $(dir $(SRC)))
$(TARGETS):	$(OBJ) 
	$(COMPILER) -o $@ $(OBJ) $(LINK)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c; $(COMPILER) $(CFLAGS) $(INCLUDES) -c $< -o $@
.PHONY:
database: ;make clean; bear make
directories:; mkdir -p $(OBJDIRS)
clean: ; rm -f $(TARGETS) $(OBJ) 
depend: ; ./gcc-makedepend $(addprefix -p ,$(OBJDIRS)) $(CFLAGS) $(INCLUDES) $(SRC)
# DO NOT DELETE
obj/linbcg.o: src/nr/linbcg.c inc/nr/nrutil.h
obj/nr/linbcg.o: src/nr/linbcg.c inc/nr/nrutil.h
obj/asolve.o: src/nr/asolve.c
obj/nr/asolve.o: src/nr/asolve.c
obj/atimes.o: src/nr/atimes.c
obj/nr/atimes.o: src/nr/atimes.c
obj/dsprstx.o: src/nr/dsprstx.c
obj/nr/dsprstx.o: src/nr/dsprstx.c
obj/nrutil.o: src/nr/nrutil.c
obj/nr/nrutil.o: src/nr/nrutil.c
obj/snrm.o: src/nr/snrm.c
obj/nr/snrm.o: src/nr/snrm.c
obj/dsprsax.o: src/nr/dsprsax.c
obj/nr/dsprsax.o: src/nr/dsprsax.c
obj/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
obj/nr/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
obj/laplace_matrix.o: src/laplace_matrix.c inc/laplace_matrix.h
obj/nr/laplace_matrix.o: src/laplace_matrix.c inc/laplace_matrix.h
obj/main.o: src/main.c inc/nr/nr.h inc/laplace_matrix.h inc/gnuplot_i.h
obj/nr/main.o: src/main.c inc/nr/nr.h inc/laplace_matrix.h inc/gnuplot_i.h
