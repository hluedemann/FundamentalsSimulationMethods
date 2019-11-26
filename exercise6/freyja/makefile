#SRC := main.c ns_differential.c integrator.c rk_classic.c lagrange_interpolation.c
CFLAGS = -O3 -Wall -g
COMPILER = clang 
LINK := -g -lm
INCLUDE_OPTIONS := 
TARGETS := diffusion

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
obj/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
obj/tridiag.o: src/tridiag.c inc/tridiag.h
obj/heat.o: src/heat.c inc/heat.h inc/tridiag.h
obj/main.o: src/main.c inc/tridiag.h inc/heat.h inc/gnuplot_i.h
