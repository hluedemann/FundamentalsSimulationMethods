#SRC := main.c ns_differential.c integrator.c rk_classic.c lagrange_interpolation.c
CFLAGS = -O3 -Wall
COMPILER = clang 
LINK := -lm 
INCLUDE_OPTIONS := 
TARGETS := double_pendulum

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
obj/main.o: src/main.c inc/integrator/rk2.h inc/integrator/integrator.h \
 inc/integrator/rk_classic.h inc/double_pendulum.h inc/pendulum_energy.h \
 inc/double_pendulum.h inc/gnuplot_i.h
obj/integrator/main.o: src/main.c inc/integrator/rk2.h inc/integrator/integrator.h \
 inc/integrator/rk_classic.h inc/double_pendulum.h inc/pendulum_energy.h \
 inc/double_pendulum.h inc/gnuplot_i.h
obj/pendulum_energy.o: src/pendulum_energy.c inc/pendulum_energy.h \
 inc/double_pendulum.h inc/double_pendulum.h
obj/integrator/pendulum_energy.o: src/pendulum_energy.c inc/pendulum_energy.h \
 inc/double_pendulum.h inc/double_pendulum.h
obj/double_pendulum.o: src/double_pendulum.c inc/double_pendulum.h
obj/integrator/double_pendulum.o: src/double_pendulum.c inc/double_pendulum.h
obj/dopri45.o: src/integrator/dopri45.c inc/integrator/dopri45.h \
 inc/integrator/integrator.h
obj/integrator/dopri45.o: src/integrator/dopri45.c inc/integrator/dopri45.h \
 inc/integrator/integrator.h
obj/rk2.o: src/integrator/rk2.c inc/integrator/rk2.h \
 inc/integrator/integrator.h
obj/integrator/rk2.o: src/integrator/rk2.c inc/integrator/rk2.h \
 inc/integrator/integrator.h
obj/integrator.o: src/integrator/integrator.c inc/integrator/integrator.h
obj/integrator/integrator.o: src/integrator/integrator.c inc/integrator/integrator.h
obj/rk_classic.o: src/integrator/rk_classic.c inc/integrator/rk_classic.h \
 inc/integrator/integrator.h
obj/integrator/rk_classic.o: src/integrator/rk_classic.c inc/integrator/rk_classic.h \
 inc/integrator/integrator.h
obj/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
obj/integrator/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
