CFLAGS = -O3 -Wall -g 
COMPILER = clang 
LINK := -g -lm
INCLUDE_OPTIONS := -I./inc
TARGETS := rejection

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = inc






INCLUDES := $(shell find $(INC_DIR) -type d | sed s/^/-I/) $(INCLUDE_OPTIONS) 
SRC := $(shell find $(SRC_DIR) -name '*.c')
.PHONY: all database directories clean depend
all: depend directories $(TARGETS)
OBJ := $(subst $(SRC_DIR),$(OBJ_DIR),$(patsubst %.c,%.o,$(SRC)))
OBJDIRS := $(sort $(dir $(OBJ)))
SRCDIRS := $(sort $(dir $(SRC)))
$(TARGETS):	$(OBJ) 
	$(COMPILER) -o $@ $(OBJ) $(LINK)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c; $(COMPILER) $(CFLAGS) $(INCLUDES) -c $< -o $@
database: ;make clean; bear make
directories:; mkdir -p $(OBJDIRS)
clean: ; rm -f $(TARGETS) $(OBJ) 
depend: ; ./gcc-makedepend $(addprefix -p ,$(OBJDIRS)) $(CFLAGS) $(INCLUDES) $(SRC)
# DO NOT DELETE
obj/gnuplot_i.o: src/gnuplot_i.c inc/gnuplot_i.h
obj/distr.o: src/distr.c inc/distr.h
obj/monte_carlo.o: src/monte_carlo.c inc/monte_carlo.h inc/random_gen.h
obj/main.o: src/main.c inc/monte_carlo.h inc/distr.h inc/gnuplot_i.h
obj/random_gen.o: src/random_gen.c inc/random_gen.h
