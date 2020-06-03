CPP = gcc
OBJDIR = ./obj
SRCDIR = ./src
INCL_DIR = ./include

OBJ = 	edge.o\
	node.o\
	rnet.o\
	paramfitter.o\
	cvode_call.o\
	Matrix.o\
	Vector.o

MEXFILES = evalclient.cpp\
	Vector.cpp

EVALFILES = evals

CVODE_INCL_DIR = ./include/cvode
CVODE_SRCDIR =	./src/cvode
CVODE_OBJDIR =	./obj/cvode
CVODE_OBJ = 	band.o\
		cvbandpre.o\
		cvband.o\
		cvdense.o\
		cvdiag.o\
		cvode.o\
		cvspgmr.o\
		dense.o\
		iterativ.o\
		llnlmath.o\
		nvector.o\
		spgmr.o

## update the following line for your system
MEX_INCL_DIR = /usr/local/MatlabR2009b/extern/include


COMP_FLAGS = -Wall -O3 -g


CFLAGS = $(COMP_FLAGS) -I$(INCL_DIR) -I$(CVODE_INCL_DIR) -DCHECK_RESULTS -DCVODE_INTEGRATION -DGATED_GROWTH -DADAPTIVE_GG -DOPTIMIZE_NODE_PARAM


LIBS =  -lm -lstdc++

LDFLAGS = 

evalserver:	$(OBJDIR)/evalserver.o $(addprefix $(OBJDIR)/,$(OBJ)) $(addprefix $(CVODE_OBJDIR)/,$(CVODE_OBJ))
	$(CPP) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^

checkres:	$(OBJDIR)/checkres.o $(addprefix $(OBJDIR)/,$(OBJ)) $(addprefix $(CVODE_OBJDIR)/,$(CVODE_OBJ))
	$(CPP) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^

getDimension: $(SRCDIR)/getDimension.cpp $(addprefix $(SRCDIR)/,$(MEXFILES))
	mex $(CFLAGS) -I$(MEX_INCL_DIR) -I$(INCL_DIR) $^ -o $@	

requestEval: $(SRCDIR)/requestEval.cpp $(addprefix $(SRCDIR)/,$(MEXFILES))
	mex $(CFLAGS) -I$(MEX_INCL_DIR) -I$(INCL_DIR) $^ -o $@	

$(CVODE_OBJDIR)/%.o:	$(CVODE_SRCDIR)/%.c $(CVODE_INCL_DIR)/%.h
	$(CPP) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp $(INCL_DIR)/%.h $(INCL_DIR)/macros.h
	$(CPP) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp $(INCL_DIR)/macros.h
	$(CPP) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(CVODE_OBJDIR)/*.o
