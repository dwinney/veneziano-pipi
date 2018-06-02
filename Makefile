CXX 			= $(shell root-config --cxx)
LD 				= $(shell root-config --ld)

CPPFLAGS 	:= $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS 	:= -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

OBJ_DIR 	= src/obj

vpath %.cpp src
vpath %.h src
vpath %.o src/obj

objects 	= $(addprefix $(OBJ_DIR)/, pipi-amp.o venez-amp.o cgamma.o gauleg.o main.o partial-waves.o)

$(objects) : $(OBJ_DIR)/%.o :  %.cpp veneziano.h
						g++ -c $< -o $@

rho : $(objects)
			g++ -o rho $(objects)

.PHONY : clean directories clean-out

directories :
						mkdir -p ./src/obj ./output/plots

clean 			: clean-exe clean-out

clean-exe			:
						rm -f rho $(OBJ_DIR)/*.o

clean-out 	:
						rm -f ./output/*.dat ./output/plots/*.pdf
