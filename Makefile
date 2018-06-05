CXX 			= $(shell root-config --cxx)
LD 				= $(shell root-config --ld)

CPPFLAGS 	:= $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS 	:= -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

OBJ_DIR 	= src/obj

vpath %.cpp src
vpath %.h src
vpath %.o src/obj

objects 	= $(addprefix $(OBJ_DIR)/, pipi-amp.o venez-amp.o cgamma.o gauleg.o main.o partial-waves.o)
directories	 = ./src/obj ./output

$(objects) : 	$(OBJ_DIR)/%.o :  %.cpp veneziano.h
						g++ -c $< -o $@

$(directories) :
									@echo "Creating folder: $@" && \
    							mkdir -p $@

rho : 	$(directories)	$(objects)
				g++ -o rho $(objects)

.PHONY : clean clean-out

spotless :
				rm -rf ./output ./src/obj rho

clean 			:	clean-exe clean-out

clean-exe			:
						rm -rf rho ./src/obj

clean-out 	:
						rm -rf ./output
