CXX 			= $(shell root-config --cxx)
LD 				= $(shell root-config --ld)

CPPFLAGS 	:= $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS 	:= -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

OBJ_DIR 	= src/obj

vpath %.cpp src
vpath %.h src
vpath %.o src/obj

mains 		=	$(addprefix $(OBJ_DIR)/, main.o phase-shift.o)
objects 	= $(addprefix $(OBJ_DIR)/, pipi-amp.o venez-amp.o cgamma.o gauleg.o partial-waves.o)
directories	 = ./src/obj ./output

$(objects)	  : 	$(OBJ_DIR)/%.o :  %.cpp veneziano.h
						g++ -c $< -o $@

$(mains) 			:	$(OBJ_DIR)/%.o	 : 	%.cpp veneziano.h $(objects)
						g++ -c $< -o $@

$(directories) :
						@echo "Creating folder: $@" && \
    				mkdir -p $@

rho 					 : 	$(directories)	$(objects) $(OBJ_DIR)/main.o
				g++ -o rho $(objects) $(OBJ_DIR)/main.o

phase-shift		 :	$(directories)	$(objects) $(OBJ_DIR)/phase-shift.o
				g++ -o phase-shift $(objects) $(OBJ_DIR)/phase-shift.o


.PHONY		 : clean clean-out spotless

spotless	 :
						rm -rf ./output ./src/obj rho

clean 		 :
						clean-exe clean-out

clean-exe	 :
						rm -rf rho phase-shift ./src/obj

clean-out	 :
						rm -rf ./output
