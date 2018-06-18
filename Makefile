CXX 			= $(shell root-config --cxx)
LD 				= $(shell root-config --ld)

CPPFLAGS 	:= $(shell root-config --libs --cflags) $(STDINCDIR)
LDFLAGS 	:= -Xlinker -rpath . $(shell root-config --glibs) -lMinuit $(STDLIBDIR)

OBJ_DIR 	= src/obj

vpath %.cpp src
vpath %.h src
vpath %.o src/obj

ROOT_objects 		=	$(addprefix $(OBJ_DIR)/, plotting.o fitting.o)
objects 	= $(addprefix $(OBJ_DIR)/, main.o VENEZ.o math-misc.o GKPY.o)
directories	 = ./src/obj ./output

pipi		 :	$(directories)	$(objects) $(ROOT_objects)
					$(LD) -o pipi $(objects) $(ROOT_objects) $(LDFLAGS)

$(OBJ_DIR)/pipi.h.gch : pipi.h
						g++ -o $@ -c $<

$(objects)	  : 	$(OBJ_DIR)/%.o :  %.cpp $(OBJ_DIR)/pipi.h.gch
						g++ -c $< -o $@

$(directories) :
						@echo "Creating folder: $@" && \
    				mkdir -p $@

$(ROOT_objects)	: $(OBJ_DIR)/%.o : %.cpp $(OBJ_DIR)/pipi.h.gch
							$(CXX) $(CPPFLAGS) -o $@ -c $<


.PHONY		 : clean spotless

spotless	 :
						rm -rf ./output ./src/obj fit pipi

clean 		 :
						rm -rf pipi ./src/obj
