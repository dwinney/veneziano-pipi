CXX 			= $(shell root-config --cxx)
LD 				= $(shell root-config --ld)

CPPFLAGS 	:= $(shell root-config --libs --cflags) $(STDINCDIR)
LDFLAGS 	:= -Xlinker -rpath . $(shell root-config --glibs) -lMinuit $(STDLIBDIR)

OBJ_DIR 	= src/obj

vpath %.cpp src
vpath %.h src
vpath %.o src/obj

mains 		=	$(addprefix $(OBJ_DIR)/, fitting.o)
objects 	= $(addprefix $(OBJ_DIR)/, main.o pipi-amp.o venez-amp.o cgamma.o gauleg.o GKPRY.o)
directories	 = ./src/obj ./output

pipi		 :	$(directories)	$(objects) $(OBJ_DIR)/plotting.o
					$(LD) -o pipi $(objects) $(OBJ_DIR)/plotting.o $(LDFLAGS)

$(OBJ_DIR)/pipi.h.gch : pipi.h
						g++ -o $@ -c $<

$(OBJ_DIR)/plotting.o : plotting.cpp $(OBJ_DIR)/pipi.gch
						$(CXX) $(CPPFLAGS) -o $@ -c $<

$(objects)	  : 	$(OBJ_DIR)/%.o :  %.cpp $(OBJ_DIR)/pipi.h.gch
						g++ -c $< -o $@

$(directories) :
						@echo "Creating folder: $@" && \
    				mkdir -p $@

$(mains)	: $(OBJ_DIR)/%.o : %.cpp $(OBJ_DIR)/pipi.h.gch
							$(CXX) $(CPPFLAGS) -o $@ -c $<


fit				: $(directories)	$(objects) $(OBJ_DIR)/fitting.o
						$(LD) -o fit $(objects) $(OBJ_DIR)/fitting.o $(LDFLAGS)


.PHONY		 : clean clean-out spotless clean-exe

spotless	 :
						rm -rf ./output ./src/obj fit pipi

clean 		 :
						clean-exe clean-out

clean-exe	 :
						rm -rf fit pipi fit ./src/obj

clean-out	 :
						rm -rf ./output
