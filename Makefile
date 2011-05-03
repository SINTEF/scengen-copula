# project parameters
BASENAME = cop-gen_heur
OBJDIR = .obj
EXEDIR = bin
SRCS = $(wildcard *.cpp)
OBJS = $(addprefix $(OBJDIR)/, $(SRCS:.cpp=.o))
EXENAME = $(EXEDIR)/$(BASENAME)

# general compiler settings (do I need -fexceptions?)
CXX = g++
CXXFLAGS = -O3 -DNDEBUG -Wall -fexceptions
CXXLINKFLAGS = -s
INCL = 
LIBS = -lQuantLib

# Both boost, tclap, quantlib and glpk should be in the default include path,
# if they were installed properly. Since I don't have admin rights,
# I have installed them to my home, so I need to point to it here:
#INCL += -I/home/kaut/local/include
#CXXLINKFLAGS += -L /home/kaut/local/lib

# default target
all: $(EXENAME)

# rule for compiling:
$(OBJDIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<

# rule for building
$(EXENAME) : $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)


# ADDITIONAL RULES

clean:
	@rm -rf $(EXEDIR)/$(EXENAME)
	@rm -rf $(OBJDIR)/*.o
