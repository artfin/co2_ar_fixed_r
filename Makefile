.PHONY: clean all

all: main
	@echo Compilation finished

CCXX = mpic++

EIGEN := -I"/usr/local/include/eigen3/"
INCLUDEHEP := -I"/home/artfin/Downloads/hep-mc-0.5/include" 
LIBS := -lgfortran -lgsl -lgslcblas 

CXXFLAGS := -std=c++11 -g -O2 -lm -Wall -Wextra
BUILDDIR := ./build/

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
	$(CCXX) $(CXXFLAGS) $(LIBS) $(EIGEN) $(INCLUDEHEP) $^ -o $@

$(BUILDDIR)%.o: %.cpp
	@echo ">> (mpic++) Compiling $<...";
	@$(CCXX) $(CXXFLAGS) $(LIBS) $(EIGEN) $(INCLUDEHEP) -c -MD $< -o $@

include $(wildcard *.d)

clean:
	@rm -f main
	@rm -f $(BUILDDIR)*.o
	@rm -f $(BUILDDIR)*.d

