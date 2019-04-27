# make macros: $@ = file name of target, $< = name of first dependency, $^ = names of all dependencies (separated by spaces, no duplicates)
# the include and library paths below are intended to accomodate compilation on OSX


CC = g++
code_repo_path=./code_repository/

CFLAGS  = -I${code_repo_path} -I/usr/local/Homebrew/Cellar/gsl/2.5/include/  -ggdb -O2
LDFLAGS = -L/usr/local/Homebrew/Cellar/gsl/2.5/lib/  -lgsl -lgslcblas -lm

DEPS = GA_GC_N.h ${code_repo_path}GA_absolute_standards.h ${code_repo_path}bren_lib.h
OBJ  = GA_GC_N.obj GA_GC_N_progs.obj

%.obj: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

GA_GC_N.x: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

ALL_OBJS = GA_GC_N.obj GA_GC_N_progs.obj

clean:
	$(RM) $(ALL_OBJS)

debug:
	g++ -ggdb GA_GC_N.cpp GA_GC_N_progs.cpp -o GA_GC_N.x $(CFLAGS) $(LDFLAGS)

voidcorr:
	g++ -ggdb GA_voidsize_corr.cpp -o GA_voidsize_corr.x $(CFLAGS) $(LDFLAGS)

