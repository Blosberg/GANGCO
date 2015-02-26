#IDIR =/usr/include/
CC=g++
code_repo_path=/home/t30/ger/ga79moz/grad_research_phd/project/code_repository/

CFLAGS  =  -I/sw/include -I${code_repo_path} -ggdb
LDFLAGS =  -L/sw/lib -lgsl -lgslcblas -lm

DEPS = GA_GC_N.h ${code_repo_path}GA_absolute_standards.h ${code_repo_path}bren_lib.h 
OBJ  = GA_GC_N.obj  GA_GC_N_progs.obj

%.obj: %.cpp $(DEPS)
	$(CC) -c -o  $@ $< $(CFLAGS) 

GA_GC_N.x: $(OBJ)
	$(CC) -ggdb -g -o  $@ $^ $(LDFLAGS) 

ALL_OBJS = GA_GC_N.obj   GA_GC_N_progs.obj

clean:
	$(RM) $(ALL_OBJS)

debug:
	g++ -ggdb GA_GC_N.cpp  GA_GC_N_progs.cpp -o GA_GC_N.x $(CFLAGS) $(LDFLAGS) 

voidcorr:
	g++ -ggdb GA_voidsize_corr.cpp -o GA_voidsize_corr.x $(CFLAGS) $(LDFLAGS) 
