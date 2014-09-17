#IDIR =/usr/include/
CC=g++
code_repo_path=/home/t30/ger/ga79moz/grad_research_phd/project/code_repository/

CFLAGS  =  -I/sw/include -I/home/t30/ger/ga79moz/grad_research_phd/project/code_repository/ -ggdb
LDFLAGS =  -L/sw/lib -lgsl -lgslcblas -lm

DEPS = GA_GC_NTF.h ${code_repo_path}GA_absolute_standards.h ${code_repo_path}bren_lib.h ${code_repo_path}GA_NTF_standards.h
OBJ  = GA_GC_NTF.obj  GA_GC_NTF_progs.obj

%.obj: %.cpp $(DEPS)
	$(CC) -c -o  $@ $< $(CFLAGS) 

GA_GC_NTF.x: $(OBJ)
	$(CC) -ggdb -g -o  $@ $^ $(LDFLAGS) 

ALL_OBJS = GA_GC_NTF.obj   GA_GC_NTF_progs.obj

clean:
	$(RM) $(ALL_OBJS)

debug:
	g++ -ggdb GA_GC_NTF.cpp  GA_GC_NTF_progs.cpp -o GA_GC_NTF.x $(CFLAGS) $(LDFLAGS) 



