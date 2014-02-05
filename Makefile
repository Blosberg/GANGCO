#IDIR =/usr/include/
CC=g++
CFLAGS  =  -I/sw/include -I/home/b/Brendan.Osberg/grad_research_phd/project/code_repository/ -ggdb
LDFLAGS =  -L/sw/lib -lgsl -lgslcblas -lm

DEPS = GA_GC_NTF.h /home/b/Brendan.Osberg/grad_research_phd/project/code_repository/GA_absolute_standards.h /home/b/Brendan.Osberg/grad_research_phd/project/code_repository/bren_lib.h /home/b/Brendan.Osberg/grad_research_phd/project/code_repository/GA_NTF_standards.h
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


