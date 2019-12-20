EXEC   = HsmlAndProjectYB_Tau.so

OBJS   = main.o tree.o peano.o allvars.o

#OPTIONS += -DPERIODIC


INCL   = proto.h allvars.h

CFLAGS =  -O3 $(OPTIONS) -g -fPIC -fopenmp -Wall 
LNKCMD =  ld -shared -L/usr/lib -L/usr/local/lib -lgomp


#SYSTYPE="Regatta"


ifeq ($(SYSTYPE),"Regatta")
CFLAGS =  -O3   $(OPTIONS)
LNKCMD =   gcc   -bM:SRE -bnoentry -bE:export.txt
endif




#LIBS   =  -lm

CC     =  gcc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)
