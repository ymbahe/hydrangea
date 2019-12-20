EXEC   = ckat.so

OBJS   = ckat.o

#OPTIONS += -DPERIODIC

INCL   = 

CFLAGS =  -O2 $(OPTIONS) -g -fPIC -Wall #-fopenmp
LNKCMD =  ld -shared  -L/usr/lib -L/usr/local/lib #-l:libgomp.so.1


CC     =  cc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)
