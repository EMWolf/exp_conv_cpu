CC       = pgcc

CCFLAGS  = -lm -I../common

OMPFLAGS = -fast -mp -Minfo

SERIALFLAGS = 


BIN =  exp_conv_serial exp_conv_omp


all: $(BIN)


exp_conv_serial: exp_conv.c
	
	$(CC) $(CCFLAGS) $(SERIALFLAGS) -o $@ $<


exp_conv_omp: exp_conv_omp.c
	
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ $<
		


clean:

	$(RM) $(BIN)