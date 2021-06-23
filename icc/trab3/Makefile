# GRR20197152 Guilherme Costa Pateiro
# Universidade Federal do Parana , Introducao a computacao cientifica
# ultima edicao: 11/06/2021 02:20
    CC     = gcc -g -std=c11
    CFLAGS = -Wall
    LFLAGS = -lm

      PROG = EDOP 
      OBJS = utils.o  SistemasLineares.o Myfunc.o

.PHONY: limpa faxina clean purge all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

