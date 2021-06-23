# GRR20197152 Guilherme Costa Pateiro
# Universidade Federal do Parana , Introducao a computacao cientifica
# ultima edicao: 18/05/2021 18:44
CFLAGS = -Wall -lm

all:
	gcc floattest.c -o floattest $(CFLAGS) 

run: all
	./testes/run.sh

clean:
	rm -f *.out
	rm -f resp.txt 

purge:	clean
	rm -f floattest