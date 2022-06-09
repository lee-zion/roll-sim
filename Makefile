algorithm.o: src/algorithm.h src/algorithm.c
		g++ -c src/algorithm.c

main.o: src/algorithm.h src/main.c
		g++ -c src/main.c

main: algorithm.o main.o
		g++ algorithm.o main.o -o main

.PHONY: clean
clean:
		rm -f algorithm.o main.o main