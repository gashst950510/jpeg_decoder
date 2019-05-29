all: clean main
main: main.c 
	gcc -O2 -std=c99 main.c -o main -lm
clean:
	rm -f main
