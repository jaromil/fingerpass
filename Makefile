all:
	gcc -o fingerpass fingerpass.c -lfprint -lm -O2

debug:
	gcc -o fingerpass fingerpass.c -lfprint -lm -ggdb -O0 -DDEBUG -Wall
