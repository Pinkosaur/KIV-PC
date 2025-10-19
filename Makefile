CC = gcc
CFLAGS = -Wall -pedantic -ansi -Wextra -O2
BIN = calc.exe
OBJ = calc.o

%.o: src/%.c
	$(CC) -c $(CFLAGS) $< -o build/$@

$(BIN): $(OBJ)
	$(CC) build/$^ -o $@
