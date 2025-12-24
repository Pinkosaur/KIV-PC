CC = gcc
CFLAGS = -Wall -pedantic -ansi -Wextra -O2 -Iinclude
BIN = calc.exe

SRCDIR = src
BUILDDIR = build

SRCS = mp_int.c mp_print.c parser.c exp.c fact.c main.c error.c
OBJS = $(patsubst %.c,$(BUILDDIR)/%.o,$(SRCS))

.PHONY: all clean dirs

all: dirs $(BIN)

dirs:
	@mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(BIN)

clean:
	@ echo 'removing binaries...'
	@rm -rf $(BUILDDIR) $(BIN)
