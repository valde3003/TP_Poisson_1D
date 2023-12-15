CC = gcc
CFLAGS = -I include -Wall
SRCDIR = src
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = $(BINDIR)/run_target

.PHONY: all build_dirs clean

all: build_dirs $(EXECUTABLE)

build_dirs:
	@mkdir -p $(BINDIR)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(EXECUTABLE)

include ambre.mk
