
CC = gcc

CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
          -Wshadow -mtune=native -O3

HEADERS = NTT_params.h tools.h gen_table.h naive_mult.h ntt_c.h
SRCs = tools.c gen_table.c naive_mult.c ntt_c.c

all: $(HEADERS) $(SRCs) main.c
	$(CC) -o main $(CFLAGS) $(SRCs) main.c

.PHONY: clean
clean:
	rm -f main

