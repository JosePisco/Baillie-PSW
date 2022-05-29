CC = gcc
CFLAGS = -lm -Wall -Wextra -Werror -std=c99 -pedantic

OBJS = crypto_utils.o miller_rabin.o baillie_psw.o

BIN = baillie_psw

all: $(OBJS)
	$(CC) $^ -o $(BIN) $(CFLAGS)


clean:
	$(RM) $(OBJS) $(BIN) *~

.PHONY: clean
