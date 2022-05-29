CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c99 -pedantic -lcrypto

OBJS = crypto_utils.o miller_rabin.o baillie_psw.o tests_bpsw.o main.o

BIN = baillie_psw

all: $(OBJS)
	$(CC) $^ -o $(BIN) $(CFLAGS)


clean:
	$(RM) $(OBJS) $(BIN) *~

.PHONY: clean
