CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c99 -pedantic
LDFLAGS = -lcrypto

CDEBUG = -fsanitize=address
LDEBUG = -lasan

CRYPTO_OBJS = crypto_utils.o miller_rabin.o baillie_psw.o tests_bpsw.o
OBJS = main.o

BIN = baillie_psw

all: $(CRYPTO_OBJS) $(OBJS)
	$(CC) $^ -o $(BIN) $(CFLAGS) $(LDFLAGS)

debug: $(CRYPTO_OBJS) $(OBJS)
	$(CC) $^ -o $(BIN) $(CFLAGS) $(CDEBUG) $(LDFLAGS) $(LDEBUG)

clean:
	$(RM) $(CRYPTO_OBJS) $(OBJS) $(BIN) *~

.PHONY: clean
