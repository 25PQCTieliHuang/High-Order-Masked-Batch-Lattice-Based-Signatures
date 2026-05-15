#	Makefile

XBIN	?=	xtest
CC		?=	gcc -f
CFLAGS	+=	-Iinc $(RACCF)
CFLAGS	+=	-Wall -Wextra -O3 -march=native -fstack-usage


MODE	?=	batch

ifeq ($(MODE),batch)
    CFLAGS += -DBATCH_MODE
endif

# -----------------

CSRC	+= 	$(wildcard *.c util/*.c)
OBJS	= 	$(CSRC:.c=.o)
SUFILES	= 	$(CSRC:.c=.su)
LDLIBS	+=

$(XBIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(XBIN) $(OBJS) $(LDLIBS)

%.o:	%.c
	$(CC) $(CFLAGS) -c $< -o $@

obj-clean:
	$(RM) -f $(XBIN) $(OBJS) $(SUFILES) nist/*.o nist/*.su

clean:	obj-clean
	$(RM) -f bench_* xtest
	$(RM) -rf kat