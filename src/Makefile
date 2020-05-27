PROG = solver
CC   = /opt/cuda/bin/nvcc
SRCS =  mmio.cpp \
	solver.cpp 

OBJS = $(SRCS:.cpp=.o)

.SUFFIXES: .cpp .o

LDFLAGS = -lcusolver -lm

all: $(PROG)
$(PROG) : $(OBJS)
	$(CC) -g main.cpp $(OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -g -c -I. $<

clean: ;
	rm -f $(PROG) $(OBJS)

test: $(PROG)
	../test/test.sh

