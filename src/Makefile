PROG = solver
SRCS =  mmio.cpp \
	maxerror.cpp \
	solver.cpp 

OBJS = $(SRCS:.cpp=.o)

.SUFFIXES: .cpp .o

LDFLAGS = -lcusolver -lm

all: $(PROG)
$(PROG) : $(OBJS)
	/opt/cuda/bin/nvcc main.cpp $(OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	/opt/cuda/bin/nvcc -c -I. $<

clean: ;
	rm -f $(PROG) $(OBJS)

test: $(PROG)
	./test.sh

