cuLU: src/main.cpp 
	/opt/cuda/bin/nvcc  -c -I. src/mmio.cpp
	/opt/cuda/bin/nvcc  -c -I. src/maxerror.cpp
	/opt/cuda/bin/nvcc  -c -I. src/cuLU.cpp
	/opt/cuda/bin/nvcc  -c -I. src/main.cpp
	/opt/cuda/bin/nvcc  cuLU.o mmio.o maxerror.o main.o -lcusolver -lm -o cuLU

clean: ;
	rm -f cuLU cuLU.o mmio.o maxerror.o main.o

