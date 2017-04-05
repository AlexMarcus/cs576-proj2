all: prog2
prog2: prog2.cu
	nvcc prog2.cu -o prog2
clean:
	rm -f *~ prog2
