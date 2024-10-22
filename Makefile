all:
	g++ hits9171.cpp -o main.o && ./main.o 7 1 input.txt

clean:
	rm -rf *.o