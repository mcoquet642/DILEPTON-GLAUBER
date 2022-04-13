all:
	gcc -Wall -Wno-stringop-truncation -std=c++11  src/Main.cpp -o Main.exe -I/usr/include -L/usr/lib -fopenmp -lstdc++ -lm -ldl -lgsl -lgslcblas
