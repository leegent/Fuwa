fuwa:	fuwa.cpp
	g++ fuwa.cpp -o fuwa -O3 -Wall -lz -Isamtools-0.1.19 -lbam -Lsamtools-0.1.19 -fopenmp -static
