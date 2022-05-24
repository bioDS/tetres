CC=gcc
default:
	$(CC) -fPIC -Wall -c -g geweke_parallel.c -fopenmp
	$(CC) -shared -o geweke_parallel.so geweke_parallel.o -lgomp