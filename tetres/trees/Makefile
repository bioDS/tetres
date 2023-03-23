CC=gcc
default:
	$(CC) -fPIC -Wall -c -g findpath.c -fopenmp
	$(CC) -shared -o findpath.so findpath.o -lgomp