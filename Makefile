default:
	gcc -fPIC -Wall -c -g tree.c
	gcc -shared -o tree.so tree.o