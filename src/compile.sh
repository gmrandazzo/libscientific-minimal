gcc -c *.c
ar rcs libscientific.a *.o
gcc -c tests/testpca.c  -I.
gcc testpca.o -lscientific -lm -L. -o testpca


