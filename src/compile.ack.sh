cc="/home/marco/ACK/bin/ack"

$cc -o testpca_pc86 *.c  tests/testpca.c -I.
rm *.o

#ar rcs libscientific.a *.o
#gcc -c tests/testpca.c  -I.
#cc testpca.o -lscientific -lm -L. -o testpca


