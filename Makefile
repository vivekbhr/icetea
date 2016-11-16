all:
	$(CC) -g -O2 trimFastq.c -o trimFastq -lz

clean:
	rm -f *.o trimFastq
