all:
	$(CC) -g -O2 trimFETISH.c -o trimFETISH -lz

clean:
	rm -f *.o trimFETISH
