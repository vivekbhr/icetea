all:
	$(CC) -g -O2 trimFastq.c -o trimFastq -lz
	$(CC) -g -O2 trimFastq_rampage.c -o trimFastq_rampage -lz

clean:
	rm -f *.o trimFastq
	rm -f *.o trimFastq_rampage
