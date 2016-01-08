#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
//
#include "trimFETISH.h"

//copy the index sequence (position 6 to 12)
void copy_index(char d[], char s[]) {
   int strt = 0;
   int idx = 5;
   int bcstart = 11;

   while (idx < bcstart) {
      d[strt] = s[idx];
      idx++;
      strt++;
   }
   d[strt] = '\0';
}

//copy pcr barcode (pos 1 to 5 + pos 12 to 13)
void copy_barcode(char d[], char s[]) {
   int bcstrt = 0;
   int bcend = 5;
   int newbcstrt = 11;
   //copy first five bases
   while (bcstrt < bcend) {
      d[bcstrt] = s[bcstrt];
      bcstrt++;
   }
   //copy last two bases
   while (newbcstrt < 13) {
      d[bcstrt] = s[newbcstrt];
      bcstrt++;
      newbcstrt++;
   }
   d[bcstrt] = '\0';
}
//copy replicate demultiplexing barcode  (pos 3 to 4)
void copy_repcode(char d[], char s[]) {
   int idx = 0;
   int strt = 2;
   //copy first five bases
   while (strt < 4) {
      d[idx] = s[strt];
      idx++;
      strt++;
   }
   d[idx] = '\0';
}
//trim the seq off the barcode and index (13 bases removed)
void trim_seq(char d[], char s[]) {
   int idx = 13;
   int end = strlen(s);
   int i = 0;
   //copy first five bases
   while (idx < end) {
      d[i] = s[idx];
      i++;
      idx++;
   }
   d[i] = '\0';
}


void make_header(char * output, char s[]) {
  char index[50];
  char barcode[20];
  char repcode[5];
  output[0] = '\0';
  copy_index(index,s);
  copy_barcode(barcode,s);
  copy_repcode(repcode,s);
  strcat(output,index);
  strcat(output,":");
  strcat(output,barcode);
  strcat(output,":");
  strcat(output,repcode);
}

// put them together
KSEQ_INIT(gzFile, gzread);

int main(int argc, char *argv[]) //main_trimFETISH
{
	gzFile fp;
	FILE * fpout;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in_R2.fastq.gz> <out_R2.fastq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	fpout = gzopen(argv[2], "w+");

	seq = kseq_init(fp);
	char * sequence;
  char * name;
  char * qual;
  char * comment;
  char header[100];
  char trimseq[50];
  char trimqual[50];

	while ((l = kseq_read(seq)) >= 0) {
		sequence = seq->seq.s;
    // other components
    name = seq->name.s;
    if (seq->qual.l) qual = seq->qual.s;//quality
    if (seq->comment.l) comment = seq->comment.s;//comment
    //now print
    make_header(header,sequence);
    gzrintf(fpout,"@%s#%s %s\n", name,header,comment);
    trim_seq(trimseq,sequence);
    gzrintf(fpout,"%s\n", trimseq);//trimmed seq:
    gzrintf(fpout,"%s\n","+" );
    trim_seq(trimqual,qual);
    gzrintf(fpout,"%s\n", trimqual);//trimmed qual:

	}

	kseq_destroy(seq);
	gzclose(fp);
  gzclose(fpout);

  return 0;

}
