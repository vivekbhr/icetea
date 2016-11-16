#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
//
#include "trimFastq.h"

// THIS PROGRAM WILL TAKE TWO FASTA FILES (R1 & R2) AND TRIM THE INDEX OFF R2 SEQUENCES TO ADD IT ON BOTH R1 AND R2 SEQ HEADERS

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
  gzFile fpt;
	gzFile fpout;
  gzFile fptout;

	kseq_t *seq;
  kseq_t *seqt;

	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in_R1.fastq.gz> <in_R2.fastq.gz> <out_R1.fastq> <out_R2.fastq>\n", argv[0]);
		return 1;
	}

	fp = gzopen(argv[2], "r"); // file pointer for R2
  fpt = gzopen(argv[1], "r"); // file pointer for R1

	fpout = gzopen(argv[4], "wb"); // output file pointer for R2
  fptout = gzopen(argv[3], "wb"); // output file pointer for R1

	seq = kseq_init(fp);
  seqt = kseq_init(fpt);
// save data for R2
	char * sequence;
  char * name;
  char * qual;
  char * comment;
// for R1
  char * tsequence;
  char * tname;
  char * tqual;
  char * tcomment;

  char header[100];
  char trimseq[50];
  char trimqual[50];

	while ((l = kseq_read(seq)) >= 0) {
    kseq_read(seqt);
		sequence = seq->seq.s;
    // other components
    name = seq->name.s;
    if (seq->qual.l) qual = seq->qual.s;//quality
    if (seq->comment.l) comment = seq->comment.s;//comment
    // for file R1
    tsequence = seqt->seq.s;
    tname = seqt->name.s;
    if (seqt->qual.l) tqual = seqt->qual.s;//quality
    if (seqt->comment.l) tcomment = seqt->comment.s;//comment

    //now print data for both R1 and R2

    make_header(header,sequence); // header only taken from R2
    gzprintf(fpout,"@%s#%s %s\n", name,header,comment);
    gzprintf(fptout,"@%s#%s %s\n", tname,header,tcomment);

    trim_seq(trimseq,sequence);
    gzprintf(fpout,"%s\n", trimseq);//trimmed seq for R2
    gzprintf(fptout,"%s\n", tsequence);//normal seq for R1

    gzprintf(fpout,"%s\n","+" );
    gzprintf(fptout,"%s\n","+" );

    trim_seq(trimqual,qual);
    gzprintf(fpout,"%s\n", trimqual);//trimmed qual for R2
    gzprintf(fptout,"%s\n", tqual);//normal qual for R1:
	}

	kseq_destroy(seq);
  kseq_destroy(seqt);

	gzclose(fp);
  gzclose(fpt);

  gzclose(fpout);
  gzclose(fptout);

  return 0;

}
