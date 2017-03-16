#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
//
#include "trimFastq_rampage.h"

// THIS PROGRAM WILL TAKE TWO FASTA FILES (R1 & R2) AND TRIM THE INDEX OFF R1 and R2 SEQUENCES TO ADD IT ON BOTH R1 AND R2 SEQ HEADERS

//copy the sequence barcode (position 1 to pos). pos = 15 for R2 and pos = 6 for R1
void copy_barcode(char d[], char s[], int pos) {
   int strt = 0;
   int idx = 0;

   while (strt < pos) {
      d[strt] = s[idx];
      idx++;
      strt++;
   }
   d[strt] = '\0';
}

//trim the seq off the barcode and index (given no of bases removed)
void trim_seq(char * trimmed, char s[], int pos) { // idx = 6 for R1 and 15 for R2
   int end = strlen(s);
   //int idx = 6;
   int i = 0;
   //copy given no of bases
   while (pos < end) {
      trimmed[i] = s[pos];
      i++;
      pos++;
   }
   trimmed[i] = '\0';
}

void make_header(char * output, char s[], char readtype[]) {
  output[0] = '\0';
  if (readtype == "r1") {
      char index[7]; // to save the index
      copy_barcode(index,s, 6);// copy pos 1 to 6
      strcat(output,index);
  }
  else {
      char barcode[16]; // to save the RTPCR barcode
      copy_barcode(barcode,s, 15);//copy pos 1 to 15
      strcat(output,barcode);
  }

}

// put them together
KSEQ_INIT(gzFile, gzread);

int main(int argc, char *argv[]) //main_trimFETISH
{
	gzFile fp_r2;
  gzFile fp_r1;
	gzFile fpout_r2;
  gzFile fpout_r1;

	kseq_t *seq_r2;
  kseq_t *seq_r1;


	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in_R1.fastq.gz> <in_R2.fastq.gz> <out_R1.fastq> <out_R2.fastq>\n", argv[0]);
		return 1;
	}

  fp_r2 = gzopen(argv[2], "r"); // file pointer for R2
  fp_r1 = gzopen(argv[1], "r"); // file pointer for R1

	fpout_r2 = gzopen(argv[4], "wb"); // output file pointer for R2
  fpout_r1 = gzopen(argv[3], "wb"); // output file pointer for R1

	seq_r2 = kseq_init(fp_r2);
  seq_r1 = kseq_init(fp_r1);
// save data for R2
  char * r2_sequence;
  char * r2_name;
  char * r2_qual;
  char * r2_comment;
// for R1
  char * r1_sequence;
  char * r1_name;
  char * r1_qual;
  char * r1_comment;

// to save the trimmed header
  char r1_header[100];
  char r2_header[100];

  char r1_trimseq[200];
  char r2_trimseq[200];

  char r1_trimqual[200];
  char r2_trimqual[200];

  int l;
	while ((l = kseq_read(seq_r2)) >= 0) {
    kseq_read(seq_r1);
	   r2_sequence = seq_r2->seq.s;
    // other components
    r2_name = seq_r2->name.s;
    if (seq_r2->qual.l) r2_qual = seq_r2->qual.s;//quality
    if (seq_r2->comment.l) r2_comment = seq_r2->comment.s;//comment

    // for file R1
    r1_sequence = seq_r1->seq.s;
    r1_name = seq_r1->name.s;
    if (seq_r1->qual.l) r1_qual = seq_r1->qual.s;//quality
    if (seq_r1->comment.l) r1_comment = seq_r1->comment.s;//comment

    //now print data for both R1 and R2

    make_header(r1_header,r1_sequence,"r1"); // header taken from R1
    make_header(r2_header,r2_sequence,"r2"); // header taken from R2

    gzprintf(fpout_r2,"@%s#%s %s\n", r2_name,r2_header,r2_comment);//doesn't work with GEO output
    gzprintf(fpout_r1,"@%s#%s %s\n", r1_name,r1_header,r1_comment);

    trim_seq(r2_trimseq,r2_sequence, 15);
    trim_seq(r1_trimseq,r1_sequence, 6);

    gzprintf(fpout_r2,"%s\n", r2_trimseq);//trimmed seq for R2
    gzprintf(fpout_r1,"%s\n", r1_trimseq);//trimmed seq for R1

    gzprintf(fpout_r2,"%s\n","+" );
    gzprintf(fpout_r1,"%s\n","+" );

    trim_seq(r2_trimqual,r2_qual, 15);
    trim_seq(r1_trimqual,r1_qual, 6);

    gzprintf(fpout_r2,"%s\n", r2_trimqual);//trimmed qual for R2
    gzprintf(fpout_r1,"%s\n", r1_trimqual);//normal qual for R1:
	}

	kseq_destroy(seq_r2);
  kseq_destroy(seq_r1);

	gzclose(fp_r2);
  gzclose(fp_r1);

  gzclose(fpout_r2);
  gzclose(fpout_r1);

  return 0;

}
