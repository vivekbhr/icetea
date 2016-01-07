#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string.h>
//
#include <Rversion.h>
#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif
#include "trimFETISH.h"
#include <R.h>

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

// put them together (the MAIN function)
KSEQ_INIT(gzFile, gzread);

void trimFETISH(char **infile, char **outfile)
{
  kseq_t *seq;
	int l;
  gzFile fp;
  FILE * fpout;
  fp = gzopen(*infile, "r");
	fpout = fopen(*outfile, "w+");
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
    fprintf(fpout,"@%s#%s %s\n", name,header,comment);
    trim_seq(trimseq,sequence);
    fprintf(fpout,"%s\n", trimseq);//trimmed seq:
    fprintf(fpout,"%s\n","+" );
    trim_seq(trimqual,qual);
    fprintf(fpout,"%s\n", trimqual);//trimmed qual:

	}

	kseq_destroy(seq);
	gzclose(fp);
  fclose(fpout);

}
