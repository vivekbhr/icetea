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
#include "trimFastq.h"
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

void trimFastq(char **inR1, char **inR2, char **outR1, char **outR2)
{
  gzFile fp;
  gzFile fpt;
  gzFile fpout;
  gzFile fptout;

  kseq_t *seq;
  kseq_t *seqt;

  fp = gzopen(*inR2, "r"); // file pointer for R2
  fpt = gzopen(*inR1, "r"); // file pointer for R1

  fpout = gzopen(*outR2, "wb"); // output file pointer for R2
  fptout = gzopen(*outR1, "wb"); // output file pointer for R1

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

  int l;
while ((l = kseq_read(seq)) >= 0) {
    kseq_read(seqt);
    // save data for R2
    sequence = seq->seq.s;
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

// close the files
  kseq_destroy(seq);
  kseq_destroy(seqt);

  gzclose(fp);
  gzclose(fpt);

  gzclose(fpout);
  gzclose(fptout);


}
