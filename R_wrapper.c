#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <Rversion.h>
#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif
#include "trimFETISH.h"
int main_trimFETISH(int argc, char *argv[]);

void R_trimFETISH(char infile, char outfile)
{
  n = 2;
  c_argv = (char **) calloc(n,sizeof(char *));
  for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
  strcpy(c_argv[0],strtok(r_argv,","));
  for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

  main_trimFETISH(c_argv);

  for(i=0;i<n;i++) free(c_argv[i]);
  free(c_argv);
}
