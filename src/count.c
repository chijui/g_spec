#include "count.h"

/***********************
* Counting subroutines *
************************/

int count_lines(char* fname) {
  int lines=0, length=0;
  char temp, *c;
  FILE* fp = fopen(fname, "r");
  if(fp==NULL) {
    printf("Error opening file %s: please check input", fname);
    return -1;
  }
  do {
    temp = fgetc(fp);
    if( temp == EOF ) return 0;
    length++;
  } while( temp!='\n' );
  lines++;
  c= (char *) malloc(2*length);    

  while(fgets(c,2*length,fp) != NULL){
    lines++;
  }
  fclose(fp);
  free(c);
  return lines;
}

// MER 04/06/2016
// Modified string tokenizer to include additional tokens beyond '\t'. 
// Also set num to start at zero rather than -1. Previously files whose
// lines did not include a terminating \t character were undercounted. 
int count_entries(char* fname) {
  int num = 0, length=0;
  char *c, temp, *ch,*d="\t, \n";
  FILE* fp = fopen(fname, "r");
  if(fp==NULL) {
    printf("Error openting file %s: please check input", fname);
    return -1;
  }
  do{
    temp = fgetc(fp);
    if( temp == EOF ) return 0;
    length++;
  } while( temp!='\n' );
  rewind(fp);
  c = (char *) malloc (2*length);
  fgets(c,2*length,fp);
  ch=strtok(c,d);
  while(ch){
    ch=strtok(NULL,d); 
    num++;
  }
  fclose(fp);
  free(c);
  return num; 
}

// cjfeng 07/05/2017
int count_fr(int *fr, int frame, int readframe, int nread, int twin, int nbuffer, int freq) {
  // Step through all recently read frames to see whether frame
  // is suitable as the end point of a calculation, or frame is
  // good to dump the currenct spectra. 
  // The "frequency" of starting a new realization or dumping 
  // spectra is given by freq.
  // The time window twin is either window or win2d, depending 
  // on whether do2d = 0 or 1.
  
  for(frame; frame<readframe; frame++) {
    // We check if at least window or win2d frames have been read
    // and if frame is a multiple of freq.
    if( ((frame % freq)==0) && (frame >= twin - 1) ) {
      *fr = (frame-twin+1) % nbuffer;
      break;
    }
  }
  if(frame == readframe) *fr = -2;
  return 0;
}
