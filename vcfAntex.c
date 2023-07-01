// Adds TNC annotation to the info field if a VCF SBS record.


#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/faidx.h>
#include "hts_internal.h"


inline void revc(char *base, char *comp)
{
   if(*base =='A')
      *comp = 'T';
   else if(*base =='T')
      *comp = 'A';
   else if(*base =='G')
      *comp = 'C';
   else if(*base =='C')
      *comp = 'G';
}

int main(int argc,char* argv[])
{

   FILE *fp = fopen(argv[1], "r");
   if(fp == NULL) {
      fprintf(stderr, "\nCan't open VCF file %s.\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }


   faidx_t *fai = fai_load(argv[2]);  // Reference sequence
   if(fai==NULL) {
      fprintf(stderr,"Can't load faidx: %s\n", argv[2]);
      exit(EXIT_FAILURE);
   }


   // First read the SBS signature file and work out how much burden we need
   // to spike in to each context. Setup targetCtxBurdenArray which records this.
   char *line = NULL;
   size_t len = 0;
   long locus = 0;
   ssize_t characters;

   while((characters = getline(&line, &len, fp)) != -1) {
      char *pch, *contigPtr = NULL, *refAllelePtr, *altAllelePtr = NULL;
      int field=0;

      // Skip blank lines.
      if(characters == 0)
         continue;

      // Skip blank lines.
      if(line[0] == '\n') {
         printf("%s",line);
         continue;
      }

      // Skip metadata/comments.
      if(line[0] == '#') {
         printf("%s",line);
         continue;
      }

      pch = strtok(line,"\t");
      if(pch != NULL) {
         field++;
         // Strip newline if its there...
         pch[strcspn(pch, "\n")] = 0;

         // "CHROM" field in the VCF record.
         // We're at the first field in the VCF record which is the
         // contig name.
         contigPtr=pch;
         printf("%s",pch);
      }

      while((pch = strtok(NULL, "\t")) != NULL) {
         //  pch = strtok(NULL, "\t");
         field++;

         // Strip newline if its there...
         pch[strcspn(pch, "\n")] = 0;
         printf("\t%s",pch);

         // "POS" field from the VCF record.
         if(field==2) {
            char* strEnd;
            // Remember locus is 1 based.
            locus = strtol(pch, &strEnd, 10)-1;
            if(pch == strEnd) {
               fprintf(stderr, "Error: Unrecognised SBS locus, \"%s\".\n", pch);
               exit(1);
            }
         }


         // "REF" field from the VCF record.
         if(field==4) {
            refAllelePtr=pch;

         }

         // "ALT" field from the VCF record.
         if(field==5) {
            altAllelePtr=pch;

         }

         // "INFO" field from the VCF record.
         // This is where we add our annotation.
         if(field==8) {
            // Append out TNC annotation to the end of this field.

            //#########################


            // Only print context for SBSs..
            // TODO!!!! We can add others later if anyone is interested in them..
            if(strlen(refAllelePtr) != 1  || strlen(altAllelePtr) != 1) {
               continue;
            }

            hts_pos_t faiRefLen;
            char *faiRef = faidx_fetch_seq64(fai, contigPtr, locus-1, locus+1, &faiRefLen);
            if(faiRefLen < 3) {
               fprintf(stderr,"Failed to fetch the sequence, trying next one..\n"); // no worries, move onto next one if so...
               continue;
            }


            if(!altAllelePtr) {
               fprintf(stderr,"Format error, %s..\n",argv[1]);
               exit(EXIT_FAILURE);
            }

            // Rem. this is strand agnostic..
            // Always print out context in the standard format.
            // Middle ref. base must be a 'C' or a 'T'

            if(faiRef[1]!='C' && faiRef[1]!='T') {
               char revcTnc[4] = {'\0','\0','\0','\0'};
               char revcAlt = '\0';

               // OK we need to print this in reverse complement.
               revc(&faiRef[2],&revcTnc[0]);
               revc(&faiRef[1],&revcTnc[1]);
               revc(&faiRef[0],&revcTnc[2]);
               revc(altAllelePtr,&revcAlt);
               printf(";TNC=%c[%c>%c]%c",revcTnc[0],revcTnc[1],revcAlt,revcTnc[2]);
            }
            else {
               printf(";TNC=%c[%c>%s]%c",faiRef[0],faiRef[1],altAllelePtr,faiRef[2]);
            }



            //##################

         }
      }
      printf("\n");



   }
}
