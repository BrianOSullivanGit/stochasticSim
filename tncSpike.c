// This program places the specified target burden at loci that match the required SBS signature profile
// (as specified by the profile file, arg 2).
// It takes as its input a random (or otherwise, upto you) set of loci from the stdin into which it spikes the required profile.
// Note the total of burden fractions across all contexts (ie., the sum of column 2) must equal 1.


#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/faidx.h>
#include "hts_internal.h"

#define BUF_SIZE 1024
#define NUM_TNC_TYPES 96

typedef struct {
   char context[4]; // The target base before substitution plus one adjacent base, either side, forward genomic strand.
   char substitution; // The base that will replace the target base, forward genomic strand.
   char contextPrime[4]; // The context complement, reverse genomic strand.
   char substitutionPrime; // The substitution complement, reverse genomic strand.
   long targetBurden; // Target burden, ie., how many of these substitutions you require.
} TargetBurdenElement;



int main(int argc,char* argv[])
{
   if(!(argc >= 3))
      fprintf(stderr, "\nMissing SBS signature profile (second arg).\n");

   long totalBurden=atol(argv[1]);

   FILE *fp = fopen(argv[2], "r");
   if(fp == NULL) {
      fprintf(stderr, "\nCan't open SBS signature file %s.\n\n",argv[2]);
      exit(EXIT_FAILURE);
   }


   faidx_t *fai = fai_load(argv[3]);  // Reference sequence
   if(fai==NULL) {
      fprintf(stderr,"Can't load faidx: %s\n", argv[3]);
      exit(-1);
   }

   FILE *fwdFp = fopen("fwd.targetLoci.txt", "w");
   if(fp == NULL) {
      fprintf(stderr, "\nCan't open forward strand orientation spike-in output file \"fwd.spike\".\n\n");
      exit(EXIT_FAILURE);
   }

   FILE *revFp = fopen("rev.targetLoci.txt", "w");
   if(fp == NULL) {
      fprintf(stderr, "\nCan't open reverse strand orientation spike-in output file \"rev.spike\".\n\n");
      exit(EXIT_FAILURE);
   }


   TargetBurdenElement targetCtxBurdenArray[NUM_TNC_TYPES];

   int ctxIdx=0;



   // First read the SBS signature file and work out how much burden we need
   // to spike in to each context. Setup targetCtxBurdenArray which records this.
   // Remember we need to record both the SBS signature entry in the file
   // (which is w.r.t. the forward genomic strand) and it's reverse complement
   // (which is on the reverse genomic strand).
   // Both provide opportunities to spike in the required mutation to BAM
   // files that have been pre-seperated into inserts aligning to the forward and
   // reverse genomic strands respectively.
   char *line = NULL;
   size_t len = 0;
   ssize_t characters;
   ssize_t tncArraySize=NUM_TNC_TYPES;

   while((characters = getline(&line, &len, fp)) != -1) {
      char * pch;

      // Skip blank lines.
      if(characters == 0)
         continue;

      // Skip blank lines.
      if(line[0] == '\n')
         continue;

      // Skip metadata/comments.
      if(line[0] == '#')
         continue;

      // Parse the record.

      // TNC.
      pch = strtok(line,"\t");
      if(pch != NULL) {
         //  printf(">>>%s<<<\n",pch);


         targetCtxBurdenArray[ctxIdx].context[0]=pch[0];
         targetCtxBurdenArray[ctxIdx].context[0]=pch[0];
         targetCtxBurdenArray[ctxIdx].context[1]=pch[2];
         targetCtxBurdenArray[ctxIdx].context[2]=pch[6];
         targetCtxBurdenArray[ctxIdx].context[3]='\0';

         targetCtxBurdenArray[ctxIdx].substitution=pch[4];


         // Reverse complement on the other strand.

         // Complement.
         for(int i=0; targetCtxBurdenArray[ctxIdx].context[i]!='\0'; i++) {
            if(targetCtxBurdenArray[ctxIdx].context[i]=='A')
               targetCtxBurdenArray[ctxIdx].contextPrime[i]='T';
            else if(targetCtxBurdenArray[ctxIdx].context[i]=='T')
               targetCtxBurdenArray[ctxIdx].contextPrime[i]='A';
            else if(targetCtxBurdenArray[ctxIdx].context[i]=='G')
               targetCtxBurdenArray[ctxIdx].contextPrime[i]='C';
            else if(targetCtxBurdenArray[ctxIdx].context[i]=='C')
               targetCtxBurdenArray[ctxIdx].contextPrime[i]='G';
         }
         targetCtxBurdenArray[ctxIdx].contextPrime[3]='\0';

         // Now reverse the order of the sequence for the reverse strand..
         // Its only tri-nucleotide so just swap outer bases.
         char tmpC = targetCtxBurdenArray[ctxIdx].contextPrime[0];
         targetCtxBurdenArray[ctxIdx].contextPrime[0] = targetCtxBurdenArray[ctxIdx].contextPrime[2];
         targetCtxBurdenArray[ctxIdx].contextPrime[2] = tmpC;

         // Complement of substitution base.
         if(targetCtxBurdenArray[ctxIdx].substitution=='A')
            targetCtxBurdenArray[ctxIdx].substitutionPrime='T';
         else if(targetCtxBurdenArray[ctxIdx].substitution=='T')
            targetCtxBurdenArray[ctxIdx].substitutionPrime='A';
         else if(targetCtxBurdenArray[ctxIdx].substitution=='G')
            targetCtxBurdenArray[ctxIdx].substitutionPrime='C';
         else if(targetCtxBurdenArray[ctxIdx].substitution=='C')
            targetCtxBurdenArray[ctxIdx].substitutionPrime='G';

         // printf("HERE1111: %s > %c.\n",targetCtxBurdenArray[ctxIdx].context,targetCtxBurdenArray[ctxIdx].substitution);
         // printf("HERE2222: %s > %c.\n",targetCtxBurdenArray[ctxIdx].contextPrime,targetCtxBurdenArray[ctxIdx].substitutionPrime);

      }

      // Burden fraction in this context.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         char* strEnd;
         float f = strtof(pch, &strEnd);
         if(pch == strEnd) {
            fprintf(stderr, "Unrecognised context burden fraction, \"%s\" (have you removed the header line?)\n", pch);
            exit(1);
         }
         targetCtxBurdenArray[ctxIdx].targetBurden=(totalBurden*f);
         ctxIdx++;
      }
   }

   tncArraySize=ctxIdx+1;

   // Now read the set of random loci from the stdin and spike in as required.
   hts_pos_t beg;
   hts_pos_t end;
   char *endPtr= NULL;
   int overflow = 0;

   while((characters = getline(&line, &len, stdin)) != -1) {

      // Exit as soon as we hit blank lines or anything that isn't a genomic range.
      char* contigPtr = strtok(line, "\t");
      if(contigPtr == NULL) {
         continue;
      }
      // Strip newline if its there...
      contigPtr[strcspn(contigPtr, "\n")] = 0;

      char* sbsTargetPtr = strtok(NULL, "\t");
      if(sbsTargetPtr == NULL) {
         continue;
      }
      // Strip newline if its there...
      sbsTargetPtr[strcspn(sbsTargetPtr, "\n")] = 0;

      beg = hts_str2uint(sbsTargetPtr, &endPtr, 63, &overflow);
      if(overflow) {
         fprintf(stderr, "\nPosition value '%s' is too large\n", argv[3]);
         exit(-2);
      }
      else {
         // Get the TNC range.
         beg -= 1;
         end = beg+3;
      }


      hts_pos_t fai_ref_len;
      char *fai_ref = faidx_fetch_seq64(fai, contigPtr, beg, end, &fai_ref_len);
      if(fai_ref_len < 3) {
         fprintf(stderr,"Failed to fetch the sequence, trying next one..\n"); // no worries, move onto next one if so...
         continue;
      }
      else {
         // Check if a mutation within this context is required?
         int burdenCompleteCount=0;
         for(int i=0; i<tncArraySize; i++) {
            // If we've finished creating mutation records for this context then continue..
            if(!targetCtxBurdenArray[i].targetBurden) {
               if(i==0)
                  burdenCompleteCount=0;
               else
                  burdenCompleteCount++;

               if(burdenCompleteCount==95) {
                  // We're done!! Don't hang around, exit..
                  fprintf(stderr, "Target burden hit...Bye!\n");
                  fflush(fwdFp);
                  fflush(revFp);
                  fclose(fwdFp);
                  fclose(revFp);
                  exit(0);
               }

               continue;
            }

            if(!strncmp(targetCtxBurdenArray[i].context, fai_ref, 3)) {

               //    if((beg+2)==924769 || (beg+1)==924769 || (beg)==924769)
               //    {
               // Context found on forward genomic strand.
               //       fprintf(stderr, "%s\t%ld\ti=%d\tURHERE Found \"%c%c%c\" on forward strand with \"%s\" with substitution %c>%c (%ld SBSs remaining).\n",contigPtr,beg+2,i, fai_ref[0],fai_ref[1],fai_ref[2],targetCtxBurdenArray[i].context, fai_ref[1],targetCtxBurdenArray[i].substitution,targetCtxBurdenArray[i].targetBurden);

               //      }


               //printf("%c%c%c\n",fai_ref[0],fai_ref[1],fai_ref[2]);

               // Print out spike-in record for the config and adjust remaing burden total.
               // Adjust the target pointer. With htslib we are working off zero base.
               // However the spike-in config must be in one base.
               fprintf(fwdFp,"%s\t%ld\t%c\n",contigPtr,beg+2,targetCtxBurdenArray[i].substitution);
               // fprintf(stderr, "Found \"%c%c%c\" on forward strand with \"%s (%s)\" with substitution %c>%c (%ld SBSs remaining).\n", fai_ref[0],fai_ref[1],fai_ref[2],targetCtxBurdenArray[i].context, targetCtxBurdenArray[i].context, fai_ref[1],targetCtxBurdenArray[i].substitution,targetCtxBurdenArray[i].targetBurden);

               targetCtxBurdenArray[i].targetBurden--;

               break;

            }

            else if(!strncmp(targetCtxBurdenArray[i].contextPrime, fai_ref, 3)) {
               // Context found on reverse genomic strand.
               // fprintf(stderr, "Found \"%c%c%c\" on reverse strand with \"%s (%s)\" with substitution %c>%c (%ld SBSs remaining).\n", fai_ref[0],fai_ref[1],fai_ref[2],targetCtxBurdenArray[i].context, targetCtxBurdenArray[i].contextPrime, fai_ref[1],targetCtxBurdenArray[i].substitutionPrime,targetCtxBurdenArray[i].targetBurden);

               // Print out spike-in record for the config and adjust remaing burden total.
               // Adjust the target pointer. With htslib we are working off zero base.
               // However the spike-in config must be in one base.
               fprintf(revFp,"%s\t%ld\t%c\n",contigPtr,beg+2,targetCtxBurdenArray[i].substitutionPrime);
               //printf("i=%d, REF=>\"%c%c%c\"\t%s\t%ld\t%s -> %c\n",i,fai_ref[0],fai_ref[1],fai_ref[2], contigPtr,beg+1,&(targetCtxBurdenArray[i].contextPrime[0]), targetCtxBurdenArray[i].substitutionPrime);
               targetCtxBurdenArray[i].targetBurden--;

               break;
            }

         }

      }


   }
   // Did we cover all the burden required by the SBS signature profile input file?
   for(int i=0; i<tncArraySize; i++) {
      if(targetCtxBurdenArray[i].targetBurden != 0) {
         fprintf(stderr, "Failed to set required burden for \"%s\" (%ld SBSs remaining)\nConsider increasing the number of random loci as input.\n", targetCtxBurdenArray[i].context,targetCtxBurdenArray[i].targetBurden);
         fflush(NULL);
         fclose(fwdFp);
         fclose(revFp);
         exit(1);
      }
   }
   fflush(NULL);
   fclose(fwdFp);
   fclose(revFp);
}
