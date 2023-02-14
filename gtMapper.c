#define  _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


// Correlate a list of target genomic coordinates (for example, candidate variants from a caller VCF)
// with corresponding GT records. Print a null record if for entries when no corelation exists.
//
// Usage:
// gtMapper <GT VCF file> <Tab sep. file containing target list>

int main(int argc,char* argv[])
{


   // Open the liftover file..

   FILE * vcfFp = fopen(argv[1], "r");
   if(vcfFp == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[1]);
      return -1;
   }

   FILE * targetFp = fopen(argv[2], "r");
   if(targetFp == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[2]);
      return -1;
   }


   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;

   char * targetLineP = NULL;
   size_t targetLen = 0;
   ssize_t targetRead;
   char * targetChrPtr = NULL;
   int targetContigNum = 0;
   char * targetPosPtr = NULL;
   long targetPos;
   char * pch;

   //Outer loop => caller VCF we want to map to the ground truth.
   while((read = getline(&lineP, &len, vcfFp)) != -1) {
      //printf("%s\n",lineP);
      char * chrPtr = NULL;
      int contigNum = 0;
      char * posPtr = NULL;
      long pos;
      // char * pch;

      // Skip blank lines.
      if(lineP[0] == '\n') {
         continue;
      }

      // Skip metadata.
      if(lineP[0] == '#') {
         continue;
      }


      // Parse the VCF record.

      // CHROM.
      // Strip off 'chr' prefix if it exists and convert to contig number. X maps to 23, Y to 24.
      chrPtr = strtok(lineP,"\t");
      if(chrPtr == NULL) {
         continue;
      }

      char *p;

      if(chrPtr[0]=='c' && chrPtr[1]=='h' && chrPtr[2]=='r') {
         p = chrPtr+3;
      }
      else {
         p = chrPtr;
      }

      if(*p == 'X')
         contigNum = 23;
      else if(*p == 'Y')
         contigNum = 24;
      else
         contigNum = atoi(p);


      // POS.
      posPtr = strtok(NULL,"\t");
      if(posPtr == NULL) {
         continue;
      }
      pos = atol(posPtr);

      // Remaining fields.
      pch = strtok(NULL,"");

      // printf(">>>>GT: chrPtr = %s, posPtr = %s, pch = %s\n",chrPtr, posPtr, pch);

      // Inner loop => Search for the target record corresponding to this ground truth record.
      // Exits from this loop are,
      //
      // continue => move onto the next line in the target file.
      // break; => move onto next GT record.
      // targetChrPtr = NULL;break; => move onto next target and next GT record.
      // targetChrPtr = NULL;continue; => move onto next target.
      while(1) {
         if(targetChrPtr == NULL) {
            // Read in the next target if we're not in the process of still handling one already.
            read = getline(&targetLineP, &targetLen, targetFp);

            if(read == -1) {
               // No more targets left, we're done here...
               exit(0);
            }

            // Skip blank lines.
            if(targetLineP[0] == '\n') {
               continue;
            }

            // Skip metadata.
            if(targetLineP[0] == '#') {
               continue;
            }

            // Parse the VCF record.

            // CHROM.
            targetChrPtr = strtok(targetLineP,"\t");
            if(targetChrPtr == NULL) {
               continue;
            }

            char *p;

            if(targetChrPtr[0]=='c' && targetChrPtr[1]=='h' && targetChrPtr[2]=='r') {
               p = targetChrPtr+3;
            }
            else {
               p = targetChrPtr;
            }

            if(*p == 'X')
               targetContigNum = 23;
            else if(*p == 'Y')
               targetContigNum = 24;
            else
               targetContigNum = atoi(p);

            // POS.
            targetPosPtr = strtok(NULL,"\t");
            if(targetPosPtr == NULL) {
               continue;
            }
            targetPos = atol(targetPosPtr);

         }
         else {
            // printf(">>%s<<>>>%ld<<<\n",targetChrPtr,targetPos);
         }

         // Have we found our target location?
         if(contigNum == targetContigNum && targetPos == pos) {
            // printf("FOUND!!: targetChrPtr = %s, targetPos = %ld, chrPtr = %s, pos = %ld, %s\n",targetChrPtr, targetPos,chrPtr, pos, pch);
            printf("%s\t%s\t%s",chrPtr, posPtr, pch);
            // If so continue onto our next target and next record in the GT VCF..
            targetChrPtr = NULL;
            break;
         }
         // Have we moved past our target on this contig?
         else if(contigNum == targetContigNum &&  pos > targetPos) {
            //printf("PASSED OUT!!: targetChrPtr = %s, targetPos = %ld, chrPtr = %s, pos = %ld, %s\n",targetChrPtr, targetPos,chrPtr, pos, pch);
            // Print out null record.
            printf("%s\t%ld\t.\t.\t.\t.\t.\t.\t.\t.\n",targetChrPtr, targetPos);
            // If so continue onto our next target..
            targetChrPtr = NULL;
            continue;
         }
         // Have we moved past our target onto a new contig? TODO!!! will have to be greater / less than target contig....
         else if(contigNum > targetContigNum) {
            // Passed out target contig.
            //printf("NEW CONTIG: contigNum = %d, targetContigNum = %d,  targetChrPtr = %s, targetPos = %ld, chrPtr = %s, pos = %ld, %s\n", contigNum, targetContigNum, targetChrPtr, targetPos,chrPtr, pos, pch);
            // Print out null record.
            printf("%s\t%ld\t.\t.\t.\t.\t.\t.\t.\t.\n",targetChrPtr, targetPos);
            // Continue onto our next target..
            targetChrPtr = NULL;
            continue;
         }

         // Have we not reached the target contig in the GT yet? TODO!!!!
         // Target contig not yet reached
         if(targetContigNum > contigNum) {
            //printf("NOT YET, DIFFERENT CONTIG!!: targetChrPtr = %s, targetPos = %ld, chrPtr = %s, pos = %ld, %s\n",targetChrPtr, targetPos,chrPtr, pos, pch);
            // Move onto the next record in the ground truth.
            break;
         }

         else {
            // Otherwise assume we havent reached our target in the ground truth yet..
            //printf("IGNORE: targetChrPtr = %s, targetPos = %ld, chrPtr = %s, pos = %ld, %s\n",targetChrPtr, targetPos,chrPtr, pos, pch);
            // Move onto the next record in the ground truth.
            break;
         }
      }
   }



   if(targetLineP) {
      free(targetLineP);
   }

   // These files can get big (particulary liftover file).
   // Free as we go along to avoid running out of memory.
   if(lineP) {
      free(lineP);
   }

   //printf("Target not found (end of file)...\n");
   exit(0);

}
