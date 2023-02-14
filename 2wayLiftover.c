#define  _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

// chr refMarker altMarker step

// step = both

// We need to map from a current locus, to a liftover locus.
// The current locus can be in either the standard reference or the personalised reference.
// If the current locus is in the standard reference then the liftover is the personal and vice versa.
typedef enum { BOTH=0, REFERENCE, DONOR} StepType;

const char* stepTypeNames[] = {"BOTH", "REFERENCE", "DONOR" };


//typedef enum { BOTH=0, PERSONAL, STANDARD } StepType;
typedef enum { PERSONAL, STANDARD } CurrentType;

int lift(int currentTargetLocus, int currentMarker, int liftoverMarker, bool stepliftTarget)
{
   // If we're at the first in range just return the liftoverMarker.
   // This will also cover the case where stepType == LIFTOVER
   // Remember if liftover is not stepping in this range, then you will always return liftoverMarker regardless.
   if(!stepliftTarget || currentTargetLocus == currentMarker)
      return(liftoverMarker);
   // if step == both, then add the delta to the liftover marker to get your liftover target.
   else {
      // Subtract the current marker (either standard or personal) from your current target. This is the delta.
      int delta=currentTargetLocus-currentMarker;
      // printf("\nURHERE delta=%d, currentTargetLocus=%d, currentMarker=%d\n",delta, currentTargetLocus, currentMarker);
      return(liftoverMarker+delta);
   }
}


typedef struct {
   char *contigPtr; // Pointer to chromosome
   long refPosition; // position in reference genome
   long donorPosition; // position in donor genome
   StepType stepType;
} LiftoverRecord;

void readNextCondensedLiftRecord(LiftoverRecord *target, FILE * liftFp, char *skipContig)
{
   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   while((read = getline(&lineP, &len, liftFp)) != -1) {

      // Skip blank lines.
      if(lineP[0] == '\n') {
         if(lineP) {
            free(lineP);
         }
         continue;
      }

      // Skip metadata.
      if(lineP[0] == '#') {
         if(lineP) {
            free(lineP);
         }
         continue;
      }

      // Parse the VCF record.

      // CHROM.
      pch = strtok(lineP,"\t");
      if(pch != NULL) {
         // Have we been told to skip all targets on this contig?
         if(skipContig && !strcmp(skipContig,pch)) {
            continue;
         }
         strcpy(target->contigPtr,pch);
      }

      // REF.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->refPosition = atol(pch);
      }

      // DONOR.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->donorPosition = atol(pch);
      }

      // STEP TYPE, 0=both, 1=reference, 2=donor.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->stepType = atoi(pch);
      }

      // These files can get big (particulary liftover file).
      // Free as we go along to avoid running out of memory.
      if(lineP) {
         free(lineP);
      }

      return;
   }
   target->contigPtr = NULL;
}


int main(int argc,char* argv[])
{


   // Open the liftover file..

   FILE * condensedLiftFileB = fopen(argv[3], "r");
   if(condensedLiftFileB == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[3]);
      return -1;
   }

   LiftoverRecord *thisRecordPtr;
   LiftoverRecord *previousRecordPtr;
   LiftoverRecord *tmpRecordPtr;


   LiftoverRecord liftTargetA;
   // Set up a location to store the contig name/chromosome.
   // TODO!!!! check for buffer overflow!!!!
#define MAX_LINE_LEN 1024
   char liftContigA[MAX_LINE_LEN] ="XXX";
   liftTargetA.contigPtr = liftContigA;


   LiftoverRecord liftTargetB;
   char liftContigB[MAX_LINE_LEN] ="XXX";
   liftTargetB.contigPtr = liftContigB;


   StepType thisStep = DONOR;
   StepType previousStep = DONOR;

   thisRecordPtr = &liftTargetA;
   previousRecordPtr = &liftTargetB;

   long liftLocus;
   char *locusPtr = argv[2];
   char *currentTargetContigPtr = argv[2];


   // Parse out target contig and locus from command line arg.
   while(*locusPtr != '\0') {
      if(*locusPtr == ':') {
         *locusPtr='\0';
         locusPtr++;
         break;
      }
      locusPtr++;
   }

   // Format check.
   if(locusPtr == (currentTargetContigPtr + strlen(currentTargetContigPtr))) {
      fprintf(stderr,"%s: Required format of genomic coordinate is <contig>:<locus>, ie., chr3:12345678 .\n", argv[0]);
      exit(1);
   }

   int currentTargetLocus = atoi(locusPtr);

   // What direction have we been asked to translate the target genomic coordinate from,
   // ie., REFERENCE -> DONOR or DONOR -> REFERENCE

   // If its REFERENCE -> DONOR
   if(!strcmp(argv[1],"r")) {
      while(1) {

         readNextCondensedLiftRecord(thisRecordPtr, condensedLiftFileB, NULL);

         // Continue until we get to the contig of interest.
         if(!thisRecordPtr->contigPtr || strcmp(thisRecordPtr->contigPtr,currentTargetContigPtr)) {
            // OK, this record is not on the contig we're interested in.
            // Check if we have we gone through all records relating to the contig of interest?
            if(!strcmp(previousRecordPtr->contigPtr,currentTargetContigPtr)) {
               // The previous record was the last record relating to this contig.
               // Is currentTargetLocus beyond the end of the current contig?
               if(previousRecordPtr->stepType == DONOR) {
                  printf("Target not found...\n");
                  exit(0);
               }

               // Otherwise, do the liftover..
               liftLocus = lift(currentTargetLocus, previousRecordPtr->refPosition, previousRecordPtr->donorPosition, previousRecordPtr->stepType != REFERENCE);

               // Print result.
               printf("%s\t%ld\n",previousRecordPtr->contigPtr,liftLocus);
               exit(0);

            }
            else {
               // We just havent reached the contig we're interested in yet.
               // If we are at the end of the file, let user know it has not been found.
               // Otherwise keep going until we get to the records relating to this contig.
               if(!thisRecordPtr->contigPtr)
                  break;
               else
                  continue;
            }
         }


         if(currentTargetLocus < thisRecordPtr->refPosition) {
            // Translating from REFERENCE -> DONOR
            liftLocus = lift(currentTargetLocus, previousRecordPtr->refPosition, previousRecordPtr->donorPosition, previousRecordPtr->stepType != REFERENCE);

            // Print result.
            printf("%s\t%ld\n",previousRecordPtr->contigPtr,liftLocus);
            exit(0);
         }
         else if(currentTargetLocus == thisRecordPtr->refPosition) {
            liftLocus=thisRecordPtr->donorPosition;
            // Print result.
            printf("%s\t%ld\n",thisRecordPtr->contigPtr,liftLocus);
            exit(0);
         }

         // Rotate the record and step pointers.
         tmpRecordPtr = previousRecordPtr;
         previousRecordPtr = thisRecordPtr;
         thisRecordPtr = tmpRecordPtr;
         previousStep = thisStep;

      }

   }

   // Otherwise, if its DONOR -> REFERENCE
   else if(!strcmp(argv[1],"d")) {
      while(1) {

         readNextCondensedLiftRecord(thisRecordPtr, condensedLiftFileB, NULL);

         // Continue until we get to the contig of interest.
         if(!thisRecordPtr->contigPtr || strcmp(thisRecordPtr->contigPtr,currentTargetContigPtr)) {
            // OK, this record is not on the contig we're interested in.
            // Check if we have we gone through all records relating to the contig of interest?
            if(!strcmp(previousRecordPtr->contigPtr,currentTargetContigPtr)) {
               // The previous record was the last record relating to this contig.
               // Is currentTargetLocus beyond the end of the current contig?
               if(previousRecordPtr->stepType == REFERENCE) {
                  printf("Target not found...\n");
                  exit(0);
               }

               // Otherwise, do the liftover..
               liftLocus = lift(currentTargetLocus, previousRecordPtr->donorPosition, previousRecordPtr->refPosition, previousRecordPtr->stepType != DONOR);

               // Print result.
               printf("%s\t%ld\n",previousRecordPtr->contigPtr,liftLocus);
               exit(0);

            }
            else {
               // We just havent reached the contig we're interested in yet.
               // If we are at the end of the file, let user know it has not been found.
               // Otherwise keep going until we get to the records relating to this contig.
               if(!thisRecordPtr->contigPtr)
                  break;
               else
                  continue;
            }
         }

         if(currentTargetLocus<thisRecordPtr->donorPosition) {
            // Translating from DONOR -> REFERENCE
            liftLocus = lift(currentTargetLocus, previousRecordPtr->donorPosition, previousRecordPtr->refPosition, previousRecordPtr->stepType != DONOR);

            // Print result.
            printf("%s\t%ld\n",previousRecordPtr->contigPtr,liftLocus);
            exit(0);
         }
         else if(currentTargetLocus == thisRecordPtr->donorPosition) {
            liftLocus=thisRecordPtr->refPosition;
            // Print result.
            printf("%s\t%ld\n",thisRecordPtr->contigPtr,liftLocus);
            exit(0);
         }

         // Rotate the record and step pointers.
         tmpRecordPtr = previousRecordPtr;
         previousRecordPtr = thisRecordPtr;
         thisRecordPtr = tmpRecordPtr;
         previousStep = thisStep;


         // printf("\nHERE22222!!!%s\t%ld\t%ld\t%s\n",previousRecordPtr->contigPtr,previousRecordPtr->refPosition, previousRecordPtr->donorPosition,stepTypeNames[previousRecordPtr->stepType]);

      }
   }

   printf("Target not found (end of file)...\n");
   exit(0);

}


// ./createSomaticVariantDistr.R 0.000005 0.001 .8 60 75828711 75828391

// ./createSomaticVariantDistr.R 0.000005 0.001 .8 70 75828711 75828391

// ./createSomaticVariantDistr.R 0.000009 0.002 .8 65 75828711 75828391

// ./createSomaticVariantDistr.R 0.000005 0.002 .7 60 75828711 75828391





