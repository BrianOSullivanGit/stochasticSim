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

int lift(int currentTarget, int currentMarker, int liftoverMarker, bool stepliftTargetB)
{

   // If we're at the first in range just return the liftoverMarker.
   // This will also cover the case where stepType == LIFTOVER
   // Remember if liftover is not stepping in this range, then you will always return liftoverMarker regardless.
   if(currentTarget == currentMarker || !stepliftTargetB)
      return(liftoverMarker);
   // if step == both, then add this delta to the liftover marker to get your liftover target.
   else {
      // Subtract the current marker (either standard or personal) from your current target. This is the delta.
      int delta=currentTarget-currentMarker;
      return(liftoverMarker+delta);
   }

}


typedef struct {
   char *contigPtr; // Pointer to chromosome
   long refPosition; // position in reference genome
   long donorPosition; // position in donor genome
} LiftoverRecord;

void readNextLiftoverRecord(LiftoverRecord *target, FILE * liftFp, char *skipContig)
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

   FILE * liftFileB = fopen(argv[1], "r");
   if(liftFileB == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[1]);
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

   readNextLiftoverRecord(previousRecordPtr, liftFileB, NULL);




   while(thisRecordPtr->contigPtr) {


      readNextLiftoverRecord(thisRecordPtr, liftFileB, NULL);
      
      if(!thisRecordPtr->contigPtr)
      break;

      if(thisRecordPtr->refPosition == previousRecordPtr->refPosition + 1 && thisRecordPtr->donorPosition == previousRecordPtr->donorPosition + 1)
         thisStep = BOTH;
      else if(thisRecordPtr->refPosition == previousRecordPtr->refPosition && thisRecordPtr->donorPosition == previousRecordPtr->donorPosition + 1)
         thisStep = DONOR;
      else
         thisStep = REFERENCE;

      // printf("Step = %s\n",stepTypeNames[thisStep]);
      if(previousStep != thisStep || strcmp(thisRecordPtr->contigPtr,previousRecordPtr->contigPtr))
         printf("%s\t%ld\t%ld\t%d\n",previousRecordPtr->contigPtr,previousRecordPtr->refPosition, previousRecordPtr->donorPosition,(int)thisStep);

      // Rotate the record and step pointers.
      tmpRecordPtr = previousRecordPtr;
      previousRecordPtr = thisRecordPtr;
      thisRecordPtr = tmpRecordPtr;
      previousStep = thisStep;

      // Which is stepping.. ref, personal both..?
      // stepType

      // Convert section...

      // Start with stepType = both... and ref = 0, personal = 0

      // Read line.

      // Has the stepType (or chromosome) changed?

      // If so print out the old step range and set up the next one.

      // Otherwise keep going.


      // Read next line
      //   printf("%s\t%ld\t%ld\n",liftTargetB.contigPtr,liftTargetB.refPosition, liftTargetB.donorPosition);
      //   readNextLiftoverRecord(&liftTargetB, liftFileB, argv[1]);
   }

   exit(0);



   StepType stepType = BOTH;
   CurrentType currentType = PERSONAL;

   /////////////////
   // or


   // if(currentType==PERSONAL)
   //{
   // currentMarker = liftRec.personal;
   // liftoverMarker = liftRec.standard;
   // else
   // currentMarker = liftRec.standard;
   // liftoverMarker = liftRec.personal;

   ////////////////

   int currentMarker = 39;


   int liftoverMarker = 56;

   int currentTarget = atoi(argv[1]);

   int liftoverTarget;

   // printf("currentTarget = %d, liftoverTarget = %d\n",currentTarget, liftoverTarget);
   liftoverMarker=7;
   for(int i=5; i<20; i++) {
      currentTarget=i;
      currentMarker=i-4;


      liftoverTarget =  lift(currentTarget,currentMarker,liftoverMarker, true);
      printf("%d\t%d\n",currentTarget, liftoverTarget);
   }




   // Keep reading until you read past your target.

   // Now examine the record before
}


