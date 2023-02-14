// Masks out (sets to 'N') regions of in the reference that do not overlap with the target BED.
// The modified reference is then used with a read simulator like ART to simulated reads
// from the target (exome, gene panel etc.) only.

#define  _GNU_SOURCE

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LINE_LEN 1024

typedef struct {
   char *contigPtr; // Pointer to chromosome of next target
   char base; // Pointer to the base that is going to replace the reference at the target.
   long start; // Range start
   long end; // Range end
} NextTargetSnp;

// Get the next target from the BED.
// All regions in the reference, except those that fall within a target
// will be masked out with 'N's.
void getNextTargetRange(NextTargetSnp *target, FILE * bedp)
{
   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   target->start = 0;
   target->end = 0;

   while((read = getline(&lineP, &len, bedp)) != -1) {
      // Skip metadata.
      if(lineP[0] == '#' || lineP[0] == '\n')
         continue;
      else
         break;
   }

   // Return unchanged if we're at the end of the BED file.
   if(read<0)
      return;

   // Parse the BED record.

   // CHROM.
   pch = strtok(lineP,"\t");
   if(pch != NULL) {
      // printf("%s\n",pch);
      strcpy(target->contigPtr,pch);
   }

   // START.
   pch = strtok(NULL,"\t");
   if(pch != NULL) {
      target->start = atol(pch);
   }

   // END.
   // Convert the end location to also be zero based for consistency.
   // (BED is funny like that, start of range is zero based,
   // end of range is one based).
   pch = strtok(NULL,"\t");
   if(pch != NULL) {
      target->end = atol(pch)-1;
   }

}

int main(int argc,char* argv[])
{

   FILE * refFile;
   FILE * bedFile;
   char * line = NULL;
   size_t len = 0;
   ssize_t read;

   // TODO!!!! check for buffer overflow.
   char targetContig[MAX_LINE_LEN];

   NextTargetSnp nextTargetRegion;

   // Read target loci in 1 based format from input VCF containing (germline) SNPs only.
   // If phasing, remember to subset germline VCF into two seperately phased VCFs first.

   // Set up a location to store the contig name/chromosome.
   nextTargetRegion.contigPtr = targetContig;

   char * pch;

   if(argc != 3) {
      fprintf(stderr, "\n   Usage: %s <personalised reference> <target BED>  >  <masked personalised reference>\n\n",argv[0]);
      exit(0);
   }
   refFile = fopen(argv[1], "r");
   if(refFile == NULL) {
      fprintf(stderr, "%s: Can't open reference %s..\n\n",argv[0],argv[1]);
      exit(1);
   }

   bedFile = fopen(argv[2], "r");
   if(bedFile == NULL) {
      fprintf(stderr, "%s: Can't open BED %s..\n\n",argv[0],argv[2]);
      exit(2);
   }

   long currentPos = 0;

   // TODO!!!! check for buffer overflow.
   char currentContigStr[MAX_LINE_LEN];

   getNextTargetRange(&nextTargetRegion, bedFile);

   // Go through the reference.
   while((read = getline(&line, &len, refFile)) != -1) {

      // We've moved onto a new contig.
      // Parce out its name and set it to the current contig.
      if(line[0] == '>') {
         printf("%s",line);

         for(int i=0; line[i]!='\0'; i++) {
            if(line[i] == ' ' || line[i] == '\t' ||  line[i] == '\n') {
               line[i]='\0';
               break;
            }
         }

         strcpy(currentContigStr,(line+1));

         currentPos = 0;

         continue;
      }

      // Is there a section from the reference we need to print out from this line?
      // First check if the ranges overlap
      // ie, (StartLine <= EndBed) and (EndLine >= StartBed)

      if(currentPos <= nextTargetRegion.end && currentPos+read-1 >= nextTargetRegion.start &&
            !strcmp(currentContigStr,nextTargetRegion.contigPtr)) {

         // OK, there is an overlap withing this line with the last range read from the bed.
         // Go through it one base at a time & process the overlap.
         for(int i=0; i<read-1; i++) {
            // printf("\nURHERE222b!!! nextTargetRegion.contigPtr = %s, nextTargetRegion.start=%ld, nextTargetRegion.end=%ld, currentContigStr = %s, currentPos=%ld, read=%ld, i=%d\n",nextTargetRegion.contigPtr, nextTargetRegion.start, nextTargetRegion.end, currentContigStr, currentPos, read, i);
            // If we're within the region specified print out the reference base.
            // Otherwise print "N"
            if(currentPos >= nextTargetRegion.start &&
                  currentPos <= nextTargetRegion.end &&
                  !strcmp(currentContigStr,nextTargetRegion.contigPtr)) {
               printf("%c", line[i]);
            }
            else
              printf("N");

            // Have we delt with all of the current target region?
            // If so load up the next one from the bed.
            // If not, leave it as is, it will be completed with subsequent iterations.
            if(currentPos > nextTargetRegion.end && !strcmp(currentContigStr,nextTargetRegion.contigPtr)) {
               getNextTargetRange(&nextTargetRegion, bedFile);
            }
            currentPos++;
         }
         printf("\n");
      }

      else {
         // Nothing on this line falls within the target region.
         // Mask it out (with 'N').
         for(int i=0; i<read-1; i++)
            printf("N");
         printf("\n");
         currentPos += read-1;
      }

      if(strcmp(currentContigStr,nextTargetRegion.contigPtr)) {
         // If the current target is now on a new contig, skip on until we get to that contig in the reference.
         // fprintf(stderr, "%s Moving on, target range, %s:%ld-%ld, current contig %s.\n",argv[0], nextTargetRegion.contigPtr, nextTargetRegion.start, nextTargetRegion.end,currentContigStr);
         continue;
      }

      // Note: this was copied out of the if section above in case the
      // the bed file contained duplicate or overlapping regions
      // (even though the bed should really be sorted and merged with bedtools previously).
      // Have we delt with all of this region?
      // If so load up the next one from the bed.
      // If not, leave it as is, it will be completed with subsequent iterations.
      if(currentPos > nextTargetRegion.end) {
         getNextTargetRange(&nextTargetRegion, bedFile);
      }
   }
   return 0;
}

