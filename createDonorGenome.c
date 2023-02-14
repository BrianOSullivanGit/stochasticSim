// This binary creates one half of a diploid reference corresponding to a donor genotype.
// Simulate doner genome from a reference by combining it with a germline VCF.
// A set of SNPs are inserted into the reference so it then corresponds with a given individuals set of haplotypes.

// This program assumes the input VCF only contains the autosomes plus X and Y.
// You should remove (grep/whatever) all other contigs.
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <sys/types.h>

#define MAX_LINE_LEN 1024
#define FA_LINE_LENGTH 70

typedef struct {
   char *contigPtr; // Pointer to chromosome of next target
   char *ref; // Pointer to ref. allele that is going to replaced.
   char *alt; // Pointer to alt. allele that is going to replace the reference at the target.
   long position; // Target locus
} NextTargetSnp;

void getNextTargetSnp(char *progName, NextTargetSnp *target, FILE * vcfp, char *haploT, char *skipContig)
{
   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   char *contigPtr; // Pointer to chromosome of next target
   char *ref; // Pointer to ref. allele that is going to replaced.
   char *alt; // Pointer to alt. allele that is going to replace the reference at the target.
   long position; // Target locus

   // TODO!!!! free lineP as you are done with each line!!!!!

   while((read = getline(&lineP, &len, vcfp)) != -1) {
      // Skip blank lines.
      if(lineP[0] == '\n')
         continue;

      // Skip metadata.
      if(lineP[0] == '#')
         continue;

      // Parse the VCF record.

      // CHROM.
      pch = strtok(lineP,"\t");
      if(pch != NULL) {
         // Have we been told to skip all targets on this contig?
         if(skipContig && !strcmp(skipContig,pch)) {
            continue;
         }
         contigPtr = pch;
      }

      // POS.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         // Rem. target loci in 1 based format from input VCF containing (germline) SNPs only.
         // Convert to zero base for ease of use here.
         position = atol(pch)-1;
      }

      // ID (skip it).
      pch = strtok(NULL,"\t");

      // REF (skip it, we've got that already..).
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         ref = pch;
      }

      // ALT.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         alt = pch;
      }

      // QUAL (skip it).
      pch = strtok(NULL,"\t");

      // FILTER
      pch = strtok(NULL,"\t");

      // Skip variant if it is not annotated as PASS.
      if(strcmp("PASS",pch))
         continue;

      // INFO (skip it).
      pch = strtok(NULL,"\t");

      // FORMAT (skip it).
      pch = strtok(NULL,"\t");

      // DONOR (sample id field, parse out phasing/GT).
      // We assume if this is a multi-field the genotype (GT) is the first field.
      pch = strtok(NULL,"\t");

      // Only take first, GT field.
      // TODO!!!! Update this to handle case where GT field is not the first field.
      int i=0;
      while(pch[i] != '\0') {
         // Truncate if more than one field.
         if(pch[i] == ':' || pch[i] == '\n')
            pch[i] = '\0';
         i++;
      }

      // Continue if it is not a haplotype of interest.
      // (rem., homzygous loci (1|1) are always of interest)
      if(strcmp(haploT,pch) && strcmp("1|1",pch))
         continue;

      // Check for VCF inconsistency.
      // If the genomic position of the next target is before the previous target on the
      // same contig then we print a warning and skip it.
      // This may occur if the VCF is not sorted in coordinate order.
      // It may also occur if we encounter a SNP in a region that has already been removed
      // by a previous deletion record (ie., see HG00110, chr1:932613 and chr1:932617 .
      // You should avoid putting such incostiencies in your VCF.
      // However if they do end up in there then it will result in the record being skipped.

      //if(target->contigPtr[0] != '\0')
      //{
      //fprintf(stderr, "Previous VCF record at %s:%ld\n", target->contigPtr, target->position);
      //fprintf(stderr, "Current VCF record at  %s:%ld\n", contigPtr, position);
      //}

      if((target->contigPtr[0] != '\0') && !strcmp(target->contigPtr,contigPtr) && (position < target->position + strlen(target->ref))) {
       //  fprintf(stderr, "%s: Skipping inconsistent/redundant VCF record at %s:%ld\n", progName, contigPtr, position+1);
         continue;
      }
      else {
         strcpy(target->contigPtr,contigPtr);
         target->position = position;
         target->alt = alt;
         target->ref = ref;

         return;
      }
   }
}

int main(int argc,char* argv[])
{

   FILE * refFile;
   FILE * vcfFile;
   FILE * liftFp;
   char * line = NULL;
   size_t len = 0;
   ssize_t read;


   char targetContig[MAX_LINE_LEN];
   targetContig[0] = '\0';

   NextTargetSnp nextTarget;

   // Set up a location to store the contig name/chromosome.
   nextTarget.contigPtr = targetContig;

   char * pch;

   if(argc <= 4) {
      fprintf(stderr, "%s: Usage: createDonorGenome <reference fasta> <donor VCF> <genotype, ie. \"1|0\"> <sex chromosome to skip>\n", argv[0]);
      exit(-1);
   }
   refFile = fopen(argv[1], "r");
   if(refFile == NULL) {
      fprintf(stderr, "%s: Can't open %s..\n\n",argv[0],argv[1]);
      return -1;
   }

   vcfFile = fopen(argv[2], "r");
   if(vcfFile == NULL) {
      fprintf(stderr, "%s: Can't open %s..\n\n",argv[0], argv[2]);
      return -1;
   }


   // Liftover file, to map from reference genome to donor genome.
   if(argc>5) {
      liftFp = fopen(argv[5], "w");
      if(liftFp == NULL) {
         fprintf(stderr, "%s: Can't open %s..\n\n",argv[0],argv[5]);
         return -1;
      }
   }
   else {
      liftFp = fopen("liftover.txt", "w");
      if(liftFp == NULL) {
         fprintf(stderr, "%s: Can't open liftover.txt..\n\n",argv[0]);
         return -1;
      }
   }

   long currentPos = 0;
   long newRefPos = 0;

   char currentContigStr[MAX_LINE_LEN];

   currentContigStr[0] = '\0';

   getNextTargetSnp(argv[0], &nextTarget, vcfFile, argv[3], argv[4]);

   int charCntInLine=0;

   // Bases to be remove from the reference that overflow to the next line..
   int overflowBases=0;

   // Go through the reference.
   // Perhaps should not really have gone the route of line at a time IO but its done now.. TODO!!!! revisit..
   while((read = getline(&line, &len, refFile)) != -1) {
      // We've moved onto a new contig.
      // Parce out its name and set it to the current contig.
      if(line[0] == '>') {
         int i=0;
         for(; line[i]!='\0'; i++) {
            if(line[i] == ' ' || line[i] == '\t' ||  line[i] == '\n') {
               //    line[i]='\0';
               break;
            }
         }

         // Close off previous contig if there was one.
         if(currentContigStr[0] != '\0') printf("\n");

         // printf("%s\n",(pch+1));
         // strcpy(currentContigStr,(line+1));
         strncpy(currentContigStr,(line+1),i-1);
         currentContigStr[i-1] = '\0';

         // Have we been told to skip this contig?
         // When creating a diploid reference to match the doner we will need to leave out
         // one of the sex chromosomes in each reference of the diploid set.
         // For example the first reference of the diploid pair for a doner male will contain
         // the autosome plus the Y chromosome. The second reference, the autosomes plus
         // the X chromosome.
         // A female donor will also contain 2 references per set, each containing the autosomes
         // and X chromosome only.
         // The caller lets us know in argv[4] the requirements depending on which reference
         // pair they are constructing here.
         if(argc>4) {
            if(!strcmp(argv[4],currentContigStr)) {
               // fprintf(stderr,"\nSkipping %s\n",currentContigStr);
               continue;
            }
         }


         // Print out the contig header.
         // TODO!!!! we should really zero out the M5 field etc. as we are modifying the contents..
         // It is no longer valid.
         //printf("%s",line);
         printf(">%s\n",currentContigStr);
         currentPos = 0;
         newRefPos = 0;
         charCntInLine=0;
         continue;
      }


      if(argc>4) {
         if(!strcmp(argv[4],currentContigStr)) {
            // fprintf(stderr,"\nSkipping %s\n",currentContigStr);
            continue;
         }
      }


      // Go through the line one base at a time and modify as required.
      // Note: the way its structured at moment reference must end in a newline.
      // TODO!!!! ensure that or change accordingly...
      for(int i=0; i<read-1; i++) {
         // Did any bases we were asked to skip overflow onto the next line?
         if(overflowBases) {
            // TODO!!! what if overflow is > line length! Potential bug here..
            // Skip bases as required..

            // Fix with check if overflowBases > read-1
            // TODO!!!! neaten this up.
            // If the overflow goes past the current line..
            if(overflowBases >= read-1) {

               //fprintf(stderr, "1>>> %s\t%ld\t%ld\t%s\t%s\n", currentContigStr, 1+currentPos, 1+newRefPos,nextTarget.ref,nextTarget.alt);

               i+=read-1;
               currentPos+=read-1;
               overflowBases -= read-1;


               if(!overflowBases) {
                  //fprintf(stderr, "1.5 >>> %s\t%ld\t%ld\t%s\t%s\n", currentContigStr, 1+currentPos, 1+newRefPos,nextTarget.ref,nextTarget.alt);
                  // Now we've skipped the required bases, put back in the alt base.
                  printf("%c", nextTarget.alt[0]);

                  // TODO!!!! Bug!!!! did we forget this???? Fixed..
                  currentPos++;


                  // Format location of newline.
                  if(++charCntInLine == FA_LINE_LENGTH) {
                     // Print out newline.
                     printf("\n");
                     charCntInLine=0;
                  }

                  // TODO!!!! BUG???? should there be a getNextTargetSnp here?
                  getNextTargetSnp(argv[0], &nextTarget, vcfFile, argv[3], argv[4]);
               }

               break;
            }

            //fprintf(stderr, "2>>> %s\t%ld\t%ld\t%s\t%s\n", currentContigStr, 1+currentPos, 1+newRefPos,nextTarget.ref,nextTarget.alt);

            i+=overflowBases-1;
            currentPos+=overflowBases-1;
            overflowBases = 0;

            // Now we've skipped the required bases, put back in the alt base.
            printf("%c", nextTarget.alt[0]);

            // TODO!!!! Bug!!!! did we forget this???? Fixed..
            currentPos++;

            // Format location of newline.
            if(++charCntInLine == FA_LINE_LENGTH) {
               // Print out newline.
               printf("\n");
               charCntInLine=0;
            }

            // TODO!!!! BUG???? should there be a getNextTargetSnp here? Fixed..
            getNextTargetSnp(argv[0], &nextTarget, vcfFile, argv[3], argv[4]);
         }

         if(currentPos == nextTarget.position && !strcmp(nextTarget.contigPtr,currentContigStr)) {

            // If the base at this position is currently unknown, skip putting in the alt. allele..
            if(line[i]!='N') {

               int refLen=strlen(nextTarget.ref);
               int altLen=strlen(nextTarget.alt);

               // What kind of a substitution are we dealing with?
               if(refLen == 1) {
                  // One or more bases are replacing a reference base.
                  int j=0;
                  while(nextTarget.alt[j] != '\0') {
                     printf("%c", nextTarget.alt[j]);
                     j++;

                     // TODO!!! Check to see ref bases and actual bases in reference correspond..

                     // Format location of newline.
                     if(++charCntInLine == FA_LINE_LENGTH) {
                        // Print out newline.
                        printf("\n");
                        charCntInLine=0;
                     }
                  }

                  // Update liftover file. Rem., output is, like VCF, 1 based.
                  for(int x=0; x<altLen; x++)
                     fprintf(liftFp, "%s\t%ld\t%ld\n", currentContigStr, 1+currentPos, 1+newRefPos++);

                  currentPos++;
               }
               else {
                  // One or more bases are removed from the reference.
                  // Check if the bases to be removed overflows to next line?
                  // NOTE!!!! Be very careful with this one..
                  // The last character on the line is usually a newline character..
                  // Also these are array indices so run from 0 to size-1
                  // hence the number 2 to push past these.
                  // TODO!!!! strip newline to make this consistant!!!!
                  overflowBases=0;
                  int basesLeftOnThisLine=(read-2-i);
                  if(refLen > basesLeftOnThisLine) {
                     // At this point, since we know the start of this reference skip,
                     // we can update the liftover file before we complete the deletion.
                     //fprintf(stderr, "0>>> %s\t%ld\t%ld\t%s\t%s\n", currentContigStr, 1+currentPos, 1+newRefPos,nextTarget.ref,nextTarget.alt);
                     for(int x=0; x<refLen; x++)
                        fprintf(liftFp, "%s\t%ld\t%ld\t,%s,%d,%ld partA\n", currentContigStr, 1+currentPos+x, 1+newRefPos, nextTarget.ref, refLen, strlen(nextTarget.ref));

                     // Update position in the output fasta so it points to next base.
                     newRefPos++;

                     // Move on past what's here.
                     // We'll pick this up again the on the next line
                     // when we deal with the overflow.
                     overflowBases=refLen-basesLeftOnThisLine;
                     i+=basesLeftOnThisLine;
                     currentPos+=basesLeftOnThisLine;
                  }
                  else {
                     // Skip bases as required..
                     i+=refLen-1;

                     // Update liftover file.
                     // TODO!!!! prob. should take the "currentPos++" out of the printf..
                     // Make it explicit, put it no a line of its own (?).
                     for(int x=0; x<refLen; x++)
                        fprintf(liftFp, "%s\t%ld\t%ld\t,%s,%d,%ld partB\n", currentContigStr, 1+currentPos++, 1+newRefPos, nextTarget.ref, refLen, strlen(nextTarget.ref));

                     // Now we've skipped the required bases, put back in the alt base.
                     printf("%c", nextTarget.alt[0]);

                     // Update position in the output fasta so it points to next base.
                     newRefPos++;

                     // Format location of newline.
                     if(++charCntInLine == FA_LINE_LENGTH) {
                        // Print out newline.
                        printf("\n");
                        charCntInLine=0;
                     }

                  }

               }

               // BUG???? Is next line in twice???? This should be up top in overflow section..
               //  getNextTargetSnp(argv[0], &nextTarget, vcfFile, argv[3], argv[4]);

            }
            else {
               fprintf(stderr,"%s: Skipping 'N' base at %s:%ld.\n",argv[0],currentContigStr,currentPos);
               // Do not modify the reference here as we don't know the ref. base at this location..
               // TODO!!!! Should we just modify it anyways??
               printf("%c", line[i]);

               // Update liftover map.
               // Nothing skipped or inserted for this base so we've got a 1 to 1 map here..
               fprintf(liftFp, "%s\t%ld\t%ld\n", currentContigStr, 1+currentPos, 1+newRefPos);


               if(++charCntInLine == FA_LINE_LENGTH) {
                  // Print out newline as we skiped it in for loop above.
                  printf("\n");
                  charCntInLine=0;
               }
               currentPos++;
               newRefPos++;
            }

            // If we are not dealing with skipping bases that have overflowed onto
            // lines that will follow then we are done with this target.
            // Load up next target base..
            // If there are no more targets left (ie., we're at the end of the VCF)
            // the getNextTargetSnp will return a null target & the rest of the .fa will
            // be printed out unaltered.
            if(!overflowBases) {
               getNextTargetSnp(argv[0], &nextTarget, vcfFile, argv[3], argv[4]);
            }

         }
         else {
            // Nothing to modify in the reference here..
            printf("%c", line[i]);

            // Update liftover map.
            // Nothing skipped or inserted for this base so we've got a 1 to 1 map here..
            fprintf(liftFp, "%s\t%ld\t%ld\n", currentContigStr, 1+currentPos, 1+newRefPos);


            if(++charCntInLine == FA_LINE_LENGTH) {
               // Print out newline as we skiped it in for loop above.
               printf("\n");
               charCntInLine=0;
            }
            currentPos++;
            newRefPos++;
         }


      }


   }

   printf("\n");

   // Close Liftover file.
   // We will use this file to map from reference genome to donor genome & vice versa.
   fprintf(stderr, "%s: Done.\n",argv[0]);
   fclose(liftFp);

   return 0;
}

