// Once you modify a reference with an inserrtion or deletion everything downstream
// of that modification is not in the same genomic location anymore.
// This means, if you create a personalised genome and are interested in for example, the exome
// part of that genome you must map your standard bed file containing the ranges from
// within the standard reference that comprise the exome to a new bedfile that corresponds
// to the personalised reference you have created.
// Using a liftover file that was created while creating the personalised genome,
// that is exactly what this code does.

// Maps locus in donor genome to the equivalent in reference.

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <sys/types.h>

#define MAX_LINE_LEN 1024
#define FA_LINE_LENGTH 70

typedef struct {
   char *contigPtr; // Pointer to chromosome
   long refPosition; // position in reference genome
   long donorPosition; // position in donor genome
} LiftoverRecord;


typedef struct {
   char *contigPtr; // Pointer to chromosome
   long rangeStart; // Range start is 0 based in BED
   long rangeEnd; // Range start is 1 based in BED
   char *miscFields; // Pointer to any remaining fields
} BedRecord;

void getBedRecord(BedRecord *target, FILE * bedFp, char *skipContig)
{
   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   while((read = getline(&lineP, &len, bedFp)) != -1) {

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

      // Parse the BED record.

      // CHROM.
      pch = strtok(lineP,"\t");
      if(pch != NULL) {
         // Have we been told to skip all targets on this contig?
         if(skipContig && !strcmp(skipContig,pch)) {
            continue;
         }
         strcpy(target->contigPtr,pch);
      }

      // Range start.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->rangeStart = atol(pch)+1; // Range start is 0 based in BED, convert to 1 based to match liftover file.
      }

      // Range end.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->rangeEnd = atol(pch); // Range end is 1 based in BED
      }


      // Misc fields.
      pch = strtok(NULL,"\t");
      target->miscFields = pch; // Any misc fields if they're there, NULL otherwise.

      // These files can get big (particulary liftover file).
      // Free as we go along to avoid running out of memory.
      // (ok we use static but no need to make this re-entrant, it will do)
      if(lineP) {
         free(lineP);
      }

      return;
   }
   target->contigPtr = NULL;
}


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
   if(argc != 3 &&  argc != 4) {
      fprintf(stderr, "%s: Usage: %s <liftover file> <bed file> <chromosome to skip, optional>\n", argv[0], argv[0]);
      exit(-1);
   }

   LiftoverRecord liftTarget;
   // Set up a location to store the contig name/chromosome.
   // TODO!!!! check for buffer overflow!!!!
   char liftContig[MAX_LINE_LEN];
   liftTarget.contigPtr = liftContig;

   BedRecord bedRecord;
   // Set up a location to store the contig name/chromosome.
   // TODO!!!! check for buffer overflow!!!!
   char bedContig[MAX_LINE_LEN];
   bedRecord.contigPtr = bedContig;

   FILE * liftFileB = fopen(argv[1], "r");
   if(liftFileB == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[1]);
      return -1;
   }

   readNextLiftoverRecord(&liftTarget, liftFileB, argv[3]);

   FILE * bedFile = fopen(argv[2], "r");
   if(liftFileB == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[2]);
      return -1;
   }

   getBedRecord(&bedRecord, bedFile, argv[3]);


   while(bedRecord.contigPtr) {

      // If our BED record is on another contig to the one we're on at the moment then fast forward to there.
      if(liftTarget.contigPtr && strcmp(liftTarget.contigPtr,bedRecord.contigPtr)) {
         while(liftTarget.contigPtr &&  strcmp(liftTarget.contigPtr,bedRecord.contigPtr)) {
            // TODO!!!!
            // Just print out these contigs, as is..
            // They are likely unlocalised/unplaced contigs for mapping purposes that will not change in the liftover..
            // NO remove them!!!! we don't want to simulate reads from here!!!!
            // REMOVE THEM BEFORE YOU GET TO THIS PROGRAM...
            readNextLiftoverRecord(&liftTarget, liftFileB, argv[3]);
         }
      }


      //fprintf(stderr, "\n%s:%ld-%ld\n", bedRecord.contigPtr, bedRecord.rangeStart, bedRecord.rangeEnd);
      //fprintf(stderr, "\n%s:%ld-%ld\n", liftTarget.contigPtr, liftTarget.refPosition, liftTarget.donorPosition);
      while(liftTarget.contigPtr && !strcmp(liftTarget.contigPtr,bedRecord.contigPtr) && (liftTarget.refPosition != bedRecord.rangeStart)) {
         readNextLiftoverRecord(&liftTarget, liftFileB, argv[3]);
      }

      if(!liftTarget.contigPtr) {
         fprintf(stderr,"\nEnd of liftover file...(?). Last BED record was: %s:%ld-%ld\n", bedRecord.contigPtr, bedRecord.rangeStart, bedRecord.rangeEnd);
         exit(0);
      }

      // Remember we have to convert the first number in the range to zero base as per BED format.
      printf("%s\t%ld\t",liftTarget.contigPtr,liftTarget.donorPosition-1);

      while(liftTarget.contigPtr && !strcmp(liftTarget.contigPtr,bedRecord.contigPtr) && (liftTarget.refPosition != bedRecord.rangeEnd)) {
         readNextLiftoverRecord(&liftTarget, liftFileB, argv[3]);
      }

      printf("%ld\n",liftTarget.donorPosition);

      getBedRecord(&bedRecord, bedFile, argv[3]);
      readNextLiftoverRecord(&liftTarget, liftFileB, argv[3]);

   }


}
// BLKWF13


// Notes.
// Get your bed from,

// Whole gene.
// wget 'https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=273356791_npiNAcOyK98nyqxfnTBDfT9ue7zC&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED'  -O whole_gene_exome_hg38.bed.gz

// Coding exons.
// wget 'https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=273356936_r2aHSDfHA27pmamBP4EClOX8Ac6s&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED' -O coding_exons_hg38.bed.gz

// 5' UTR
// https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=273356936_r2aHSDfHA27pmamBP4EClOX8Ac6s&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=utr5&fbDownBases=200&hgta_doGetBed=get+BED

// 5' UTR
// https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=273356936_r2aHSDfHA27pmamBP4EClOX8Ac6s&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=utr5&fbDownBases=200&hgta_doGetBed=get+BED


// Merge ranges that overlap with,
// bedtools merge -i <(zcat ./whole_gene_exome_hg38.bed.gz|cut -f1-3) > whole_gene_exome_hg38.merged.bed

// Now you're good to go to liftover the ranges to your personalised genome.

