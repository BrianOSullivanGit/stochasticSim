
#define  _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "htslib/sam.h"
#include "samtools.h"
#include "sam_opts.h"
#include <htslib/faidx.h>

#include <stdbool.h>

// This programs takes 1 BAMs as input,
// Adjust as required (here and below).
#define NUM_BAMS 1
#define MAX_PILEUP_SIZE 10000
#define MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP MAX_PILEUP_SIZE

// Set this to the longest insertion sequence (ie., number of bases) you will require.
// Change this value and recompile as required.
#define MAX_INSERTION_LENGTH 3000

// According to SAM spec..
#define MAX_QNAME_LENGTH 254

// Can't find a max RNAME length (ie., contig/chromosome name in BAM) in the spec anywhere
// but this should cover anything (I think!)..
#define MAX_RNAME_LENGTH 254

// Initial values for instance of a reference fa structure mplp_ref_t
#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

// Macro to mutate the base.
#define bam_seqi_set(s,i,b) ((s)[(i)>>1]=((s)[(i)>>1]&(0xf<<((i&1)<<2)))|((b)<<((~i&1)<<2)))

typedef struct {
   char qName[MAX_QNAME_LENGTH];
} AlignmentRecord;


// Struct to hold info about reference (for caching purposes)..
typedef struct {
   char *ref[2];
   int ref_id[2];
   hts_pos_t ref_len[2];
} mplp_ref_t;

// Put all BAM file related stuff together in a structure.
// Below, declare an instance of this structure for every BAM file we open to keep track of stuff..
// Btw., It's declared in a structure like this so we can pack all this info together and
// pass it to functions like read_bam below nice and neatly in the first arg of the function call.
// I guess there would be other ways of doing this without defining a struct, but it would probably
// be messy and a nightmare to keep track of, particulary it you were dealing with more than one BAM.
typedef struct {     // auxiliary data structure
   samFile *fpIn;     // file handle for input
   samFile *fpOut;     // file handle for input
   sam_hdr_t *hdr;  // the file header
   hts_itr_t *iter; // NULL if a region not specified
   char *sampName;  // the sample name taken from SM field in read group.

   // I've put reference stuff in here, seemed as good a place as any..
   mplp_ref_t *ref;  // Reference struct for caching (below).
   faidx_t *fai;  // Indexed reference file.

   // Any other business..
   int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;


// This record defines the required indel or base substitution(s) to be performed.
// It is parsed from a record in the config file.




bool qNameFound(AlignmentRecord *alignmentRecord, char *qNamePtr, int size)

{
   int i=0;
   while(i<size) {
      if(!strncmp(alignmentRecord[i].qName, qNamePtr, MAX_QNAME_LENGTH)) return true;
      i++;
   }
   return false;
}

// Returns the simulated result of a loaded coin toss.
// rand() returns a random number between 0 and RAND_MAX.
// Multiply RAND_MAX by whatever you want your coin probability to be.
// Is the random number (returned by rand()) less than this?
bool coinToss(double probability)
{
   return rand() < (probability*((double)RAND_MAX + 1.0));
}



// This function handles reading from the reference and caching to speed it up.
// I've hacked this function and pulled it in from samtools mpileup (bam_plcmd.c).
// You could literally do this in 3 or 4 lines back in main().
// However this function does some basic caching. I am not prepared to invest any time
// to write any improved caching of my own so I stole it.
// I had to hack it (in the aux_t struct etc) a bit to get it to fit but it works.
// I would expect it to be as fast as samtools mpileup.
// Btw., if you loaded the reference into ram before running then caching or not shouldn't matter,
// a blunt 4 or 5 lines back in main() should work just as fast (I think..)
// Btw.,
// It looks fairly straight forward to exclude BAQ from the the pileup if you want to
// (which, unless you know what it is doing and why you need it, you want to..).
// Just don't include the call to sam_prob_realn(), the function that kicks it all off.
//
static int mplp_get_ref(aux_t *ma, int tid, char **ref, hts_pos_t *ref_len)
{
   mplp_ref_t *r = ma->ref;

   if(!r || !ma->fai) {
      *ref = NULL;
      return 0;
   }

   // Do we need to reference count this so multiple mplp_aux_t can
   // track which references are in use?
   // For now we just cache the last two. Sufficient?
   if(tid == r->ref_id[0]) {
      *ref = r->ref[0];
      *ref_len = r->ref_len[0];
      return 1;
   }
   if(tid == r->ref_id[1]) {
      // Last, swap over
      int tmp_id;
      hts_pos_t tmp_len;
      tmp_id  = r->ref_id[0];
      r->ref_id[0]  = r->ref_id[1];
      r->ref_id[1]  = tmp_id;
      tmp_len = r->ref_len[0];
      r->ref_len[0] = r->ref_len[1];
      r->ref_len[1] = tmp_len;

      char *tc;
      tc = r->ref[0];
      r->ref[0] = r->ref[1];
      r->ref[1] = tc;
      *ref = r->ref[0];
      *ref_len = r->ref_len[0];
      return 1;
   }

   // New, so migrate to old and load new
   free(r->ref[1]);
   r->ref[1]     = r->ref[0];
   r->ref_id[1]  = r->ref_id[0];
   r->ref_len[1] = r->ref_len[0];

   r->ref_id[0] = tid;
   r->ref[0] = faidx_fetch_seq64(ma->fai,
                                 sam_hdr_tid2name(ma->hdr, r->ref_id[0]),
                                 0,
                                 HTS_POS_MAX,
                                 &r->ref_len[0]);

   if(!r->ref[0]) {
      r->ref[0] = NULL;
      r->ref_id[0] = -1;
      r->ref_len[0] = 0;
      *ref = NULL;
      return 0;
   }

   *ref = r->ref[0];
   *ref_len = r->ref_len[0];
   return 1;
}




// This function reads a BAM alignment from one BAM file.
// Every time we move up the pileup onto a new base, this function is called and the read
// to which that base belongs is examined to make sure we want to include it.
// You will see below that a pointer to this function is included when initiating the
// pileup ( bam_mplp_init(NUM_BAMS, read_bam, (void**)data) ) so htslib knows what function
// we want called to do our read level filtering etc..
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
   aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
   int ret;
   while(1) {
      ret = aux->iter? sam_itr_next(aux->fpIn, aux->iter, b) : sam_read1(aux->fpIn, aux->hdr, b);
      if(ret<0) break;

      // Debug Brian.
      // Only take full matches for now.
      // TODO!!!! Remove this!!!!
      // if(!(b->core.n_cigar == 1 && bam_cigar_op(bam_get_cigar(b)[0])==BAM_CMATCH)) continue;

      // Read must not be unmapped, secondary, (possibly) failed QC or dupicate.
      // Strictly paired reads only.
      if(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
      if((int)b->core.qual < aux->min_mapQ) continue;
      if(aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len) continue;
      // Filter out orphan reads (so it matches mpileup output).
      // Rem.,
      // #define BAM_FPAIRED        1
      // #define BAM_FPROPER_PAIR   2
      // So should also be able to do this externally with with,
      // samtools view -f 0x1 -F 0x2 aa.sam -o aa.filtered.sam ???
      // Brian.
      if((b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) continue;

      break;
   }
   return ret;
}

void writeOutAlignment(aux_t* dataPtr, bam1_t * b)
{
   int r = sam_write1(dataPtr->fpOut, dataPtr->hdr, b);

   if(r < 0) {
      fprintf(stderr, "Couldn't write out ...\n");
      exit(1);
   }
}

// Generates random number between min and max, inclusive.
// https://cboard.cprogramming.com/c-programming/145187-how-pick-random-number-between-x-y.html
int randomNum(int min, int max)
{
   int range, result, cutoff;

   if(min >= max) {
      fprintf(stderr, "Invalid parameters passed to randomNum()...\n");
      exit(1);
   }

   range = max-min+1;
   cutoff = (RAND_MAX / range) * range;

   // Rejection method, to be statistically unbiased.
   do {
      result = rand();
   }
   while(result >= cutoff);

   return result % range + min;
}

char selectMutantAllele(char wildType)
{
   static char bases[] = "GCAT";

   // Choose a base other than what's at this locus in the pileup at the moment.
   int i;

   // bases is a string so rem. to account for the '\0' at the end..
   int numOfBases=(sizeof(bases)/sizeof(bases[0])-1);

   do {
      i=randomNum(0,numOfBases-1);
   }
   while(bases[i] == wildType);

   // Return the base that replaced the previous one at this locus.
   return(bases[i]);
}




// The purpose of this function is to identify reads at the pileup that we want to mutate.
// Maybe put the idx in data also and then we could get rid of 3 args in this function????
int selectMutatedQnames(char *targetRegion,
                        samFile *in,
                        hts_idx_t *idx,
                        bam_hdr_t *header,
                        AlignmentRecord *mutatedAlignmentRecord,
                        int *mutatedAlignmentCountPtr,
                        int *totalNumberOfMutatedAlignments,
                        float mutantAlleleFreq)
{
   static AlignmentRecord alignmentRecord[MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP] = {{{0}}};
   int alignmentCount=0;

   static AlignmentRecord uniqueAlignmentRecord[MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP] = {{{0}}};
   int uniqueAlignmentCount=0;

   *mutatedAlignmentCountPtr = 0;
   *totalNumberOfMutatedAlignments = 0;

   hts_itr_t *iter = NULL;

   iter  = sam_itr_querys(idx, header, targetRegion); // parse a region in the format like `chr4:500-2000'
   if(iter==NULL) {
      fprintf(stderr, "Couldn't create iterator for target %s...",targetRegion);
      exit(-3);
   }

   // Initialise your alignment
   // (ie., get a location where you're going to store the read/alignment returned by the iterator)
   bam1_t *b = bam_init1();

   // Debug
   printf("STARTED....\n");

   // Iterate through all reads within this region.
   // This will return >= 0 so long are there are still leads left in the pileup.
   while(sam_itr_next(in, iter, b) >= 0) {

      //printf("\nFMALIST<<<<%s>>>>\n",bam_get_qname(b));

      // Filter out unwanted reads.
      // You may also want to include !BAM_FPROPER_PAIR below also / plus other stuff..
      // (although all sim. data, before it is realigned, will be proper read pairs)
      // If any of these are detected inn the alignment flags we skip onto the next one.
      if(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;

      if(alignmentCount >= MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP) {
         fprintf(stderr, "Number of alignments in pileup exceeds maximum allowable (%d). Consider increasing MAX_NUM_ALIGNMENTS_IN_PILEUP in the source and recompiling this binary.\n", MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP);
      }

      // Copy alignment qName.
      strcpy(alignmentRecord[alignmentCount++].qName,bam_get_qname(b));

      // Debug
      // printf("\nMoving onto next read...\n");
   }



   // Compile a unique list of qNames (where overlapping reads appear only once..).

   int i=0;
   while(i<alignmentCount) {

      // Check for alignment and add the alignment name if its mate is not in there already.
      if(!qNameFound(uniqueAlignmentRecord, alignmentRecord[i].qName, uniqueAlignmentCount))
         strcpy(uniqueAlignmentRecord[uniqueAlignmentCount++].qName,alignmentRecord[i].qName);
      else {
         // Debug
         // printf("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%s OVERLAP FOUND!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n",alignmentRecord[i].qName);
      }

      i++;
   }

   // Now select which of these alignments are to be mutated.
   i=0;
   while(i<uniqueAlignmentCount) {
      // Toss a loaded coin with p=malignant allele frequency.
      // If it's malignant DNA, note it in the mutatedAlignmentRecord.

      if(coinToss(mutantAlleleFreq)) {
         strcpy(mutatedAlignmentRecord[*mutatedAlignmentCountPtr].qName,uniqueAlignmentRecord[i].qName);
         (*mutatedAlignmentCountPtr)++;
      }
      i++;
   }


   // Finally find the total number of alignments, including mates in overlap, that we need to mutate.
   // We use this number to tell when we are done spiking in this mutation.
   i=0;
   *totalNumberOfMutatedAlignments=0;
   while(i<*mutatedAlignmentCountPtr) {
      int j=0;
      while(j<alignmentCount) {
         if(!strcmp(mutatedAlignmentRecord[i].qName,alignmentRecord[j].qName)) {
            (*totalNumberOfMutatedAlignments)++;
            // Debug
            // printf("%s FMALIST\n",alignmentRecord[j].qName);
         }
         j++;
      }
      i++;
   }

   // Cleanup.
   hts_itr_destroy(iter);
   bam_destroy1(b);

   // Return the coverage in this region.
   return uniqueAlignmentCount;
}



typedef struct {
   char *chrom;    // The target chromosome name as per BAM header.
   long targetLocus; // The target locus.
   char *insertionSequence; // The DNA sequence to be inserted (null if this is a deletion only)
   long deletionLength; // Length of the deletion (zero if this is an insertion only).
   float mutantAlleleFreq; // The true allele frequency of mutant allele.
} IndelSubRecord;


void getNextTarget(IndelSubRecord *target, FILE * configFileP)
{
   static char *recordPtr = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   // Set targetLocus to zero to indicate no more targets if there is
   // nothing more to read from the config file.
   target->targetLocus = 0;

   while((read = getline(&recordPtr, &len, configFileP)) != -1) {
      // Skip metadata.
      if(recordPtr[0] == '#')
         continue;
      if(!strlen(recordPtr))
         continue;

      // Parse the record.

      // CHROM.
      pch = strtok(recordPtr,"\t");
      if(pch != NULL) {
         int i = snprintf(target->chrom, MAX_RNAME_LENGTH, "%s", pch);

         if(i != strlen(pch)) {
            fprintf(stderr,"Chromosome name exceeds max allowed length (%d>%d).\nConsider increasing MAX_RNAME_LENGTH in the source and recompiling this binary.\n",i,MAX_RNAME_LENGTH);
            exit(-1);
         }
      }
      else continue;

      // TARGET LOCUS.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->targetLocus = atol(pch);
      }
      else continue;

      // INSERTION SEQUENCE.
      pch = strtok(NULL,"\t");

      if(pch != NULL) {
         // printf("%s\n",pch);
         if(pch[0] == '-') {
            // A '-' indicates there's nothing to insert, just terminate the insertion string.
            target->insertionSequence[0]='\0';
         }
         else {
            int i = snprintf(target->insertionSequence, MAX_INSERTION_LENGTH, "%s", pch);

            if(i != strlen(pch)) {
               fprintf(stderr,"Insertion sequence exceeds max allowed length (%d>%d).\nConsider increasing MAX_INSERTION_LENGTH in the source and recompiling this binary.\n",i,MAX_INSERTION_LENGTH);
               exit(-1);
            }


         }
      }
      else continue;

      // DELETION LENGTH.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->deletionLength = atol(pch);
      }
      else continue;

      // MUTANT AF; True underlying allele frequency of mutant allele.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->mutantAlleleFreq = atof(pch);
         // printf("\nSETTING!!!! %s, %ld, %c, %f\n", target->contigPtr, target->locus, target->base, target->mutFreq);
         return;
      }
      else continue;
   }

   //free(recordPtr);
   // If we've got to the end of the config file,
   // set the config to an empty string to flag the caller to ignore.
   if(read == -1) {
      target->chrom[0] = '\0';
      return;
   }

}


void setUpNextTarget(FILE *configFilePtr,
                     IndelSubRecord *indelSubRecordPtr,
                     long *readReconstructOffsetPtr,
                     char *targetRegionPtr)
{

   // Get the genomic coordinate of the first somatic target.
   getNextTarget(indelSubRecordPtr, configFilePtr);

   if(!indelSubRecordPtr->targetLocus) {
      strcpy(targetRegionPtr,"");
      *readReconstructOffsetPtr=0;
      printf(">>>>>>>>>>>>>>>NO MORE TARGETS<<<<<<<<<<<<<<<<<<<<<<<\n");
      return;
   }

   // Set up the mutation target region.
   // This identifies the reads that cover, or partly cover the deleted section of DNA.
   // Also set up the read reconstruction offset used to reconstruct the remaining
   // section of the read after the mutation has been inserted.

   if(strlen(indelSubRecordPtr->insertionSequence) && !indelSubRecordPtr->deletionLength) {
      sprintf(targetRegionPtr,"%s:%ld-%ld",indelSubRecordPtr->chrom, indelSubRecordPtr->targetLocus, indelSubRecordPtr->targetLocus);
      *readReconstructOffsetPtr = -strlen(indelSubRecordPtr->insertionSequence);
      // Debug
      // printf("MARK URHERE A...*readReconstructOffsetPtr=%ld, indelSubRecordPtr->insertionSequence = >%s<\n", *readReconstructOffsetPtr, indelSubRecordPtr->insertionSequence);
   }
   else if(strlen(indelSubRecordPtr->insertionSequence) && indelSubRecordPtr->deletionLength) {
      sprintf(targetRegionPtr,"%s:%ld-%ld",indelSubRecordPtr->chrom, indelSubRecordPtr->targetLocus, indelSubRecordPtr->targetLocus+indelSubRecordPtr->deletionLength);
      *readReconstructOffsetPtr = (indelSubRecordPtr->deletionLength-strlen(indelSubRecordPtr->insertionSequence));
      // Debug
      // printf("MARK URHERE B...*readReconstructOffsetPtr=%ld, indelSubRecordPtr->insertionSequence = >%s<\n", *readReconstructOffsetPtr, indelSubRecordPtr->insertionSequence);
   }
   else if(!strlen(indelSubRecordPtr->insertionSequence) && indelSubRecordPtr->deletionLength) {
      sprintf(targetRegionPtr,"%s:%ld-%ld",indelSubRecordPtr->chrom, indelSubRecordPtr->targetLocus, indelSubRecordPtr->targetLocus+indelSubRecordPtr->deletionLength);
      *readReconstructOffsetPtr = indelSubRecordPtr->deletionLength;
      // Debug
      // printf("MARK URHERE C...*readReconstructOffsetPtr=%ld, indelSubRecordPtr->insertionSequence = >%s<\n", *readReconstructOffsetPtr, indelSubRecordPtr->insertionSequence);
   }
   else {
      fprintf(stderr,"Error in config file at %s:%ld\n",indelSubRecordPtr->chrom, indelSubRecordPtr->targetLocus);
      exit(-1);
   }
}

// Here we go...
int main(int argc,char* argv[])
{

   // Strip path.
   // TODO!!!! Should really hardcode the app name, not take from the binary name.
   char *commandName;
   commandName = strrchr(argv[0],'/');

   if(!commandName)
      commandName=argv[0];
   else
      commandName++;

   if(argc!= 3) {
      fprintf(stderr, "\n Usage: %s <donor BAM> <donor reference> <genomic rearrangement config file> <seed>\n\n",commandName);
      exit(0);
   }


   // NB: Rem. to use static (as we don't want an array of this size ending up on the stack..)
   static AlignmentRecord mutatedAlignmentRecord[MAX_NUM_ALIGNMENTS_IN_TARGET_PILEUP] = {{{0}}};
   int mutatedAlignmentCount = 0;
   int totalNumberOfMutatedAlignments = 0;
   int numberOfAlignmentsMutatedSoFar = 0;



   // Read in mutation record

   IndelSubRecord indelSubRecord;

   char tmpChrom[MAX_RNAME_LENGTH];
   char tmpInsert[MAX_INSERTION_LENGTH];
   indelSubRecord.chrom=tmpChrom;
   indelSubRecord.insertionSequence = tmpInsert;

   // Load the config file containing mutated loci and their allele frequencies.
   FILE * configFileP = NULL;
   configFileP = fopen(argv[3], "r");
   if(configFileP == NULL) {
      fprintf(stderr, "\nCan't open mutation config file, %s..\n\n",argv[3]);
      return -1;
   }

   // For this simple example we're taking that the region to iterate through is supplied to the second
   // argument to this binary when its called.
   // (The first arguement argv[1], is the BAM file etc..)
   char targetRegion[2*MAX_RNAME_LENGTH];
   long readReconstructOffset;
   char *insertionPointer = indelSubRecord.insertionSequence;


   /* Intializes random number generator */
   srand((unsigned)atoi(argv[4]));

   // Open two sets of file descriptors for the input BAM.

   // The read iterator (sam_itr_querys) uses the first descriptor set to select the
   // alignments within the target region that will undergo mutation.
   // Subsequently, utilizing the second file descriptor set, we proceed to traverse the pileup,
   // altering the read contents one locus at a time (bam_mplp_init) to incorporate the mutation.


   // Mutant reads list, 1st file descriptor set
   samFile *mutReadsFd = NULL;
   mutReadsFd = sam_open(argv[1], "r");
   if(mutReadsFd == 0) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   // Load the header in that BAM and its index.
   bam_hdr_t *mutReadsHdr = NULL;
   mutReadsHdr = sam_hdr_read(mutReadsFd);
   if(mutReadsHdr == 0) {
      fprintf(stderr, "\nCan't open header in %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   hts_idx_t *mutReadsIdx = sam_index_load(mutReadsFd,  argv[1]);
   if(mutReadsIdx == 0) {
      fprintf(stderr, "\nCan't load index for %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   // Pileup across mutation region, 2nd file descriptor set

   // First arg. is the SAM/BAM file of interest.
   // We are using the first command line arg to give the BAM file (plus path if required).
   samFile *mutRegionFd = NULL;
   mutRegionFd = sam_open(argv[1], "r");
   if(mutRegionFd == 0) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   // After that load the header in that BAM and its index.
   bam_hdr_t *mutRegionHdr = NULL;
   mutRegionHdr = sam_hdr_read(mutRegionFd);
   if(mutRegionHdr == 0) {
      fprintf(stderr, "\nCan't open header in %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   hts_idx_t *mutRegionIdx = NULL;
   mutRegionIdx = sam_index_load(mutRegionFd,  argv[1]);
   if(mutRegionIdx == 0) {
      fprintf(stderr, "\nCan't load index for %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }


   /////////////////////////////////////////////////////////////////////////////////



   mutatedAlignmentCount = 0;
   totalNumberOfMutatedAlignments = 0;
   numberOfAlignmentsMutatedSoFar = 0;
   // TODO!!!! load up next target mutation

   // Load the next target that has coverage within the target region,
   // or skip if we are out of targets (indelSubRecord.targetLocus == 0).
   // TODO!!!! we need to handle the coverage = 0 case here......
   int coverageAtMutatedRegion = 0;

   // Setup first target mutation.
   // Go through the list of targets until you find one with coverage...
   while(coverageAtMutatedRegion == 0) {
      numberOfAlignmentsMutatedSoFar = 0;

      // Get next target mutation from config file.
      setUpNextTarget(configFileP,
                      &indelSubRecord,
                      &readReconstructOffset,
                      targetRegion);

      // While there is a remaining target, check the coverage at this target.
      if(indelSubRecord.targetLocus) {
         // Check the coverage at this target, selecting mutant alignments if coverage > 0
         coverageAtMutatedRegion = selectMutatedQnames(targetRegion,
                                   mutReadsFd,
                                   mutReadsIdx,
                                   mutReadsHdr,
                                   mutatedAlignmentRecord,
                                   &mutatedAlignmentCount,
                                   &totalNumberOfMutatedAlignments,
                                   indelSubRecord.mutantAlleleFreq);
      }

      //Debug
      if(!indelSubRecord.targetLocus)
         fprintf(stderr,"No targets with coverage..\n");

   }

   ////////////////////////////////////////////////////////////////////////////////

   // Create a read iterator.
   // This is basically a function you will use to cycle through the pileup from 5' to 3',
   // returning the value of the base recorded in each read as it hits it.
   // The procedure is basically you create an iterator and then you call it within a while loop.
   // With each call it will return the read, you do stuff with the contents of that read and then call it
   // again to get the next one. When there are no more reads the iterator will return a value less than zero
   // and you're out of your loop.





   //####################################################

   // Set min mapping quality.
   // Any reads below this are dropped.
   int mapQ = 30;

   int tid, pos, *n_plp, min_len = 0;
   int previousPos=0;
   int ret;

   bam_mplp_t mplp;
   const bam_pileup1_t **plp;
   aux_t **data;

   // Init reference.
   mplp_ref_t mp_ref = MPLP_REF_INIT;

   // Allocate data pointer for our instance of our aux_t structure described above.
   // Could have written this simpler and avoided alloc but I will leave it like this as
   // it will be useful if you ever deal with a number of BAMs (that is decided only at run time).
   data = calloc(NUM_BAMS, sizeof(aux_t*)); // data[i] for the i-th input

   // ################ Initialise BAM structures ######################
   data[0] = calloc(1, sizeof(aux_t));
   data[0]->fpIn = mutRegionFd; // the input BAM we opened earlier

   data[0]->min_mapQ = mapQ;                    // set the mapQ filter
   data[0]->min_len  = min_len;                 // set the qlen filter


   data[0]->hdr = mutRegionHdr;    // the BAM header we opened earlier

   // Create the output SAM filename and record this in data...
   char *outFnamePtr = malloc(strlen(argv[1])+strlen(".indel_spike.sam"));
   strcpy(outFnamePtr,argv[1]);
   char *tmpP;
   tmpP=outFnamePtr;
   tmpP+=strlen(argv[1])-4;
   strcpy(tmpP,".indel_spike.sam");
   printf("Output file containing spike-in store as %s.\n",outFnamePtr);

   sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

   data[0]->fpOut = sam_open_format(outFnamePtr, "w", &ga.out);// open BAM output.


   // Prepare the output BAM.
   // Write out header.
   // TODO!!!! add a @PG for this app to the header!!!!
   // See sam_hdr_add_pg() in sam_view.c
   if(sam_hdr_write(data[0]->fpOut, data[0]->hdr) != 0) {
      fprintf(stderr, "failed to write the tumour SAM header\n");
   }


   // Read sample name.
   // TODO!!!! should we free smField when done? prob not worth
   kstring_t smField = { 0, 0, NULL };

   // Init sample name. This will be changed to whats in the read group if it's there.
   char noName[] = "SAMPLE";

   // Init to default in case there is no name found.
   data[0]->sampName = noName;
   // In case more than one @RG with a given SM:..
   int nRg = sam_hdr_count_lines(data[0]->hdr, "RG");
   if(nRg < 0) {
      fprintf(stderr, "Couldn't parse BAM header");
      exit(-6);
   }

   // Iterate through read groups.
   for(int i = 0; i < nRg; i++) {
      if(sam_hdr_find_tag_pos(data[0]->hdr, "RG", 0, "SM", &smField) <= -1) {
         fprintf(stderr, "Couldn't read sample name from header...\n");
         exit(-1);
      }
   }

   // Record name.
   if(smField.s)
      data[0]->sampName = smField.s;


   //Debug,, URHERE commenting out for test....
   //   if((data[0]->iter=sam_itr_querys(mutRegionIdx, data[0]->hdr, "chr1")) == 0) {
   //      fprintf(stderr, "Failed to parse region 'NULL'\n");
   //      exit(EXIT_FAILURE);
   //   }

   // Go through the entire BAM.
   data[0]->iter=NULL;

   // Load the reference (get it from argv[2] for now)
   data[0]->fai = fai_load(argv[2]);
   if(data[0]->fai==NULL) {
      fprintf(stderr,"Could not load faidx: %s\n", argv[2]);
      exit(EXIT_FAILURE);
   }

   // Set this up so we can cache reference chunks (in mplp_get_ref() above).
   data[0]->ref = &mp_ref;

   // ################ end initialisation ######################


   // the core multi-pileup loop
   mplp = bam_mplp_init(NUM_BAMS, read_bam, (void**)data); // initialization

   // Put in overlapping read pair detection to rework base quality based on mate,
   // (more like samtools mpileup does it).
   // For each overlapp base in pair this resets the base quality of the
   // lower-quality effectively discarding it from calling and recalculates the other base qual.
   // If you want this, leave this line in, otherwise comment it out...
   // bam_mplp_init_overlaps(mplp);

   // Set max depth to 10000..(ie., whatever MAX_PILEUP_SIZE has been defined as above).
   bam_mplp_set_maxcnt(mplp,MAX_PILEUP_SIZE);


   n_plp = calloc(NUM_BAMS, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
   plp = calloc(NUM_BAMS, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)



   // ##################################################################################################
   // To move on to the next covered position.
   // Locations with zero depths are skipped.
   // It returns when it hits something in one or more of the BAMs.

   while((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) {
      // tid is the current target ID (header @SQ lines, contig/chromosome etc.).
      if(tid >= sam_hdr_nref(data[0]->hdr)) continue;      // diff number of @SQ lines per file?

      // Align previous marker with where coverage starts.
      if(!previousPos)
         previousPos=pos-1;

      // Debug
      // printf("MARK A, pos=%d, mutant bases left >%s<\n",pos,insertionPointer);

      hts_pos_t ref_len;
      char *ref;
      mplp_get_ref(data[0], tid, &ref, &ref_len);

      // For each BAM we have..
      for(int i = 0; i < NUM_BAMS; ++i) {
         // Run through the pileup in this BAM.
         for(int j = 0; j < n_plp[i]; ++j) {
            // So we are at the pileup in the ith BAM.
            // How far into that pileup are we?
            // We need to keep track of that ourselves with j above.
            // Again, the size of the pileup (number of bases stacking at this locus)
            // is stored for us in n_plp[i] by htslib.
            // A little C aside here,
            // Remember we had to pass the addresses (with the '&' symbol) of where we wanted htslib to
            // store tid and pos for us. However we didn't need to use the '&' symbol with n_plp and plp
            // as these are already defined as pointers which we have assigned to memory which
            // we allocated for them.
            const bam_pileup1_t *p = plp[i] + j;

            // Do we need to mutate this alignment?
            // Is the qname in the mutate list and have we passed the mutation start locus?


            // Have we hit the target we want to mutate?
            // Rem., we need to convert the target locus to zero base before checking below.
            // Mutate at this locus?

            if(indelSubRecord.targetLocus && (pos+1 >= indelSubRecord.targetLocus)) {
               // And on this chromosome?
               if(!strcmp(indelSubRecord.chrom,sam_hdr_tid2name(data[0]->hdr, tid))) {

                  // Get what the reference base is at this position using the mplp_get_ref() function
                  hts_pos_t ref_len;
                  char *ref;
                  mplp_get_ref(data[0], tid, &ref, &ref_len);

                  // Is there anything we need to mutate in this alignment?
                  if(qNameFound(mutatedAlignmentRecord, bam_get_qname(p->b), mutatedAlignmentCount)) {

                     // Are there are still some bases we need to insert here?
                     if(*insertionPointer!='\0') {
                        // Check first to see if the base at this position in the
                        // read contains a sequence error.

                        if(seq_nt16_str[bam_seqi(bam_get_seq(p->b),p->qpos)] != ref[pos]) {

                           // Does this sequence error correspond to the base we were trying to insert?
                           // If so change the base to something else.
                           if(seq_nt16_str[bam_seqi(bam_get_seq(p->b),p->qpos)] == *insertionPointer) {
                              // Pick a base other than the mutant allele.
                              char differentBase = selectMutantAllele(*insertionPointer);
                              // Update the read base to something other than the mutant allele.
                              bam_seqi_set(bam_get_seq(p->b),
                                           p->qpos,
                                           seq_nt16_table[(int)differentBase]);
                           }
                        }
                        else {
                           // Update the base to the mutant allele.
                           bam_seqi_set(bam_get_seq(p->b),
                                        p->qpos,
                                        seq_nt16_table[(int)*insertionPointer]);
                        }
                     }
                     else {
                        // If there were any bases to insert, and if we hit this, we know they are now in place.
                        // If the genomic alteration resulted in a net insertion or deletion of bases
                        // then we need to replace the remainder of the contents of the read from the reference,
                        // allowing for any bases that were deleted or inserted.
                        // Otherwise we can leave the remaing bases in the read as they are.

                        // Is there a net insertion or deletion of bases?
                        if(readReconstructOffset) {
                           // Check first to see if the base at this position in the
                           // read contains a sequence error.
                           if(seq_nt16_str[bam_seqi(bam_get_seq(p->b),p->qpos)] != ref[pos]) {
                              // Does this sequence error correspond to the base we were trying to insert?
                              // If so change the base to something else.
                              if(seq_nt16_str[bam_seqi(bam_get_seq(p->b),p->qpos)] ==
                                    ref[pos+readReconstructOffset]) {
                                 // Pick a base other than the mutant allele.
                                 char differentBase = selectMutantAllele(ref[pos+readReconstructOffset]);
                                 // Update the read base to a different allele.
                                 bam_seqi_set(bam_get_seq(p->b),
                                              p->qpos,
                                              seq_nt16_table[(int)differentBase]);
                              }
                           }

                           else {
                              // Replace with contents from reference, accounting for offset.
                              bam_seqi_set(bam_get_seq(p->b),
                                           p->qpos,
                                           seq_nt16_table[(int)ref[pos+readReconstructOffset]]);
                           }
                        }
                     }
                  }
               }
            }


            int readEndLocus = p->b->core.pos + bam_cigar2rlen(p->b->core.n_cigar, bam_get_cigar(p->b));

            // Are we at this alignment's last appearance in a pileup?
            if(pos == readEndLocus - 1) {
               // If so time to write it out.
               writeOutAlignment(data[0], p->b);

               // If we've just written out a mutated alignment,
               // decrement our total and reset mutatedAlignmentCount, load up next target etc.
               // when we're done.

               // Check to see if we are within a valid target region and the alignment in question
               // was on the mutate list.
               if(indelSubRecord.targetLocus &&
                     (pos+1 >= indelSubRecord.targetLocus) &&
                     !strcmp(indelSubRecord.chrom,sam_hdr_tid2name(data[0]->hdr, tid)) &&
                     qNameFound(mutatedAlignmentRecord, bam_get_qname(p->b), mutatedAlignmentCount)) {

                  numberOfAlignmentsMutatedSoFar++;

                  // Debug
                  //printf("%s mutated alignment written out..\n",bam_get_qname(p->b));
                  // Debug
                  // printf("TRACK URHERE B: *insertionPointer= >%s<\n",insertionPointer);

                  // If we have modified all the reads containing this mutation then reset our
                  // counts and load up the next target mutation.
                  if(numberOfAlignmentsMutatedSoFar == totalNumberOfMutatedAlignments) {

                     // Debug
                     // printf("numberOfAlignmentsMutatedSoFar = %d,  totalNumberOfMutatedAlignments = %d  LAST ONE!!!!!!!!!!!!!!\n",numberOfAlignmentsMutatedSoFar, totalNumberOfMutatedAlignments);
                     // We are finished with this target.
                     // Reset everything and load up next target.
                     numberOfAlignmentsMutatedSoFar = 0;
                     mutatedAlignmentCount = 0;
                     totalNumberOfMutatedAlignments = 0;
                     insertionPointer = indelSubRecord.insertionSequence;

                     int coverageAtMutatedRegion = 0;

                     // Look through the list of remaining targets until you find the first one with coverage...
                     while(coverageAtMutatedRegion == 0 && indelSubRecord.targetLocus) {

                        // Setup next target mutation.
                        setUpNextTarget(configFileP,
                                        &indelSubRecord,
                                        &readReconstructOffset,
                                        targetRegion);

                        // If there are remaining targets, check the coverage at the next target.
                        if(indelSubRecord.targetLocus) {
                           // Check the coverage at this target, selecting mutant alignments if coverage > 0
                           coverageAtMutatedRegion = selectMutatedQnames(targetRegion,
                                                     mutReadsFd,
                                                     mutReadsIdx,
                                                     mutReadsHdr,
                                                     mutatedAlignmentRecord,
                                                     &mutatedAlignmentCount,
                                                     &totalNumberOfMutatedAlignments,
                                                     indelSubRecord.mutantAlleleFreq);
                        }


                        //if(!indelSubRecord.targetLocus)
                        //   printf(">>>>>>>>>>>>>>>AT MARK, END OF CONFIG<<<<<<<<<<<<<<<<<<<<<<<\n");
                        //else
                        //   printf(">>>>>>>>>>>>>>>AT MARK, MORE IN CONFIG<<<<<<<<<<<<<<<<<<<<<<<%d,%ld,%s\n",pos,indelSubRecord.targetLocus,indelSubRecord.insertionSequence);

                     }



                     if(!indelSubRecord.targetLocus)
                        printf(">>>>>>>>>>>>>>>>>>>>DONE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
                  }
               }
            }
         }
         // Debug
         // printf("numberOfAlignmentsMutatedSoFar = %d,  totalNumberOfMutatedAlignments = %d\n",numberOfAlignmentsMutatedSoFar, totalNumberOfMutatedAlignments);

      }

      // Increment the insertion pointer if we have just handled a locus that required an insertion.
      // The difference, pos-previousPos will normally equal 1 if there is coverage across
      // the region of interest.
      // However working it out from (pos-previousPos) allows is to skip parts of the insertion
      // if there is no coverage to support it.

      // Make sure we have something to insert and are at or past the
      // target locus (chromosome & offset) before we increment.
      if(*insertionPointer!='\0' &&
            pos+1 >= indelSubRecord.targetLocus &&
            !strcmp(indelSubRecord.chrom,sam_hdr_tid2name(data[0]->hdr, tid))) {
         insertionPointer += (pos-previousPos);
      }


      // Update our previous reference position marker so we can tell how many reference
      // positions have been consumed by the next iteration.
      previousPos=pos;

   }

   printf("CLEANING UP..\n");

   // Clean up.
   for(int i = 0; i < NUM_BAMS && data[i]; ++i) {
      sam_hdr_destroy(data[i]->hdr);
      if(data[i]->fpIn) sam_close(data[i]->fpIn);
      if(data[i]->fpOut) sam_close(data[i]->fpOut);
      hts_itr_destroy(data[i]->iter);
      free(data[i]);
   }
   free(data);
   free(mp_ref.ref[0]);
   free(mp_ref.ref[1]);

   sam_close(mutReadsFd);

   free(outFnamePtr);

}


