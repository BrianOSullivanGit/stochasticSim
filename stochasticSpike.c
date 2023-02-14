#define VERSION "0.01"

// TODO!!!! you prob want rbeta to generate your config to test this out..
#define bam_seqi_set(s,i,b) ((s)[(i)>>1]=((s)[(i)>>1]&(0xf<<((i&1)<<2)))|((b)<<((~i&1)<<2)))


// This code extracts allele breakdown, coverage and base quality scores from each locus in a pileup.
// The info is extracted as it traverses the pileup one locus at a time (from 5' to 3' or left to right) and at each stop traversing (upwards) the bases stacking up on it.

// It's compiled with, on my machine,

// gcc -g -O2 -Wall -I. -I../../htslib -I./lz4 -L../../htslib -o somaticSim somaticSim.c libst.a ../../htslib/libhts.a -lz -lm -lbz2 -llzma -lcurl -lncurses -lm -lz  -lpthread

// On syd that's probably something like,
// gcc  -std=c99 -g -O2 -Wall -I. -Ihtslib-1.11 -Ilz4  -o simpleMplpExample simpleMplpExample.c libst.a htslib-1.11/libhts.a -lz -lm -lm -lz  -lpthread
// (I havent tried it there but any problems give me a shout)

// You run it with something like,
// ./baseQualsForAdib  /home/sully/BIO_INFORMATICS/DATA/TCGA/SAMPLES/BAMS/C509.TCGA-78-7537-01A-11D-2063-08.1_gdc_realn_hg38exome.bam chr1:21706893-21706904 /home/sully/BIO_INFORMATICS/REFERENCE/GRCh38.d1.vd1.fa

#define  _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "htslib/sam.h"
#include "samtools.h"
#include "sam_opts.h"
#include <htslib/faidx.h>

#include <stdbool.h>

#include <math.h>

// This programs takes 1 BAMs as input, if you want more,
// Adjust as required (here and below).
#define NUM_BAMS 1
#define MAX_PILEUP_SIZE 10000

enum FilterField {NONE=0, PASS, MASKED, MASKED_OVL, UNDETECTED, NO_COVERAGE};
const char* filterFieldNames[] = {"NONE", "PASS", "MASKED", "MASKED_OVL", "UNDETECTED", "NO_COVERAGE"};


// A place to store the pileup indices.
// Used to generate a random selection of the for mutating.
int pileupIndices[MAX_PILEUP_SIZE];

// Initial values for instance of a reference fa structure mplp_ref_t
#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

// Struct to hold info about reference (for caching purposes, see below)..
typedef struct {
   char *ref[2];
   int ref_id[2];
   hts_pos_t ref_len[2];
} mplp_ref_t;

// Struct to hold info about error alleles.
typedef struct {
   char allele;
   int count;
} ErrorAlleleDepths;


// Put all BAM file related stuff together in a structure.
// Below, declare an instance of this structure for every BAM file we open to keep track of stuff..
// Btw., It's declared in a structure like this so we can pack all this info together and
// pass it to functions like read_bam below nice and neatly in the first arg of the function call.
// I guess there would be other ways of doing this without defining a struct, but it would probably
// be messy and a nightmare to keep track of, particulary it you were dealing with more than one BAM.
typedef struct {     // auxiliary data structure
   samFile *fpIn;     // file handle for input
   samFile *fpOut;     // file handle for output
   FILE *cfgFile; // File containing target somatic mutations to spike in.
   FILE *vcfOutFile; // File containing ground truth of mutations that were spiked in.
   sam_hdr_t *hdr;  // the file header
   char *sampName;  // the sample name taken from SM field in read group.
   hts_itr_t *iter; // NULL if a region not specified

   // I've put reference stuff in here, seemed as good a place as any..
   mplp_ref_t *ref;  // Reference struct for caching (below).
   faidx_t *fai;  // Indexed reference file.

   // Any other business..
   int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

typedef struct {
   char *contigPtr; // Pointer to chromosome of next target
   int c_tid;    // The target identifier of that chromosome in the BAM header.
   long locus; // SNP locus
   char base; // Pointer to the new base that is going to replace the current base at the target.
   float mutFreq; // The true underlying allele frequency of mutant allele.
} NextTargetSnp;



void getNextTarget(NextTargetSnp *target, FILE * configFileP, sam_hdr_t *hdr)
{
   char * lineP = NULL;
   size_t len = 0;
   ssize_t read;
   char * pch;

   target->locus = 0;

   while((read = getline(&lineP, &len, configFileP)) != -1) {
      // Skip metadata.
      if(lineP[0] == '#')
         continue;
      if(!strlen(lineP))
         continue;

      // Parse the record.

      // CHROM.
      pch = strtok(lineP,"\t");
      if(pch != NULL) {
         // printf("%s\n",pch);
         strcpy(target->contigPtr,pch);
         // Pull out this chromsomes target identifier.
         target->c_tid = sam_hdr_name2tid(hdr, target->contigPtr);

      }
      else continue;

      // LOCUS.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->locus = atol(pch)-1; // Change to zero based.
      }
      else continue;

      // BASE.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->base = *pch;
      }
      else continue;

      // MUTANT AF; True underlying allele frequency of mutant allele.
      pch = strtok(NULL,"\t");
      if(pch != NULL) {
         target->mutFreq = atof(pch);
         // printf("\nSETTING!!!! %s, %ld, %c, %f\n", target->contigPtr, target->locus, target->base, target->mutFreq);
         return;
      }
      else continue;
   }

   // If we've got to the end of the config file,
   // set the config to an empty string to flag the caller to ignore.
   if(read == -1) {
      target->contigPtr[0] = '\0';
      return;
   }

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
// Btw.,  BAQ is not included (deliberatly).
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

// Shuffle the elements of the array pointed to by ptr.
// Length of the array is given by length.
void shuffle(int *ptr, int last)
{
   // DEBUG!!!! REMOVE THIS LINE!!!!
   //return;
   int currentLastElement=last;

   int r,tmp;

   while(currentLastElement) {
      // Pick a random index between array start (0) and its last element.
      r = randomNum(0,currentLastElement);

      // Swap the last element with the randomly selected element.
      tmp=ptr[r];
      ptr[r]=ptr[currentLastElement];
      ptr[currentLastElement]=tmp;

      // Decrement currentLastElement and go again until you reach the array start.
      currentLastElement--;
   }
}

// Returns the simulated result of a loaded coin toss.
// rand() returns a random number between 0 and RAND_MAX.
// Multiply RAND_MAX by whatever you want your coin probability to be.
// Is the random number (returned by rand()) less than this?
bool coinToss(double probability)
{
   return rand() < (probability*((double)RAND_MAX + 1.0));
}


char selectMutantAllele(char wildType)
{
   static char bases[] = "GCAT";

   //printf("URHERE_A wildType=%c\n",wildType);

   // Choose a base other than what's at this locus in the pileup at the moment.
   int i;

   // bases is a string so rem. to account for the '\0' at the end..
   int numOfBases=(sizeof(bases)/sizeof(bases[0])-1);

   do {

      i=randomNum(0,numOfBases-1);
   }
   while(bases[i] == wildType);

   // printf("URHERE_A i=%d size=%d mutant=%c\n",i,numOfBases,bases[i]);

   // Return the base that replaced the previous one at this locus.
   return(bases[i]);
}


int findOverlappingMate(const bam_pileup1_t *p, int currentLocationInPileup, const bam_pileup1_t * bottomOfPileup, int pileupSize)
{
   // Is this base part of a read pair overlap?
   // If so then we also need to locate (above us inn the pileup) and mutate the
   // corresponding base in its mate and mark that it has been processed
   // (set it's pilupIndices entry to -1).

   // If the insert length associated with this bases alignment is less than
   // the combined length of the alignment and its mate then we need to look for the
   // overlap.
   // Assume read and mate have same length (otherwise get MC tag and work it out
   // from that, but that would be overkill).


   // #################################################
   // NOTE!!! This assumes duplicates have already been marked and removed.

   // TODO!!!! should we put the following shortcut in to speed things up???
   //      if(abs(p->b->core.isize) >
   //         2*bam_cigar2rlen(p->b->core.n_cigar, bam_get_cigar(p->b))) return -1;

   // no overlap possible, unless some wild cigar
   //    if ( (p->b.core->mtid >= 0 && p->b->core.tid != p->b->core.mtid)
   //         || (llabs(p->b->core.isize) >= 2*p->b->core.l_qseq
   //         && p->b->core.mpos >= p->end) // for those wild cigars
   //       ) return -1;


   //  l_qseq is calculated from the total length of an alignment block on reading or from CIGAR
   // Take a short cut and return if no overlap possible
   // (next 10 lines lifted from stats.c)
#define READ_ORDER_FIRST 1
#define READ_ORDER_LAST 2
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)

   uint32_t order = (IS_READ1(p->b) ? READ_ORDER_FIRST : 0) + (IS_READ2(p->b) ? READ_ORDER_LAST : 0);
   if(!(p->b->core.flag & BAM_FPAIRED) ||
         (p->b->core.flag & BAM_FMUNMAP) ||
         (llabs(p->b->core.isize) >= 2*p->b->core.l_qseq) ||
         (order != READ_ORDER_FIRST && order != READ_ORDER_LAST)) {
      return -1;
   }


   // Run through the pileup to check if this read's mate is in there also.
   // TODO!!!! logic in this loop could be made a bit more succinct.
   for(int x = currentLocationInPileup+1; x < pileupSize; ++x) {

      const bam_pileup1_t *mate = bottomOfPileup + x;

      // Should we prefix this with something like,
      // if(llabs(node->b.core.isize) >= 2*node->b.core.l_qseq) continue;
      // to speed it up??
      // Possible mate here? if so look in more detail..
      if(p->b->core.isize == -(mate->b->core.isize) &&
            p->b->core.mpos == mate->b->core.pos &&
            p->b->core.pos == mate->b->core.mpos
        ) {
         // And just in case...
         if((p->b->core.flag&BAM_FPAIRED) && !(p->b->core.flag&BAM_FPROPER_PAIR)) continue;
         if((mate->b->core.flag&BAM_FPAIRED) && !(mate->b->core.flag&BAM_FPROPER_PAIR)) continue;

         long iStart;
         long iEnd;

         // Try to confirm the first and second segment in this template.
         // If you can't, its not a read pair overlap, move on.
         if(p->b->core.flag & BAM_FREAD1 && mate->b->core.flag & BAM_FREAD2) {
            iStart = p->b->core.pos;
            iEnd = mate->b->core.mpos + bam_cigar2rlen(mate->b->core.n_cigar, bam_get_cigar(mate->b));
            if(p->b->core.isize == (iEnd-iStart)) {
               // TODO!!!! remove this after more testing...
               if(strcmp(bam_get_qname(p->b), bam_get_qname(mate->b)))
                  printf("\nERROR: YOU MISTOOK ONE MATE 1: %s(%d) : %s(%d), p->b->core.isize = %ld\n", bam_get_qname(p->b),p->b->core.flag, bam_get_qname(mate->b),mate->b->core.flag, p->b->core.isize);
               // We've found it's overlapping mate.
               return x;
            }

         }
         else if(p->b->core.flag & BAM_FREAD2 && mate->b->core.flag & BAM_FREAD1) {
            iStart = mate->b->core.pos;
            iEnd = p->b->core.mpos + bam_cigar2rlen(p->b->core.n_cigar, bam_get_cigar(p->b));
            if(mate->b->core.isize == (iEnd-iStart)) {
               // TODO!!!! remove this after more testing...
               if(strcmp(bam_get_qname(p->b), bam_get_qname(mate->b)))
                  printf("\nERROR: YOU MISTOOK ONE MATE 2: %s : %s\n", bam_get_qname(p->b), bam_get_qname(mate->b));

               // We've found it's overlapping mate.
               return x;
            }
         }

         else {
            // TODO!!!! remove this after more testing...
            if(!strcmp(bam_get_qname(p->b), bam_get_qname(mate->b)))
               printf("\nERROR: YOU MISSED ONE MATE 3: %s : %s\n", bam_get_qname(p->b), bam_get_qname(mate->b));
         }
      }
   }
   return -1;

}

// Return base and BQ of read, and if its in read pair overlap at this locus,
// the base and BQ of ist mat at this locus also.
void getBaseWithRPOcheck(int pileupSize,
                         int currentLocationInPileup,
                         const bam_pileup1_t *bottomOfPileup,
                         char *readBase,
                         char *mateBase,
                         int *readBQ,
                         int *mateBQ,
                         int *mateIdx)
{
   *readBase = 0;
   *mateBase = 0;
   *readBQ = 0;
   *mateBQ = 0;
   *mateIdx = 0;

   int mIdx;

   // Pointer to current pileup location.
   const bam_pileup1_t *p = bottomOfPileup + currentLocationInPileup;

   // Read Base
   *readBase = seq_nt16_str[bam_seqi(bam_get_seq(p->b),p->qpos)];

   // Read BQ
   *readBQ = bam_get_qual(p->b)[p->qpos];

   // Is this reads mate somewhere above us in the pileup, overlapping at this locus?
   // If not we can move on and return, but if it is we need to handle this read pair overlap.
   mIdx=findOverlappingMate(p,currentLocationInPileup,bottomOfPileup,pileupSize);

   // Are we in read pair overlap?
   if(mIdx>=0) {

      if(pileupIndices[mIdx] >= 0) {

         // Mates position in the pileup, if in read pair overlap.
         *mateIdx = mIdx;

         const bam_pileup1_t *mate = bottomOfPileup + mIdx;

         // Mate Base
         *mateBase = seq_nt16_str[bam_seqi(bam_get_seq(mate->b),mate->qpos)];

         // Mate BQ
         *mateBQ = bam_get_qual(mate->b)[mate->qpos];
      }
   }
}



void incErrorAlleleDepths(ErrorAlleleDepths *errorAlleleDepthsPtr, char errorAllele, int totalNumPossibleAlleles)

{
   int i=0;
   while(i<totalNumPossibleAlleles) {

      if(errorAlleleDepthsPtr[i].allele == errorAllele) {
         errorAlleleDepthsPtr[i].count++;
         return;
      }
      i++;
   }
}


int getErrorAlleleDepth(ErrorAlleleDepths *errorAlleleDepthsPtr, char errorAllele, int totalNumPossibleAlleles)

{
   int i=0;
   while(i<totalNumPossibleAlleles) {
      if(errorAlleleDepthsPtr[i].allele == errorAllele) {
         return errorAlleleDepthsPtr[i].count;
      }
      i++;
   }
   // return error.
   return(-1);
}


bool errorAlleleInPileup(ErrorAlleleDepths *errorAlleleDepthsPtr, int totalNumPossibleAlleles)

{
   int i=0;
   while(i<totalNumPossibleAlleles) {
      if(errorAlleleDepthsPtr[i].count>0) return true;
      i++;
   }
   return false;
}

int totalDepthOfErrorAlleles(ErrorAlleleDepths *errorAlleleDepthsPtr, int totalNumPossibleAlleles)

{
   int i=0,total=0;
   while(i<totalNumPossibleAlleles) {
      total += errorAlleleDepthsPtr[i].count;
      i++;
   }
   return total;
}


void resetErrorAlleleDepths(ErrorAlleleDepths *errorAlleleDepthsPtr, int totalNumPossibleAlleles)

{
   int i=0;
   while(i<totalNumPossibleAlleles) {
      errorAlleleDepthsPtr[i].count=0;
      i++;
   }

}

void sortDescendErrorAlleleDepths(ErrorAlleleDepths *errorAlleleDepthsPtr, int *array, int size)

{
   for(int step = 0; step < size - 1; ++step) {

      bool swapped = false;

      for(int i = 0; i < size - step - 1; ++i) {
         if(errorAlleleDepthsPtr[array[i]].count < errorAlleleDepthsPtr[array[i + 1]].count) {

            int temp = array[i];
            array[i] = array[i + 1];
            array[i + 1] = temp;

            swapped = true;
         }
      }

      if(!swapped) {
         break;
      }

   }
}


void attemptToMutateBase(int pileupSize,
                         int currentLocationInPileup,
                         const bam_pileup1_t *bottomOfPileup,
                         char refBase,
                         char mutantAllele,
                         float mutFreq, // The true underlying allele frequency of mutant allele.
                         int *pileupIndices,
                         enum FilterField *filterField,
                         int *mutAlleleCnt,
                         int *refAlleleCnt,
                         int *errorAlleleCnt,
                         ErrorAlleleDepths *errorAlleleDepths, int numErrorAlleleTypes)
{


   //########################
   // Has this base in this read, in this pile up has been handled
   // before as part of an overlapping read pair?
   if(pileupIndices[currentLocationInPileup] < 0) return;



   // Does this base already contain a non-ref allele due to sequencing error?

   char readBase,mateBase;
   int readBQ,mateBQ, mateIdx;

   getBaseWithRPOcheck(pileupSize, currentLocationInPileup, bottomOfPileup, &readBase,&mateBase,&readBQ,&mateBQ,&mateIdx);

   //################################################
   // Ignore certain base calls, like 'N'.
   // This makes the out more similar to the likes of bcftools, and simplifies out validation process.
   // If anyone wants 'N' base calls in their VCF its easy to change this
   // but until it's requested, its probably better to leave it out.

   // Set the base quality of an 'N' base call to zero so it is handled as required below.
   if(mateBase == 'N')
      mateBQ = 0;
   if(readBase == 'N')
      readBQ = 0;

   //################################################
   char base = readBase;

   // Are we in read pair overlap?
   // Get a consensus value for the base if so.
   // For a much more comprehensive job see tweak_overlap_quality() in sam.c
   // This will serve our purposes here..
   if(mateBase)
      if(mateBase != readBase)
         if(mateBQ > readBQ)
            base = mateBase;


   // #####################################

   // If after all that we're still only left with an 'N' here,
   // mark that we have effectively no coverage from this read/read pair at this locus.
   if(base == 'N') {

      // Flag that this base and its mate (if in read pair overlap)
      // have now been handled.
      pileupIndices[currentLocationInPileup] = -1;
      if(mateBase) pileupIndices[mateIdx] = -1;

      return;
   }

   // #################################

   // Pointer to current pileup location.
   const bam_pileup1_t *p = bottomOfPileup + currentLocationInPileup;

   // Pointer to mate pileup location, if we're in read pair overlap here.
   const bam_pileup1_t *mate = NULL;
   if(mateBase) mate = bottomOfPileup + mateIdx;


   // See if we get a mutatant allele at this base.
   // Reads in this region come from cells with healthy donor DNA,
   // and also from malignant donor cells.
   // The probability that this read is from malignant cell DNA depends
   // on the frequency of the malignant allele.

   // Toss a loaded coin with p=malignant allele frequency.
   // If it's malignant DNA, then update with the mutatant allele as required.

   if(!coinToss(mutFreq)) {

      // OK, we didn't get any mutantant DNA in this read at this locus.
      // If there is a sequencing error note it, flag that this read has been handled and move on.

      // Does this base already contain a non-ref allele due to sequencing error?
      if(base == refBase)
         (*refAlleleCnt)++;
      else {
         // We are dealing with a sequencing error.
         // Record new alleles due to sequencing error at this locus..
         incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);

         // Flag that this error base and its mate (if in read pair overlap)
         // have already been handled.
         pileupIndices[currentLocationInPileup] = -1;
         if(mateBase) pileupIndices[mateIdx] = -1;
      }

      return;

   }

   // Otherwise, handle the mutant allele according to scenarios below.

   // Case 1: No read pair overlap, readBase is reference.
   if(!mateBase && readBase == refBase) {

      // No mate to worry about here.
      // Just update base to the mutant allele,
      // mark that it has been handled update count/flags and we're done.
      bam_seqi_set(bam_get_seq(p->b),
                   p->qpos,
                   seq_nt16_table[(int)mutantAllele]);

      // Update count.
      (*mutAlleleCnt)++;

      // Unless it has previously been set to something else,
      // set the filterField (for truth.vcf) to PASS.
      if(*filterField == UNDETECTED || *filterField == NO_COVERAGE)
         *filterField = PASS;

      // For completeness, flag we're done with this one.
      pileupIndices[currentLocationInPileup] = -1;
   }

   // Case 2: No read pair overlap, readBase is not reference.
   else if(!mateBase && readBase != refBase) {
      // Mutant allele has been masked by sequencing error.
      // Also, no mate to worry about here.


      // Unless it has previously been set to MASKED_OVL,
      // set the filterField (for truth.vcf) to MASKED.
      if(*filterField != MASKED_OVL)
         *filterField = MASKED;

      // Does the read base already contain the alt. allele we attempted to spike in?
      if(readBase == mutantAllele) {
         // If so we need to change this as the read simulator intended to put in a sequencing error here,
         // not the alt. base we attempted to spike in.

         // printf("\nDEBUG: readBase =%c, mutantAllele=%c", readBase, mutantAllele);

         // Pick a base other than the mutantAllele.
         char differentBase = selectMutantAllele(mutantAllele);

         // printf("\nDEBUG: differentBase =%c, mutantAllele=%c", differentBase, mutantAllele);

         // Update the read base to something other than the mutant allele.
         bam_seqi_set(bam_get_seq(p->b),
                      p->qpos,
                      seq_nt16_table[(int)differentBase]);
         
         // Update consensus base            
         base = differentBase;

      }

      // Record new alleles due to sequencing error at this locus..
      incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);


      // For completeness, flag we're done with this one.
      pileupIndices[currentLocationInPileup] = -1;
   }

   // Case 3: Read pair overlap, readBase and mateBase both reference.
   else if(mateBase && readBase == refBase && mateBase == refBase) {

      // First update read base to the mutant allele &
      // mark that it has been handled.
      bam_seqi_set(bam_get_seq(p->b),
                   p->qpos,
                   seq_nt16_table[(int)mutantAllele]);

      // For completeness, flag we're done with this alignment.
      pileupIndices[currentLocationInPileup] = -1;

      // Now handle the mate. Update mate base to the mutant allele &
      // mark that it has been handled.
      bam_seqi_set(bam_get_seq(mate->b),
                   mate->qpos,
                   seq_nt16_table[(int)mutantAllele]);

      // Flag we're done with this mate.
      pileupIndices[mateIdx] = -1;

      // Update count.
      (*mutAlleleCnt)++;

      // Unless it has previously been set to something else,
      // set the filterField (for truth.vcf) to PASS.
      if(*filterField == UNDETECTED || *filterField == NO_COVERAGE)
         *filterField = PASS;

   }


   // Case 4: Read pair overlap, readBase is reference, mateBase is not.
   else if(mateBase && readBase == refBase && mateBase != refBase) {

      // First update read base to the mutant allele &
      // mark that it has been handled.
      bam_seqi_set(bam_get_seq(p->b),
                   p->qpos,
                   seq_nt16_table[(int)mutantAllele]);

      // For completeness, flag we're done with this alignment.
      pileupIndices[currentLocationInPileup] = -1;

      // Mate base has been masked by a sequencing error.
      // Make sure that that sequencing error did not result in the required alt allele being inserted.
      // If it did change it to something else.
      // Otherwise leave it alone.

      // TODO!!!! URHERE!!!!
      // Does the mate base already contain the alt. allele we attempted to spike in?
      if(mateBase == mutantAllele) {
         // We need to change this as the read simulator intended to put in a sequencing error here,
         // not the alt. base we attempted to spike in.

         // Pick a base other than the mutantAllele.
         char differentBase = selectMutantAllele(mutantAllele);

         // Update the mate base to something other than the mutant allele.
         bam_seqi_set(bam_get_seq(mate->b),
                      mate->qpos,
                      seq_nt16_table[(int)differentBase]);

         // If the consensus base came from the mate then we need to update this too.
         if(base == mateBase)
            base = differentBase;
      }

      // Mark that it has been handled.
      pileupIndices[mateIdx] = -1;

      // Set VCF filter to MASKED_OVL
      *filterField = MASKED_OVL;

      // Look at the consensus base of both read and mate and decide which
      // count (errorAlleleCnt or mutAlleleCnt) we need to increment.
      if(base == refBase)
         (*mutAlleleCnt)++;
      else {
         // Record new alleles due to sequencing error at this locus..
         incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);
      }
   }

   // Case 5: Read pair overlap, mateBase is reference, readBase is not.
   else if(mateBase && readBase != refBase && mateBase == refBase) {
      // Pretty much the mirror image of Case 4.

      // First update mate base to the mutant allele &
      // mark that it has been handled.
      bam_seqi_set(bam_get_seq(mate->b),
                   mate->qpos,
                   seq_nt16_table[(int)mutantAllele]);

      // For completeness, flag we're done with this alignment.
      pileupIndices[mateIdx] = -1;

      // Read base has been masked by a sequencing error.
      // Make sure that that sequencing error did not result in the mutant alt allele being inserted.
      // If it did change it to something else.
      // Otherwise leave it alone.

      // Does the read base already contain the alt. allele we attempted to spike in?
      if(readBase == mutantAllele) {
         // We need to change this as the read simulator intended to put in a sequencing error here,
         // not the alt. base we attempted to spike in.

         // Pick a base other than the mutantAllele.
         char differentBase = selectMutantAllele(mutantAllele);

         // Update the read base to something other than the mutant allele.
         bam_seqi_set(bam_get_seq(p->b),
                      p->qpos,
                      seq_nt16_table[(int)differentBase]);

         // If the consensus base came from the read base then we need to update this too.
         if(base == readBase)
            base = differentBase;
      }

      // Flag read base has been handled.
      pileupIndices[currentLocationInPileup] = -1;

      // Set VCF filter to MASKED_OVL
      *filterField = MASKED_OVL;

      // Look at the consensus base of both read and mate and decide which
      // count (errorAlleleCnt or mutAlleleCnt) we need to increment.
      if(base == refBase)
         (*mutAlleleCnt)++;
      else {
         // Record new alleles due to sequencing error at this locus..
         incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);
      }

   }

   // Case 6: Read pair overlap, both readBase and mateBase are not reference.
   else if(mateBase && readBase != refBase && mateBase != refBase) {
      // Not much else to do other than flag the error.

      // Set VCF filter to MASKED_OVL
      *filterField = MASKED_OVL;

      // This mutant allele has been masked by a sequencing error.

      // Make sure that that sequencing error did not result in the required alt allele being inserted.
      // If it did change it to something else.
      // Otherwise leave it alone.

      // Does the read base already contain the alt. allele we attempted to spike in?
      if(readBase == mutantAllele) {
         // We need to change this as the read simulator intended to put in a sequencing error here,
         // not the alt. base we attempted to spike in.

         // Pick a base other than the mutantAllele.
         char differentBase = selectMutantAllele(mutantAllele);

         // Update the read base to something other than the mutant allele.
         bam_seqi_set(bam_get_seq(p->b),
                      p->qpos,
                      seq_nt16_table[(int)differentBase]);

         // If the consensus base came from the read base then we need to update this too.
         if(base == readBase)
            base = differentBase;
      }

      // Does the mate base already contain the alt. allele we attempted to spike in?
      if(mateBase == mutantAllele) {
         // We need to change this as the read simulator intended to put in a sequencing error here,
         // not the alt. base we attempted to spike in.

         // Pick a base other than the mutantAllele.
         char differentBase = selectMutantAllele(mutantAllele);

         // Update the mate base to something other than the mutant allele.
         bam_seqi_set(bam_get_seq(mate->b),
                      mate->qpos,
                      seq_nt16_table[(int)differentBase]);

         // If the consensus base came from the read base then we need to update this too.
         if(base == mateBase)
            base = differentBase;
      }

      // Record new alleles due to sequencing error at this locus..
      incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);

      // Flag that both have been handled.
      pileupIndices[mateIdx] = -1;
      pileupIndices[currentLocationInPileup] = -1;


   }

   // "And you may ask yourself, how did i get here?"
   else {
      fprintf(stderr, "How did I get here? %c:%c:%c", refBase,readBase,mateBase ? mateBase : '?');
   }
}





// Here we go...
int main(int argc,char* argv[])
{
   // int beg, end;
   // Set min mapping quality.
   // Any reads below this are dropped.
   int mapQ = 30;

   int tid, pos, *n_plp, min_len = 0;
   int ret;

   bam_mplp_t mplp;
   const bam_pileup1_t **plp;
   aux_t **data;

   // Init sample name. This will be changed to whats in the read group if it's there.
   char noName[] = "SAMPLE";

   // Strip path.
   // TODO!!!! I should really hardcode the app name, not take from the binary name.
   char *commandName;
   commandName = strrchr(argv[0],'/');

   if(!commandName)
      commandName=argv[0];
   else
      commandName++;

   if(argc!= 6) {
      fprintf(stderr, "\n Usage: %s <donor BAM> <donor reference> <somatic mutation config file> <seed> <output SAM filename>\n\n",commandName);
      exit(0);
   }

   // Open output vcf vile.
   FILE *vcfOutFile = fopen("truth.vcf", "w");


   /* Intializes random number generator */
   srand((unsigned)atoi(argv[4]));

   // Init reference.
   mplp_ref_t mp_ref = MPLP_REF_INIT;

   // Allocate data pointer for our instance of our aux_t structure described above.
   // Could have written this simpler and avoided alloc but I will leave it like this as
   // it will be useful if you ever deal with a number of BAMs (that is decided only at run time).
   data = calloc(NUM_BAMS, sizeof(aux_t*)); // data[i] for the i-th input

   // ################ Initialise BAM structures ######################
   data[0] = calloc(1, sizeof(aux_t));
   data[0]->fpIn =  sam_open(argv[1], "r"); // open BAM input

   if(data[0]->fpIn == NULL) {
      fprintf(stderr, "Couldn't open bam...\n");
      exit(1);
   }

   data[0]->min_mapQ = mapQ;                    // set the mapQ filter
   data[0]->min_len  = min_len;                 // set the qlen filter


   data[0]->hdr = sam_hdr_read(data[0]->fpIn);    // read the BAM header
   if(data[0]->hdr == NULL) {
      fprintf(stderr, "Couldn't read header...\n");
      exit(1);
   }

   sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
   // data[0]->fpOut = sam_open_format("-", "w", &ga.out);// open BAM output (stdout in this case, TODO!!!! change.)
   // Put sam output file name in command line arg instead.
   data[0]->fpOut = sam_open_format(argv[5], "w", &ga.out);// open BAM output.

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

   // Initialise VCF file with sample name etc.
   fprintf(vcfOutFile,"##fileformat=VCFv4.2\n");
   fprintf(vcfOutFile,"##FILTER=<ID=PASS,Description=\"All filters passed\">\n");

   fprintf(vcfOutFile,"##%sVersion=%s\n",commandName, VERSION);
   fprintf(vcfOutFile,"##%sCommand=",commandName);

   // Prints out args this command was called with.
   for(int i=1; i<argc; i++)
      if(i==1)
         fprintf(vcfOutFile,"%s",argv[i]);
      else
         fprintf(vcfOutFile," %s",argv[i]);

   fprintf(vcfOutFile,"\n");

   fprintf(vcfOutFile,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",data[0]->sampName);


   hts_idx_t *idx = sam_index_load(data[0]->fpIn,  argv[1]);
   if(idx == 0) {
      fprintf(stderr, "\nCan't load index for %s..\n\n",argv[1]);
      exit(EXIT_FAILURE);
   }

   // Load the reference.
   data[0]->fai = fai_load(argv[2]);
   if(data[0]->fai==NULL) {
      fprintf(stderr,"Could not load faidx: %s\n", argv[2]);
      exit(EXIT_FAILURE);
   }

   // Load the config file comtaining mutated loci and their allele frequencies.
   data[0]->cfgFile = fopen(argv[3], "r");
   if(data[0]->cfgFile == NULL) {
      fprintf(stderr, "\nCan't open %s..\n\n",argv[3]);
      return -1;
   }

   // TODO!!!! check for buffer overflow.
#define MAX_LINE_LEN 1024
   char targetContig[MAX_LINE_LEN] = "EMPTY";

   NextTargetSnp nextTargetRegion;


   // Set up a location to store the contig name/chromosome.
   nextTargetRegion.contigPtr = targetContig;

   // Get the genomic coordinate of the first somatic target.
   getNextTarget(&nextTargetRegion, data[0]->cfgFile, data[0]->hdr);


   // On second thoughs let's allow use to pass in a region here rather than going
   // through full BAM.
   // We'll pass in the region in the last (6th) arg to this binary..
   // We'll need an index in that case so, if we're specifying a region..

   // data[0]->iter = NULL;   // No region specified, we're working through the entire BAM


   // If we don't pass a region as a command line argument take the entire BAM as pileup region.
   if(argc==6) {
      data[0]->iter=NULL;
   }
   else {

      if((data[0]->iter=sam_itr_querys(idx, data[0]->hdr, argv[6])) == 0) {
         fprintf(stderr, "Failed to parse region '%s'\n",argv[6]);
         exit(EXIT_FAILURE);
      }
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
   //  bam_mplp_init_overlaps(mplp);

   // Set max depth to 10000..(ie., whatever MAX_PILEUP_SIZE has been defined as above).
   bam_mplp_set_maxcnt(mplp,MAX_PILEUP_SIZE);


   n_plp = calloc(NUM_BAMS, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
   plp = calloc(NUM_BAMS, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)


   // Info. stats., print to stdout at the end...
   long numberOfLociCovered = 0;
   long totalFoldCoverage = 0;
   long maxDepth = 0;
   long alignmentCount = 0;



   // ##################################################################################################
   // To move on to the next covered position.
   // Locations with zero depths are skipped.
   // It returns when it hits something in one or more of the BAMs.
   // When it hits the next locus with coverage you move on to the inner loop
   // and traverse up the pileup.

   while((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) {
      // tid is the current target ID (header @SQ lines, contig/chromosome etc.).
      if(tid >= sam_hdr_nref(data[0]->hdr)) continue;      // diff number of @SQ lines per file?

      // There are 2 for loops here.
      // The first one says "for every BAM file you have included..." (we have just one here for an example).
      // The next is "run through the pileup at this locus in this BAM.."
      // The inner loop is basically where all the action happens.
      // Everything else in this soure file is pretty much only paperwork to set that loop up
      // and get us to that point.
      //
      // It is likely, you will come back another day, needing to do a completely different job in the pileup
      // and the only lines you will need to change will be within this while loop.
      //
      // "n_plp" is an array.
      // The size of this array (ie., number of elements it contains) ==
      //                     the number of BAMs you want to run through (we only have one here).
      //
      // Each element in this array conatins an int. which gives you the size of the pileup at this locus in BAMi,
      // (or n_plp[i] equals the size of the pileup in BAMi at position pos)

      // "plp" is also an array.
      // It contains a pointer to a struct telling you all about stuff relating to a position in the pileup.
      // j is basically an offset which you increment as you work through the pileup one base
      // at a time upwards (ie through each stacked read). It ensures you're always looking at the
      // right pileup pointer when you want info about whats happing where you are at the
      // current position in the pileup.

      // Remember above when you called the iterator you created before with
      // bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) (in the while loop above)?
      // Well, in that call you also passed it the addresses of two ints (&tid, &pos).
      // Through the magic of the htslib API, every time we move onto the next locus and a new pileup at
      // that locus (ie., with each iteration of the outer while loop), htslib fills in the values of these
      // int variables to be the target id (tid, thats hts speak for chromosome, basically) and
      // your location in that chromosome (pos). Inside this loop you use those two ints to
      // orientate yourself as to where you are in the genome & then run through the pile at that locus
      // below and do stuff.
      // Remember each iteration of the outer while only returns at a locus where there is coverage.
      // It skips any interviening zero coverage regions. If you need to know about these you'll
      // have to write code to print out the regions of zero coverage yourself (I don't bother here,,
      // see mpileup or samtools depth source for how they go about it).

      // Start off by printing out where we are in the genome at this point.
      // After this is printed, the for loop(s) below will print out the pileup(s)
      // (plural if there's more than one BAM...which there is not in this simple example
      // but you know what I mean..)
      // Btw., I've left it as zero based, you can adjust accordingly...

      // Get what the reference base is at this position using the mplp_get_ref() function
      // that I hacked and pinched from bam_plcmd.c (samtools mpileup).
      hts_pos_t ref_len;
      char *ref;
      mplp_get_ref(data[0], tid, &ref, &ref_len);
      int depthAtThisLocus = 0;


      // Set up the VCF filter field in case we encounter a mutant allele at this locus.
      enum FilterField filterField = NONE;
      int refAlleleCnt=0;
      int errorAlleleCnt = 0;
      int mutAlleleCnt = 0;
      // Initialise mutant allele in case its required below.
      // Btw., ref[pos] is the wt allele (/doner allele)
      // that was at this at this locus prior to mutation.
      // The function returns mutantAllele, a random pick of one of the 3 other bases.
      // TODO!!!! take the one from the config file if there is one there.
      // Otherwise go sith a random pick.
      char mutantAllele = selectMutantAllele(ref[pos]);

      // Stores the counts of the various non-ref alleles encountered as we traverse up the pileup.
      ErrorAlleleDepths errorAlleleDepths[] = {{'G',0},{'C',0},{'A',0},{'T',0}};
      int numErrorAlleleTypes=sizeof(errorAlleleDepths)/sizeof(errorAlleleDepths[0]);
      int tmpIdx[sizeof(errorAlleleDepths)/sizeof(errorAlleleDepths[0])] = {0,1,2,3};

      // TODO!!!! relace this duplicate ref with numErrorAlleleTypes
      int sizeE=sizeof(errorAlleleDepths)/sizeof(errorAlleleDepths[0]);


      // For each BAM we have..
      // TODO!!!! remove multiple BAM handling..
      // It seemed like a good idea at the time but we will likely never use it..
      for(int i = 0; i < NUM_BAMS; ++i) {

         // Run through the pileup in this BAM.

         // Keep a tally..
         numberOfLociCovered++;
         // Keep track as to what the max depth of any given locus is.
         if(n_plp[i] > maxDepth)
            maxDepth=n_plp[i];

         // Check the next target from the current config record.
         // Do we need to mutate anything at this locus?

         bool mutateHere=false;

         // Mutate at this locus?
         if(pos==nextTargetRegion.locus) {
            // And on this chromosome?
            if(!strcmp(nextTargetRegion.contigPtr,sam_hdr_tid2name(data[0]->hdr, tid))) {
               // Flag it so we know to mutate as we work through the pileup.
               mutateHere = true;
               // There is an underlying somatic mutation at this locus.
               // It is as yet undetected.
               // Depending on coverage at this locus and chance, that may change
               // (as dealt with in the logic below).
               filterField = UNDETECTED;
            }
         }


         // Record indices of each entry in the pileup.
         // We will use these to track read pair overlaps as we work our way up the pileup.
         for(int j = 0; j < n_plp[i]; ++j) {
            pileupIndices[j] = j;
         }

         // Now work your way up the pileup.
         for(int j = 0; j < n_plp[i]; ++j) {

            const bam_pileup1_t *p = plp[i] + j;

            totalFoldCoverage++;

            // Skip dels or ref skips, we don't want to end up spiking something in there.
            // Also, pileupIndices[j] == -1 indicates this base in the pile up has been handled
            // before as part of an overlapping read pair.
            // Do not attempt to mutate anything. Move on.
            // Also, base quality must be non-zero.
            // Remember, even if we changed this to filter on a min BQ, read pair
            // overlap handling can still result in a base with BQ = 0.
            // Therefore we always need to check.
            // TODO!!!! We can remove the BQ==0 check now as we have implemented our read pair overlap handling?
            if(p->is_del || p->is_refskip || !bam_get_qual(p->b)[p->qpos] || (pileupIndices[j] == -1)) {
               // Just write out the alignment as is and continue.
               int readEndLocus = p->b->core.pos + bam_cigar2rlen(p->b->core.n_cigar, bam_get_cigar(p->b));

               // Are we at this alignment's last appearance in a pileup?
               if(pos == readEndLocus-1) {
                  // If so increment alignment counter.
                  alignmentCount++;

                  // Now write it out.
                  writeOutAlignment(data[0], p->b);

                  // For completeness, mark alignment has been handled.
                  pileupIndices[j] = -1;

               }
               continue;
            }

            // Do we attempt to mutate (as flagged above)?
            if(mutateHere) {
               // See if we have a mutation here..
               attemptToMutateBase(n_plp[i], j,
                                   plp[i], ref[pos], mutantAllele,nextTargetRegion.mutFreq,
                                   pileupIndices, &filterField, &mutAlleleCnt, &refAlleleCnt, &errorAlleleCnt, errorAlleleDepths, numErrorAlleleTypes);
            }
            else {
               // OK, no mutant allele here, just record any sequencing error
               // and move on.

               char readBase,mateBase;
               int readBQ,mateBQ, mateIdx;

               // Are we in read pair overlap at this locus?
               // If so get the details plus the consensus for the base at this locus.
               getBaseWithRPOcheck(n_plp[i], j, plp[i], &readBase,&mateBase,&readBQ,&mateBQ,&mateIdx);

               //################################################
               // Ignore certain base calls, like 'N'.
               // This makes the output more similar to the likes of bcftools, and simplifies our
               // validation process.
               // If anyone wants 'N' base calls in their VCF its easy to change this
               // but until it's requested, its probably better to leave it out.

               // Set the base quality of an 'N' base call to zero so it is handled as required below.
               if(mateBase == 'N')
                  mateBQ = 0;
               if(readBase == 'N')
                  readBQ = 0;

               //################################################

               char base = readBase;

               // Are we in read pair overlap?
               // Get a consensus value for the base if so.
               // For a much more comprehensive job see tweak_overlap_quality() in sam.c
               // This will serve our purposes here..
               if(mateBase)
                  if(mateBase != readBase)
                     if(mateBQ > readBQ)
                        base = mateBase;

               // If after all that we're still only left with an 'N' then we have
               // effectively no coverage from this read/read pair at this locus.
               // In this case just mark it has been handled and move on.
               if(base == 'N') {
                  pileupIndices[j] = -1;
                  if(mateBase) pileupIndices[mateIdx] = -1;
               }

               // Does this base already contain a non-ref allele due to sequencing error?
               else if(base == ref[pos])
                  // No, just the reference allele.
                  refAlleleCnt++;
               else {
                  // We are dealing with a sequencing error.
                  incErrorAlleleDepths(errorAlleleDepths, base, numErrorAlleleTypes);

                  // Flag that this error base and its mate (if in read pair overlap)
                  // have already been handled.
                  pileupIndices[j] = -1;
                  if(mateBase) pileupIndices[mateIdx] = -1;

               }

            }


            // See if we need to write any SAM records out and then move on.

            int readEndLocus = p->b->core.pos + bam_cigar2rlen(p->b->core.n_cigar, bam_get_cigar(p->b));

            // Are we at this alignment's last appearance in a pileup?
            if(pos == readEndLocus-1) {
               // If so increment alignment counter.
               alignmentCount++;

               // Now write it out...TODO!!!!
               // No, this belongs in the outer loop..??
               writeOutAlignment(data[0], p->b);

               // Btw., you can validate this section & alignmentCount above with,
               // samtools view  -q 30 -f 0x02 -F 0x400 /home/sully/BIO_INFORMATICS/DATA/TCGA/SAMPLES/BAMS/C509.TCGA-78-7537-01A-11D-2063-08.1_gdc_realn_hg38exome.bam chr22:0-25000000 | wc -l
               //
               // samtools mpileup -q 30 -Q 0 --no-BAQ  --incl-flags 0x02 --excl-flags 0x400 /home/sully/BIO_INFORMATICS/DATA/TCGA/SAMPLES/BAMS/C509.TCGA-78-7537-01A-11D-2063-08.1_gdc_realn_hg38exome.bam chr22:0-25000000

               // Or better still,

               // bcftools mpileup -f /home/sully/BIO_INFORMATICS/REFERENCE/GRCh38.d1.vd1.fa  --annotate FORMAT/AD -q 30 -Q 0 --no-BAQ  --incl-flags 0x02 --excl-flags 0x400 /home/sully/BIO_INFORMATICS/DATA/TCGA/SAMPLES/BAMS/C509.TCGA-78-7537-01A-11D-2063-08.1_gdc_realn_hg38exome.bam -r "chr22:15528239-15528239"| grep 15528239


               // ./somaticSim  /home/sully/BIO_INFORMATICS/DATA/TCGA/SAMPLES/BAMS/C509.TCGA-78-7537-01A-11D-2063-08.1_gdc_realn_hg38exome.bam /home/sully/BIO_INFORMATICS/REFERENCE/GRCh38.d1.vd1.fa tst.bak.cfg 934 out.sam chr22:0-25000000

               // egrep -v ERR ./truth.vcf > subtruth.vcf

               // bcftools mpileup -f /home/sully/BIO_INFORMATICS/REFERENCE/GRCh38.d1.vd1.fa  --annotate FORMAT/AD -q 30 -Q 0 --no-BAQ  --incl-flags 0x02 --excl-flags 0x400 out.bam -R ./subtruth.vcf|egrep -v '^#'

            }
         }
         // Keep track of the number reads we've processed.
         depthAtThisLocus++;

      }


      // We're finished with the pileup at this locus.
      // Close it off by writing the details out to the output VCF file
      // and move onto next locus.


      // Write the details of any mutant allele encountered to the ground truth, output VCF.
      // Note, ground truth vcf is at present built to handle one bam only as input..


      if(filterField != NONE && filterField != NO_COVERAGE) {
         int totalErrDepth = totalDepthOfErrorAlleles(errorAlleleDepths, sizeE);
         int totalDepth = refAlleleCnt+mutAlleleCnt+totalErrDepth;

         // Sort (ie., put most significant error alleles first).
         sortDescendErrorAlleleDepths(errorAlleleDepths,
                                      tmpIdx,
                                      sizeE);

         fprintf(vcfOutFile,"%s\t%d\t%s\t%c\t%c",
                 sam_hdr_tid2name(data[0]->hdr, tid),  // CHROM
                 pos+1, // POS (changed to 1 based)
                 ".", // ID
                 ref[pos], // REF
                 mutantAllele); // ALT

         if(errorAlleleInPileup(errorAlleleDepths, sizeE)) {
            fprintf(vcfOutFile,",");
            for(int i=0; i<sizeE; i++) { // ALT (error alleles)
               if(errorAlleleDepths[tmpIdx[i]].count == 0)
                  break;

               fprintf(vcfOutFile,"%c",errorAlleleDepths[tmpIdx[i]].allele);
               if((i+1)==sizeE || errorAlleleDepths[tmpIdx[i+1]].count == 0)
                  break;
               else
                  fprintf(vcfOutFile,",");
            }
         }

         // TODO!!!! the AF field below need to be changed.
         // We need to define a new info field,
         // ##INFO=<ID=TMAF,Number=A,Type=Float,Description="The true, underlying mutant Allele Frequency within the tumour, set as a parameter of the simulated distribution of mutant alleles.">
         //
         // Also, we should perhaps split up the allele depth, defining,
         // ##INFO=<ID=TMAD,Number=A,Type=Float,Description="The true, underlying mutant Allele depth.">
         // to discern mutant alleles due to error and thosee sue to spikein?

         fprintf(vcfOutFile,"\t%s\t%s\tDP=%d;AF=%.6g\tAD\t%d,%d",
                 ".", // QUAL
                 filterFieldNames[filterField], // FILTER
                 totalDepth,nextTargetRegion.mutFreq, // INFO
                 refAlleleCnt,mutAlleleCnt  // FORMAT
                );

         //               printf("\nURHERE!!!!: %c:%d, %c:%d, %c:%d, %c:%d\n",errorAlleleDepths[0].allele,errorAlleleDepths[0].count,
         //errorAlleleDepths[1].allele,errorAlleleDepths[1].count,
         //errorAlleleDepths[2].allele,errorAlleleDepths[2].count,
         //errorAlleleDepths[3].allele,errorAlleleDepths[3].count);

         if(errorAlleleInPileup(errorAlleleDepths, sizeE)) {
            fprintf(vcfOutFile,",");
            for(int i=0; i<sizeE; i++) { // FORMAT, AD
               if(errorAlleleDepths[tmpIdx[i]].count == 0)
                  break;

               fprintf(vcfOutFile,"%d",errorAlleleDepths[tmpIdx[i]].count);
               if((i+1)==sizeE || errorAlleleDepths[tmpIdx[i+1]].count == 0)
                  break;
               else
                  fprintf(vcfOutFile,",");
            }
         }

         fprintf(vcfOutFile,"\n");
      }

      else if(filterField == NO_COVERAGE) {

         if(pos!=nextTargetRegion.locus) {
            fprintf(vcfOutFile,"%s\t%ld\t%s\t%c\t%c\t%c\t%s\t%s\t%s\n",
                    nextTargetRegion.contigPtr,  // CHROM
                    nextTargetRegion.locus + 1, // POS (change to 1 based)
                    ".", // ID
                    '.', // REF
                    '.', // ALT
                    '.', // QUAL
                    filterFieldNames[NO_COVERAGE], // FILTER
                    ".", // INFO
                    "." // FORMAT
                   );
         }
      }


      // Else record errors if that's all there was at this locus.
      // Otherwise we've nothing for the VCF at this locus
      // (we do not bother printing out reference entries).
      else {
         if(errorAlleleInPileup(errorAlleleDepths, sizeE)) {

            fprintf(vcfOutFile,"%s\t%d\t%s\t%c\t",
                    sam_hdr_tid2name(data[0]->hdr, tid),  // CHROM
                    pos+1, // POS (changed to 1 based)
                    ".", // ID
                    ref[pos]); // REF

            // Sort (ie., put most significant error alleles first).
            sortDescendErrorAlleleDepths(errorAlleleDepths,
                                         tmpIdx,
                                         sizeE);

            for(int i=0; i<sizeE; i++) { // ALT
               if(errorAlleleDepths[tmpIdx[i]].count == 0)
                  break;

               fprintf(vcfOutFile,"%c",errorAlleleDepths[tmpIdx[i]].allele);
               if((i+1)==sizeE || errorAlleleDepths[tmpIdx[i+1]].count == 0)
                  break;
               else
                  fprintf(vcfOutFile,",");
            }


            int totalErrDepth = totalDepthOfErrorAlleles(errorAlleleDepths, sizeE);
            int totalDepth = refAlleleCnt+totalErrDepth;

            // NOTE: Because deletions and reference skips are skipped as we work through the pileup
            // (refer to "if(p->is_del || p->is_refskip...." at top of pileup look),
            // no entries in this VCF will appear for indels due to sequencing error.
            // I think we can live with that but if not it is then probably better to split this
            // up into a 2 stage process, the first records all sequencing errors in a seperate VCF,
            // the second spikes in and records all the details of how that went in a second VCF.
            // TODO!!!! Another possibility might be to put the code that prints out a SEQ_ERR entry in the VCF
            // for a indel as a result of sequencing error in the "if(p->is_del || p->is_refskip...."
            // statment above?
            // As it stands at present, if its not an indel in the germline VCF, then it is a sequencing
            // error (as we only spike in SSBs here, no indels).

            // TODO!!!! As a comparison, from looking at their code,
            // I think bamsurgeon will spike into refskips and non refs regardless, check this out.

            fprintf(vcfOutFile,"\t%s\t%s\tDP=%d;AF=%.6g\tAD\t%d",
                    ".", // QUAL
                    "SEQ_ERROR", // FILTER
                    totalDepth,(float)totalErrDepth/totalDepth,  // INFO
                    refAlleleCnt);  // FORMAT, AD


            fprintf(vcfOutFile,",");
            for(int i=0; i<sizeE; i++) { // FORMAT, AD
               if(errorAlleleDepths[tmpIdx[i]].count == 0)
                  break;

               fprintf(vcfOutFile,"%d",errorAlleleDepths[tmpIdx[i]].count);
               if((i+1)==sizeE || errorAlleleDepths[tmpIdx[i+1]].count == 0)
                  break;
               else
                  fprintf(vcfOutFile,",");
            }
            fprintf(vcfOutFile,"\n");
         }

         // ##############################################################




      }
      // One last thing before we move on to the next locus.
      // Check to see if we have passed out current somatic target.
      // This loop only returns for BAM regions with coverage.
      // We need to check in case our somatic target is within a region without coverage.
      // TODO!!!! we should really only do a '>=' check in one place here and load the next target..
      // TODO!!!! we need to ochange this so it checks if we've moved onoto the next chromosome..
      // Have we passed or already processed this locus?


      // Are we at end of config file (ie., no more somatic targets left)?
      // nextTargetRegion.contigPtr[0] will be zero in that case.
      // In that case just finish out the BAM, writing out the rest  of
      // the alignments as we come to them.
      if(nextTargetRegion.contigPtr[0] == '\0') continue;

      // Have we just finished processing a target in this iteration?
      // If so load up the next one.
      // TODO!!!! Update this to use either the tid / c_tid (preferably) or a strcmp call to check if you're on a target contig.
      // Don't mix (convert to usingtid/c_tid)
      if(pos==nextTargetRegion.locus &&
            !strcmp(nextTargetRegion.contigPtr,sam_hdr_tid2name(data[0]->hdr, tid))) {
         // OK, we've just processed that target, load the next one.

         // Load the genomic coordinates of the next somatic target.
         getNextTarget(&nextTargetRegion, data[0]->cfgFile, data[0]->hdr);

      }
      else
         // Have we skipped past the target due to no coverage at that locus??
         // TODO!!!!, new chromosome check!!!! this should be up top!!!!

         // BUG!!!!! TODO!!!!
         // This block should be within a while loop, not an if statement???
         // In case you've skipped more than one target due to no coverage etc...
         if((tid==nextTargetRegion.c_tid && pos>nextTargetRegion.locus) ||
               (tid>nextTargetRegion.c_tid)) {
            // OK, we ended up skipping that target as there was no coverage there.
            // Printout a notification and move onto the next one.
            if(pos!=nextTargetRegion.locus) {
               fprintf(vcfOutFile,"%s\t%ld\t%s\t%c\t%c\t%c\t%s\t%s\t%s\n",
                       nextTargetRegion.contigPtr,  // CHROM
                       nextTargetRegion.locus + 1, // POS (change to 1 based)
                       ".", // ID
                       '.', // REF
                       '.', // ALT
                       '.', // QUAL
                       filterFieldNames[NO_COVERAGE], // FILTER
                       ".", // INFO
                       "." // FORMAT
                      );
            }

            // Load the genomic coordinates of the next somatic target.
            getNextTarget(&nextTargetRegion, data[0]->cfgFile, data[0]->hdr);
         }

      // Otherwise, just keep going onto the next pileup.

   }

   // No chrY BUG FIX.
   // OK, we've gone through the pileup, but what happens if we've still targets left
   // in the cfg file? Mark anything remaining in the config file that
   // (obviously at this stage) can not make it into the BAM as NO_COVERAGE.

   while(nextTargetRegion.contigPtr[0] != '\0') {

      fprintf(vcfOutFile,"%s\t%ld\t%s\t%c\t%c\t%c\t%s\t%s\t%s\n",
              nextTargetRegion.contigPtr,  // CHROM
              nextTargetRegion.locus + 1, // POS (change to 1 based)
              ".", // ID
              '.', // REF
              '.', // ALT
              '.', // QUAL
              filterFieldNames[NO_COVERAGE], // FILTER
              ".", // INFO
              "." // FORMAT
             );

      getNextTarget(&nextTargetRegion, data[0]->cfgFile, data[0]->hdr);

   }


   // Clean up.
   // Remember, if we have more than one BAM then we have more than one BAM to clean up after..
   // (I know in this example we've only 1 but for completeness I've left as a loop).
   // Loop through all the BAMs we opened, freeing in turn, all the memory that
   // we (or rather htslib, on our behalf) acquired to manage them when we were running through them.
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
   fclose(vcfOutFile);


   // Print out some overall stats.
   printf("\nDONE...\nalignmentCount (#reads) = %ld,\nnumberOfLociCovered = %ld\ntotalFoldCoverage = %ld\nmaxDepth = %ld, Avg. coverage = %ld\n",alignmentCount,numberOfLociCovered,totalFoldCoverage,maxDepth,(long)(totalFoldCoverage/numberOfLociCovered));


}


