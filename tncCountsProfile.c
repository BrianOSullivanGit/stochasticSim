// Takes a reference bed file as input and gives you the tri-nucleotide context breakdown from
// the reference sequences within the ranges specified in the bed.

// This could be made *a lot* more efficient by caching sectons of the reference but it will do for now...


// Maybe next program that checks the bam should be one that outputs a bed file of reqions covered by the bam at greater than X coverage and Y mapping quality...

#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AAAidx 0
#define AACidx 1
#define AAGidx 2
#define AATidx 3
#define ACAidx 4
#define ACCidx 5
#define ACGidx 6
#define ACTidx 7
#define AGAidx 8
#define AGCidx 9
#define AGGidx 10
#define AGTidx 11
#define ATAidx 12
#define ATCidx 13
#define ATGidx 14
#define ATTidx 15
#define CAAidx 16
#define CACidx 17
#define CAGidx 18
#define CATidx 19
#define CCAidx 20
#define CCCidx 21
#define CCGidx 22
#define CCTidx 23
#define CGAidx 24
#define CGCidx 25
#define CGGidx 26
#define CGTidx 27
#define CTAidx 28
#define CTCidx 29
#define CTGidx 30
#define CTTidx 31
#define GAAidx 32
#define GACidx 33
#define GAGidx 34
#define GATidx 35
#define GCAidx 36
#define GCCidx 37
#define GCGidx 38
#define GCTidx 39
#define GGAidx 40
#define GGCidx 41
#define GGGidx 42
#define GGTidx 43
#define GTAidx 44
#define GTCidx 45
#define GTGidx 46
#define GTTidx 47
#define TAAidx 48
#define TACidx 49
#define TAGidx 50
#define TATidx 51
#define TCAidx 52
#define TCCidx 53
#define TCGidx 54
#define TCTidx 55
#define TGAidx 56
#define TGCidx 57
#define TGGidx 58
#define TGTidx 59
#define TTAidx 60
#define TTCidx 61
#define TTGidx 62
#define TTTidx 63

// faidx_fetch_seq64(fai, c_name, p_beg_i, p_end_i, &len64);
// char *faidx_fetch_seq64(const faidx_t *fai, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len)

// FYI.....
//static inline uint64_t hts_str2uint(const char *in, char **end, int bits,
//htslib-1.11/textutils_internal.h-                                      int *failed) {
/** @param[in]  in     Input string
    @param[out] end    Returned end pointer
    @param[in]  bits   Bits available for the converted value
    @param[out] failed Location of overflow flag
    @return String value converted to a uint64_t

Converts an unsigned decimal string to a uint64_t.  The string should
consist of an optional '+' sign followed by one or more of the digits 0
to 9.  The output value will be limited to fit in the given number of bits.
If the value is too big, the largest possible value will be returned
and *failed will be set to 1.

The address of the first character following the converted number will
be stored in *end.

Both end and failed must be non-NULL.
 */
// Visualise with frange function in  vcfView...

// Increment triplet context count that includes the triples and its reverse complement.
void incCtx(char *refStr, long *ctxCntPtr)
{
   if(!strncmp("AAA", refStr, 3)) {
      ctxCntPtr[AAAidx]++;
      return;
   }
   else if(!strncmp("AAC", refStr, 3)) {
      ctxCntPtr[AACidx]++;
      return;
   }
   else if(!strncmp("AAG", refStr, 3)) {
      ctxCntPtr[AAGidx]++;
      return;
   }
   else if(!strncmp("AAT", refStr, 3)) {
      ctxCntPtr[AATidx]++;
      return;
   }
   else if(!strncmp("ACA", refStr, 3)) {
      ctxCntPtr[ACAidx]++;
      return;
   }
   else if(!strncmp("ACC", refStr, 3)) {
      ctxCntPtr[ACCidx]++;
      return;
   }
   else if(!strncmp("ACG", refStr, 3)) {
      ctxCntPtr[ACGidx]++;
      return;
   }
   else if(!strncmp("ACT", refStr, 3)) {
      ctxCntPtr[ACTidx]++;
      return;
   }
   else if(!strncmp("AGA", refStr, 3)) {
      ctxCntPtr[AGAidx]++;
      return;
   }
   else if(!strncmp("AGC", refStr, 3)) {
      ctxCntPtr[AGCidx]++;
      return;
   }
   else if(!strncmp("AGG", refStr, 3)) {
      ctxCntPtr[AGGidx]++;
      return;
   }
   else if(!strncmp("AGT", refStr, 3)) {
      ctxCntPtr[AGTidx]++;
      return;
   }
   else if(!strncmp("ATA", refStr, 3)) {
      ctxCntPtr[ATAidx]++;
      return;
   }
   else if(!strncmp("ATC", refStr, 3)) {
      ctxCntPtr[ATCidx]++;
      return;
   }
   else if(!strncmp("ATG", refStr, 3)) {
      ctxCntPtr[ATGidx]++;
      return;
   }
   else if(!strncmp("ATT", refStr, 3)) {
      ctxCntPtr[ATTidx]++;
      return;
   }
   else if(!strncmp("CAA", refStr, 3)) {
      ctxCntPtr[CAAidx]++;
      return;
   }
   else if(!strncmp("CAC", refStr, 3)) {
      ctxCntPtr[CACidx]++;
      return;
   }
   else if(!strncmp("CAG", refStr, 3)) {
      ctxCntPtr[CAGidx]++;
      return;
   }
   else if(!strncmp("CAT", refStr, 3)) {
      ctxCntPtr[CATidx]++;
      return;
   }
   else if(!strncmp("CCA", refStr, 3)) {
      ctxCntPtr[CCAidx]++;
      return;
   }
   else if(!strncmp("CCC", refStr, 3)) {
      ctxCntPtr[CCCidx]++;
      return;
   }
   else if(!strncmp("CCG", refStr, 3)) {
      ctxCntPtr[CCGidx]++;
      return;
   }
   else if(!strncmp("CCT", refStr, 3)) {
      ctxCntPtr[CCTidx]++;
      return;
   }
   else if(!strncmp("CGA", refStr, 3)) {
      ctxCntPtr[CGAidx]++;
      return;
   }
   else if(!strncmp("CGC", refStr, 3)) {
      ctxCntPtr[CGCidx]++;
      return;
   }
   else if(!strncmp("CGG", refStr, 3)) {
      ctxCntPtr[CGGidx]++;
      return;
   }
   else if(!strncmp("CGT", refStr, 3)) {
      ctxCntPtr[CGTidx]++;
      return;
   }
   else if(!strncmp("CTA", refStr, 3)) {
      ctxCntPtr[CTAidx]++;
      return;
   }
   else if(!strncmp("CTC", refStr, 3)) {
      ctxCntPtr[CTCidx]++;
      return;
   }
   else if(!strncmp("CTG", refStr, 3)) {
      ctxCntPtr[CTGidx]++;
      return;
   }
   else if(!strncmp("CTT", refStr, 3)) {
      ctxCntPtr[CTTidx]++;
      return;
   }
   else if(!strncmp("GAA", refStr, 3)) {
      ctxCntPtr[GAAidx]++;
      return;
   }
   else if(!strncmp("GAC", refStr, 3)) {
      ctxCntPtr[GACidx]++;
      return;
   }
   else if(!strncmp("GAG", refStr, 3)) {
      ctxCntPtr[GAGidx]++;
      return;
   }
   else if(!strncmp("GAT", refStr, 3)) {
      ctxCntPtr[GATidx]++;
      return;
   }
   else if(!strncmp("GCA", refStr, 3)) {
      ctxCntPtr[GCAidx]++;
      return;
   }
   else if(!strncmp("GCC", refStr, 3)) {
      ctxCntPtr[GCCidx]++;
      return;
   }
   else if(!strncmp("GCG", refStr, 3)) {
      ctxCntPtr[GCGidx]++;
      return;
   }
   else if(!strncmp("GCT", refStr, 3)) {
      ctxCntPtr[GCTidx]++;
      return;
   }
   else if(!strncmp("GGA", refStr, 3)) {
      ctxCntPtr[GGAidx]++;
      return;
   }
   else if(!strncmp("GGC", refStr, 3)) {
      ctxCntPtr[GGCidx]++;
      return;
   }
   else if(!strncmp("GGG", refStr, 3)) {
      ctxCntPtr[GGGidx]++;
      return;
   }
   else if(!strncmp("GGT", refStr, 3)) {
      ctxCntPtr[GGTidx]++;
      return;
   }
   else if(!strncmp("GTA", refStr, 3)) {
      ctxCntPtr[GTAidx]++;
      return;
   }
   else if(!strncmp("GTC", refStr, 3)) {
      ctxCntPtr[GTCidx]++;
      return;
   }
   else if(!strncmp("GTG", refStr, 3)) {
      ctxCntPtr[GTGidx]++;
      return;
   }
   else if(!strncmp("GTT", refStr, 3)) {
      ctxCntPtr[GTTidx]++;
      return;
   }
   else if(!strncmp("TAA", refStr, 3)) {
      ctxCntPtr[TAAidx]++;
      return;
   }
   else if(!strncmp("TAC", refStr, 3)) {
      ctxCntPtr[TACidx]++;
      return;
   }
   else if(!strncmp("TAG", refStr, 3)) {
      ctxCntPtr[TAGidx]++;
      return;
   }
   else if(!strncmp("TAT", refStr, 3)) {
      ctxCntPtr[TATidx]++;
      return;
   }
   else if(!strncmp("TCA", refStr, 3)) {
      ctxCntPtr[TCAidx]++;
      return;
   }
   else if(!strncmp("TCC", refStr, 3)) {
      ctxCntPtr[TCCidx]++;
      return;
   }
   else if(!strncmp("TCG", refStr, 3)) {
      ctxCntPtr[TCGidx]++;
      return;
   }
   else if(!strncmp("TCT", refStr, 3)) {
      ctxCntPtr[TCTidx]++;
      return;
   }
   else if(!strncmp("TGA", refStr, 3)) {
      ctxCntPtr[TGAidx]++;
      return;
   }
   else if(!strncmp("TGC", refStr, 3)) {
      ctxCntPtr[TGCidx]++;
      return;
   }
   else if(!strncmp("TGG", refStr, 3)) {
      ctxCntPtr[TGGidx]++;
      return;
   }
   else if(!strncmp("TGT", refStr, 3)) {
      ctxCntPtr[TGTidx]++;
      return;
   }
   else if(!strncmp("TTA", refStr, 3)) {
      ctxCntPtr[TTAidx]++;
      return;
   }
   else if(!strncmp("TTC", refStr, 3)) {
      ctxCntPtr[TTCidx]++;
      return;
   }
   else if(!strncmp("TTG", refStr, 3)) {
      ctxCntPtr[TTGidx]++;
      return;
   }
   else if(!strncmp("TTT", refStr, 3)) {
      ctxCntPtr[TTTidx]++;
      return;
   }
}


int main(int argc,char* argv[])
{


   long ctxCnt[32*2];

   memset(ctxCnt,0,sizeof(ctxCnt));


   FILE * fp;
   char * line = NULL;
   size_t len = 0;
   ssize_t read;

   fp = fopen(argv[1], "r");
   if(fp == NULL)
      exit(EXIT_FAILURE);

   char tncBuff[4] = {'\0','\0','\0','\0'};
   char *tncBuffEnd=&tncBuff[3];


   // Total number of contexts checked within  reference.
   long totalNumOfContextsChecked=0;

   while((read = getline(&line, &len, fp)) != -1) {

      // Strip newline if its there...
      line[strcspn(line, "\n")] = 0;

      // Skip any blank lines.
      // Also skip any line with a contig name (ie., starting with a '>' character).
      if(line[0] == '\0' || line[0] == '>') continue;

      // If there is just a string of N's don't bother with it..
      int i=0;
      for(; line[i]!='\0'; i++) {
         if(line[i]=='G' || line[i]=='C' || line[i]=='A' || line[i]=='T') {
            break;
         }
      }
      if(line[i]=='\0') continue;

      for(i=0; i<read;) {
         char *tncChPtr=tncBuff;

         // See if there is anything in the buffer remaining from the last attempt.
         while(tncChPtr!= tncBuffEnd && *tncChPtr!='\0') tncChPtr++;

         // We need 3 contiguous bases to check for the occurance of a context..
         while(tncChPtr != tncBuffEnd && i<read) {
            //printf("i=%d, len=%ld\n",i,read);
            *tncChPtr = line[i];
            tncChPtr++;
            i++;
         }

         // We havent got 3 contiguous bases yet, and we're at the end of this line, read next line.
         if(tncChPtr!= tncBuffEnd) {
            break;
         }
         else {

            // Only count bases (ie., no N's) in our estimation.
            if((tncBuff[0]=='G' || tncBuff[0]=='C' || tncBuff[0]=='A' || tncBuff[0]=='T') &&
                  (tncBuff[1]=='G' || tncBuff[1]=='C' || tncBuff[1]=='A' || tncBuff[1]=='T') &&
                  (tncBuff[2]=='G' || tncBuff[2]=='C' || tncBuff[2]=='A' || tncBuff[2]=='T')) {

               //printf("checking %s\n",tncBuff);
               // Check for context and update count.
               incCtx(tncBuff,ctxCnt);
               totalNumOfContextsChecked++;
            }

            // Move the bases along and see what the next character brings
            tncBuff[0] = tncBuff[1];
            tncBuff[1] = tncBuff[2];
            tncBuff[2]='\0';

         }
      }
   }


   // Unstranded counts...
   // (as TNC signature does not differentiate on which DNA strand (forward or reverse) the substitution initially occurred).
   printf("ACA\t%ld\n", (long)(ctxCnt[ACAidx] + ctxCnt[TGTidx]));
   printf("ACC\t%ld\n", (long)(ctxCnt[ACCidx] + ctxCnt[GGTidx]));
   printf("ACG\t%ld\n", (long)(ctxCnt[ACGidx] + ctxCnt[CGTidx]));
   printf("ACT\t%ld\n", (long)(ctxCnt[ACTidx] + ctxCnt[AGTidx]));
   printf("ATA\t%ld\n", (long)(ctxCnt[ATAidx] + ctxCnt[TATidx]));
   printf("ATC\t%ld\n", (long)(ctxCnt[ATCidx] + ctxCnt[GATidx]));
   printf("ATG\t%ld\n", (long)(ctxCnt[ATGidx] + ctxCnt[CATidx]));
   printf("ATT\t%ld\n", (long)(ctxCnt[ATTidx] + ctxCnt[AATidx]));
   printf("CCA\t%ld\n", (long)(ctxCnt[CCAidx] + ctxCnt[TGGidx]));
   printf("CCC\t%ld\n", (long)(ctxCnt[CCCidx] + ctxCnt[GGGidx]));
   printf("CCG\t%ld\n", (long)(ctxCnt[CCGidx] + ctxCnt[CGGidx]));
   printf("CCT\t%ld\n", (long)(ctxCnt[CCTidx] + ctxCnt[AGGidx]));
   printf("CTA\t%ld\n", (long)(ctxCnt[CTAidx] + ctxCnt[TAGidx]));
   printf("CTC\t%ld\n", (long)(ctxCnt[CTCidx] + ctxCnt[GAGidx]));
   printf("CTG\t%ld\n", (long)(ctxCnt[CTGidx] + ctxCnt[CAGidx]));
   printf("CTT\t%ld\n", (long)(ctxCnt[CTTidx] + ctxCnt[AAGidx]));
   printf("GCA\t%ld\n", (long)(ctxCnt[GCAidx] + ctxCnt[TGCidx]));
   printf("GCC\t%ld\n", (long)(ctxCnt[GCCidx] + ctxCnt[GGCidx]));
   printf("GCG\t%ld\n", (long)(ctxCnt[GCGidx] + ctxCnt[CGCidx]));
   printf("GCT\t%ld\n", (long)(ctxCnt[GCTidx] + ctxCnt[AGCidx]));
   printf("GTA\t%ld\n", (long)(ctxCnt[GTAidx] + ctxCnt[TACidx]));
   printf("GTC\t%ld\n", (long)(ctxCnt[GTCidx] + ctxCnt[GACidx]));
   printf("GTG\t%ld\n", (long)(ctxCnt[GTGidx] + ctxCnt[CACidx]));
   printf("GTT\t%ld\n", (long)(ctxCnt[GTTidx] + ctxCnt[AACidx]));
   printf("TCA\t%ld\n", (long)(ctxCnt[TCAidx] + ctxCnt[TGAidx]));
   printf("TCC\t%ld\n", (long)(ctxCnt[TCCidx] + ctxCnt[GGAidx]));
   printf("TCG\t%ld\n", (long)(ctxCnt[TCGidx] + ctxCnt[CGAidx]));
   printf("TCT\t%ld\n", (long)(ctxCnt[TCTidx] + ctxCnt[AGAidx]));
   printf("TTA\t%ld\n", (long)(ctxCnt[TTAidx] + ctxCnt[TAAidx]));
   printf("TTC\t%ld\n", (long)(ctxCnt[TTCidx] + ctxCnt[GAAidx]));
   printf("TTG\t%ld\n", (long)(ctxCnt[TTGidx] + ctxCnt[CAAidx]));
   printf("TTT\t%ld\n", (long)(ctxCnt[TTTidx] + ctxCnt[AAAidx]));

}
