HTSLIB_static_LIBS := $(shell bash -c "grep HTSLIB_static_LIBS ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/htslib_static.mk | sed -e 's/.* = //1'")

all: stochasticSpike stochasticIndel tncSpike vcfAntex createDonorGenome liftover targetRef condenseLift gtMapper 2wayLiftover tncCountsProfile

stochasticSpike: Makefile
	cc -std=c99 -Wall -g -O2 -I${SAMTOOLS_BUILD_PATH} -I${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION} -I${SAMTOOLS_BUILD_PATH}/lz4  -o ./bin/stochasticSpike stochasticSpike.c ${SAMTOOLS_BUILD_PATH}/libst.a ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/libhts.a ${HTSLIB_static_LIBS}

stochasticIndel: Makefile
	cc -std=c99 -Wall -g -O2 -I${SAMTOOLS_BUILD_PATH} -I${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION} -I${SAMTOOLS_BUILD_PATH}/lz4  -o ./bin/stochasticIndel stochasticIndel.c ${SAMTOOLS_BUILD_PATH}/libst.a ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/libhts.a ${HTSLIB_static_LIBS}

tncSpike: Makefile
	cc -std=c99 -Wall -g -O2 -I${SAMTOOLS_BUILD_PATH} -I${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION} -I${SAMTOOLS_BUILD_PATH}/lz4  -o ./bin/tncSpike tncSpike.c ${SAMTOOLS_BUILD_PATH}/libst.a ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/libhts.a ${HTSLIB_static_LIBS}

vcfAntex: Makefile
	cc -std=c99 -Wall -g -O2 -I${SAMTOOLS_BUILD_PATH} -I${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION} -I${SAMTOOLS_BUILD_PATH}/lz4  -o ./bin/vcfAntex vcfAntex.c ${SAMTOOLS_BUILD_PATH}/libst.a ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/libhts.a ${HTSLIB_static_LIBS}

createDonorGenome: Makefile	
	cc -std=c99 createDonorGenome.c -o ./bin/createDonorGenome

liftover: Makefile	
	cc -std=c99 liftover.c -o ./bin/liftover
	
targetRef: Makefile	
	cc -std=c99 targetRef.c -o ./bin/targetRef

condenseLift: Makefile	
	cc -std=c99 condenseLift.c -o ./bin/condenseLift
	
gtMapper: Makefile	
	cc -std=c99 gtMapper.c -o ./bin/gtMapper
	
2wayLiftover: Makefile	
	cc -std=c99 2wayLiftover.c -o ./bin/2wayLiftover

tncCountsProfile: Makefile	
	cc -std=c99 tncCountsProfile.c -o ./bin/tncCountsProfile
clean:
	rm ./bin/stochasticSpike ./bin/stochasticIndel ./bin/createDonorGenome ./bin/liftover ./bin/targetRef ./bin/condenseLift ./bin/gtMapper ./bin/2wayLiftover ./bin/tncSpike ./bin/vcfAntex ./bin/tncCountsProfile

