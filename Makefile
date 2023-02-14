all: stochasticSpike createDonorGenome liftover targetRef condenseLift gtMapper 2wayLiftover

stochasticSpike: Makefile
	cc -std=c99 -Wall -g -O2 -I${SAMTOOLS_BUILD_PATH} -I${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION} -I${SAMTOOLS_BUILD_PATH}/lz4  -o ./bin/stochasticSpike stochasticSpike.c ${SAMTOOLS_BUILD_PATH}/libst.a ${SAMTOOLS_BUILD_PATH}/htslib-${HTSLIB_VERSION}/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl -lcrypto

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

clean:
	rm ./bin/stochasticSpike ./bin/createDonorGenome ./bin/liftover ./bin/targetRef ./bin/condenseLift ./bin/gtMapper ./bin/2wayLiftover
