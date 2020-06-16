all: trimPhred.c parseAlignment.cpp countMutations.cpp parsed_reads.h ref_seq.h string_funcs.h
	gcc -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 trimPhred.c -o trimPhred -O3
	g++ parseAlignment.cpp -o parseAlignment -O3
	g++ countMutations.cpp -o countMutations -O3

