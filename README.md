# Rank-Select-Queries-with-Adjusted-Anchoring

https://github.com/simongog/sdsl-lite sdsl-lite framwork is used in the implementation. 
Thus, you need to have sdsl-lite installed to run the RSAA.

You can compile with 

g++ main.cpp -O3 -lsdsl -o RSAA -I /sdsl-lite-include-file-directory/ -L /the-directory-including-the-sdsl-lib/ -std=c++17 -march=native 

Your CPU should be supporting the popcnt,pdep, tzcnt instrcutions.

The program accepts two types of input 

1) Benchmark the R/S schemes on randomly generated bitmaps via the command 

RSAA 0 percentage-of-the-density-of-the-to-be-generated-random-bitmap size-of-the-random-bitmap-in-megabits s_parameter 
e.g., ./RSAA 0 40 10 128 generates a random bitmap of length 10x1024x1024 bits long where 40% of the bits are randomly set to 1, and then runs the benchmark on this bitmap with s=128

2) Benchmark the schemes on the number files GOV2, URL, 5GRAM, DNA from the links provided in https://github.com/aboffa/Learned-Rank-Select-ALENEX21 via the command

RSAA 1 inputNumberFile s_parameter
e.g., ./RSAA 1 /RelatedDirectoryPath/GOV2/10M-/11_gov.bin.docs 16 runs the benchmark on the /GOV2/10M-/11_gov.bin.docs file with s=16

You can feed the files GOV2(files under each directory in the mentioned link) URL(_1/_2/_3) 5GRAM(_1/_2/_3) DNA(_1/_2/_3)   directly to the RSAA program.

