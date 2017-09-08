#To compile:
#module load valgrind/3.10
#module load zlib/1.2.8
gcc ./sources/*.c -lm -o ./bin/weight4tfa -Wall -O3 -g -O0 -lz -g

cd ./Examples

valgrind --leak-check=full --track-origins=yes -v ../bin/weight4tfa -h > ../weight4tfa_help.txt 2> valgrind00.txt
valgrind --leak-check=full --track-origins=yes -v ../bin/weight4tfa -i 3Kallchr.tfa.gz -o 3Kallchr.tfa_weights.gz -g ./3Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 3000,2000,3000 -G 2  2> valgrind01.txt
valgrind --leak-check=full --track-origins=yes -v ../bin/weight4tfa -i 100Kallchr.tfa.gz -o 100Kallchr.tfa_weights.gz -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -n chr10,chr12,chr14 -l 50000,100000,75000 -G 2  2> valgrind02.txt
