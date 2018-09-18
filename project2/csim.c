//516030910259 Xinpeng Liu
#include "cachelab.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <memory.h>
typedef struct {
	int valid;
	long tag;
	int t;
} cacheLine;

cacheLine *cache;// Cache
int s,e,b;// S, E, B
int verboseFlag=0;// Whether verbose mode is used.
int hit=0,miss=0,eviction=0
int nowt=0;// Time for the access.
char temp[25];
char* filePath=NULL;// Path of the trace file.

void printHelp();
void accessCache(char type, long unsigned int addr, int size);

int main(int argc,char* argv[]){
	int tmp=0,in=0;
	while ((tmp=getopt(argc,argv,"s:E:b:t:hv"))!=-1){
		in=1;
		switch(tmp){
			case 's':
				s=(int) pow(2,atoi(optarg));break;
			case 'E':
				e=atoi(optarg);break;
			case 'b':
				b=(int) pow(2,atoi(optarg));break;
			case 't':
				filePath=optarg;break;
			case 'v':
				verboseFlag=1;break;
			default:
				printHelp();return 1;break;
				//Error type 1: invalid arguments and -h. Output Help doc.
		}
	}//Handling arguments

	if (in==0 && tmp==-1){
		printHelp();
		return 1;
	}
	//Error type 1: invalid arguments. Output Help doc.

	FILE *file=fopen(filePath,"r");
	if (file==NULL){
		printf("File not found.");
		return 2;
	}
	//Error type 2: File not found.
	cache = (cacheLine *) malloc(s*e*sizeof(cacheLine));
	if (cache==NULL){
		printf("Fail to allocate cache.");
		return 3;
	}
	//Error type 3: Fail to allocate cache.
	cacheLine* ptr;
	for (int i=0;i<s*e;i++) {
		ptr=(cache+i);
		ptr->valid=0;
	}
	//Initialize cache.
	char type;
	int size;
	long unsigned int addr;
	while (!feof(file)){
		int tr=fscanf(file, " %c %lx,%x",&type, &addr, &size);
		if (tr!=3) continue;
		if (type=='I') continue;
		//Invalid insructions.
		accessCache(type,addr,size);
		nowt++;
	}
	free(cache);
	cache=NULL;
	printSummary(hit, miss, eviction);
	return 0;
}

void printHelp(){
	printf("Usage: ./csim-wrc [-hv] -s <s> -E <E> -b <b> -t <tracefile>\n");
	printf("• -h: Optional help flag that prints usage info\n");
    printf("• -v: Optional verbose flag that displays trace info\n");
    printf("• -s <s>: Number of set index bits (S = 2^s is the number of sets)\n");
    printf("• -E <E>: Associativity (number of lines per set)\n");
    printf("• -b <b>: Number of block bits (B = 2^b is the block size)\n");
    printf("• -t <tracefile>: Name of the valgrind trace to replay\n");
}

void accessCache(char type, long unsigned int addr, int size){
	int sIndex=(int) ((addr/b)%s);
	int tag=(int) (addr/(b*s));
	int insert=-1;// Address used by situation 2.
	int replace=-1;// Address used by situation 3.
	int find=-1;// Address used by situation 1.
	int i,mint=nowt;
	cacheLine* ptr;
	for (i=sIndex*e;i<(sIndex+1)*e;i++){
		ptr=(cache+i);
		if (ptr->tag==tag && ptr->valid) {
			find=i;
			break;
			// Situation 1
		}
		else if (!ptr->valid) {insert=i;}// Situation 2
		else if (ptr->t<mint){
			replace=i;mint=ptr->t;
			// Situation 3
		}
	}
	if (verboseFlag) printf(" %c %lx,%x ",type,addr,size);
	// Verbose mode.
	if (find!=-1){
		ptr=cache+find;
		ptr->t=nowt;
		hit++;
		if (verboseFlag) printf("hit");
		// Situation 1
	}
	else if (insert!=-1){
		ptr=cache+insert;
		ptr->valid=1;
		ptr->tag=tag;
		ptr->t=nowt;
		miss++;
		if (verboseFlag) printf("miss");
		// Situation 2
	}
	else {
		ptr=cache+replace;
		ptr->tag=tag;
		ptr->t=nowt;
		miss++;eviction++;
		if (verboseFlag) printf("miss eviction");
		// Situation 3
	}
	if (type=='M') {
		hit++;
		if (verboseFlag) printf(" hit");
		// M instructions.
	}
	if (verboseFlag) printf("\n");
}
