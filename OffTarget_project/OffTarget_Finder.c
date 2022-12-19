#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ALLOCATE_MEMORY_FOR_TEXT_FILE 300000000 //300000000
#define MAX_OFF_TARGETS 7000
#define MAX_LINE_SIZE 100
#define FILE_PATH "/Users/amichaim/CLionProjects/OffTarget/OffTarget_project/inout"
#define MAX_DISTANCE 5 //Maximum distance between two matches - Number of R vectors
#define ALPHABET_SIZE 4
#define CHAR_TO_MASK(char) (char == 'A' ? 0 : char == 'C' ? 1 : char == 'G' ? 2 : char == 'T' ? 3 : 4)

//struct
typedef struct {
    int inx;
    int distance;
    char Reverse;
} OffTarget;

//Prototypes
int FileLength(FILE *file);
int  SetPatternBitMaskVectors(unsigned long *PatternBitMaskVectors);
void InitRdVectors(unsigned long *RdVectors, int PatternLength);
void BitapCalc(unsigned long PmVector, unsigned long *RdVectors);
int CheckForMatch(unsigned long *RdVectors, int inx, OffTarget *offTargetList, int NumOfOffTargets);
int addOffTargetToList(OffTarget *offTargetList, int index, int distance, int offTargetListSize);
void sortOffTargetLintByInx(OffTarget *offTargetList, int offTargetListSize);
void printOffTargets(FILE *OutputFile, OffTarget *offTargetList, int NumOfOffTargets);

int main(){
    clock_t start = clock();
    FILE *TextFile = fopen(FILE_PATH"/text.txt", "r");
    FILE *OutputFile = fopen(FILE_PATH"/output.txt", "w");
    int inx = 0, PatternLength, TextLength;
    int NumOfOffTargets = 0;
    OffTarget *offTargetList = (OffTarget *)malloc(MAX_OFF_TARGETS*sizeof(OffTarget));
    unsigned long PatternBitMaskVectors[ALPHABET_SIZE];
    unsigned long RdVectors[MAX_DISTANCE+1];
    //  initialize the PatternMasks & RdVectors
    PatternLength = SetPatternBitMaskVectors(PatternBitMaskVectors);
    InitRdVectors(RdVectors, PatternLength);
    if (ALLOCATE_MEMORY_FOR_TEXT_FILE){
        char *Text = (char *)malloc(ALLOCATE_MEMORY_FOR_TEXT_FILE*sizeof(char));
        char lineBuffer[MAX_LINE_SIZE];
        while (fscanf(TextFile,"%[^\n]c",lineBuffer) != EOF){ //read until new line
            fscanf(TextFile,"%*c"); //skip the new line '\n'
            strcpy(Text+inx,lineBuffer);
            inx = (int)strlen(lineBuffer)+inx;
        }
        clock_t mid = clock();
        double elapsed0 = (double)(mid - start) / CLOCKS_PER_SEC;
        printf("Finish reading, time: %f seconds\n", elapsed0);
        fclose(TextFile);
        while (inx != 0){
            BitapCalc(PatternBitMaskVectors[CHAR_TO_MASK(Text[inx])], RdVectors);
            NumOfOffTargets = CheckForMatch(RdVectors, inx, offTargetList, NumOfOffTargets);
            inx--;
        }
        free(Text);
    } else {
        char charBuffer;
        fseek(TextFile, -1L, SEEK_END);
        TextLength = (int) ftell(TextFile) + 1;
        //TextLength = FileLength(TextFile);
        //fseek(TextFile, -1L, SEEK_END);
        do {
            fscanf(TextFile, "%c", &charBuffer);
            if (charBuffer == 'A' || charBuffer == 'C' || charBuffer == 'G' || charBuffer == 'T') {
                BitapCalc(PatternBitMaskVectors[CHAR_TO_MASK(charBuffer)], RdVectors);
                NumOfOffTargets = CheckForMatch(RdVectors, TextLength - inx - 1, offTargetList, NumOfOffTargets);
                inx++;
            }
        } while (fseek(TextFile, -2L, SEEK_CUR) == 0);
        fclose(TextFile);
    }
    clock_t end = clock();
    double elapsed1 = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Finish running time: %f seconds\n", elapsed1);
    sortOffTargetLintByInx(offTargetList,NumOfOffTargets);
    printOffTargets(OutputFile, offTargetList, NumOfOffTargets);
    free(offTargetList);
    fclose(OutputFile);
}
int FileLength(FILE *file){ // no bug but slower function
    int length = 0;
    char charBuffer;
    fseek(file, 0L, SEEK_SET);
    while(fscanf(file,"%c",&charBuffer)!=EOF){
        if (charBuffer == 'A' || charBuffer == 'C' || charBuffer == 'G' || charBuffer == 'T' || charBuffer == 'N') length++;
    }
    fseek(file, 0L, SEEK_SET);
    return length;
}
/*int FileLength(FILE *file){ //There is a bug in this function!! +1 error for each \n..
    fseek(file, 0L, SEEK_END);
    int length = (int)ftell(file);
    fseek(file, 0L, SEEK_SET);
    return length;
}*/
int SetPatternBitMaskVectors(unsigned long *PatternBitMaskVectors){
    FILE *PatternFile = fopen(FILE_PATH"/pattern.txt", "r");
    unsigned long MaskVector = (unsigned long)pow(2,sizeof(unsigned long)*4-1);
    int PatternLength = FileLength(PatternFile);
    char letter;
    int i, initPMValue = (int)pow(2,(int)sizeof(unsigned long)*4-PatternLength)-1;
    for (i=0; i<ALPHABET_SIZE; i++){
        PatternBitMaskVectors[i] = initPMValue;
    }
    //PatternBitMaskVectors = {PatternLength, PatternLength, PatternLength, PatternLength};
    while (fscanf(PatternFile,"%c",&letter) != EOF){
        if (letter == 'A'){
            PatternBitMaskVectors[CHAR_TO_MASK('A')] = PatternBitMaskVectors[CHAR_TO_MASK('A')] | MaskVector;
        } else if (letter == 'C'){
            PatternBitMaskVectors[CHAR_TO_MASK('C')] = PatternBitMaskVectors[CHAR_TO_MASK('C')] | MaskVector;
        } else if (letter == 'G'){
            PatternBitMaskVectors[CHAR_TO_MASK('G')] = PatternBitMaskVectors[CHAR_TO_MASK('G')] | MaskVector;
        } else if (letter == 'T'){
            PatternBitMaskVectors[CHAR_TO_MASK('T')] = PatternBitMaskVectors[CHAR_TO_MASK('T')] | MaskVector;
        };
        MaskVector = MaskVector >> 1;
    }
    for (i=0;i<ALPHABET_SIZE;i++){
        PatternBitMaskVectors[i] = ~PatternBitMaskVectors[i];
    }
    fclose(PatternFile);
    return PatternLength;
}
void InitRdVectors(unsigned long *RdVectors, int PatternLength){
    int i;
    for (i=0;i<=MAX_DISTANCE;i++){
        RdVectors[i] = -1;
        RdVectors[i] = RdVectors[i] << (sizeof(unsigned long)*4 - PatternLength);
    }
}
void BitapCalc(unsigned long PmVector, unsigned long *RdVectors){
    int d;
    unsigned long RdMinus1Vector;
    unsigned long RdVectorInsertion, RdVectorDeletion, RdVectorMatch, RdVectorMismatch;
    RdMinus1Vector = RdVectors[0];
    RdVectors[0] = (RdVectors[0]<<1 | PmVector);
    for (d=1;d<=MAX_DISTANCE;d++){
        RdVectorMatch = (RdVectors[d]<<1 | PmVector);
        RdVectorMismatch = RdMinus1Vector<<1;
        RdVectorDeletion = RdMinus1Vector;
        RdVectorInsertion = RdVectors[d-1]<<1;
        RdMinus1Vector = RdVectors[d];  // save RdMinus1Vector for next iteration
        RdVectors[d] = RdVectorMatch & RdVectorMismatch & RdVectorInsertion & RdVectorDeletion;
    }
}
int CheckForMatch(unsigned long *RdVectors, int inx, OffTarget *offTargetList, int NumOfOffTargets){
    int d;
    unsigned long CheckMsbIsZero = (unsigned long)pow(2, sizeof(unsigned long)*4-1);
    for (d = 0; d <=MAX_DISTANCE; d++) {
        if ((RdVectors[d] & CheckMsbIsZero) == 0) {
            NumOfOffTargets = addOffTargetToList(offTargetList, inx, d, NumOfOffTargets);
            return NumOfOffTargets;
        }
    }
    return NumOfOffTargets;
}
int addOffTargetToList(OffTarget *offTargetList, int index, int distance, int offTargetListSize){
    int i;
    for (i=0;i<offTargetListSize; i++){
        if (offTargetList[i].inx == index){//} && offTargetList[i].Reverse == Reverse){
            if (offTargetList[i].distance > distance){
                offTargetList[i].distance = distance;
            }
            return offTargetListSize;
        }
    }
    if (i==offTargetListSize){
        offTargetList[i].inx = index;
        offTargetList[i].distance = distance;
        offTargetList[i].Reverse = '+'; //= Reverse;
        offTargetListSize++;
    }
    return offTargetListSize;
}
void sortOffTargetLintByInx(OffTarget *offTargetList, int offTargetListSize){
    int i,j;
    OffTarget temp;
    for (i=0;i<offTargetListSize; i++){
        for (j=i+1;j<offTargetListSize; j++){
            if (offTargetList[i].inx > offTargetList[j].inx){
                temp = offTargetList[i];
                offTargetList[i] = offTargetList[j];
                offTargetList[j] = temp;
            }
        }
    }
}
void printOffTargets(FILE *OutputFile, OffTarget *offTargetList, int NumOfOffTargets){
    fprintf(OutputFile,"OffFinder v1.0\n");
    //printf("OffFinder v1.0\n");
    fprintf(OutputFile,"Index\tDistance\tReverse\n");
    //printf("Index\tDistance\tReverse\n");
    for (int i = 0; i < NumOfOffTargets; ++i) {
        fprintf(OutputFile, "Inx : %10d\tdistance : %d\tReverse : %c\n", offTargetList[i].inx,
                    offTargetList[i].distance, offTargetList[i].Reverse);
        //printf("Inx : %10d\tdistance : %d\tReverse : %c\n", offTargetList[i].inx, offTargetList[i].distance,
        //           offTargetList[i].Reverse);
    }
}