#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// run parameters
#define FILE_PATH "/Users/amichaim/CLionProjects/OffTarget/OffTarget_project/inout"
#define ALLOCATE_MEMORY_FOR_TEXT_FILE 300000000 //300000000
#define MAX_OFF_TARGETS 10000
#define MAX_LINE_SIZE 100
enum RunMode {forward, reverse, both};

// pattern parameters
#define MAX_PATTERN_LENGTH 23


// distance parameters
#define MAX_MISMATCH 1
#define MAX_BALCH 1
#define MAX_DISTANCE 2

//structs
typedef struct {
    int inx;
    int distance;
    int balch;
    int mismatch;
    char TextTarget[MAX_PATTERN_LENGTH+2];    // leave place for \0
    char alignmentCode[MAX_PATTERN_LENGTH+2]; // leave place for \0
} OffTarget;

//Prototypes
int ReadTextFile(char *Text);
int ReadPatternFile(char *Pattern);
OffTarget *NWAligner(char *patternPtr, int patternLength, char *textWindowPtr, int textWindowLength, int maxMismatch, int maxBalch, int *alignmentInx, OffTarget *OffTargetPtr);
int OffFinderRunLoop(char *patternPtr, int PatternLength, char *Text, int TextLength, OffTarget *offTargetList);
void sortOffTargetLintByInx(OffTarget *offTargetList, int offTargetListSize);
void printOffTargets(OffTarget *offTargetList, int NumOfOffTargets);


int main(){
    clock_t start = clock();
    // Text variables
    char *Text = (char *) malloc(ALLOCATE_MEMORY_FOR_TEXT_FILE * sizeof(char));
    char *patternPtr = (char *) malloc(MAX_PATTERN_LENGTH * sizeof(char));
    int PatternLength, TextLength;
    // Pattern variables & initialization
    //  initialize the Rd & OffTarget Parameters
    int NumOfOffTargets = 0;
    OffTarget *offTargetList = (OffTarget *)malloc(MAX_OFF_TARGETS*sizeof(OffTarget));
    TextLength = ReadTextFile(Text);
    PatternLength =  ReadPatternFile(patternPtr);

    clock_t mid = clock();
    double elapsed0 = (double) (mid - start) / CLOCKS_PER_SEC;
    printf("init phase done, time for init stage: %f seconds\n", elapsed0);

    NumOfOffTargets = OffFinderRunLoop(patternPtr, PatternLength, Text, TextLength, offTargetList);

    clock_t end = clock();
    double elapsed1 = (double)(end - mid) / CLOCKS_PER_SEC;
    printf("run phase done, time for run stage: %f seconds\n", elapsed1);
    sortOffTargetLintByInx(offTargetList,NumOfOffTargets);
    printOffTargets(offTargetList, NumOfOffTargets);
    free(Text);
    free(patternPtr);
    free(offTargetList);
    return 0;
}

OffTarget *NWAligner(char *patternPtr, int patternLength, char *textWindowPtr, int textWindowLength, int maxMismatch, int maxBalch, int *alignmentInx, OffTarget *OffTargetPtr){ //, int mismatchIsPrio, int deletionIsPrio){
    while (patternLength>0 && textWindowLength>0){
        if (patternPtr[patternLength-1] == textWindowPtr[textWindowLength-1]){ // match
            //patternLength--;
            //textWindowLength--;
            if (NWAligner(patternPtr, patternLength-1, textWindowPtr, textWindowLength-1, maxMismatch,maxBalch, alignmentInx, OffTargetPtr) != NULL) {  // insertion
                OffTargetPtr->alignmentCode[(*alignmentInx)++] = 'M';
                return OffTargetPtr;
            }
        } else {
            if (maxBalch>0){  // balch
                if (NWAligner(patternPtr, patternLength-1, textWindowPtr, textWindowLength, maxMismatch,maxBalch-1, alignmentInx, OffTargetPtr) != NULL) {  // insertion
                    OffTargetPtr->alignmentCode[(*alignmentInx)++] = 'I';
                    return OffTargetPtr;
                }
                if (NWAligner(patternPtr, patternLength, textWindowPtr, textWindowLength-1, maxMismatch, maxBalch-1, alignmentInx, OffTargetPtr) != NULL) { // deletion
                    OffTargetPtr->alignmentCode[(*alignmentInx)++] = 'D';
                    return OffTargetPtr;
                }
            }
            if (maxMismatch>0){  // mismatch
                if (NWAligner(patternPtr, patternLength-1, textWindowPtr, textWindowLength-1, maxMismatch-1, maxBalch, alignmentInx, OffTargetPtr) != NULL) {  // deletion
                    OffTargetPtr->alignmentCode[(*alignmentInx)++] = 'S';
                    return OffTargetPtr;
                }
            }
        }
        return NULL;        // no match
    }
    OffTargetPtr->mismatch = MAX_MISMATCH - maxMismatch;
    OffTargetPtr->balch = MAX_BALCH - maxBalch;
    OffTargetPtr->distance = OffTargetPtr->balch + OffTargetPtr->mismatch;
    *alignmentInx = 0;
    return OffTargetPtr;
}

int OffFinderRunLoop(char *patternPtr, int PatternLength, char *Text, int TextLength, OffTarget *offTargetList){
    int NumOfOffTargets = 0, TextInx = TextLength - PatternLength;
    int *alignmentInx = (int *)malloc(sizeof(int));
    OffTarget *tempOffTarget = (OffTarget *)malloc(sizeof(OffTarget));
    while (TextInx >= 0){
        if (NWAligner(patternPtr, PatternLength, (Text+TextInx), PatternLength, MAX_MISMATCH, MAX_BALCH, alignmentInx, tempOffTarget) != NULL) {
            tempOffTarget->inx = TextInx;
            tempOffTarget->alignmentCode[(*alignmentInx)++] = '\0';
            offTargetList[NumOfOffTargets] = *tempOffTarget;
            NumOfOffTargets++;
        }
        TextInx--;
    }
    free(tempOffTarget);
    free(alignmentInx);
    return NumOfOffTargets;
}

int ReadTextFile(char *Text) {
    int TextInx = 0;
    FILE *TextFile = fopen(FILE_PATH"/text.txt", "r");
    char lineBuffer[MAX_LINE_SIZE];
    while (fscanf(TextFile, "%[^\n]c", lineBuffer) != EOF) { //read until new line
        fscanf(TextFile, "%*c"); //skip the new line '\n'
        strcpy(Text + TextInx, lineBuffer);
        TextInx = TextInx + (int) strlen(lineBuffer);
    }
    fclose(TextFile);
    return TextInx;
}

int ReadPatternFile(char *Pattern) {
    //int Inx = 0;
    FILE *PatternFile = fopen(FILE_PATH"/pattern.txt", "r");
    //int PatternLength = FileLength(PatternFile);
    fscanf(PatternFile, "%[^\n]c", Pattern);
    int PatternLength = (int) strlen(Pattern);
    //while (fscanf(TextFile, "%[^\n]c", lineBuffer) != EOF) { //read until new line
    //    fscanf(TextFile, "%*c"); //skip the new line '\n'
    //    strcpy(Text + Inx, lineBuffer);
    //    Inx = Inx + (int) strlen(lineBuffer);
    //}
    fclose(PatternFile);
    return PatternLength;
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

void printOffTargets(OffTarget *offTargetList, int NumOfOffTargets){
    FILE *OutputFile = fopen(FILE_PATH"/naive_output.txt", "w");
    fprintf(OutputFile,"naive OffFinder v1.0  (c)\n");
    //printf("OffFinder v3.0\n");
    //fprintf(OutputFile,"Index\tDistance\tReverse\n");
    //printf("Index\tDistance\tMismatch\tBalch\tReverse\n");
    for (int i = 0; i < NumOfOffTargets; ++i) {
        fprintf(OutputFile, "Inx : %10d\tTarget: %24s\talignment: %24s\tdistance : %d\tMismatch : %d\tBalch : %d\n", offTargetList[i].inx,
                offTargetList[i].TextTarget, offTargetList[i].alignmentCode, offTargetList[i].distance, offTargetList[i].mismatch,  offTargetList[i].balch);
        //printf("Inx : %10d\t Target: %24s\t alignment: %24s\t distance : %d\tMismatch : %d\tBalch : %d\tReverse : %c\n", offTargetList[i].inx, offTargetList[i].TextTarget,
        //       offTargetList[i].alignmentCode, offTargetList[i].distance, offTargetList[i].mismatch,  offTargetList[i].balch, offTargetList[i].Reverse);
    }
    fclose(OutputFile);
}