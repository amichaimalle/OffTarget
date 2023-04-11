#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//--------------------------------------------
//        Definition of Constants
//--------------------------------------------
// run parameters
#define FILE_PATH "/Users/amichaim/CLionProjects/OffTarget/OffTarget_project"
//#define FILE_PATH "inout"
#define ALLOCATE_MEMORY_FOR_TEXT_FILE 300000000 //300000000
#define MAX_OFF_TARGETS 1000000
#define MAX_LINE_SIZE 100
enum RunMode {forward, reverse, both};

// pattern parameters
#define ALPHABET_SIZE 4
#define CHAR_TO_MASK(char) (char == 'A' ? 0 : char == 'C' ? 1 : char == 'G' ? 2 : char == 'T' ? 3 : 4)
#define MAX_PATTERN_LENGTH 23

// distance parameters
#define MAX_MISMATCH 4
#define MAX_BALCH 1
#define MAX_DISTANCE 5

//structs
typedef struct {
    int inx;
    int distance;
    int balch;
    int mismatch;
    char Reverse;
    char TextTarget[MAX_PATTERN_LENGTH+2];    // leave place for \0
    char alignmentCode[MAX_PATTERN_LENGTH+2]; // leave place for \0
} OffTarget;

typedef struct {
    //unsigned long RdVectors;
    unsigned long RdMatchVectors[MAX_DISTANCE+1];
    unsigned long RdMismatchVectors[MAX_DISTANCE+1];
    unsigned long RdInsertionVectors[MAX_DISTANCE+1];
    unsigned long RdDeletionVectors[MAX_DISTANCE+1];
} RdMatrix;

//--------------------------------------------
//        Prototypes
//--------------------------------------------
//
int FileLength(FILE *file);
int ReadTextFile(char *Text ,int *TextInx);
int  SetPatternBitMaskVectors(unsigned long *PatternBitMaskVectors);
void InitRdVectorsAndMatrix(unsigned long *RdVectors, RdMatrix **RdMatrixs, int PatternLength);
int OffFinderRunLoop(enum RunMode runMode,int PatternLength, unsigned long *PatternBitMaskVectors, int TextInx, char *Text, unsigned long *RdVectors, RdMatrix **RdMatrixs, OffTarget *offTargetList);
void BitapCalc(unsigned long PmVector, unsigned long *RdVectors, RdMatrix *RdMatrix);
OffTarget *CheckForMatch(unsigned long *RdVectors, int inx, RdMatrix **RdMatrixs, int PatternLength, OffTarget *offTarget);
OffTarget *DPTargetTB(int InitInx, int PatternLength, RdMatrix **RdMatrixs, int Errors, unsigned long CheckInxBitZero, char *AlignmentCode, int PatternInx, int TextInx, int AlignmentInx, int Mismatch, int Balch, int EnDeletion, int EnInsertion, OffTarget *offTarget);
void insertOffTarget(OffTarget *offTargetList, OffTarget *offTarget, int offTargetListIndex);
void sortOffTargetLintByInx(OffTarget *offTargetList, int offTargetListSize);
int addOffTargetToList(OffTarget *offTargetList, int offTargetListSize, OffTarget *offTarget);
void printOffTargets(OffTarget *offTargetList, int NumOfOffTargets);
void FreeRdMatrixs(RdMatrix **RdMatrixs, int PatternLength);


int main(){
    enum RunMode runMode = both;
    clock_t start = clock();
    // Text variables
    char *Text = (char *) malloc(ALLOCATE_MEMORY_FOR_TEXT_FILE * sizeof(char));
    int TextInx = 0;
    // Pattern variables & initialization
    unsigned long PatternBitMaskVectors[ALPHABET_SIZE];
    unsigned long RdVectors[MAX_DISTANCE+1];
    int PatternLength = SetPatternBitMaskVectors(PatternBitMaskVectors);
    //  initialize the Rd & OffTarget Parameters
    RdMatrix **RdMatrixs = (RdMatrix **)malloc((PatternLength+1)*sizeof(RdMatrix *));
    InitRdVectorsAndMatrix(RdVectors, RdMatrixs, PatternLength);
    int NumOfOffTargets = 0;
    OffTarget *offTargetList = (OffTarget *)malloc(MAX_OFF_TARGETS*sizeof(OffTarget));
    ReadTextFile(Text ,&TextInx);
    clock_t mid = clock();
    double elapsed0 = (double) (mid - start) / CLOCKS_PER_SEC;
    printf("init phase done, time for init stage: %f seconds\n", elapsed0);

    NumOfOffTargets = OffFinderRunLoop(runMode, PatternLength, PatternBitMaskVectors, TextInx, Text, RdVectors, RdMatrixs, offTargetList);
    clock_t end = clock();
    double elapsed1 = (double)(end - mid) / CLOCKS_PER_SEC;
    printf("run phase done, time for run stage: %f seconds\n", elapsed1);

    printf("number of target found - %d\n", NumOfOffTargets);

    sortOffTargetLintByInx(offTargetList,NumOfOffTargets);
    printOffTargets(offTargetList, NumOfOffTargets);
    FreeRdMatrixs(RdMatrixs, PatternLength);
    free(Text);
    free(offTargetList);
    return 0;
}

int ReadTextFile(char *Text ,int *TextInx) {
    FILE *TextFile = fopen(FILE_PATH"/chr1.txt", "r");
    char lineBuffer[MAX_LINE_SIZE];
    while (fscanf(TextFile, "%[^\n]c", lineBuffer) != EOF) { //read until new line
        fscanf(TextFile, "%*c"); //skip the new line '\n'
        strcpy(Text + *TextInx, lineBuffer);
        *TextInx = *TextInx + (int) strlen(lineBuffer);
    }
    fclose(TextFile);
    return 0;
}

int FileLength(FILE *file){ // no bug but slower function - can be faster if there is no \n in file.
    int length = 0;
    char charBuffer;
    fseek(file, 0L, SEEK_SET);
    while(fscanf(file,"%c",&charBuffer)!=EOF){
        if (charBuffer == 'A' || charBuffer == 'C' || charBuffer == 'G' || charBuffer == 'T' || charBuffer == 'N') length++;
    }
    fseek(file, 0L, SEEK_SET);
    return length;
}

int SetPatternBitMaskVectors(unsigned long *PatternBitMaskVectors){
    FILE *PatternFile = fopen(FILE_PATH"/guides.txt", "r");
    int i, PatternLength = FileLength(PatternFile);
    char letter;
    unsigned long MaskVector = (unsigned long)pow(2,sizeof(unsigned long)*8-1);
    unsigned long initPMValue = (unsigned long)(pow(2,(int)sizeof(unsigned long)*8-PatternLength))-1;
    for (i=0; i<ALPHABET_SIZE; i++){
        PatternBitMaskVectors[i] = initPMValue;
    }
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

void InitRdVectorsAndMatrix(unsigned long *RdVectors, RdMatrix **RdMatrixs, int PatternLength){
    int i;
    for (i=0;i<=MAX_DISTANCE;i++){
        RdVectors[i] = -1;
        RdVectors[i] = RdVectors[i] << (sizeof(unsigned long)*8 - PatternLength);
    }
    for (i=0;i<PatternLength+1;i++){
        RdMatrixs[i] = (RdMatrix*)malloc(sizeof(RdMatrix));
        RdMatrixs[i]->RdMatchVectors[0] = -1;
        RdMatrixs[i]->RdMismatchVectors[0] = -1;
        RdMatrixs[i]->RdInsertionVectors[0] = -1;
        RdMatrixs[i]->RdDeletionVectors[0] = -1;
    }
}

int OffFinderRunLoop(enum RunMode runMode,int PatternLength, unsigned long *PatternBitMaskVectors, int TextInx, char *Text, unsigned long *RdVectors, RdMatrix **RdMatrixs, OffTarget *offTargetList){ // TextInx need to be at the end of the text file
    int NumOfOffTargets = 0;
    OffTarget *tempOffTarget = (OffTarget *)malloc(sizeof(OffTarget));
    if ((runMode == forward) || (runMode == both)) {
        while (TextInx >= 0) { // Reverse = "+" (forward)
            if (Text[TextInx] == 'A' || Text[TextInx] == 'C' || Text[TextInx] == 'G' || Text[TextInx] == 'T') {
                BitapCalc(PatternBitMaskVectors[CHAR_TO_MASK(Text[TextInx])], RdVectors, RdMatrixs[TextInx%(PatternLength+1)]);
                if (CheckForMatch(RdVectors, TextInx, RdMatrixs, PatternLength, tempOffTarget) != NULL) {
                    for (int j = 0; j < (PatternLength+1); j++) {  // add the relevant text to the off target
                        tempOffTarget->TextTarget[j] = Text[TextInx + j];
                    }
                    NumOfOffTargets = addOffTargetToList(offTargetList, NumOfOffTargets, tempOffTarget);
                }
            }
            TextInx--;
        }
    }
    /*if (runMode == reverse) | (runMode == both) {
        while (inx > TextLength) { // Reverse = "-" (reverse)
            if (Text[inx] == 'A' || Text[inx] == 'C' || Text[inx] == 'G' || Text[inx] == 'T') {
                BitapCalc(PatternBitMaskVectors[CHAR_TO_MASK(Text[inx])], RdVectors, RdMatrixs[PatternLength-1 - (TextInx%PatternLength)]); check the ordet of matrix saved
                tempOffTarget = CheckForMatch(RdVectors, inx, RdMatrixs, PatternLength);
                if (tempOffTarget != NULL) {
                    for (int j = 0; j < 24; j++) {  // add the relevant text to the off target
                        tempOffTarget->TextTarget[j] = Text[inx + j];
                    }
                    NumOfOffTargets = addOffTargetToList(offTargetList, NumOfOffTargets, tempOffTarget);
                }
            }
            inx++;
        }
    }*/
    free(tempOffTarget);
    return NumOfOffTargets;
}

void BitapCalc(unsigned long PmVector, unsigned long *RdVectors, RdMatrix *RdMatrix){
    int d;
    unsigned long RdMinus1Vector;
    unsigned long RdVectorInsertion, RdVectorDeletion, RdVectorMatch, RdVectorMismatch;
    RdMinus1Vector = RdVectors[0];
    RdVectors[0] = (RdVectors[0]<<1 | PmVector);
    RdMatrix->RdMatchVectors[0] = RdVectors[0];
    for (d=1;d<=MAX_DISTANCE;d++){
        RdVectorMatch = (RdVectors[d]<<1 | PmVector);
        RdVectorMismatch = RdMinus1Vector<<1;
        RdVectorDeletion = RdMinus1Vector;
        RdVectorInsertion = RdVectors[d-1]<<1;
        RdMatrix->RdMatchVectors[d] = RdVectorMatch;
        RdMatrix->RdMismatchVectors[d] = RdVectorMismatch;
        RdMatrix->RdDeletionVectors[d] = RdVectorDeletion;
        RdMatrix->RdInsertionVectors[d] = RdVectorInsertion;
        RdMinus1Vector = RdVectors[d];  // save RdMinus1Vector for next iteration
        RdVectors[d] = RdVectorMatch & RdVectorMismatch & RdVectorInsertion & RdVectorDeletion;
    }
}

OffTarget *CheckForMatch(unsigned long *RdVectors, int inx, RdMatrix **RdMatrixs, int PatternLength, OffTarget *offTarget){
    int d;
    unsigned long CheckMsbIsZero = (unsigned long)pow(2, sizeof(unsigned long)*8-1);
    char AlignmentCode[MAX_PATTERN_LENGTH+2];
    for (d = 0; d <=MAX_DISTANCE; d++) {
        if ((RdVectors[d] & CheckMsbIsZero) == 0) {
            if(DPTargetTB(inx,PatternLength, RdMatrixs, d, CheckMsbIsZero, AlignmentCode, 0,
                          PatternLength, 0, 0, 0, 1, 1, offTarget) != NULL) {
                offTarget->inx = inx;
                return offTarget;
            }
        }
    }
    return NULL;
}

OffTarget *DPTargetTB(int InitInx, int PatternLength, RdMatrix **RdMatrixs, int Errors, unsigned long CheckInxBitZero, char *AlignmentCode, int PatternInx, int TextInx, int AlignmentInx, int Mismatch, int Balch, int EnDeletion, int EnInsertion, OffTarget *offTarget){
    int CurError = Errors;
    int MatrixInx = InitInx%(PatternLength+1);    // initialize MatrixInx to the 1st matrix for this alignment
    while ((PatternInx < PatternLength) & (TextInx > 0)){
        if ((RdMatrixs[MatrixInx]->RdMatchVectors[CurError] & CheckInxBitZero) == 0) {  //Match
            PatternInx++;
            TextInx--;
            CheckInxBitZero = CheckInxBitZero >> 1;
            MatrixInx = (MatrixInx+1)%(PatternLength+1);
            AlignmentCode[AlignmentInx] = 'M';
            AlignmentInx++;
        } else if (((RdMatrixs[MatrixInx]->RdDeletionVectors[CurError] & CheckInxBitZero) == 0) && (Balch<MAX_BALCH) && (EnDeletion == 1)){ //Deletion
            if(DPTargetTB(MatrixInx, PatternLength, RdMatrixs, CurError, CheckInxBitZero, AlignmentCode, PatternInx,
                          TextInx, AlignmentInx, Mismatch, Balch, 0, 0, offTarget) != NULL) return offTarget;
            if(DPTargetTB(MatrixInx, PatternLength, RdMatrixs, CurError, CheckInxBitZero, AlignmentCode, PatternInx,
                          TextInx, AlignmentInx, Mismatch, Balch, 0, 1, offTarget) != NULL) return offTarget;
            CurError--;
            Balch++;
            PatternInx++;
            MatrixInx = (MatrixInx+1)%(PatternLength+1);
            AlignmentCode[AlignmentInx] = 'D';
            AlignmentInx++;
        } else if (((RdMatrixs[MatrixInx]->RdInsertionVectors[CurError] & CheckInxBitZero) == 0) && (Balch<MAX_BALCH) && (EnInsertion == 1)){ //Insertion
            if (DPTargetTB(MatrixInx, PatternLength, RdMatrixs, CurError, CheckInxBitZero, AlignmentCode, PatternInx,
                           TextInx, AlignmentInx, Mismatch, Balch, 0, 0, offTarget) != NULL) return offTarget;
            CurError--;
            Balch++;
            TextInx--;
            CheckInxBitZero = CheckInxBitZero >> 1;
            AlignmentCode[AlignmentInx] = 'I';
            AlignmentInx++;
            EnDeletion = 1;
        } else if (((RdMatrixs[MatrixInx]->RdMismatchVectors[CurError] & CheckInxBitZero) == 0) && (Mismatch<MAX_MISMATCH)) { //Mismatch/Substitution
            CurError--;
            Mismatch++;
            PatternInx++;
            TextInx--;
            CheckInxBitZero = CheckInxBitZero >> 1;
            MatrixInx = (MatrixInx+1)%(PatternLength+1);
            AlignmentCode[AlignmentInx] = 'S';
            AlignmentInx++;
            EnDeletion = 1;
            EnInsertion = 1;
        } else {
            return NULL;
        }
    }
    while (AlignmentInx < MAX_PATTERN_LENGTH+1) {
        AlignmentCode[AlignmentInx] = 'M';
        AlignmentInx++;
    }
    AlignmentCode[AlignmentInx] = '\0';
    offTarget->distance = Balch + Mismatch;
    offTarget->balch = Balch;
    offTarget->mismatch = Mismatch;
    strcpy(offTarget->alignmentCode, AlignmentCode);
    return offTarget;
}

int addOffTargetToList(OffTarget *offTargetList, int offTargetListSize, OffTarget *offTarget){
    int i;
    for (i=0;i<offTargetListSize; i++){
        if (offTargetList[i].inx == offTarget->inx){//) && (offTargetList[i].Reverse == offTarget->Reverse)){
            if (offTargetList[i].distance > offTarget->distance){
                insertOffTarget(offTargetList, offTarget, i);
            }
            return offTargetListSize;
        }
    }
    if (i==offTargetListSize){
        insertOffTarget(offTargetList, offTarget, offTargetListSize);
        offTargetListSize++;
    }
    return offTargetListSize;
}

void insertOffTarget(OffTarget *offTargetList, OffTarget *offTarget, int offTargetListIndex){
    strcpy(offTargetList[offTargetListIndex].TextTarget, offTarget->TextTarget);
    strcpy(offTargetList[offTargetListIndex].alignmentCode, offTarget->alignmentCode);
    offTargetList[offTargetListIndex].TextTarget[MAX_PATTERN_LENGTH+1] = '\0';
    offTargetList[offTargetListIndex].inx = offTarget->inx;
    offTargetList[offTargetListIndex].distance = offTarget->distance;
    offTargetList[offTargetListIndex].balch = offTarget->balch;
    offTargetList[offTargetListIndex].mismatch = offTarget->mismatch;
    offTargetList[offTargetListIndex].Reverse = '+';
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
    FILE *OutputFile = fopen(FILE_PATH"/output.txt", "w");
    fprintf(OutputFile,"OffFinder v3.0  (c)\n");
    //printf("OffFinder v3.0\n");
    //fprintf(OutputFile,"Index\tDistance\tReverse\n");
    //printf("Index\tDistance\tMismatch\tBalch\tReverse\n");
    for (int i = 0; i < NumOfOffTargets; ++i) {
        fprintf(OutputFile, "Inx : %10d\tTarget: %24s\talignment: %24s\tdistance : %d\tMismatch : %d\tBalch : %d\tReverse : %c\n", offTargetList[i].inx,
                offTargetList[i].TextTarget, offTargetList[i].alignmentCode, offTargetList[i].distance, offTargetList[i].mismatch,  offTargetList[i].balch, offTargetList[i].Reverse);
        //printf("Inx : %10d\t Target: %24s\t alignment: %24s\t distance : %d\tMismatch : %d\tBalch : %d\tReverse : %c\n", offTargetList[i].inx, offTargetList[i].TextTarget,
        //       offTargetList[i].alignmentCode, offTargetList[i].distance, offTargetList[i].mismatch,  offTargetList[i].balch, offTargetList[i].Reverse);
    }
    fclose(OutputFile);
}

void FreeRdMatrixs(RdMatrix **RdMatrixs, int PatternLength){
    for (int i = 0; i < PatternLength; ++i) {
        free(RdMatrixs[i]);
    }
    free(RdMatrixs);
}