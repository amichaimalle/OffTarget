#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>

#define CasOffinderPath "/Users/amichaim/Desktop/cas-offinder-3.0.0b3"
#define OffinderPath "/Users/amichaim/CLionProjects/OffTarget/OffTarget_project/inout"
#define MAX_TARGETS 10000

#define MAX_INX 10000000000

int main() {
    FILE *CasTargetFile = fopen(CasOffinderPath"/output2.txt", "r");
    //FILE *CasTargetFile = fopen(OffinderPath"/naive_output.txt", "r");
    FILE *OffTargetFile = fopen(OffinderPath"/output.txt", "r");
    int CasTargetInxList[MAX_TARGETS];
    int OffTargetInxList[MAX_TARGETS];
    int NonMatchList[MAX_TARGETS];
    int Inx, CasTargetListSize = 0, OffTargetListSize = 0, MatchCounter = 0, NonMatchCounter = 0;
    int i = 0, j = 0;
    char lineBuffer[300];
    fgets(lineBuffer, 300, OffTargetFile);  // skip the header line
    while (fgets(lineBuffer, 300, OffTargetFile) != NULL) {
        Inx = atoi(lineBuffer + 5);
        if (Inx > MAX_INX) {break;}
        OffTargetInxList[OffTargetListSize] = Inx;
        OffTargetListSize++;
    }
    fgets(lineBuffer, 300, CasTargetFile);  // skip the header line
    while (fgets(lineBuffer, 300, CasTargetFile) != NULL) {
        Inx = atoi(lineBuffer + 5);
        if (Inx > MAX_INX) {break;}
        CasTargetInxList[CasTargetListSize] = Inx;
        CasTargetListSize++;
    }
    for (i = 0; i < OffTargetListSize; ++i) {
        for (j = 0; j < CasTargetListSize; ++j) {
            if (OffTargetInxList[i] == CasTargetInxList[j]) {
                MatchCounter++;
                break;
            }
        }
        if (j == CasTargetListSize){
            NonMatchList[NonMatchCounter] = OffTargetInxList[i];
            NonMatchCounter++;
        }
    }
    printf("---------- Match Report ----------\n");
    printf("Offinder targets Found: %d\n", OffTargetListSize);
    printf("CasOffinder targets Found %d\n", CasTargetListSize);
    printf("Match: %d\tNon-Match: %d\n", MatchCounter,NonMatchCounter);
    printf("---------- Offinder Non-Match Targets ----------\n");
    for (i = 0; i < NonMatchCounter; ++i) {
        printf("%d\n", NonMatchList[i]);
    }
    printf("---------- Cas Non-Match Targets ---------------\n");
    for (i = 0; i < CasTargetListSize; ++i) {
        for (j = 0; j < OffTargetListSize; ++j) {
            if (CasTargetInxList[i] == OffTargetInxList[j]) {
                break;
            }
        }
        if (j == OffTargetListSize){
            printf("%d\n", CasTargetInxList[i]);
        }
    }
    fclose(CasTargetFile);
    fclose(OffTargetFile);
}
