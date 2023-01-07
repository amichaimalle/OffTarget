#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <time.h>

#define MAX_DISTANCE 6
#define MAX_MISMATCH 5
#define MAX_BALCH 1
#define MAX_PATTERN_LENGTH 23
#define CasOffinderPath "/Users/amichaim/Desktop/cas-offinder-3.0.0b3"

typedef struct {
    int inx;
    int mismatch;
    int balch;
    int distance;
    char Reverse;
    char TextTarget[MAX_PATTERN_LENGTH+1];
} OffTarget;

int checkLine(char* lineBuffer);
int addOffTarget(OffTarget *offTargetList, int offTargetListSize, char* lineBuffer, int BalchFilter, int MismatchFilter);
void sortOffTargetLineByInx(OffTarget *offTargetList, int offTargetListSize);

int main(int argc, char* argv[]) {
    int BalchFilter = (int)atoi(argv[1]);
    int MismatchFilter = (int)atoi(argv[2]);
    char status;
    status = *argv[3];
    if (status == 'r') {
        char *command = "./cas-offinder input.txt G output.txt";
        system(command);
    }
    FILE *OutputFile = fopen(CasOffinderPath"/output.txt", "r");
    FILE *FilteredOutputFile = fopen(CasOffinderPath"/output2.txt", "w");
    OffTarget *offTargetList = (OffTarget *)malloc(800*sizeof(OffTarget));
    int offTargetListSize = 0;
    char lineBuffer[300];
    fgets(lineBuffer,300,OutputFile);
    fputs(lineBuffer,FilteredOutputFile); // print the header line
    while (fgets(lineBuffer,300,OutputFile) != NULL){
        if (checkLine(lineBuffer)){
            //fputs(lineBuffer,FilteredOutputFile); - print only goot off-targets to the output2 file
            //puts(lineBuffer);
            offTargetListSize = addOffTarget(offTargetList, offTargetListSize, lineBuffer, BalchFilter, MismatchFilter);
        }
    }
    sortOffTargetLineByInx(offTargetList,offTargetListSize);
    for (int i = 0; i < offTargetListSize; ++i) {
        if (offTargetList[i].Reverse == '+') {
            fprintf(FilteredOutputFile, "Inx : %10d\t Target: %23s\tdistance : %d\tMismatch : %d\tBalch : %d\tReverse : %c\n", offTargetList[i].inx,
                    offTargetList[i].TextTarget, offTargetList[i].distance, offTargetList[i].mismatch,  offTargetList[i].balch, offTargetList[i].Reverse);
            //printf("Inx : %10d\tdistance : %d\tReverse : %c\n", offTargetList[i].inx, offTargetList[i].distance,
            //       offTargetList[i].Reverse);
        }
    }
    fclose(OutputFile);
    fclose(FilteredOutputFile);
    return 0;
}

int checkLine(char* lineBuffer){
    int i, lineLength = (int)strlen(lineBuffer);
    int tabCounter = 0;
    //read line from the 3 \t to the 4 \t if there is 'n' in this part return 0 else 1
    for (i=0;i<lineLength; i++){
        if (lineBuffer[i] == '\t'){
            tabCounter++;
        }
        if (tabCounter == 3){
            if (lineBuffer[i] == 'n'){
                return 0;
            }
        }
        if (tabCounter == 4){
            return 1;
        }
    }
    return 0;
}

int addOffTarget(OffTarget *offTargetList, int offTargetListSize, char* lineBuffer, int BalchFilter, int MismatchFilter){
    int i, lineLength = (int)strlen(lineBuffer);
    int Flag = 0;
    int tabCounter = 0;
    int TargetInx = 0;
    int distance = 0, mismatch = 0, balch = 0;
    char TextTarget[MAX_PATTERN_LENGTH+1];
    int TextTargetInx = 0;
    char Reverse;
    for (i=0;i<lineLength; i++) {
        if (lineBuffer[i] == '\t') {
            tabCounter++;
        }else if ((tabCounter == 2) && (TextTargetInx < MAX_PATTERN_LENGTH)) {
            TextTarget[TextTargetInx] = lineBuffer[i];
            TextTargetInx++;
        } else if (tabCounter == 5 && Flag == 0) {
            TargetInx = atoi(lineBuffer + i);
            Flag = 1;
        } else if (tabCounter == 6) {
            Reverse = lineBuffer[i];
        } else if (tabCounter == 7) {
            mismatch = atoi(lineBuffer + i);
            distance = mismatch;
        } else if (tabCounter == 8) {
            balch = balch + atoi(lineBuffer + i);
            distance = distance + balch;
            tabCounter ++; // for end of line behevior
        }
    }
    TextTarget[MAX_PATTERN_LENGTH] = '\0';
    if ((mismatch > MismatchFilter) || (balch > BalchFilter)){return offTargetListSize;}
    //if ((mismatch > MAX_MISMATCH) || (balch > MAX_BALCH)){return offTargetListSize;}
    for (i=0;i<offTargetListSize; i++){
        if (offTargetList[i].inx == TargetInx && offTargetList[i].Reverse == Reverse){
            if (offTargetList[i].distance > distance){
                offTargetList[i].distance = distance;
                offTargetList[i].mismatch = mismatch;
                offTargetList[i].balch = balch;
                strcpy(offTargetList[i].TextTarget, TextTarget);
            }
            return offTargetListSize;
        }
    }
    if (i==offTargetListSize){
        offTargetList[i].inx = TargetInx;
        offTargetList[i].distance = distance;
        offTargetList[i].mismatch = mismatch;
        offTargetList[i].balch = balch;
        offTargetList[i].Reverse = Reverse;
        strcpy(offTargetList[i].TextTarget, TextTarget);
        offTargetListSize++;
    }
    return offTargetListSize;
}

void sortOffTargetLineByInx(OffTarget *offTargetList, int offTargetListSize){
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


