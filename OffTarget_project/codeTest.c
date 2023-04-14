//
// Created by Amichai Malle on 20/03/2023.
//
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>

#define MAX_BULGE 2

int main() {
    char b = 0x11;
   for (int i=0; i<4; i++){
         printf("%d\n",b);
         b = b << 1;
   }
    return 0;
}

/*int main(){
    FILE *fp = fopen("/Users/amichaim/CLionProjects/OffTarget/OffTarget_project/test.txt", "r");
    if (fp == NULL){
        printf("Error opening file");
        exit(1);
    }
    int Total = 0;
    char **lineBuffer = (char **) malloc(sizeof(char *));
    char *buffer;
    int a;
    //fscanf(fp,"%s",lineBuffer);
    //printf("%s\n",lineBuffer);
    //----------------------------------------------
    int length = 0;
    char charBuffer;
    //fscanf(fp,"%c",&charBuffer);
    while (fscanf(fp,"%c",&charBuffer) != EOF){ // read all guide
    //do { // read all guide
        if (charBuffer == '\n' || charBuffer == ' ') {
            a = length;
            buffer = (char *) malloc((a+1) * sizeof(char));
            lineBuffer = (char **) realloc(lineBuffer, (++Total) * sizeof(char *));
            lineBuffer[Total-1] = buffer;
            fseek(fp, -(a+1), SEEK_CUR);
            fscanf(fp,"%s",buffer);
            printf("length = %d\n",a);
            fscanf(fp,"%[^a-zA-Z]c",&charBuffer);
            length = 0;
        } else {
            length++;
        }
    //} while (fscanf(fp,"%c",&charBuffer) != EOF); // read all guide
    }
    a = length;
    buffer = (char *) malloc((a+1) * sizeof(char));
    lineBuffer = (char **) realloc(lineBuffer, (++Total) * sizeof(char *));
    lineBuffer[Total-1] = buffer;
    fseek(fp, -(a), SEEK_CUR);
    fscanf(fp,"%s",buffer);
    printf("length = %d\n",a);
    printf("total lines = %d\n",Total);
    for (int i = 0; i < Total; i++){
        printf("%s\n",lineBuffer[i]);
    }
    return 0;
}*/