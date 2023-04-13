#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h> // for lower case only
#include <time.h> // for time measurement

//--------------------------------------------
//        Definition of Constants
//--------------------------------------------
// run parameters
#define MAX_LINE_SIZE 100

// pattern parameters
#define ALPHABET_SIZE 4
#define CHAR_TO_MASK(char) (char == 'A' ? 0 : char == 'C' ? 1 : char == 'G' ? 2 : char == 'T' ? 3 : 4)
#define COMPLEMENT_NUCLEOTIDE(char) (char == 'A' ? 'T' : char == 'C' ? 'G' : char == 'G' ? 'C' : char == 'T' ? 'A' : char)

//structs
typedef struct {
    int ChromosomePosition;
    int ChromosomeNum;
    int Distance;
    int Bulge;
    int Mismatch;
    char Strand;
    char *Guide;
    char *SiteAlignment;
    char *GuideAlignment;
    struct OffTarget *NextOffTarget;
} OffTarget;

typedef struct {
    unsigned long *MatchVectors;
    unsigned long *MismatchVectors;
    unsigned long *InsertionVectors;
    unsigned long *DeletionVectors;
} BitapMatrix;

typedef struct {
    char *Read;  // guide RNA - for result printing
    int Length;
    char Strand; // '+' - forward, '-' - reverse
    unsigned long BitMaskVectors[ALPHABET_SIZE];
    unsigned long *DistanceVectors;
    BitapMatrix **LastBitapMatrices;
    int *EncodeAlignment;
} Guide;

typedef struct {
    char *Text;
    int TextInx;
    int ChromosomeNum;
} ChromosomeInfo;


//--------------------------------------------
//        Global variables
//--------------------------------------------
int NumOfChromosomes = 0;
int NumOfGuides = 0;
int NumOfOffTargets = 0;
unsigned long MsbMaskVectorConst = 0;
// distance parameters
int max_mismatch = -1;
int max_bulge = -1;
int max_distance = -1;
// file & memory parameters
char *output_file_path = "./output.txt";
char *chromosome_file_path = "./chr1.txt";
char *guide_file_path = "./guides.txt";
long int max_chromosome_size = 300000000;
//int max_line_size = 100;


//--------------------------------------------
//        Prototypes
//--------------------------------------------
// Genome Read
FILE *OpenNextChromosomeFile();
int ReadChromosome(ChromosomeInfo *Chromosome);
// Guide Initialization
Guide **InitializeGuidesInfo();
Guide **AllocateAndReadGuide(FILE *GuideFile);
void AllocateAndSetReverseGuide(Guide **guideLst);
void InitializeBitMaskVectors(Guide *GuideInfo);
void InitializeBitVectorsAndMatrix(Guide *GuideInfo);
// OffTarget Finder
void OffFinderMainLoop(Guide **guideLst, ChromosomeInfo *Chromosome, OffTarget **OffTargetHead);
void DistanceVectorCalc(unsigned long BitMaskVector, unsigned long *DistanceVectors, BitapMatrix *BitapMatrix);
int CheckForAlignment(Guide *GuideInfo, OffTarget *offTarget, int ChromosomeInx);
// OffTarget TraceBack
int TargetTraceBack(Guide *GuideInfo, OffTarget *offTarget, int MatrixInx, unsigned long MaskBitVector, int GuideInx, int SiteInx, int CurDistance, int CurMismatch, int CurBulge, int AlignmentInx);
void AddOffTargetToList(Guide *GuideInfo, ChromosomeInfo *Chromosome, OffTarget **OffTargetHead, OffTarget *offTargetToAdd);
void DecodeAlignment(Guide *GuideInfo, ChromosomeInfo *Chromosome, OffTarget *offTarget);
// Post Processing
void PrintOffTargets(OffTarget *OffTargetHead);
void FreeAllMemory(ChromosomeInfo *Chromosome, Guide **guideLst, OffTarget *OffTargetHead);

//--------------------------------------------
//              Main
//--------------------------------------------

int main(int argc, char* argv[]){

    for (int i = 1; i < argc; i+=2) {
        if (strcmp(argv[i], "-m") == 0 && i+1 < argc) {
            max_mismatch = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-b") == 0 && i+1 < argc) {
            max_bulge = atoi(argv[i+1]);
        } else if (strcmp(argv[i], "-output") == 0 && i+1 < argc) {
            output_file_path = argv[i+1];
        } else if (strcmp(argv[i], "-chr") == 0 && i+1 < argc) {
            chromosome_file_path = argv[i+1];
        } else if (strcmp(argv[i], "-guide") == 0 && i+1 < argc) {
            guide_file_path = argv[i+1];
        } else if (strcmp(argv[i], "-max_chr") == 0 && i+1 < argc) {
            max_chromosome_size = atol(argv[i+1]);
        } else {
            printf("Unknown Flag or Value was not given: %s\n", argv[i]);
            exit(1);
        }
    }
    if (max_mismatch == -1 || max_bulge == -1) {
        printf("Error: max_mismatch or max_bulge was not provide\n");
        exit(1);
    } else {
        max_distance = max_mismatch + max_bulge;
    }

    clock_t start = clock();

    clock_t chr_start;
    double elapsed;

    MsbMaskVectorConst = (unsigned long)pow(2,sizeof(unsigned long)*8-1);; // 10000000 vector - to reduce calculation
    ChromosomeInfo *Chromosome = (ChromosomeInfo *) malloc(sizeof(ChromosomeInfo));
    Guide **guideLst = InitializeGuidesInfo();
    OffTarget *OffTargetHead = NULL;

    while(ReadChromosome(Chromosome) == 0){ // run foe each chromosome
        chr_start = clock();
        printf("Start search in Chromosome %d\n", Chromosome->ChromosomeNum);
        OffFinderMainLoop(guideLst, Chromosome, &OffTargetHead);
        elapsed = (double) (clock() - chr_start) / CLOCKS_PER_SEC;
        printf("Chromosome %d search time: %f seconds\n\n",Chromosome->ChromosomeNum, elapsed);
    }

    elapsed = (double) (clock() - start) / CLOCKS_PER_SEC;
    printf("Search is done. total search time: %f seconds\n", elapsed);

    PrintOffTargets(OffTargetHead);
    FreeAllMemory(Chromosome, guideLst, OffTargetHead);
    return 0;
}

//--------------------------------------------
//        Genome Read functions
//--------------------------------------------

FILE *OpenNextChromosomeFile() {
    char *FilePath = (char *)malloc(sizeof(char)* strlen(chromosome_file_path));
    strcpy(FilePath, chromosome_file_path);
    FilePath[strlen(FilePath)-5] = '1' + NumOfChromosomes;
    FILE *ChrFile = fopen(FilePath, "r");
    if (ChrFile == NULL && NumOfChromosomes == 0) {
        printf("Error: 1st Chromosome file not found\n");
        printf("File Name provide was %s\n",chromosome_file_path);
        exit(1);
    }
    free(FilePath);
    if (ChrFile == NULL) {
        return NULL;
    }
    NumOfChromosomes++;
    return ChrFile;
}

int ReadChromosome(ChromosomeInfo *Chromosome) {
    FILE *ChrFile = OpenNextChromosomeFile();
    if (ChrFile == NULL) {
        return -1;
    }
    Chromosome->Text = (char *) malloc(max_chromosome_size * sizeof(char));
    Chromosome->TextInx = 0;
    Chromosome->ChromosomeNum = NumOfChromosomes;
    char lineBuffer[MAX_LINE_SIZE];
    while (fscanf(ChrFile, "%[^\n]c", lineBuffer) != EOF) { //read until new line
        fscanf(ChrFile, "%*c"); //skip '\n'
        strcpy(Chromosome->Text + Chromosome->TextInx, lineBuffer);
        Chromosome->TextInx = Chromosome->TextInx + (int) strlen(lineBuffer);
    }
    fclose(ChrFile);
    return 0;
}

//--------------------------------------------
//    Initialize Guide Info functions
//--------------------------------------------

Guide **InitializeGuidesInfo(){
    FILE *GuideFile = fopen(guide_file_path, "r");
    if (GuideFile == NULL) {
        printf("Error: Guide file not found\n");
        printf("File Name provide was %s\n",guide_file_path);
        exit(1);
    }
    Guide **GuideLst = AllocateAndReadGuide(GuideFile);
    if (NumOfGuides == 0){
        printf("No guide was found");
        exit(1);
    }
    fclose(GuideFile);
    AllocateAndSetReverseGuide(GuideLst);
    for (int i=0; i<NumOfGuides; i++){
        InitializeBitMaskVectors(GuideLst[i]);
        InitializeBitVectorsAndMatrix(GuideLst[i]);
        GuideLst[i]->EncodeAlignment = (int *)malloc(sizeof(int)*(GuideLst[i]->Length+max_bulge+1));
    }
    return GuideLst;
}

Guide **AllocateAndReadGuide(FILE *GuideFile) {
    int length = 0;
    char GuideBuffer[MAX_LINE_SIZE]; // consider bigger size
    Guide *GuideInfo;
    while (fscanf(GuideFile, "%[^\n]c", GuideBuffer) != EOF) { //read until new line
        fscanf(GuideFile, "%*c"); //skip '\n'
        if (GuideBuffer[0] != '\n') {
            NumOfGuides++;
        }
    }
    rewind(GuideFile);
    Guide **GuideList = (Guide **) malloc((NumOfGuides * 2) * sizeof(Guide *)); // each guide have reverse complement guide
    for (int i=0; i<NumOfGuides; i++){
        GuideInfo = (Guide *) malloc(sizeof(Guide));
        fscanf(GuideFile, "%s", GuideBuffer);
        fscanf(GuideFile, "%*c"); //skip '\n'
        length = (int) strlen(GuideBuffer);
        GuideInfo->Read = (char *) malloc((length + 1) * sizeof(char));
        strcpy(GuideInfo->Read, GuideBuffer);
        GuideInfo->Length = length;
        GuideInfo->Strand = '+';
        GuideList[i*2] = GuideInfo; // forward guide in all even places
    }
    return GuideList;
}

void AllocateAndSetReverseGuide(Guide **GuideLst){
    int NumOfReverseGuides = NumOfGuides;
    Guide *ReverseGuideInfo;
    for (int i=0; i<NumOfReverseGuides; i++){
        ReverseGuideInfo = (Guide *) malloc(sizeof(Guide));
        ReverseGuideInfo->Length = GuideLst[i*2]->Length;
        ReverseGuideInfo->Strand = '-';
        ReverseGuideInfo->Read = (char *) malloc((ReverseGuideInfo->Length + 1) * sizeof(char));
        for (int j=0; j<ReverseGuideInfo->Length; j++){
            ReverseGuideInfo->Read[j] = COMPLEMENT_NUCLEOTIDE(GuideLst[i*2]->Read[ReverseGuideInfo->Length - j - 1]);
        }
        NumOfGuides++;
        GuideLst[i*2 + 1] = ReverseGuideInfo; // reverse guide in all odd places
    }
}

void InitializeBitMaskVectors(Guide *GuideInfo){ //TODO: reverse all bit! msb is the 0th bit
    unsigned long MaskVector = MsbMaskVectorConst; // 10000000 vector
    unsigned long initPMValue = (MsbMaskVectorConst >> (GuideInfo->Length-1))-1;// 00000111 vector
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        GuideInfo->BitMaskVectors[i] = initPMValue;
    }
    for (int i=0; i<GuideInfo->Length; i++){
        if (GuideInfo->Read[i] == 'A'){
            GuideInfo->BitMaskVectors[CHAR_TO_MASK('A')] = GuideInfo->BitMaskVectors[CHAR_TO_MASK('A')] | MaskVector;
        } else if (GuideInfo->Read[i] == 'T'){
            GuideInfo->BitMaskVectors[CHAR_TO_MASK('T')] = GuideInfo->BitMaskVectors[CHAR_TO_MASK('T')] | MaskVector;
        } else if (GuideInfo->Read[i] == 'G'){
            GuideInfo->BitMaskVectors[CHAR_TO_MASK('G')] = GuideInfo->BitMaskVectors[CHAR_TO_MASK('G')] | MaskVector;
        } else if (GuideInfo->Read[i] == 'C'){
            GuideInfo->BitMaskVectors[CHAR_TO_MASK('C')] = GuideInfo->BitMaskVectors[CHAR_TO_MASK('C')] | MaskVector;
        };
        MaskVector = MaskVector >> 1;
    }
    for (int i=0;i<ALPHABET_SIZE;i++){
        GuideInfo->BitMaskVectors[i] = ~GuideInfo->BitMaskVectors[i]; // forward
    }
}

void InitializeBitVectorsAndMatrix(Guide *GuideInfo){
    GuideInfo->DistanceVectors = (unsigned long *)malloc(sizeof(unsigned long)*(max_distance+1));
    for (int i = 0; i <= max_distance; i++) {
        GuideInfo->DistanceVectors[i] = -1;
        GuideInfo->DistanceVectors[i] = GuideInfo->DistanceVectors[i] << (sizeof(unsigned long) * 8 - GuideInfo->Length);
    }
    GuideInfo->LastBitapMatrices = (BitapMatrix **)malloc((GuideInfo->Length+1)*sizeof(BitapMatrix *));
    for (int i = 0; i < GuideInfo->Length + 1; i++) {
        GuideInfo->LastBitapMatrices[i] = (BitapMatrix *) malloc(sizeof(BitapMatrix));
        GuideInfo->LastBitapMatrices[i]->MatchVectors = (unsigned long *)malloc(sizeof(unsigned long)*(max_distance+1));
        GuideInfo->LastBitapMatrices[i]->MatchVectors[0] = -1;
        GuideInfo->LastBitapMatrices[i]->MismatchVectors = (unsigned long *)malloc(sizeof(unsigned long)*(max_distance+1));
        GuideInfo->LastBitapMatrices[i]->MismatchVectors[0] = -1;
        GuideInfo->LastBitapMatrices[i]->InsertionVectors = (unsigned long *)malloc(sizeof(unsigned long)*(max_distance+1));
        GuideInfo->LastBitapMatrices[i]->InsertionVectors[0] = -1;
        GuideInfo->LastBitapMatrices[i]->DeletionVectors = (unsigned long *)malloc(sizeof(unsigned long)*(max_distance+1));
        GuideInfo->LastBitapMatrices[i]->DeletionVectors[0] = -1;
    }
}

//--------------------------------------------
//    Edit Distance Calculation functions
//--------------------------------------------

void OffFinderMainLoop(Guide **guideLst, ChromosomeInfo *Chromosome, OffTarget **OffTargetHead){
    char nucleotide;
    unsigned long BitMaskVector;
    OffTarget *tempOffTarget = (OffTarget *)malloc(sizeof(OffTarget));
    while (Chromosome->TextInx > 0) { // run from tail to head of chromosome
        nucleotide = Chromosome->Text[Chromosome->TextInx];
        if (nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'G' || nucleotide == 'T') { // if nucleotide is valid
            for (int GuideInx = 0; GuideInx < NumOfGuides; GuideInx++) {
                BitMaskVector =  guideLst[GuideInx]->BitMaskVectors[CHAR_TO_MASK(nucleotide)];
                //TODO: efficient: 1. calc MatrixInx once for all guide 2. pass init matrix value to TB here.
                DistanceVectorCalc(BitMaskVector, guideLst[GuideInx]->DistanceVectors,
                                   guideLst[GuideInx]->LastBitapMatrices[Chromosome->TextInx%(guideLst[GuideInx]->Length+1)]);
                if (CheckForAlignment(guideLst[GuideInx], tempOffTarget, Chromosome->TextInx) ==0) {
                    AddOffTargetToList(guideLst[GuideInx], Chromosome, OffTargetHead, tempOffTarget);
                    tempOffTarget = (OffTarget *)malloc(sizeof(OffTarget));
                }
            }
        }
        Chromosome->TextInx--;
    }
    free(tempOffTarget);
}

void DistanceVectorCalc(unsigned long BitMaskVector, unsigned long *DistanceVectors, BitapMatrix *BitapMatrix){
    unsigned long InsertionVector, DeletionVector, MatchVector, MismatchVector;
    unsigned long LastDistanceVectors = DistanceVectors[0];
    DistanceVectors[0] = (DistanceVectors[0]<<1 | BitMaskVector);
    BitapMatrix->MatchVectors[0] = DistanceVectors[0];
    for (int d=1;d<=max_distance;d++){
        MatchVector = (DistanceVectors[d]<<1 | BitMaskVector);
        MismatchVector = LastDistanceVectors<<1;
        DeletionVector = LastDistanceVectors;
        InsertionVector = DistanceVectors[d-1]<<1;
        BitapMatrix->MatchVectors[d] = MatchVector;
        BitapMatrix->MismatchVectors[d] = MismatchVector;
        BitapMatrix->DeletionVectors[d] = DeletionVector;
        BitapMatrix->InsertionVectors[d] = InsertionVector;
        LastDistanceVectors = DistanceVectors[d];  // save RdMinus1Vector for next iteration
        DistanceVectors[d] = MatchVector & MismatchVector & InsertionVector & DeletionVector;
    }
}

int CheckForAlignment(Guide *GuideInfo, OffTarget *offTarget, int ChromosomeInx){
    unsigned long MsbMaskVector = MsbMaskVectorConst; // TODO: efficiency tip: change to 1!!
    for (int Distance = 0; Distance <= max_distance; Distance++) {
        if ((GuideInfo->DistanceVectors[Distance] & MsbMaskVector) == 0) { // Edit distance match threshold
            if (TargetTraceBack(GuideInfo, offTarget, (ChromosomeInx)%(GuideInfo->Length+1), MsbMaskVector,
                                0, GuideInfo->Length, Distance, 0, 0, 0)==0) {
                return 0;
            }
        }
    }
    return -1;
}

//--------------------------------------------
//           TraceBack functions
//--------------------------------------------

int TargetTraceBack(Guide *GuideInfo, OffTarget *offTarget, int MatrixInx, unsigned long MaskBitVector, int GuideInx, int SiteInx, int CurDistance, int CurMismatch, int CurBulge, int AlignmentInx){
    while ((GuideInx <= GuideInfo->Length) & (SiteInx >= 0)){
        if ((GuideInfo->LastBitapMatrices[MatrixInx]->MatchVectors[CurDistance] & MaskBitVector) == 0) { //Match
            GuideInx++;
            SiteInx--;
            MaskBitVector = MaskBitVector >> 1;
            MatrixInx=(MatrixInx+1)%(GuideInfo->Length+1);
            GuideInfo->EncodeAlignment[AlignmentInx++] = 1; // Match
        } else { // check for mismatch
            if (CurMismatch < max_mismatch & (GuideInfo->LastBitapMatrices[MatrixInx]->MismatchVectors[CurDistance] & MaskBitVector) == 0) {
                if (TargetTraceBack(GuideInfo, offTarget, (MatrixInx+1)%(GuideInfo->Length+1), MaskBitVector >> 1,
                                    GuideInx+1, SiteInx-1, CurDistance-1, CurMismatch+1, CurBulge, AlignmentInx+1)==0){
                    GuideInfo->EncodeAlignment[AlignmentInx] = 2; // Mismatch
                    return 0;
                }
            } // check for deletion
            if (CurBulge < max_bulge & (GuideInfo->LastBitapMatrices[MatrixInx]->InsertionVectors[CurDistance] & MaskBitVector) == 0) {
                if (TargetTraceBack(GuideInfo, offTarget, MatrixInx, MaskBitVector >> 1,
                                    GuideInx, SiteInx-1, CurDistance-1, CurMismatch, CurBulge+1, AlignmentInx+1)==0){
                    GuideInfo->EncodeAlignment[AlignmentInx] = 3; // Deletion
                    return 0;
                }
            } // check for insertion
            if (CurBulge < max_bulge & (GuideInfo->LastBitapMatrices[MatrixInx]->DeletionVectors[CurDistance] & MaskBitVector) == 0) {
                if (TargetTraceBack(GuideInfo, offTarget, (MatrixInx+1)%(GuideInfo->Length+1), MaskBitVector,
                                    GuideInx+1, SiteInx, CurDistance-1, CurMismatch, CurBulge+1, AlignmentInx+1)==0){
                    GuideInfo->EncodeAlignment[AlignmentInx] = 4; // Insertion
                    return 0;
                }
            }
            return -1;
        }
    }
    //update OffTarget data
    offTarget->Bulge = CurBulge;
    offTarget->Mismatch = CurMismatch;
    offTarget->Distance = CurMismatch + CurBulge;
    GuideInfo->EncodeAlignment[AlignmentInx++] = 1; // Match
    return 0;
}

void AddOffTargetToList(Guide *GuideInfo, ChromosomeInfo *Chromosome, OffTarget **OffTargetHead, OffTarget *offTargetToAdd){
    offTargetToAdd->ChromosomePosition = Chromosome->TextInx;
    offTargetToAdd->ChromosomeNum = Chromosome->ChromosomeNum;
    offTargetToAdd->Strand = GuideInfo->Strand;
    offTargetToAdd->Guide = (char *)malloc((GuideInfo->Length+1)*sizeof(char));
    strcpy(offTargetToAdd->Guide, GuideInfo->Read);
    DecodeAlignment(GuideInfo, Chromosome, offTargetToAdd);
    // TODO: efficiency - maybe write to file now?
    NumOfOffTargets++;
    offTargetToAdd->NextOffTarget = (struct OffTarget *) *OffTargetHead;
    *OffTargetHead = offTargetToAdd;
}

void DecodeAlignment(Guide *GuideInfo, ChromosomeInfo *Chromosome, OffTarget *offTarget){
    int AlignmentLength = 0;
    int GuideInx = 0;
    int SiteInx = Chromosome->TextInx;
    char SiteAlignment[GuideInfo->Length + offTarget->Bulge + 1];
    char GuideAlignment[GuideInfo->Length + offTarget->Bulge + 1];
    for (int i=0; i<GuideInfo->Length; i++){
        if (GuideInfo->EncodeAlignment[AlignmentLength]==4) i--;
        AlignmentLength++;
    }
    for (int i=0; i<AlignmentLength; i++) {
        switch (GuideInfo->EncodeAlignment[i]){
            case 1: // Match
                SiteAlignment[i] = Chromosome->Text[SiteInx++];
                GuideAlignment[i] = GuideInfo->Read[GuideInx++];
                break;
            case 2: // Mismatch
                SiteAlignment[i] = (char)tolower(Chromosome->Text[SiteInx++]);// + 32;
                GuideAlignment[i] = GuideInfo->Read[GuideInx++];
                break;
            case 3: // Deletion
                SiteAlignment[i] = '-';
                GuideAlignment[i] = GuideInfo->Read[GuideInx++];
                break;
            case 4: // Insertion
                SiteAlignment[i] = Chromosome->Text[SiteInx++];
                GuideAlignment[i] = '-';
                break;
            default:
                printf("Error in DecodeAlignment\n");
                break;
        }
        GuideInfo->EncodeAlignment[i] = 0;
    }
    SiteAlignment[AlignmentLength] = '\0';
    GuideAlignment[AlignmentLength] = '\0';
    offTarget->SiteAlignment = (char *)malloc((AlignmentLength+1)*sizeof(char)); // +1?!
    strcpy(offTarget->SiteAlignment, SiteAlignment);
    offTarget->GuideAlignment = (char *)malloc((AlignmentLength+1)*sizeof(char));
    strcpy(offTarget->GuideAlignment, GuideAlignment);
}

//--------------------------------------------
//         Post Processing functions
//--------------------------------------------

void PrintOffTargets(OffTarget *OffTargetHead){
    FILE *OutputFile = fopen(output_file_path, "w");
    if (OutputFile == NULL) {
        printf("Error opening file output file");
        exit(1);
    }
    int GuideLengthPrint = (int)strlen(OffTargetHead->Guide)+max_bulge;
    fprintf(OutputFile,"%*s\t%-*s\t%-*s\t%-*s\t%-*s\t%-*s\t%-s\t%-s\t%-s\n",
            10,"Chromosome",10,"Index",6,"Strand",GuideLengthPrint,"Guide",
            GuideLengthPrint,"SiteAlignment",GuideLengthPrint,"GuideAlignment",
            "Distance","Mismatch","Bulge");
    while (OffTargetHead != NULL) {
        fprintf(OutputFile,"Chr%-7d\t%10d\t%c\t%-*s\t%-*s\t%-*s\t%-8d\t%-8d\t%-8d\n",
                OffTargetHead->ChromosomeNum, OffTargetHead->ChromosomePosition, OffTargetHead->Strand,
                GuideLengthPrint, OffTargetHead->Guide, GuideLengthPrint, OffTargetHead->SiteAlignment,
                GuideLengthPrint, OffTargetHead->GuideAlignment,
                OffTargetHead->Distance, OffTargetHead->Mismatch, OffTargetHead->Bulge);
        OffTargetHead = (OffTarget *) OffTargetHead->NextOffTarget;
    }
    fclose(OutputFile);
    printf("\n");
    printf("OffFinder v4.0  (c)\n");
    printf("Output Summery:\n");
    printf("Number of guides (2 Strand per guide):%d\n", NumOfGuides/2);
    printf("Number of chromosomes: %d\n", NumOfChromosomes);
    printf("Number of off-targets found: %d\n\n", NumOfOffTargets);
}

void FreeAllMemory(ChromosomeInfo *Chromosome, Guide **guideLst, OffTarget *OffTargetHead){
    // free guides
    for (int i = 0; i < NumOfGuides; ++i) {
        free(guideLst[i]->EncodeAlignment);
        free(guideLst[i]->Read);
        free(guideLst[i]->DistanceVectors);
        for (int j = 0; j < guideLst[i]->Length + 1; j++) {
            free(guideLst[i]->LastBitapMatrices[j]->MatchVectors);
            free(guideLst[i]->LastBitapMatrices[j]->MismatchVectors);
            free(guideLst[i]->LastBitapMatrices[j]->InsertionVectors);
            free(guideLst[i]->LastBitapMatrices[j]->DeletionVectors);
            free(guideLst[i]->LastBitapMatrices[j]);
        }
        free(guideLst[i]);
    }
    free(guideLst);
    // free off-targets
    while (OffTargetHead != NULL) {
        OffTarget *tmp = (OffTarget *) OffTargetHead->NextOffTarget;
        free(OffTargetHead->Guide);
        free(OffTargetHead->SiteAlignment);
        free(OffTargetHead);
        OffTargetHead = tmp;
    }
    // free chromosomes
    free(Chromosome->Text);
    free(Chromosome);
}