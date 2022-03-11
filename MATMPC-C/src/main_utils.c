//
// Created by pico on 02/03/22.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "main_utils.h"
#include "typedef_utils.h"

/* utils functions */
void printStructInput(TypeInput input){
    printf("Size of input.x0 is %d \n", (int)sizeof(input));
    printDoubleArray(input.x0, (int)( sizeof(input.x0) / sizeof(input.x0[0])));
    printf("Size of input.x is %d \n", (int)sizeof(input));
    printDoubleArray(input.x, (int)( sizeof(input.x) / sizeof(input.x[0])));
}

void printDoubleArray(double array[], size_t n){
    printf("No. of elements in the array is %d \n", (int)n);
    for (int i = 0; i < n; i++) {
        printf("Element %d of the array is %f \n", i, array[i]);
    }
}

void printDoubleArrayLine(double array[], size_t n){
    printf("[");
    for (int i = 0; i < n-1; i++) {
        printf("%f, ", array[i]);
    }
    printf("%f]", array[n-1]);
}

void printIntArray(int array[], size_t n){
    printf("No. of elements in the array is %d \n", (int)n);
    for (int i = 0; i < n; i++) {
        printf("Element %d of the array is %d \n", i, array[i]);
    }
}

void setArrayDoubleValues(double array_out[], double array_in[], size_t n) {
    for (int i = 0; i < n; i++) {
        array_out[i] = array_in[i];
    }
}

void setMatrixDoubleValues(int n_row, int n_cols, double matrix_out[n_row][n_cols], double matrix_in[n_row][n_cols]) {
    for (int i = 0; i < n_row; i++) {
        for (int j = 0; j<n_cols;j++) {
            matrix_out[i][j] = matrix_in[i][j];
        }
    }
}

void setArrayIntValues(int array_out[], int array_in[], size_t n) {
    for (int i = 0; i < n; i++) {
        array_out[i] = array_in[i];
    }
}

void repeatDoubleArray(double array_in[], double array_out[], int n_in, int n_out){
    for(int i=0;i<n_out;i++) {
        int i_in = i%n_in;
        //printf("Copying element no. %d of the input array to the %d th element of the out array, which is %f \n", i_in, i, array_in[i_in]);
        array_out[i] = array_in[i_in];
        //printf("Element %d th of the out array is %f \n", i, array_out[i]);
    }
}

void identityMatrixVector(double vector[], int n){
    for (int i=0;i<n;i++){
        vector[i*(n+1)] = 1;
    }
}

void readRef(int n_ref, double ref[n_ref][6]){
    FILE *fp;
    char buff[255];
    double tmp;

    // read drd file
    fp = fopen("../csv/ref.csv", "r");
    fgets(buff, 255, (FILE*)fp);
//    printf("skipping header line: %s\n", buff );
    for (int i=0; i<n_ref; i++) {
        for (int j = 0; j < 6; j++) {
            fscanf(fp, "%s", buff);
//            printf("Reading value : %s", buff);
            tmp = atof(buff);
//            printf(" converted as number in: %f\n", tmp);
            ref[i][j] = tmp;
        }
    }
    fclose(fp);
}