//
// Created by pico on 02/03/22.
//

#ifndef MAIN_UTILS_H
#define MAIN_UTILS_H

#include <stdio.h>
#include "typedef_utils.h"

/* function declaration */
void printStructInput(TypeInput input);
void printDoubleArray(double array[], size_t n);
void printIntArray(int array[], size_t n);
void setArrayDoubleValues(double array_out[], double array_in[], size_t n);
void setArrayIntValues(int array_out[], int array_in[], size_t n);
void repeatDoubleArray(double array_in[], double array_out[], int n_in, int n_out);
void identityMatrixVector(double vector[], int n);
void setMatrixDoubleValues(int n_row, int n_cols, double matrix_out[n_row][n_cols], double matrix_in[n_row][n_cols]);
void printDoubleArrayLine(double array[], size_t n);
void readRef(int n_ref, double ref[n_ref][6]);

#endif //MAIN_UTILS_H
