/* matrix.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "numeric.h"
#include "algebra.h"

#include "memwrapper.h"

void initMatrix(matrix **m)
{
  (*m) = xmalloc(sizeof(matrix));
  (*m)->row = 0;
  (*m)->col = 0;
  (*m)->data = NULL;
}

void NewMatrix(matrix **m, size_t row_ , size_t col_)
{
  size_t i, j;
  (*m) = xmalloc(sizeof(matrix));
  (*m)->row = row_;
  (*m)->col = col_;
  (*m)->data = xmalloc(sizeof(double*)*row_);
  for(i = 0; i < row_; i++){
    (*m)->data[i] = xmalloc(sizeof(double)*col_);
    for(j = 0; j < col_; j++)
      (*m)->data[i][j] = +0.f;
  }
}

void ResizeMatrix(matrix **m, size_t row_, size_t col_)
{
  size_t i, j;
  if(m != NULL){
    if((*m)->row == row_ && (*m)->col == col_){
      MatrixSet((*m), +0.f);
    }
    else{
      if((*m)->col > 0 && (*m)->row > 0){
        for(i = 0; i < (*m)->row; i++){
          xfree((*m)->data[i]);
        }
        xfree((*m)->data);
      }

      (*m)->data = xmalloc(sizeof(double*)*row_);
      for(i = 0; i < row_; i++){
        (*m)->data[i] = xmalloc(sizeof(double)*col_);
        for(j = 0; j < col_; j++)
          (*m)->data[i][j] = +0.f;
      }

      (*m)->row = row_;
      (*m)->col = col_;
    }
  }
  else{
    NewMatrix(m, row_, col_);
  }
}

void DelMatrix(matrix **m)
{
  if(m != NULL){
    size_t i;
    for(i = 0; i < (*m)->row; i++)
      xfree((*m)->data[i]);
    xfree((*m)->data);
    xfree((*m));
  }
}

void MatrixCheck(matrix *m)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(_isnan_(m->data[i][j]) || !_isfinite_(m->data[i][j])){
        m->data[i][j] = MISSING;
      }
    }
  }
}

void FindNan(matrix *m)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      if(_isnan_(m->data[i][j])){
        printf("Nan at %d %d\n", (int)i, (int)j);
      }
    }
  }
}

void PrintMatrix(matrix *m)
{
  size_t i, j;
  printf("Matrix of row: %u; col: %u\n", (unsigned int)m->row, (unsigned int)m->col);
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++)
      printf("%8.3f\t", m->data[i][j]);
    printf("\n");
  }
}

/*if a value is in matrix return 1 else 0*/
int ValInMatrix(matrix* m, double val)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      /*if((size_t)getMatrixValue(m, i, j) == (size_t)id) */
      if(FLOAT_EQ(m->data[i][j], val, 1*10e-8))
        return 1;
      else
        continue;
    }
  }
  return 0;
}

void MatrixSet(matrix *m, double val)
{
  size_t i, j;
  if(m->row == m->col){
    for(i = 0; i < m->row; i++){
      m->data[i][i] = val;
      for(j = i+1; j < m->col; j++){
        m->data[i][j]= m->data[j][i] = val;
      }
    }
  }
  else{
    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        m->data[i][j] = val;
      }
    }
  }
}

void MatrixInitRandomInt(matrix *m, int low, int high)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      m->data[i][j] = (double)randInt(low, high);
    }
  }
}

void MatrixInitRandomFloat(matrix *m, double low, double high)
{
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      m->data[i][j] = randDouble(low, high);
    }
  }
}

void MatrixCopy(matrix *msrc, matrix **mdst)
{
  size_t i, j;
  if((*mdst)->data == NULL){
    (*mdst)->row = msrc->row;
    (*mdst)->col = msrc->col;
    (*mdst)->data = xmalloc(sizeof(double*)*msrc->row);
    for(i = 0; i < msrc->row; i++){
      (*mdst)->data[i] = xmalloc(sizeof(double)*msrc->col);
      for(j = 0; j < msrc->col; j++)
        (*mdst)->data[i][j] = +0.f;
    }
  }
  else{
    if((*mdst)->row != msrc->row || (*mdst)->col != msrc->col){
      DelMatrix(mdst);
      initMatrix(mdst);
      (*mdst)->row = msrc->row;
      (*mdst)->col = msrc->col;

      (*mdst)->data = xmalloc(sizeof(double*)*msrc->row);
      for(i = 0; i < (*mdst)->row; i++){
        (*mdst)->data[i] = xmalloc(sizeof(double)*msrc->col);
        for(j = 0; j < msrc->col; j++){
          (*mdst)->data[i][j] = +0.f;
        }
      }
    }
  }

  for(i = 0; i < msrc->row; i++){
    for(j = 0; j < msrc->col; j++){
      (*mdst)->data[i][j] = msrc->data[i][j];
    }
  }
}

void setMatrixValue(matrix *m, size_t row, size_t col, double val)
{
  if(row < m->row && col < m->col){
    if(_isnan_(val) || _isinf_(val)){
      (*m).data[row][col] = MISSING;
    }
    else{
      (*m).data[row][col] = val;
    }
  }
  else{
    fprintf(stdout,"setMatrixValue Error: row to set: %u row max: %u; column to set %u; colum max: %u out of range.\n", (unsigned int)row, (unsigned int)m->row-1, (unsigned int)col, (unsigned int)m->col-1);
    fflush(stdout);
//     abort();
  }
}

double getMatrixValue(matrix *m, size_t row, size_t col)
{
  if(row < (*m).row && col < (*m).col){
    return (*m).data[row][col];
  }
  else{
    fprintf(stdout,"getMatrixValue Error: row to get: %u row max: %u; column to get %u; colum max: %u out of range.\n", (unsigned int)row, (unsigned int)m->row-1, (unsigned int)col, (unsigned int)m->col-1);
    fflush(stdout);
//     abort();
    return NAN;
  }
}

dvector *getMatrixRow(matrix *m, size_t row)
{
  if(row < (*m).row){
    dvector *v;
    size_t j;
    NewDVector(&v, m->col);
    for(j = 0; j < m->col; j++){
      v->data[j] = m->data[row][j];
    }

    /*
    initDVector(&v);
    memcpy(v->data, (*m).data[row], (*m).col);
    v->size = (*m).col;
    */
    return v;
  }
  else{
    fprintf(stdout,"getRow Error: row %u out of range.\n", (unsigned int)row);
    fflush(stdout);
    return NULL;
  }

}

dvector *getMatrixColumn(matrix *m, size_t col)
{
  if(col < (*m).col){
    dvector *v;
    size_t i;
    NewDVector(&v, m->row);
    for(i = 0; i < m->row; i++)
      v->data[i] = m->data[i][col];
    return v;
  }
  else{
    fprintf(stdout,"getColumn Error: column %u out of range.\n", (unsigned int)col);
    fflush(stdout);
    return NULL;
  }
}

void MatrixAppendRow(matrix** m, dvector *row)
{
  size_t i, j;
  size_t rowsize =  (*m)->row + 1;
  size_t colsize;


  if((*m)->col != 0){
    if(row->size > (*m)->col)
      colsize = row->size;
    else /*if (row->size <= (*m)->col)*/
      colsize = (*m)->col;
  }
  else{
    colsize = row->size;
  }

  /*adding a new row*/
  (*m)->data = xrealloc((*m)->data, sizeof(double*)*rowsize);

  if(colsize > (*m)->col){
    /*resize the column*/
    for(i = 0; i < (*m)->row; i++){
      (*m)->data[i] = xrealloc((*m)->data[i], sizeof(double*)*colsize);
      /*initialize the new column value*/
      for(j = (*m)->col; j < colsize; j++){
        (*m)->data[i][j] = +0.f;
      }
    }
    /*allocate the last row added*/
    (*m)->data[rowsize-1] = xmalloc(sizeof(double)*colsize);

    /*copy the row value to the new row matrix*/
    for(i = 0; i < row->size; i++){
      (*m)->data[rowsize-1][i] = row->data[i];
    }
  }
  else /*if(colsize <= (*m)->col)*/{
    /*allocate the last row added*/
    (*m)->data[rowsize-1] = xmalloc(sizeof(double)*colsize);
    for(i = 0; i < (*m)->col; i++){
      if(i < row->size){
        (*m)->data[rowsize -1][i] = row->data[i];
      }
      else{
        (*m)->data[rowsize -1][i] = +0.f;
      }
    }
  }

  if(row->size > (*m)->col)
    (*m)->col = row->size;

  (*m)->row += 1;

}

void MatrixAppendCol(matrix** m, dvector *col)
{
  size_t i, j;
  size_t lastcol;
  size_t colsize;
  size_t rowsize;

  if((*m)->col != 0){
    colsize =  (*m)->col + 1;
  }
  else{
    /* redefine anyway the row size because the column have 0 size */
    colsize = 1;
  }

  if((*m)->row != 0){
    if((*m)->row < col->size){
      rowsize = col->size;
    }
    else{
      rowsize = (*m)->row;
    }
  }
  else{
    rowsize = col->size;
  }

  if((*m)->row < rowsize){
    (*m)->data = xrealloc((*m)->data, sizeof(double*)*rowsize);
  }

  for(i = 0; i < rowsize; i++){
    if(i < (*m)->row)
      (*m)->data[i] = xrealloc((*m)->data[i], sizeof(double)*colsize);
    else
      (*m)->data[i] = xmalloc(sizeof(double)*colsize);
  }

  lastcol = (*m)->col;

  if(rowsize < (*m)->row){
    for(i = 0; i < (*m)->row; i++ ){
      if(i < rowsize)
        (*m)->data[i][lastcol] = col->data[i];
    else
      (*m)->data[i][lastcol] = +0.f;
    }
  }
  else{
    if(rowsize > (*m)->row){
      for(i = 0; i < rowsize; i++ ){
        (*m)->data[i][lastcol] = col->data[i];
      }

      /*Fill with 0 value the new rows except the last column */
      for(i = (*m)->row; i < rowsize; i++)
        for(j = 0; j < colsize-1; j++)
          (*m)->data[i][j] = +0.f;
    }
    else{
      for(i = 0; i < rowsize; i++){
        (*m)->data[i][lastcol] = col->data[i];
      }
    }
  }

  (*m)->col = colsize; /* updating the column matrix size */
  (*m)->row = rowsize; /* updating the row matrix size */
}


void MatrixAppendUIRow(matrix** m, uivector *row)
{
  size_t i, j;
  size_t rowsize =  (*m)->row + 1;
  size_t colsize;


  if((*m)->col != 0){
    if(row->size > (*m)->col)
      colsize = row->size;
    else /*if (row->size <= (*m)->col)*/
      colsize = (*m)->col;
  }
  else{
    colsize = row->size;
  }

  /*adding a new row*/
  (*m)->data = xrealloc((*m)->data, sizeof(double*)*rowsize);

  if(colsize > (*m)->col){
    /*resize the column*/
    for(i = 0; i < (*m)->row; i++){
      (*m)->data[i] = xrealloc((*m)->data[i], sizeof(double*)*colsize);
      /*initialize the new column value*/
      for(j = (*m)->col; j < colsize; j++){
        (*m)->data[i][j] = +0.f;
      }
    }
    /*allocate the last row added*/
    (*m)->data[rowsize-1] = xmalloc(sizeof(double)*colsize);

    /*copy the row value to the new row matrix*/
    for(i = 0; i < row->size; i++){
      (*m)->data[rowsize-1][i] = row->data[i];
    }
  }
  else /*if(colsize <= (*m)->col)*/{
    /*allocate the last row added*/
    (*m)->data[rowsize-1] = xmalloc(sizeof(double)*colsize);
    for(i = 0; i < (*m)->col; i++){
      if(i < row->size){
        (*m)->data[rowsize -1][i] = row->data[i];
      }
      else{
        (*m)->data[rowsize -1][i] = +0.f;
      }
    }
  }

  if(row->size > (*m)->col)
    (*m)->col = row->size;

  (*m)->row += 1;

}

void MatrixAppendUICol(matrix** m, uivector *col)
{
  size_t i, j;
  size_t lastcol;
  size_t colsize;
  size_t rowsize;

  if((*m)->col != 0){
    colsize =  (*m)->col + 1;
  }
  else{
    /* redefine anyway the row size because the column have 0 size */
    colsize = 1;
  }

  if((*m)->row != 0){
    if((*m)->row < col->size){
      rowsize = col->size;
    }
    else{
      rowsize = (*m)->row;
    }
  }
  else{
    rowsize = col->size;
  }

  if((*m)->row < rowsize){
    (*m)->data = xrealloc((*m)->data, sizeof(double*)*rowsize);
  }

  for(i = 0; i < rowsize; i++){
    if(i < (*m)->row)
      (*m)->data[i] = xrealloc((*m)->data[i], sizeof(double)*colsize);
    else
      (*m)->data[i] = xmalloc(sizeof(double)*colsize);
  }

  lastcol = (*m)->col;

  if(rowsize < (*m)->row){
    for(i = 0; i < (*m)->row; i++ ){
      if(i < rowsize)
        (*m)->data[i][lastcol] = col->data[i];
    else
      (*m)->data[i][lastcol] = +0.f;
    }
  }
  else{
    if(rowsize > (*m)->row){
      for(i = 0; i < rowsize; i++ ){
        (*m)->data[i][lastcol] = col->data[i];
      }

      /*Fill with 0 value the new rows except the last column */
      for(i = (*m)->row; i < rowsize; i++)
        for(j = 0; j < colsize-1; j++)
          (*m)->data[i][j] = +0.f;
    }
    else{
      for(i = 0; i < rowsize; i++){
        (*m)->data[i][lastcol] = col->data[i];
      }
    }
  }

  (*m)->col = colsize; /* updating the column matrix size */
  (*m)->row = rowsize; /* updating the row matrix size */
}

/*
 * p[i] =   Σ mx[i][j] * v[j]
 */
void MatrixDVectorDotProduct(matrix *mx, dvector *v, dvector *p)
{
  /* (mx*vect) where t is a column vector not transposed
     the size of the "vect" vector must be equal to the number of matrix row*/
  size_t i, j;
  double res;
  if(mx->col == v->size){
    for(i = 0; i < mx->row; i++){
      for(j = 0; j < mx->col; j++){
        if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1) ||
           FLOAT_EQ(v->data[j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = mx->data[i][j] * v->data[j];
          if(_isnan_(res) || _isinf_(res)){
            continue;
          }
          else{
            p->data[i] += res;
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixDVectorDotProduct Error while calculating product (X*v)!!\n The column vector size must be equal to the matrix column size.\n");
    fflush(stdout);
    abort();
  }
}

/*
 * p[j] =   Σ v[i] * mx[i][j]
 */
void DVectorMatrixDotProduct(matrix *mx, dvector *v, dvector *p)
{
  size_t i, j;
  double res;
  /* (vect'*mx) where vect is colunm vector transposed
     the size of the "vect" vector must be equal to the number of matrix column */
  if(mx->row == v->size){
    for(j = 0; j < mx->col; j++){
      for(i = 0; i < mx->row; i++){
        if(FLOAT_EQ(v->data[i], MISSING, 1e-1) ||
           FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = v->data[i] * mx->data[i][j];
          if(_isnan_(res) || _isinf_(res)){
            continue;
          }
          else{
            p->data[j] += res;
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"DVectorMatrixDotProduct Error while calculating product of a (v'*X)!!\n The transposed column vector size must be equal to the matrix row size.\n");
    fflush(stdout);
    abort();
  }
}

void DVectorTrasposedDVectorDotProduct(dvector *v1, dvector *v2, matrix *m)
{
  size_t i, j;

  if(m->row != v1->size && m->col != v2->size)
    ResizeMatrix(&m, v1->size, v2->size);

  for(i = 0; i < v1->size; i++){
    for(j = 0; j < v2->size; j++){
      if(FLOAT_EQ(v1->data[i], MISSING, 1e-1) ||
         FLOAT_EQ(v2->data[i], MISSING, 1e-1)){
        m->data[i][j] = MISSING;
      }
      else{
        m->data[i][j] = v1->data[i]*v2->data[j];
      }
    }
  }
}

/*
v/mx = (inverse (mx') * v')'
*/
void DVectorTransposedMatrixDivision(dvector *v, matrix *mx, dvector *r)
{
  if(mx->col == v->size){
    matrix *mx_t, *mx_inv;
    NewMatrix(&mx_t, mx->col, mx->row);
    MatrixTranspose(mx, mx_t);
    initMatrix(&mx_inv);
    MatrixInversion(mx_t, &mx_inv);
    DVectorResize(&r, mx_inv->col);
    MatrixDVectorDotProduct(mx_inv, v, r);
    DelMatrix(&mx_inv);
    DelMatrix(&mx_t);
  }
  else{
    fprintf(stderr, "Unable to compute vector / matrix. vector size %d != matrix column size %d\n", (int)v->size, (int)mx->row);
    abort();
  }
}
/*
 * R = M'M
 * opp
 *
 * The product of an m x n matrix A and an n x p matrix B is an m x p matrix C where
 *
 * c[i][j] = Sum a[i][k]*b[k][j]
 */

/*
 * m_t is the transposed matrix of m. The result is a square matrix named r.
 * void MatrixDotProduct(matrix *m_t, matrix *m, matrix **r)
{
  int i, j, k;
//  (*r).row =(*r).col = m.col; square matrix
  for(i = 0; i < (*r)->row; i++){
    for(j = 0; j < m_t->row; j++){
      for(k = 0; k < m_t->col; k++){
        (*r)->data[j][i] += m_t->data[j][k] * m->data[k][i];
      }
    }
  }
}
*/

void MatrixDotProduct(matrix *a, matrix *b, matrix *c)
{
  if(a->col == b->row){
    size_t i, j, k;
    double res;
    for(i = 0; i < a->row; i++){ /* m */
      for(j = 0; j < b->col; j++){ /* p */
        for(k = 0; k < a->col; k++){ /* n */
          if(FLOAT_EQ(a->data[i][k], MISSING, 1e-1) ||
             FLOAT_EQ(b->data[k][j], MISSING, 1e-1)){
            continue;
          }
          else{
            res = a->data[i][k] * b->data[k][j];
            if(_isnan_(res) || _isinf_(res)){
              c->data[i][j] +=  +0.f;
            }
            else{
              c->data[i][j] +=  res;
            }
          }
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixDotProduct Error!!\n The product of an m x n matrix A and an n x p matrix B is an m x p matrix C. col(a): %u != row(b) %u\n", (unsigned int)a->col, (unsigned int)b->row);
    fflush(stdout);
    abort();
  }
}


/*
 * Array product between two vector
 * X = t ⊗ p'
 *
 * X_ij = t_i*p_j
 *
 */
void RowColOuterProduct(dvector *a, dvector *b, matrix *m)
{
  size_t i, j;
  double res;
  for(i = 0; i < a->size; i++){
    for(j = 0; j < b->size; j++){
      if(FLOAT_EQ(a->data[i], MISSING, 1e-1) ||
         FLOAT_EQ(b->data[j], MISSING, 1e-1)){
        m->data[i][j] = MISSING;
      }
      else{
        res = a->data[i] * b->data[j];
        if(_isnan_(res) || _isinf_(res)){
          m->data[i][j] = +0.f;
        }
        else{
          m->data[i][j] = res;
        }
      }
    }
  }
}

void MatrixTranspose(matrix *m, matrix *r){
  size_t i, j;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++){
      r->data[j][i] = m->data[i][j];
    }
  }
}

/*Gauss-Jordan
 *
 * 1) If the first row have the first element 0, then excange this with an other row that have the first element != 0. If all the row have the first element 0 go to the step 3.
 * 2) For each row (*AT)i with the first element != 0, except the first row considered, multiply the first row for a coefficient c that must
 *
 * Spiegazione in italiano:
 * Per ogni riga Ai con primo elemento non nullo, eccetto la prima (i > 1), moltiplica la prima riga per un coefficiente scelto in maniera tale che
 * la somma tra la prima riga e Ai abbia il primo elemento nullo (quindi coefficiente = − Ai1 / A11). Sostituisci Ai con la somma appena ricavata.
 *
 *
 */

void MatrixInversion(matrix *m, matrix **m_inv)
{
  if(m->row == m->col){
    size_t i, j, k;
    double ratio, a;

    matrix *AI;
    NewMatrix(&AI, m->row, m->col*2);

    /*copy the m1 value to AI*/
    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        AI->data[i][j] = m->data[i][j];
      }
    }

    /*build the identity matrix*/
    for(i = 0; i < m->row; i++){
      for(j = m->col; j < 2*m->col; j++){
        if(i==(j-m->row)){
          AI->data[i][j] = 1.f;
        }
        else{
          AI->data[i][j] = +0.f;
        }
      }
    }


    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        if(i!=j){
          ratio = AI->data[j][i] / AI->data[i][i];
          for(k = 0; k < 2*m->col; k++){
            AI->data[j][k] -= AI->data[i][k]*ratio;
          }
        }
      }
    }

    for(i = 0; i < m->row; i++){
      a = AI->data[i][i];
      for(j = 0; j < 2*m->col; j++){
        AI->data[i][j] = AI->data[i][j]/a;
      }
    }

    if((*m_inv)->row != m->row || (*m_inv)->col != m->col){
      ResizeMatrix(m_inv, m->row, m->col);
    }

    for(i = 0; i < m->row; i++){
      for(j = 0; j < m->col; j++){
        (*m_inv)->data[i][j] = AI->data[i][m->col+j];
      }
    }

    DelMatrix(&AI);
  }
  else{
    fprintf(stdout,"Matrix Inversion Error!\n The matrix to invert must be squared!\n");
    fflush(stdout);
    abort();
  }
}

void GenIdentityMatrix(matrix **m)
{
  size_t i;
  if((*m)->row == (*m)->col){
    for(i=0; i < (*m)->row; i++){
      (*m)->data[i][i]=1;
    }
  }
}

void MeanCenteredMatrix(matrix *mx, matrix *mxc)
{
  size_t i, j, n;
  double average, res;
  for(j = 0; j < mx->col; j++){
    if(mx->row > 1){
      /*Calculate the average */
      average = +0.f;
      n = 0;
      for(i = 0; i < mx->row; i++ ){
        if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          average += mx->data[i][j];
          n++;
        }
      }

      average /= (double)n;

      for(i = 0; i < mx->row; i++){
        if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
          continue;
        }
        else{
          res = mx->data[i][j] - average;
          if(_isnan_(res) || _isinf_(res)){
            (*mxc).data[i][j] = +0.f;
          }
          else{
            (*mxc).data[i][j] = res;
          }
        }
      }
    }
    else{
      for(i = 0; i < mx->row; i++){
        (*mxc).data[i][j] = mx->data[i][j];
      }
    }
  }
}

/*
 * The equation for the Pearson product moment correlation coefficient, r, is:
 * RSQ = Sum (x_i -x_med)*(y_i - y_med) / sqrt( Sum (x_i - x_med)^2 Sum (y_i - y_med)^2 )
 */

void PearsonCorrelMatrix(matrix* mxsrc, matrix* mxdst)
{
  size_t i, j, k;
  double n, a, b, xres, yres;
  dvector *mean;
  ResizeMatrix(&mxdst, mxsrc->col, mxsrc->col);
  initDVector(&mean);
  MatrixColAverage(mxsrc, &mean);
  for(k = 0; k < mxsrc->col; k++){
    mxdst->data[k][k] = 1.f;
    for(j = k+1; j < mxsrc->col; j++){
      n = a = b = +0.f;
      for(i = 0; i < mxsrc->row; i++){
        xres = mxsrc->data[i][k] - mean->data[k];
        yres = mxsrc->data[i][j] - mean->data[j];
        n += xres * yres;
        a += square(xres);
        b += square(yres);
      }
      if((int)floor(a*b) == 0){
        mxdst->data[k][j] = mxdst->data[j][k] = +0.f;
      }
      else{
        mxdst->data[k][j] = mxdst->data[j][k] = square(n / sqrt(a*b));
      }
    }
  }
  DelDVector(&mean);
}

/*
 * The equation for the Spearman product moment correlation coefficient, r, is:
 * rho = 1 - (Sum (rank_x_i - rank_y_i)^2 / (n(n^2-1)))
 * where n are the observation
 */


void SpearmanCorrelMatrix(matrix* mxsrc, matrix* mxdst)
{
  size_t i, j, k, l;
  matrix *rankmx;
  dvector *vtosort;
  double n;
  ResizeMatrix(&mxdst, mxsrc->col, mxsrc->col);
  NewMatrix(&rankmx, mxsrc->row, 4);
  NewDVector(&vtosort, mxsrc->row);
  for(k = 0; k < mxsrc->col; k++){
    mxdst->data[k][k] = 1.f;
    for(i = 0; i < mxsrc->row; i++){
      rankmx->data[i][0] = vtosort->data[i] = mxsrc->data[i][k];
    }

    DVectorSort(vtosort);

    for(i = 0; i < vtosort->size; i++){
      for(j = 0; j < rankmx->row; j++){
        if(FLOAT_EQ(vtosort->data[i], rankmx->data[j][0], EPSILON)){
          rankmx->data[j][2] = i+1;
          break;
        }
        else{
          continue;
        }
      }
    }

    for(j = k+1; j < mxsrc->col; j++){
      for(i = 0; i < mxsrc->row; i++){
        rankmx->data[i][1] = vtosort->data[i] = mxsrc->data[i][j];
      }

      /* rank second column */
      DVectorSort(vtosort);
      for(i = 0; i < vtosort->size; i++){
        for(l = 0; l < rankmx->row; l++){
          if(FLOAT_EQ(vtosort->data[i], rankmx->data[l][1], EPSILON)){
            rankmx->data[l][3] = i+1;
            break;
          }
          else{
            continue;
          }
        }
      }

      /*calculate d^2*/
      n = +0.f;
      for(i = 0; i < rankmx->row; i++){
        n += square(rankmx->data[i][2] - rankmx->data[i][3]);
      }
      mxdst->data[k][j] = mxdst->data[j][k] = 1 - ((6*n) / (rankmx->row*((square(rankmx->row)-1))));
    }
  }
  DelMatrix(&rankmx);
  DelDVector(&vtosort);
}

void MatrixColAverage(matrix *mx, dvector **colaverage)
{
  size_t i, j, n;

  double average;
  for(j = 0; j < mx->col; j++){
    /*Calculate the average */
    average = +0.f;
    n = 0;
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += mx->data[i][j];
        n++;
      }
    }

    if(FLOAT_EQ(average, 0.f, 1e-6))
        average = 0.f;
    else
        average /= (double)n;
    DVectorAppend(colaverage, average);
  }
}

void MatrixRowAverage(matrix *mx, dvector **rowaverage)
{
  size_t i, j, n;

  double average;
  for(i = 0; i < mx->row; i++){
    /*Calculate the average */
    average = +0.f;
    n = 0;
    for(j = 0; j < mx->col; j++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += mx->data[i][j];
        n++;
      }
    }
    average /= (double)n;
    DVectorAppend(rowaverage, average);
  }
}

void MatrixColSDEV(matrix* mx, dvector** colsdev)
{
  size_t i, j, n;
  double var, average;
  for(j = 0; j < mx->col; j++){
    average=+0.f;
    n = 0;
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += mx->data[i][j];
        n++;
      }
    }

    /* average of the column j;*/
    average /= (double)n;

    var = +0.f;
    n = 0;
    for(i = 0; i< mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        var += square(mx->data[i][j] - average);
        n++;
      }
    }
    /* sample variance: is used whe the average of data is not known so you need to extimate the data average */
    var = var/(n-1);

    /* standard deviation calculation */
    DVectorAppend(colsdev, sqrt(var));
  }
}

void MatrixColRMS(matrix* mx, dvector** colrms)
{
  size_t i, j, n;
  double a;

  for(j = 0; j < mx->col; j++){
    a = +0.f;
    n = 0;
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        a += square(mx->data[i][j]);
        n++;
      }
    }

    /* average of the column j;*/
    a /= (double)n;
    /* standard deviation calculation */
    DVectorAppend(colrms, sqrt(a));
  }
}

void MatrixColVar(matrix* mx, dvector** colvar)
{
  size_t i, j, n;
  double var, average;
  for(j = 0; j < mx->col; j++){
    average = +0.f;
    n = 0;
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        average += mx->data[i][j];
        n++;
      }
    }

    /* average of the column j;*/
    average /= (double)n;

    var = +0.f;
    n = 0;
    for(i = 0; i<mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1))
        continue;
      else{
        double a = (mx->data[i][j] - average);
        var += a*a;
        n++;
      }
    }
    /* sample variance: is used whe the average of data is not known so you need to extimate the data average */
    var = var/(double)(n-1);

    /* standard deviation calculation */
    DVectorAppend(colvar, var);
  }
}

/* Calculate the matrix descriptive statistics:
 *  - Column Average
 *  - Column Median
 *  - Column Armonic Average
 *  - Column Variance Population
 *  - Column Variance Sample (Correcter Variance)
 *  - Column Standard Deviation
 *  - Column Standard Deviation Sample (Corrected Standard Deviation)
 *  - Column Max
 *  - Column Min
 *  - Coefficient of variation Population (CV)
 *  - Coefficient of variation Sample (CV)
 *  - N. zeros
 *  - N. missing values
 */
void MatrixColDescStat(matrix *mx, matrix **ds)
{
  int i, j, n;
  size_t n_zeros;
  size_t n_missing;
  double avg = 0.f;
  double median = 0.f;
  double armonic = 0.f;
  double var = 0.f;
  double min = 0.f, max = 0.f;
  dvector *v;

  /* n = mx->row if no MISSING value */
  ResizeMatrix(ds, mx->col, 13);
  NewDVector(&v, mx->row);

  for(j = 0; j < mx->col; j++){
    avg = 0.f;
    var = 0.f;
    median = 0.f;
    armonic = 0.f;

    /* get column min and max,
     * cacluate the average, the median and the armonic average
     */
    min = max = mx->data[0][j];
    n_zeros = 0;
    n_missing = 0;
    n = 0;
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], 0.f, 1e-6))
        n_zeros++;

      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
        n_missing++;
      }
      else{
        avg += mx->data[i][j];
        armonic += 1.f/mx->data[i][j];
        v->data[i] = mx->data[i][j];
        if(mx->data[i][j] > max)
          max = mx->data[i][j];

        if(mx->data[i][j] < min)
          min = mx->data[i][j];
        n++;
      }
    }
    avg /= (double)n;

    DVectorMedian(v, &median);
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
        continue;
      }
      else{
        var += square(mx->data[i][j] - avg);
      }
    }

    (*ds)->data[j][0] = avg;
    (*ds)->data[j][1] = median;
    (*ds)->data[j][2] = (double)(mx->row)/armonic;
    (*ds)->data[j][3] = var/(double)n;
    (*ds)->data[j][4] = var/(double)(n-1);
    (*ds)->data[j][5] = sqrt(var/(double)n);
    (*ds)->data[j][6] = sqrt(var/(double)(n-1));
    (*ds)->data[j][7] = (*ds)->data[j][5]/avg * 100;
    (*ds)->data[j][8] = (*ds)->data[j][6]/avg * 100;
    (*ds)->data[j][9] = min;
    (*ds)->data[j][10] = max;
    (*ds)->data[j][11] = n_zeros;
    (*ds)->data[j][12] = n_missing;
  }
  DelDVector(&v);
}

/* calculation of the covariance matrix */
void MatrixCovariance(matrix* mx, matrix** cm)
{
  size_t i, j, k;
  double sum;
  dvector *colaverage;

  ResizeMatrix(cm, mx->col, mx->col);

  initDVector(&colaverage);
  MatrixColAverage(mx, &colaverage);


  for(i = 0; i < mx->col; i++){
    for(j = 0; j < mx->col; j++){
      sum = +0.f;
      for(k =0; k < mx->row; k++){
        sum += (mx->data[k][i] - colaverage->data[i]) * (mx->data[k][j] - colaverage->data[j]);
      }
      (*cm)->data[i][j] = sum/(mx->row-1);
    }
  }

  DelDVector(&colaverage);
}

/* Transform a matrix into a logaritmic matrix */
void Matrix2LogMatrix(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++){
    for(j = 0; j < mx_in->col; j++){
      (*mx_out)->data[i][j] = log10(mx_in->data[i][j]+1);
    }
  }
}

/* Transform a matrix into a SQUARE matrix */
void Matrix2SquareMatrix(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++)
    for(j = 0; j < mx_in->col; j++)
      (*mx_out)->data[i][j] = square(mx_in->data[i][j]);
}

/* Transform a matrix into a SQRT matrix */
void Matrix2SQRTMatrix(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++)
    for(j = 0; j < mx_in->col; j++)
      (*mx_out)->data[i][j] = sqrt(mx_in->data[i][j]);
}

/* Transform a matrix into ABS matrix */
void Matrix2ABSMatrix(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++)
    for(j = 0; j < mx_in->col; j++)
      (*mx_out)->data[i][j] = fabs(mx_in->data[i][j]);
}


/* Develop an interaction factors matrix
 * Es. Use in DOE
 */
void Matrix2IntFactorsMatrix(matrix *mx_in, size_t factors, matrix **mx_out)
{
  size_t i, j, k, l, c;
  size_t nifc = 0;
  for(i = 1; i < mx_in->col; i++)
    nifc += i;
  ResizeMatrix(mx_out, mx_in->row, (size_t)nifc+(2*mx_in->col));
  for(i = 0; i < mx_in->row; i++){
    puts("#######");
    for(k = 0, c = 0; k < factors; k++){
      for(j = 0; j < mx_in->col; j++){
        double res = mx_in->data[i][j];
        for(l = 0; l < k; l++){
          res *= mx_in->data[i][j+l];
          printf("%d * %d \n", (int)j, (int)l);
        }
        (*mx_out)->data[i][c] = res;
        c++;
      }
    }
  }
}


/* Transform a matrix into a row centered scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixRowCenterScaling(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  double rowsum;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++){
    rowsum = 0.f;
    for(j = 0; j < mx_in->col; j++){
      rowsum += mx_in->data[i][j];
    }

    for(j = 0; j < mx_in->col; j++){
      (*mx_out)->data[i][j] = mx_in->data[i][j]/rowsum;
    }
  }
}

/* Transform a matrix into a SVN row scaled matrix
 * Es. Use in Spectroscopy
 */
void MatrixSVNScaling(matrix *mx_in, matrix **mx_out)
{
  size_t i, j;
  double rowaverage, rowstdev;
  ResizeMatrix(mx_out, mx_in->row, mx_in->col);
  for(i = 0; i < mx_in->row; i++){
    rowaverage = 0.f;
    for(j = 0; j < mx_in->col; j++){
      rowaverage += mx_in->data[i][j];
    }
    rowaverage /= (double)mx_in->col;

    rowstdev = 0.f;
    for(j = 0; j < mx_in->col; j++){
      rowstdev += square(mx_in->data[i][j]-rowaverage);
    }

    rowstdev /= (double)(mx_in->col-1);
    rowstdev = sqrt(rowstdev);

    for(j = 0; j < mx_in->col; j++){
      (*mx_out)->data[i][j] = (mx_in->data[i][j]-rowaverage)/rowstdev;
    }
  }
}


/*
 * ||X|| = the square root of the sum of the squares of all the elements in the matrix
 */
double Matrixnorm(matrix *mx)
{
  size_t i, j;
  double norm = +0.f, v;
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      v = mx->data[i][j];
      if(_isnan_(v) || !_isfinite_(v)){
        continue;
      }
      else{
        norm += square(v);
      }
    }
  }
  return sqrt(norm);
}

double Matrix1norm(matrix *mx)
{
  size_t i, j;
  double norm = +0.f;
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      norm += fabs(mx->data[i][j]);
    }
  }
  return norm;
}

double MatrixDeterminant(matrix *mx1)
{
  if(mx1->row == mx1->col){
    size_t i, j, k, l;
    double d = 0;
    matrix *mx2;

    if (mx1->row < 1){
      return MISSING;
    }
    else if(mx1->row == 1){
      d = mx1->data[0][0];
    }
    else if (mx1->row == 2){
      d = mx1->data[0][0] * mx1->data[1][1] - mx1->data[1][0] * mx1->data[0][1];
    }
    else{
      d = 0.;
      for(k = 0; k < mx1->row; k++){
        NewMatrix(&mx2, mx1->row-1, mx1->row-1);
        for(i = 1; i < mx1->row; i++){
          l = 0;
          for(j = 0; j < mx1->row; j++){
            if(j == k)
              continue;
            else{
              mx2->data[i-1][l] = mx1->data[i][j];
              l++;
            }
          }
        }
        d += pow(-1.0, 1.0+k+1.0) * mx1->data[0][k] * MatrixDeterminant(mx2);
        DelMatrix(&mx2);
      }
    }
    return d;
  }
  else{
    /* Error! */
    return MISSING;
  }
}

/*
double MatrixDeterminant(matrix *mx)
{
  if(mx->row == mx->col){
    return Determinant(mx->data, mx->row);
  }
  else{
    return MISSING;
  }
}
*/

void MatrixNorm(matrix *mx, matrix *nmx)
{
  if((*nmx).row == mx->row && mx->col == (*nmx).col){
    size_t i, j;
    double mod, res;

    mod = Matrixnorm(mx);

    for(i = 0; i < mx->row; i++){
      for(j = 0; j < mx->col; j++){
        res = mx->data[i][j]/mod;
        if(_isnan_(res) || _isinf_(res)){
          nmx->data[i][j] = +0.f;
        }
        else{
          nmx->data[i][j] = res;
        }
      }
    }
  }
  else{
    fprintf(stdout,"MatrixNorm Error!\n");
    fflush(stdout);
    abort();
  }
}

void MatrixColumnMinMax(matrix* mx, size_t col, double* min, double* max)
{
  if(mx->row > 0 && col < mx->col ){
    size_t i;
    double a;
    (*min) = (*max) = mx->data[0][col];
    for(i = 1; i < mx->row; i++){
      a = mx->data[i][col];
      if(FLOAT_EQ(a, MISSING, 1e-1)){
        continue;
      }
      else{
        if(a < (*min)){
          (*min) = a;
        }

        if(a > (*max)){
          (*max) = a;
        }
      }
    }
  }
  else{
    fprintf(stdout,"Get Column Max Min Error!\n");
    fflush(stdout);
    (*min) = (*max) = MISSING;
  }
}

void MatrixSort(matrix* mx, size_t col_n)
{
  size_t i, j, k;
  double temp;
  for(i = 0; i < mx->row; i++){
    for(j = i+1; j < mx->row; j++){
      if(mx->data[i][col_n] > mx->data[j][col_n]){
        for(k = 0; k < mx->col; k++){
          temp = mx->data[i][k];
          mx->data[i][k] = mx->data[j][k];
          mx->data[j][k] = temp;
        }
      }
      else{
        continue;
      }
    }
  }
}

void MatrixReverseSort(matrix* mx, size_t col_n)
{
  size_t i, j, k;
  double temp;
  for(i = 0; i < mx->row; i++){
    for(j = i+1; j < mx->row; j++){
      if(mx->data[i][col_n] < mx->data[j][col_n]){
        for(k = 0; k < mx->col; k++){
          temp = mx->data[i][k];
          mx->data[i][k] = mx->data[j][k];
          mx->data[j][k] = temp;
        }
      }
      else{
        continue;
      }
    }
  }
}

void MatrixGetMaxValueIndex(matrix* mx, size_t* row, size_t* col)
{
  size_t i, j;
  double tmp_value, best_value;

  if(col != NULL)
    (*col) = 0;

  if(row != NULL)
    (*row) = 0;

  best_value = mx->data[0][0];

  for(j = 0; j < mx->col; j++){
    for(i = 1; i < mx->row; i++){
      tmp_value = mx->data[i][j];
      if(tmp_value > best_value || FLOAT_EQ(tmp_value, best_value, EPSILON)){
        best_value = tmp_value;
        if(col != NULL)
          (*col) = j;

        if(row != NULL)
          (*row) = i;
      }
    }
  }
}


void MatrixGetMinValueIndex(matrix* mx, size_t* row, size_t* col)
{
  size_t i, j;
  double tmp_value, best_value;

  if(col != NULL)
    (*col) = 0;

  if(row != NULL)
    (*row) = 0;

  best_value = mx->data[0][0];

  for(j = 0; j < mx->col; j++){
    for(i = 1; i < mx->row; i++){
      tmp_value = mx->data[i][j];
      if(tmp_value < best_value || FLOAT_EQ(tmp_value, best_value, EPSILON)){
        best_value = tmp_value;
        if(col != NULL)
          (*col) = j;

        if(row != NULL)
          (*row) = i;
      }
    }
  }
}


int cmpfunc(const void *a, const void *b )
{
  const double *a_ = *(const double **)a;
  const double *b_ = *(const double **)b;
  return b_[0] - a_[0];
}

void conv2matrix(int m, int n, double* a, int lda, matrix **mx)
{
  size_t i, j;
  ResizeMatrix(mx, m, n);
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++)
      (*mx)->data[i][j] = a[i+j*lda];
  }
}

void QRMatrixVectNorm(matrix *mx, size_t col, dvector *nv)
{
  size_t i;
  double s = +0.f;

  if(nv->size != mx->row)
    DVectorResize(&nv, mx->row);

  for(i = 0; i < mx->row; i++){
    s += square(mx->data[i][col]);
  }
  s = sqrt(s);

  for(i = 0; i < mx->row; i++){
    nv->data[i] = mx->data[i][col]/s;
  }
}

void QRDecomposition(matrix *mx, matrix **Q, matrix **R)
{
  size_t i, j, k;
  double D, p;
  dvector *d, *v;
  matrix *a, *a1, *P, *PQ;
  initMatrix(&a);
  MatrixCopy(mx, &a);
  NewMatrix(&a1, a->row, a->col);

  NewDVector(&d, a->row);
  NewDVector(&v, a->row);


  NewMatrix(&P, v->size, v->size);
  NewMatrix(&PQ, a->row, a->row);

  for(k = 0; k < mx->col && k < mx->row - 1; k++){
    QRMatrixVectNorm(a, k, d);

    D = +0.f;
    for(i = k; i < d->size; i++){
      D += d->data[i]*d->data[i];
    }
    D = sqrt(D);

    if(d->data[k] > 0)
      D = -D;

    for(i = 0; i < k; i++)
      v->data[i] = +0.f;

    v->data[k] = sqrt((1/2.)* (1 - (d->data[k]/D)));
    p = -1*D*v->data[k];

    for(i = k+1; i < v->size; i++){
      v->data[i] = d->data[i]/(2*p);
    }

    for(i = 0; i < v->size; i++){
      for(j = 0; j < v->size; j++){
        P->data[i][j] = -2 * v->data[i] * v->data[j];
      }
    }

    for(i = 0; i < v->size; i++){
      P->data[i][i] += 1;
    }

    MatrixSet(a1, +0.f);
    MatrixDotProduct(P, a, a1);

    if(k == 0){
      MatrixCopy(P, Q);
    }
    else{
      MatrixDotProduct(P, (*Q), PQ);
      MatrixCopy(PQ, Q);
      MatrixSet(PQ, +0.f);
    }

    MatrixCopy(a1, &a);
    MatrixSet(P, +0.f);
  }

  for(i = 0; i < (*Q)->row; i++){
    if(FLOAT_EQ((*Q)->data[i][i], 0, 1e-6))
      (*Q)->data[i][i] = 1.f;
  }

  ResizeMatrix(R, mx->row, mx->col);
  MatrixDotProduct((*Q), mx, (*R));

  MatrixCopy((*Q), &a1);

  MatrixTranspose(a1, (*Q));

  DelMatrix(&P);
  DelMatrix(&PQ);
  DelMatrix(&a1);
  DelMatrix(&a);
  DelDVector(&d);
  DelDVector(&v);
}

void LUDecomposition(matrix *mx, matrix **L, matrix **U)
{
  size_t i, j, k;
  double s;
  matrix *a;
  initMatrix(&a);
  MatrixCopy(mx, &a);

  ResizeMatrix(L, a->row, a->col);
  ResizeMatrix(U, a->row, a->col);
  for(i = 0; i < a->row; i++)
    (*U)->data[i][i] = 1.f;

  for(i = 0; i < a->row; i++){
    for(j = i; j < a->row; j++){
      s = +0.f;
      for(k = 0; k < a->row; k++){
        s += (*L)->data[j][k] * (*U)->data[k][i];
      }
      (*L)->data[j][i] = a->data[i][j] - s;
    }

    for(j = i + 1; j < a->row; j++){
      s = +0.f;
      for(k = 0; k < a->row; k++){
        s += (*L)->data[i][k] * (*U)->data[k][j];
      }

      (*U)->data[i][j] = (a->data[i][j] - s)/(*L)->data[i][i];;
    }
  }

  DelMatrix(&a);
}

/*
  Householder vector vector product
  function p = Housmvp(u, x)
  % Producto p = H*x, donde H es la transformación de Householder
  %    definida por u.
  u = u(:);
  x = x(:);
  v = u/norm(u);
  p = x - 2*v*(v.'*x);
*/

void HouseholderVectorVectorProduct(dvector *u, dvector *x, dvector *p)
{
  size_t i;
  double u_norm;
  dvector *v, *c;
  matrix *vx;
  NewDVector(&v, u->size);

  u_norm = DvectorModule(u);

  for(i = 0; i < u->size; i++){
    v->data[i] = u->data[i] / u_norm;
  }

  /*vx = (v.'*x)*/
  NewMatrix(&vx, v->size, x->size);
  DVectorTrasposedDVectorDotProduct(v, x, vx);

  /*c = 2*v*(v.'*x)*/
  NewDVector(&c, u->size);
  DVectorMatrixDotProduct(vx, v, c);

  if(p->size != x->size)
    DVectorResize(&p, x->size);

  for(i = 0; i < x->size; i++){
    p->data[i] = x->data[i] - (2*c->data[i]);
  }

  DelDVector(&c);
  DelMatrix(&vx);
  DelDVector(&v);
}

/*P = H*A */
void HouseholderVectMatrixProduct(dvector *h, matrix *A, matrix *P)
{
  size_t i, j;
  double mod_h, vaj;
  dvector *v, *aj;
  initDVector(&v);
  DVectorCopy(h, &v);
  mod_h = DvectorModule(h);

  for(i = 0; i < v->size; i++){
    v->data[i] /= mod_h;
  }

  if(P->row != A->row && P->col != A->col)
    ResizeMatrix(&P, A->row, A->col);

  for(j = 0; j < P->col; j++){
    aj = getMatrixColumn(A, j);
    vaj = DVectorDVectorDotProd(v, aj);

    for(i = 0; i < P->row; i++){
      P->data[i][j] = aj->data[i] - (2*v->data[i] * vaj);
    }
    DelDVector(&aj);
  }
}

/*HOUSEHOLDER REFLECTION UNIT VECTOR u FROM THE VECTOR x*/
void HouseReflectorVect(dvector *x, dvector *u)
{
  size_t i;
  double max, tmp, su;
  max = fabs(x->data[0]);
  for(i = 1; i < x->size; i++){
    tmp = fabs(x->data[i]);
    if(max > tmp){
      continue;
    }
    else{
      max = tmp;
    }
  }

  if(u->size != x->size)
    DVectorResize(&u, x->size);

  for(i = 0; i < u->size; i++)
    u->data[i] = x->data[i] / max;

  if(FLOAT_EQ(u->data[0], 0, EPSILON))
    su = 1;
  else{
    if(u->data[0] > 0)
      su = 1;
    else
      su = -1;
  }

  /*calculating norm of vector u*/
  tmp = +0.f;
  for(i = 0; i < u->size; i++)
    tmp += u->data[i]*u->data[i];
  tmp = sqrt(tmp);

  u->data[0] = u->data[0] + su*tmp;

  /*recalculate the new norm*/
  tmp = +0.f;
  for(i = 0; i < u->size; i++)
    tmp += u->data[i]*u->data[i];
  tmp = sqrt(tmp);

  for(i = 0; i < u->size; i++)
    u->data[i] /= tmp;
}

/*build the householder matrix according the relation:
 * I - 2*u u_T
 */
void HouseholderMatrix(dvector *v, matrix *h)
{
  size_t i, j;
  dvector *u;
  initDVector(&u);
  HouseReflectorVect(v, u);

  if(h->row != u->size){
    ResizeMatrix(&h, u->size, u->size);
    for(i = 0; i < h->row; i++){
      h->data[i][i] = 1.f;
    }
  }
  else{
    for(i = 0; i < h->row; i++){
      h->data[i][i] = 1.f;
      for(j = i+1; j < h->col; j++){
        h->data[i][j] = h->data[j][i] = +0.f;
      }
    }
  }

  for(i = 0; i < h->row; i++){
    for(j = 0; j < h->col; j++){
      h->data[i][j] -= 2*u->data[i]*u->data[j];
    }
  }

  DelDVector(&u);
}


void HouseholderReduction(matrix *mx){
  size_t i, j, k, ii;
  double a, d, w, f;
  dvector *col, *v;
  initDVector(&col);
  initDVector(&v);
  for(k = 0; k < mx->col - 1; k++){
    DVectorResize(&col, mx->row-k);
    DVectorResize(&v, mx->row-k);
    a = 0.f;
    ii = 0;
    for(i = k; i < mx->row; i++){
      a += square(mx->data[i][k]);
      col->data[ii] = mx->data[i][k];
      ii++;
    }
    a = sqrt(a);
    /*
    puts("Column");
    PrintDVector(col);
    */
    if(mx->data[k][k] > 0)
      d = -a;
    else
      d = a;

    w = mx->data[k][k] - d;
    f = sqrt(-2.0*w*d);

    /*printf("f1: %f\n", f);*/

    ii = 0;
    for(i = k; i < mx->row; i++){
      if(i == k){
        mx->data[i][k] = d;
        v->data[ii] = w/f;
      }
      else{
        mx->data[i][k] = 0.f;
        v->data[ii] = col->data[ii]/f;
      }
      ii++;
    }
    /*
    puts("v");
    for(i = 0; i < v->size; i++)
      printf("%f\n", v->data[i]);

    puts("Ha");
    for(i = 0; i < mx->row; i++)
      printf("%f\n", mx->data[i][k]);
    */
    for(j = k+1; j < mx->col; j++){
      ii = 0;
      for(i = k; i < mx->row; i++){
        col->data[ii] = mx->data[i][j];
        ii++;
      }

      f = 0.f;
      for(i = 0; i < v->size; i++){
        f += v->data[i] * col->data[i];
      }
      f *= 2;

      ii = 0;
      for(i = k; i < mx->row; i++){
        mx->data[i][j] = col->data[ii] - f*v->data[ii];
        ii++;
      }
      /*
      printf("fl %f\n", f);

      puts("Hal");
      for(i = 0; i < mx->row; i++)
        printf("%f\n", mx->data[i][j]);
     */
    }
    /*
    PrintMatrix(mx);
    sleep(2);*/
  }
  DelDVector(&v);
  DelDVector(&col);
}

void Cholesky(matrix *mx){
  size_t i, j, k, n = mx->row;
  double s;
  dvector *A, *L;
  NewDVector(&A, n*n);

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++)
      A->data[i * n + j] = mx->data[i][j];
  }

  NewDVector(&L, n*n);

  for(i = 0; i < n; i++)
    for(j = 0; j < (i+1); j++){
      s = 0;
      for(k = 0; k < j; k++)
        s += L->data[i * n + k] * L->data[j * n + k];
      L->data[i * n + j] = (i == j) ? sqrt(A->data[i * n + i] - s) : (1.0 / L->data[j * n + j] * (A->data[i * n + j] - s));
    }


  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++)
      mx->data[i][j] = L->data[i * n + j];
  }

  DelDVector(&A);
  DelDVector(&L);
}
