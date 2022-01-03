/* tensor.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tensor.h"
#include "matrix.h"
#include "numeric.h"
#include "memwrapper.h"

/*  Giuseppe Marco Randazzo
 *  Sun Dec 25 20:54:11 CET 2011
 *
 * Multidimensional Matrices
 */

void initTensor(tensor **t)
{
  (*t) = xmalloc(sizeof(tensor));
  (*t)->order = 0;
  (*t)->m = NULL;
}

void NewTensor(tensor** t, size_t order_)
{
  size_t i;
  (*t) = xmalloc(sizeof(tensor));
  (*t)->order = order_;
  (*t)->m = xmalloc(sizeof(matrix*)*order_);

  for(i = 0; i < order_; i++){
    (*t)->m[i] = NULL;
  }
}

void NewTensorMatrix(tensor **t, size_t order, size_t row, size_t col)
{
  if((*t)->order != 0){
    if(order < (*t)->order){
      NewMatrix(&(*t)->m[order], row, col);
    }
    else{
      fprintf(stderr, "Error! Order not in range order. order: %u > ordersize: %u\n", (unsigned int)order, (unsigned int)(*t)->order);
    }
  }
  else{
    fprintf(stderr, "Error! Order Tensor not Defined! Please create Tensor with NewTensor(tensor **t, size_t order);\n");
    abort();
  }
}

void AddTensorMatrix(tensor **t, size_t row, size_t col)
{
   if((*t)->order > 0){
    /*if((*t)->m[(*t)->order-1]->row != row){
      fprintf(stderr, "Error while appending matrix to tensor! Object size differ: matrix: %u;  tensor: %u", (unsigned int)row, (unsigned int)(*t)->m[(*t)->order-1]->row);
      abort();
    }
    else{*/
      (*t)->order += 1;
      (*t)->m = xrealloc((*t)->m, sizeof(matrix*)*(*t)->order);
      NewMatrix(&(*t)->m[(*t)->order-1], row, col);
/*     }*/
  }
  else{
    (*t)->order = 1;
    (*t)->m = xmalloc(sizeof(matrix*)*1);
    NewMatrix(&(*t)->m[0], row, col);
  }
}

void DelTensor(tensor** t)
{
  size_t i;
  for(i = 0; i < (*t)->order; i++)
    DelMatrix(&(*t)->m[i]);
  free((*t)->m);
  free((*t));
}


void PrintTensor(tensor *t)
{
  size_t i, j, k;
  printf("Tensor - order: %u\n", (unsigned int)t->order);
  for(i = 0; i < t->order; i++){
    printf("Tensor No: %u of row: %u; col: %u\n", (unsigned int)i+1, (unsigned int)t->m[i]->row, (unsigned int)t->m[i]->col);
    for(j = 0; j < t->m[i]->row; j++){
      for(k = 0; k < t->m[i]->col; k++)
        printf("%8.3f\t", getMatrixValue(t->m[i], j, k));
      printf("\n");
    }
  }
}

void setTensorValue(tensor *t, size_t i, size_t j, size_t k, double val)
{
  if(i < (*t).order){
    if(j < (*t).m[i]->row){
      if(k < (*t).m[i]->col){
        setMatrixValue((*t).m[i], j, k, val);
      }
      else{
        fprintf(stderr, "setTensorValue Error! Wrong colum index.\n");
        fflush(stderr);
        abort();
      }
    }
    else{
      fprintf(stderr, "setTensorValue Error! Wrong row index.\n");
      fflush(stderr);
      abort();
    }
  }
  else{
    fprintf(stderr, "setTensorValue Error! Wrong order.\n");
    fflush(stderr);
    abort();
  }
}

double getTensorValue(tensor *t, size_t i, size_t j, size_t k){
    if(i < (*t).order){
    if(j < (*t).m[i]->row){
      if(k < (*t).m[i]->col){
        return getMatrixValue((*t).m[i], j, k);
      }
      else{
        fprintf(stderr, "getTensorValue Error! Wrong colum index.\n");
        fflush(stderr);
        return NAN;
      }
    }
    else{
      fprintf(stderr, "getTensorValue Error! Wrong row index.\n");
      fflush(stderr);
      return NAN;
    }
  }
  else{
    fprintf(stderr, "getTensorValue Error! Wrong order.\n");
    fflush(stderr);
    return NAN;
  }
}

void TensorAppendMatrix(tensor **tdst, matrix *msrc)
{
  if((*tdst)->order > 0){
    if((*tdst)->m[(*tdst)->order-1]->row != msrc->row){
      fprintf(stderr, "Error while appending matrix to tensor! Object size differ: matrix: %u;  tensor: %u", (unsigned int)msrc->row, (unsigned int)(*tdst)->m[(*tdst)->order-1]->row);
      abort();
    }
    else{
      (*tdst)->order += 1;
      (*tdst)->m = xrealloc((*tdst)->m, sizeof(matrix*)*(*tdst)->order);
      NewMatrix(&(*tdst)->m[(*tdst)->order-1], msrc->row, msrc->col);
      MatrixCopy(msrc, &(*tdst)->m[(*tdst)->order-1]);
    }
  }
  else{
    (*tdst)->order = 1;
    (*tdst)->m = xmalloc(sizeof(matrix*)*1);
    NewMatrix(&(*tdst)->m[0], msrc->row, msrc->col);
    MatrixCopy(msrc, &(*tdst)->m[0]);
  }
}

void TensorAppendMatrixAt(tensor **tdst, size_t order, matrix *msrc)
{
  if(order < (*tdst)->order){
    fprintf(stderr, "Module not developed. Work in progress...\n");
    fflush(stderr);
    abort();
    /*tensor *tmp;
    initTensor(&tmp);


    DelTensor(&tmp);*/
  }
  else{
    TensorAppendMatrix(tdst, msrc);
  }
}

/*
 * t is the tensor
 * n is the order number
 * column is the colum vector to append
 */
void TensorAppendColumn(tensor **t, size_t n, dvector* column)
{
  if(n < (*t)->order){
    MatrixAppendCol(&(*t)->m[n], column);
    /*
    if(column->size == (*t)->m[n]->row){
      MatrixAppendCol(&(*t)->m[n], column);
    }
    else{
      fprintf(stderr, "Error! The objects number differ %u != %u\n", (unsigned int)column->size, (unsigned int)(*t)->m[n]->row);
      fflush(stderr);
      abort();
    }
    */
  }
  else{
    if(n > (*t)->order){
      fprintf(stderr, "Error! Order number too high. %u > %u\n", (unsigned int)n, (unsigned int)(*t)->order );
    }
    else{
      fprintf(stderr, "Error! Order 0!\n");
    }
    fflush(stderr);
    abort();
  }
}

/*
 * t is the tensor
 * n is the order number
 * column is the colum vector to append
 */
void TensorAppendRow(tensor **t, size_t n, dvector* row)
{
  if(n < (*t)->order){
    if(row->size != (*t)->m[n]->col){
      MatrixAppendRow(&(*t)->m[n], row);
    }
    else{
      fprintf(stderr, "Error! The column number differ %u != %u\n", (unsigned int)row->size, (unsigned int)(*t)->m[n]->col);
      fflush(stderr);
      abort();
    }
  }
  else{
    fprintf(stderr, "Error! Order number too high. %u > %u\n", (unsigned int)n, (unsigned int)(*t)->order);
    fflush(stderr);
    abort();
  }
}

void TensorSet(tensor* t, double val)
{
  size_t i;
  for(i = 0; i < t->order; i++)
    MatrixSet(t->m[i], val);
}

void TensorCopy(tensor* asrc, tensor** adst)
{
  size_t i, j, k;
  if((*adst)->m == NULL){
    (*adst)->order = asrc->order;

    (*adst)->m = xmalloc(sizeof(matrix*)*asrc->order);

    for(k = 0; k < asrc->order; k++){
      NewMatrix(&((*adst)->m[k]), asrc->m[k]->row, asrc->m[k]->col);
    }
  }
  else{
    if(asrc->order != (*adst)->order){
      /* resize  the order */
      (*adst)->m = xrealloc((*adst)->m, sizeof(tensor*)*asrc->order);
    }

    /*chek and resize the matrix for each order if is necessary */
    for(k = 0; k < asrc->order; k++){
      if(asrc->m[k]->row != (*adst)->m[k]->row || asrc->m[k]->col != (*adst)->m[k]->col){

        (*adst)->m[k]->row = asrc->m[k]->row;
        (*adst)->m[k]->col = asrc->m[k]->col;

        (*adst)->m[k]->data = xrealloc((*adst)->m[k]->data, sizeof(double*)*asrc->m[k]->row);

        for(i = 0; i < asrc->m[k]->row; i++){
          (*adst)->m[k]->data[i] = xrealloc((*adst)->m[k]->data[i], sizeof(double)*asrc->m[k]->col);
        }
      }
    }

  }

  /*copy the data...*/
  for(k = 0; k < asrc->order; k++){
    for(i = 0; i < asrc->m[k]->row; i++){
      for(j = 0; j < asrc->m[k]->col; j++){
        setTensorValue((*adst), k, i, j, getTensorValue(asrc, k, i, j));
      }
    }
  }

}


/*
 * t is the input data tensor
 * tc is the output tensor meancentered. This must be created with NewTensor()
 */
void MeanCenteredTensor(tensor *t, tensor *tc)
{
  size_t i;
  for(i = 0; i < t->order; i++){
    MeanCenteredMatrix(t->m[i], tc->m[i]);
  }
}

void TensorColAverage(tensor *t, matrix **colaverage)
{
  size_t i;
  dvector *v;
  for(i = 0; i < t->order; i++){
    initDVector(&v);
    MatrixColAverage(t->m[i], &v);
    MatrixAppendCol(colaverage, v);
    DelDVector(&v);
  }
}

void TensorColSDEV(tensor *t, matrix **colsdev)
{
  size_t i;
  dvector *v;
  for(i = 0; i < t->order; i++){
    initDVector(&v);
    MatrixColSDEV(t->m[i], &v);
    MatrixAppendCol(colsdev, v);
    DelDVector(&v);
  }
}

void TensorTranspose(tensor* t1, tensor* t2)
{
  size_t i;
  for(i = 0; i < t1->order; i++){
    MatrixTranspose(t1->m[i], t2->m[i]);
  }
}


/*
 * k = order
 * i = row
 * j = column
 *
 * E'(j,i,k) = E(i,j,k)
 *
 * P(k,j) = Sum_i[ E'(j,i,k)*t(i) ]
 *
 */
void TransposedTensorDVectorProduct(tensor *t, dvector *v, matrix *p)
{
  size_t i, j, k;
  double res;
  for(k = 0; k < t->order; k++){
    for(i = 0; i < t->m[k]->row; i++){
      for(j = 0; j < t->m[k]->col; j++){
        res = getTensorValue(t, k, i, j)*getDVectorValue(v, j);
        if(_isnan_(res) || _isinf_(res)){
          setMatrixValue(p, k, i, getMatrixValue(p, k, i) + 0);
        }
        else{
          setMatrixValue(p, k, i, getMatrixValue(p, k, i) + res);
        }
      }
    }
  }
}


/* the output "m" is size of:
 *   colum = t->order;
 *   row = v->size
 */
void DvectorTensorDotProduct(tensor* t, dvector* v, matrix* m)
{
  /*We do not need tests because getMatrixValue and setMatrixValue and getTensorValue handle some errors*/
  size_t i, j, k;
  double res;
  /* m = v't
  *
  * k = order of tensor
  * j = column size
  * i = row size
  *
  * m[j][k] =   Σ v[i] * t[k][i][j]
  *
  */
  for(k = 0; k < t->order; k++){
    for(i = 0; i < t->m[k]->row; i++){
      for(j = 0; j < t->m[k]->col; j++){
        if(v->size == t->m[k]->row && t->m[k]->col == m->row && t->order == m->col){
          res = getDVectorValue(v, i) * getTensorValue(t, k, i, j);
          if(_isnan_(res) || _isinf_(res)){
            setMatrixValue(m, j, k, (getMatrixValue(m, j, k) + 0));
          }
          else{
            setMatrixValue(m, j, k, (getMatrixValue(m, j, k) + res));
          }
        }
        else{
          fprintf(stderr, "Error while computing DvectorTensorDotProduct.\n");
          fflush(stderr);
          abort();
        }
      }
    }
  }
}

void TensorMatrixDotProduct(tensor *t, matrix *m, dvector *v)
{
  size_t i, j, k;
  double res;
  for(k = 0; k < t->order; k++){
    if(m->col == t->order && m->row == t->m[k]->col){
      for(i = 0; i < t->m[k]->row; i++){
        for(j = 0; j < t->m[k]->col; j++){
          res = (getTensorValue(t, k, i, j)*getMatrixValue(m, j, k));
          if(_isnan_(res) || _isinf_(res)){
            setDVectorValue(v, i, (getDVectorValue(v, i) + 0));
          }
          else{
            setDVectorValue(v, i, (getDVectorValue(v, i) + res));
          }
        }
      }
    }
    else{
      fprintf(stderr, "Error while computing TensorMatrixDotProduct.");
      fflush(stderr);
      abort();
    }
  }
}


/*
 * i = row
 * j = column
 * k = order
 *
 * t(i) = Sum_j[ Sum_k[ E(i,j,k) * P(k,j)] ]
 *
 */
void TensorMatrixDotProduct2(tensor *t, matrix *m, dvector *v)
{
  size_t i, j, k;
  for(i = 0; i < v->size; i++){
    for(k = 0; k < t->order; k++){
      for(j = 0; j < t->m[k]->col; j++){
        setDVectorValue(v, i, getDVectorValue(v, i) + (getTensorValue(t, k, i, j)*getMatrixValue(m, k, j)));
      }
    }
  }
}

void KronekerProductVectorMatrix(dvector* v, matrix* m, tensor* t)
{
  size_t i, j, k;
  double res;
  for(i = 0; i < v->size; i++){
    for(j = 0; j < m->row; j++){
      for(k = 0; k < m->col; k++){
        res = getDVectorValue(v, i)*getMatrixValue(m, j, k);
        if(_isnan_(res) || _isinf_(res)){
          setTensorValue(t, k, i, j, 0.f);
        }
        else{
          setTensorValue(t, k, i, j, res);
        }
      }
    }
  }
}
