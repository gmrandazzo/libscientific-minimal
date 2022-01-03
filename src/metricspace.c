/* metricspace.c
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

#include "memwrapper.h"
#include "metricspace.h"
#include "numeric.h"
#include <math.h>

void EuclideanDistance(matrix* m1, matrix* m2, matrix **distances)
{
  if(m1->col == m2->col){

    // SINGLE THREAD IMPLEMENTATION
    size_t i, j, k;
    double dist;
    // each column is a distance that correspond to m1->row
    ResizeMatrix(distances, m2->row, m1->row);


    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        // is the same of for(j = 0; j < m2->col; j++){
        for(j = 0; j < m1->col; j++){
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = sqrt(dist);
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void SquaredEuclideanDistance(matrix *m1, matrix *m2, matrix **distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;

    ResizeMatrix(distances, m2->row, m1->row);

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){
          dist += square(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }

  }
  else{
    fprintf(stderr, "Unable to compute Euclidean Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void ManhattanDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double dist;

    ResizeMatrix(distances, m2->row, m1->row);

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        dist = 0.f;
        for(j = 0; j < m1->col; j++){
          dist += fabs(m1->data[i][j] - m2->data[k][j]);
        }
        (*distances)->data[k][i] = dist;
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Manhattan Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void CosineDistance(matrix* m1, matrix* m2, matrix** distances)
{
  if(m1->col == m2->col){
    size_t i, j, k;
    double n, d_a, d_b;

    ResizeMatrix(distances, m2->row, m1->row);

    for(i = 0; i < m1->row; i++){
      for(k = 0; k < m2->row; k++){
        n = 0.f; d_a = 0.f; d_b = 0.f;
        for(j = 0; j < m1->col; j++){
          n += m1->data[i][j] * m2->data[k][j];
          d_a += square(m1->data[i][j]);
          d_b += square(m2->data[k][j]);
        }
        (*distances)->data[k][i] = n/(sqrt(d_a)*sqrt(d_b));
      }
    }
  }
  else{
    fprintf(stderr, "Unable to compute Cosine Distance. The number of variables differ\n");
    fflush(stderr);
    abort();
  }
}

void CovarianceDistanceMap(matrix* mi, matrix **mo)
{
  size_t i, j;
  dvector *x;
  dvector *p;
  dvector *colavg;
  matrix *covmx;
  matrix *inv_cov;
  
  initDVector(&colavg);
  MatrixColAverage(mi, &colavg);
  
  NewMatrix(&covmx, mi->col, mi->col);
  MatrixCovariance(mi, &covmx);

  
  NewMatrix(&inv_cov, mi->col, mi->col);
  MatrixInversion(covmx, &inv_cov);
  DelMatrix(&covmx);

  ResizeMatrix(mo, mi->row, mi->col);

  NewDVector(&x, mi->col);
  NewDVector(&p, mi->col);

  for(i = 0; i < mi->row; i++){
    for(j = 0; j < mi->col; j++){
      x->data[j] = mi->data[i][j]-colavg->data[j];
    }

    DVectorMatrixDotProduct(inv_cov, x, p);

    for(j = 0; j < mi->col; j++){
      (*mo)->data[i][j] = p->data[j];
    }

    DVectorSet(p, 0.f);
  }

  DelDVector(&p);
  DelDVector(&x);
  DelMatrix(&inv_cov);
  DelDVector(&colavg);
}

size_t square_to_condensed_index(size_t i, size_t j, size_t n)
{
  if(i == j)
    return n+1;
  else{
    size_t ii;
    size_t jj;
    if(i < j){
      ii = j;
      jj = i;
    }
    else{
      ii = i;
      jj = j;
    }
    return n*jj - jj*(jj+1)/2 + ii - 1 - jj;
  }
}

