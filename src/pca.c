/* pca.c
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

#include <math.h>
#include <string.h>

#include "memwrapper.h"
#include "numeric.h"
#include "pca.h"
#include "vector.h"
#include "matrix.h"
#include "scientificinfo.h"

void NewPCAModel(PCAMODEL** m)
{
  (*m) = xmalloc(sizeof(PCAMODEL));
  initMatrix(&((*m)->scores));
  initMatrix(&((*m)->loadings));
  initMatrix(&((*m)->dmodx));
  initDVector(&((*m)->varexp));
  initDVector(&((*m)->colaverage));
  initDVector(&((*m)->colscaling));
}

void DelPCAModel(PCAMODEL** m)
{
  DelDVector(&((*m)->colscaling));
  DelDVector(&((*m)->colaverage));
  DelDVector(&((*m)->varexp));
  DelMatrix(&((*m)->dmodx));
  DelMatrix(&((*m)->loadings));
  DelMatrix(&((*m)->scores));
  xfree((*m));
}


void calcVarExpressed(double ss, dvector *eval, dvector **varexp) /* ss is the sum of squares, eval = eigenvalue  varexp is an object that is resized for each component */
{
  size_t i;
  /*double a;*/
  for(i = 0; i < eval->size; i++){
    DVectorAppend(varexp, (getDVectorValue(eval, i)/ss) * 100);
    #ifdef DEBUG
    printf("Variance expressed for PC %u\t %f\n", (unsigned int)i, (getDVectorValue(eval, i)/ss) * 100);
    #endif
  }

  /*The latest value is calculated from the sum of the previous explained variance
  a = 0.f;
  for(i = 0; i < (*varexp)->size; i++)
    a += getDVectorValue((*varexp), i);

  DVectorAppend(varexp, 100 - a);
  */
}

double calcObjectDistance(matrix *m)
{
  size_t i, j;
  double sum=0.f;
  for(i = 0; i < m->row; i++){
    for(j = 0; j < m->col; j++)
      sum += m->data[i][j]*m->data[i][j];
  }

  #ifdef DEBUG
  printf("Distance of Noise Matrix: %f\n", sum);
  #endif
  return sum;
}

/*
 *     NIPALS Algorithm for PCA
 * Input:
 *     E = Mean Centered Data Matrix
 * Local Variable:
 *   t, t_old = score vector with the size of the number of row of the E matrix;
 *   p = loadings vector with the size of the number of column of the E matrix;
 *
 * 1.  t := e_i                 Select a column vector e_i of the matrix E with the best value of column average and copy it to the vector t
 * 2.  p' := (t'X)/(t't) if t't > 1.0E-9         Project the matrix E onto t in order to find the corresponding loading p
 * 3.  p' := p'/|p'|            Normalize the loading vector p to length 1
 * 4.  t_old := t
 *     t := (Xp)/(p'p) if p'p > 1.0E-9          Store the score vector t into uold and project the matrix E onto p in order to find corresponding score vector t
 * 5.  d := t_old-t             In order to check for the convergence of the process calculate the difference vector d as the difference between
 *                              the previous scores and the current scores. If the difference |d| is larger than a pre-defined threshold (e.g. 10-8) then return to step 2,
 *                              else store the score and loadings to an output and follow the step 6.
 * 6.  E := E - tp'             Remove the estimated PCA component (the product of the scores and the loadings) from E.
 *                              In order to estimate the other PCA components repeat this procedure from step 1 using the matrix E -tp' as the new E
 *
 * This algorithm deal with missing values.
 * Indeed when a missing value on X input matrix is found,
 * will be skipped on the calculation of p and t.
 * See page 1 of the following document:
 * https://cran.r-project.org/web/packages/nipals/vignettes/nipals_algorithm.pdf
 *
 */
void PCA(matrix *mx, size_t scaling, size_t npc, PCAMODEL* model, ssignal *s)
{
  size_t i, j, pc;
  dvector *t;
  dvector *p;
  dvector *colvar;
  dvector *eval; /* t't */

  double min, max, mod_p, mod_t_old, mod_t_new, ss;

  matrix *E; /* data matrix of autoscaled / mean centred object */
  NewMatrix(&E, mx->row, mx->col);

  /* check if matrix have nan or infinite and convert them to MISSING */
  MatrixCheck(mx);

  /*CENTERING */
  MatrixColAverage(mx, &(model->colaverage));
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      if(FLOAT_EQ(mx->data[i][j], MISSING, 1e-1)){
        continue;
      }
      else{
        E->data[i][j] = mx->data[i][j] - model->colaverage->data[j];
      }
    }
  }

  if(scaling > 0){
    if(scaling == 1){
      MatrixColSDEV(mx, &(model->colscaling));
    }
    else if(scaling == 2){
      MatrixColRMS(mx, &(model->colscaling));
    }
    else if(scaling == 3){ /* PARETO Autoscaling */
      MatrixColSDEV(mx, &(model->colscaling));
      for(i = 0; i < model->colscaling->size; i++){
        model->colscaling->data[i] = sqrt(model->colscaling->data[i]);
      }
    }
    else if(scaling == 4){ /* Range Scaling */
      for(i = 0; i < mx->col; i++){
        MatrixColumnMinMax(mx, i, &min, &max);
        DVectorAppend(&model->colscaling, (max - min));
      }
    }
    else if(scaling == 5){ /* Level Scaling  */
      DVectorCopy(model->colaverage, &model->colscaling);
    }
    else{
      for(i = 0; i < model->colaverage->size; i++){
        DVectorAppend(&model->colscaling, 1.0);
      }
    }

    for(j = 0; j < E->col; j++){
      if(FLOAT_EQ(getDVectorValue(model->colscaling, j), 0, EPSILON)){
        for(i = 0; i< E->row; i++){
          E->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < E->row; i++){
          if(FLOAT_EQ(E->data[i][j], MISSING, 1e-1)){
            continue;
          }
          else{
            E->data[i][j] /= model->colscaling->data[j];
          }
        }
      }
    }
  }

   /* if the number of principal component selected is major of the permitted */
  if(npc > E->col)
    npc = E->col;    /* set the value to the max value */


  ss = 0.f;
  for(i = 0; i < E->row; i++){
    for(j = 0; j < E->col; j++)
      ss += square(E->data[i][j]);
  }

  NewDVector(&t, E->row);
  NewDVector(&p, E->col);
  NewDVector(&eval, npc);


  ResizeMatrix(&(model->scores), E->row, npc);
  ResizeMatrix(&(model->loadings), E->col, npc);
  ResizeMatrix(&(model->dmodx), E->row, npc);

  /*setDVectorValue(obj_d, 0, calcObjectDistance(E)); */
  for(pc = 0; pc < npc; pc++){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      #ifdef DEBUG
      printf("## Computing Component %d\n", (int)pc);
      #endif

      /* Init Step
      * calculate variance for each column for select after the vector with the greater value
      */

      /*
      * During the first iteration the variance for all the column is the same: 1
      * So we select a random column. With this can change the sign of the plot
      *
      if(pc == 0){
        srand(time(0));
        j = (size_t) rand() % E->col;
        printf("Selected %u\n", (unsigned int)j);
      }
      else{
        initDVector(&colvar);
        MatrixColVar(E, &colvar);

        Step 1: select the column vector t with the largest column variance
        j = 0;
        for(i = 1; i < E->col; i++){
          if(colvar->data[i] > colvar->data[j])
            j = i;
          else
            continue;
        }

        DelDVector(&colvar);
      }

      * In order to prevent a random sign numbering is preferable
      * select in the first calculation the first column and afther getting
      * the best one.
      */
      /*
      if(pc == 0){
        we do not know if the first variable is 0 so we have to iterate to find the first column with a value != 0
        j = 0;
        while(1){
          if(j < model->colaverage->size){
            if(FLOAT_EQ(getDVectorValue(model->colaverage, j), 0.f, PCACONVERGENCE)){
              j++;
            }
            else{
              break;
            }
          }
          else{
            j = 0;
            break;
          }
        }
      }
      else{
        initDVector(&colvar);
        MatrixColVar(E, &colvar);

        Step 1: select the column vector t with the largest column variance
        j = 0;
        for(i = 1; i < E->col; i++){
          if(getDVectorValue(colvar, i) > getDVectorValue(colvar, j))
            j = i;
          else
            continue;
        }

        DelDVector(&colvar);
      }
      */

      initDVector(&colvar);
      MatrixColVar(E, &colvar);

      /* Step 1: select the column vector t with the largest column variance */
      j = 0;
      for(i = 1; i < E->col; i++){
        if(colvar->data[i] > colvar->data[j])
          j = i;
        else
          continue;
      }

      DelDVector(&colvar);

      /* copy the vector to the score mx.t_old for computing loadings */
      for(i = 0; i < E->row; i++)
        t->data[i] = E->data[i][j];

      /* End Step 1 */

      while(1){
        /* Step 2: projection of t' in E (t'*E) */
        DVectorMatrixDotProduct(E, t, p);
        /* calc the vectors product t'*t = Sum(t[i]^2) */
        mod_t_old = DVectorDVectorDotProd(t, t);

//         if(mod_t_old > 1.0E-9){
          /* division of (t'*E)/t'*t (mx.p/mx.mod_t_old) for calculate the p' vector that represents the loadings */
        for(i = 0; i < p->size; i++)
          p->data[i] /= mod_t_old; /* now p' is the loadings vector calculated */
//         }

        /* End Step 2 */

        /* Step 3: normalizing the p loadings by overwriting the mx.p vector */
        DVectNorm(p, p);
        /* End Step 3 */

        /* Step 4: t_new = (X*p)/(p'*p) */
        DVectorSet(t, 0.f);
        MatrixDVectorDotProduct(E, p, t);

        /* p'*p = Sum(p[i]^2) */
        mod_p = DVectorDVectorDotProd(p, p);
//         if(mod_p > 1.0E-9){
        for(i = 0; i < E->row; i++)
          t->data[i] /= mod_p;
//           setDVectorValue(t, i, getDVectorValue(t, i)/mod_p);
//         }

        /* End Step 4 */

        /* Step 5: In order to check for the convergence of the process, the difference between the
        * modules of the new score vector with the oldest must be minus or equal tu 10^-8
        */
        mod_t_new = DVectorDVectorDotProd(t, t);

        #ifdef DEBUG
        printf("new eigen: %f\told eigen:%f\n", mod_t_new, mod_t_old);
        printf("new score calculated\n");
        PrintDVector(t);
        puts("....................");
        puts("loadings calculated from 'new score'");
        PrintDVector(p);
        puts("....................");
        #endif

        if(mod_t_new - mod_t_old <= PCACONVERGENCE && !isnan(mod_t_new)){
          /* copy the loadings and score to the output data matrix */
          for(i = 0; i < t->size; i++){
            model->scores->data[i][pc] = t->data[i];
          }

          for(i = 0; i < p->size; i++){
            model->loadings->data[i][pc] = p->data[i];
          }

          /*
          MatrixAppendCol(&(model->scores), t);
          MatrixAppendCol(&(model->loadings), p);*/

          /* Step 6.  E := E - tp' Remove the estimated PCA component
          * (the product of the scores and the loadings) from E. */

          for(i = 0; i < E->row; i++){
            for(j = 0; j < E->col; j++){
              E->data[i][j] -= t->data[i]*p->data[j];
              /* Calculate the dmodx
               * with the following function:
               * DmodX = SQRT(∑k (E_i_k)^2/ degree of fredom)
               * where
               * i is the index for the object
               * k is the index for the variable
               * npc is the number of principal components
               *
               * This value of DmodX is absolute. to obtain a normalized
               * DmodX it simple to divide for the variance of the residuals E
               * DmodX /= Variance(E)
               */
               model->dmodx->data[i][pc] += square(E->data[i][j]);
           }
           /*
            * the degree of freedom are calculated according to the formula
            *  (n_variables-npc)  * (nobjects/(nobject-1)
            */
            double degree_freedom;
            if(E->col < pc+1 && E->row > pc+1){
                degree_freedom = (float) (E->col-(pc+1)) / (float)(E->row/((E->row-(pc+1))-1));
            }
            else{
                degree_freedom = (pc+1);
            }
            model->dmodx->data[i][pc] = sqrt(model->dmodx->data[i][pc]/degree_freedom);

          }
          /* End Step 6 */


          eval->data[pc] = mod_t_new;
          break;
        }
        else{
          continue;
          /* End Step 5 */
        }
      }
      /*Reset all the variables*/
      DVectorSet(p, 0.f);
      DVectorSet(t, 0.f);
      mod_p = mod_t_old = mod_t_new = 0.f;
    }
  }

  calcVarExpressed(ss, eval, &(model->varexp));

  DelDVector(&eval);
  DelDVector(&t);
  DelDVector(&p);
  DelMatrix(&E);
}

/*
 * In the original PCA model
 *  X = TP' + E
 *
 * for a external dataset X*
 * X*=T*P' + E*
 *
 * and using the NIPALS algorithm, for a certain dimension a:
 *
 * t_a*=X*p_a/p_a'p_a  (look at the 4.step of NIPALS Algorithm)
 *
 */
void PCAScorePredictor(matrix *mx, PCAMODEL *model, size_t npc, matrix **pscores)
{
  size_t i ,j, k;
  matrix *E; /* data matrix of mean centred object */
  dvector *p;
  dvector *t;
  double mod_p;


  if(npc > model->loadings->col)
    npc = model->loadings->col;

  ResizeMatrix(pscores, mx->row, npc);

  NewMatrix(&E, mx->row, mx->col);

  if(model->colaverage->size > 0){
    for(j = 0; j < E->col; j++){
      for(i = 0; i < E->row; i++){
        E->data[i][j] = mx->data[i][j] - model->colaverage->data[j];
      }
    }
  }
  else{
    for(j = 0; j < mx->col; j++){
      for(i = 0; i < mx->row; i++){
        E->data[i][j] = mx->data[i][j];
      }
    }
  }

  if(model->colscaling->size > 0){
    for(j = 0; j < E->col; j++){
      if(getDVectorValue(model->colscaling, j) == 0){
        for(i = 0; i< E->row; i++){
          E->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < E->row; i++){
          E->data[i][j] /= model->colscaling->data[j];
        }
      }
    }
  }

  NewDVector(&t, E->row);
  NewDVector(&p, model->loadings->row);

  for(i = 0; i < npc; i++){
    /* Step 4: t = (X*p)/(p'*p) */
    /* copy the column of loadings to p */
    for(j = 0; j < model->loadings->row; j++)
      p->data[j] = model->loadings->data[j][i];

    mod_p = DVectorDVectorDotProd(p, p);

    MatrixDVectorDotProduct(E, p, t);

    /* p'*p = Sum(p[i]^2) */
    mod_p = DVectorDVectorDotProd(p, p);

    /*t_new->data[i] = t_new->data[i]/mod_p;*/
    for(j = 0; j < t->size; j++){
      t->data[j] /= mod_p;
      /* and set the result score to the output*/
      (*pscores)->data[j][i] = t->data[j];
    }

    /* Remove residual...*/
    /* Step 6.  E := E - tp' Remove the estimated PCA component (the product of the scores and the loadings) from E. */
    for(j = 0; j < E->row; j++){
      for(k = 0; k < E->col; k++){
        E->data[j][k] = E->data[j][k] - (t->data[j]*p->data[k]);
        //dmodx->data[j][i] += square(E->data[j][k]);
      }
      /*
      * the degree of freedom for a predicted object correspond to the number of pc
      */
      //dmodx->data[j][i] = sqrt(model->dmodx->data[j][i]/( i+1));
      /* End Step 6 */
    }

    DVectorSet(p, 0);
    DVectorSet(t, 0);
  }

  DelDVector(&p);
  DelDVector(&t);
  DelMatrix(&E);
}


/*
 *  X = t*p
 */
void PCAIndVarPredictor(matrix* t,
                        matrix* p,
                        dvector* colaverage,
                        dvector* colscaling,
                        size_t npc,
                        matrix** mx)
{

  size_t pc, i, j;

  /*Allocate the output array*/
  ResizeMatrix(mx, t->row, p->row);

  if(npc > t->col)
    npc = t->col;

  for(pc = 0; pc < npc; pc++){
    for(i = 0; i < t->row; i++){
      for(j = 0; j < p->row; j++){
        (*mx)->data[i][j] += t->data[i][pc]*p->data[j][pc];
      }
    }
  }

  if(colaverage->size > 0){
    if(colscaling->size > 0){
      for(i = 0; i < (*mx)->row; i++){
        for(j = 0; j < (*mx)->col; j++){
          (*mx)->data[i][j] *= colscaling->data[j];
          (*mx)->data[i][j] += colaverage->data[j];
        }
      }
    }
    else{
      for(i = 0; i < (*mx)->row; i++){
        for(j = 0; j < (*mx)->col; j++){
          (*mx)->data[i][j] += colaverage->data[j];
        }
      }
    }
  }
}


void PCARSquared(matrix *mx, PCAMODEL *model, size_t npc, dvector** r2)
{
  size_t pc, i, j;
  matrix *recalcy;
  matrix *recalcx;
  matrix *recalcscores;
  double n, d;

  if(npc > model->loadings->col)
    npc = model->loadings->col;

  for(pc = 1; pc <= npc; pc++){

    initMatrix(&recalcy);
    initMatrix(&recalcx);
    initMatrix(&recalcscores);

    PCAScorePredictor(mx, model, npc, &recalcscores);

    PCAIndVarPredictor(recalcscores, model->loadings, model->colaverage, model->colscaling, pc, &recalcx);

    /* r2 is cumulative because the X is all the matrix */
    n = d = 0.f;
    for(j = 0; j < recalcx->col; j++){
      for(i = 0; i < recalcx->row; i++){
        n += square(getMatrixValue(recalcx, i, j) - getMatrixValue(mx, i, j));
        d += square(getMatrixValue(mx, i, j) - getDVectorValue(model->colaverage, j));
      }
    }

    DVectorAppend(r2, 1 - (n/d));

    DelMatrix(&recalcy);
    DelMatrix(&recalcx);
    DelMatrix(&recalcscores);
  }
}


/*
 * 1) The method works dividing the dataset randomly into G groups. For example we have 5 groups: G1, G2, G3, G4, G5.
 * 2) Then the first groups G1 is taken out, computing a reduced model which is used to "predict" the values for the objects in the deleted group.
 * 3) The error in the prediction is measured in terms of a sum of squares of the prediction errors (PRESS) for this reduced models.
 * The whole procedure is repeated removing G2, G3... and accumulating all the partial prediction errors in a total PRESS.
 * The value of this error is compared with the data sum of squares (Seps) as:
 *
 *   R=PRESS/Seps
 *
 * where Seps, for dimension a, is computed as the sum of squares of the X matrix after removing the variance explained by the previous (a-1) PC's.
 *
 * The R value is calculated for every model dimensionality. When the value of R obtained is larger than 1.00, the incorporation of this PC does not improve the predictions and therefore this PC should not be included.
 *
 * A detailed description of the method can be found in: S. Wold, Cross-Validatory Estimation of the Number of Components in Factor and Principal Component Models, Technometrics 20, 397-405 (1978).
 *
 */
/* npc number of pricipal components, g = number of groups */
void PCARankValidation(matrix *mx, size_t npc, size_t scaling, size_t group, size_t iterations, dvector **r2, ssignal *s)
{
  size_t iterations_, i, j, k, n, g;
  matrix *gid;
  matrix *subX;
  PCAMODEL *subm;
  matrix *predictX;
  dvector *predictr2x;

  NewMatrix(&gid, group, (size_t)ceil(mx->row/(double)group));

  srand(group*mx->row*iterations);

  iterations_ = 0;
  while(iterations_ <  iterations){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /* Divide in group  all the Dataset */
      MatrixSet(gid, -1);

      /* step 1 generate the random groups */
      k = 0;
      for(i = 0; i <  gid->row; i++){
        for(j = 0; j <  gid->col; j++){
          do{
            n = (size_t)rand() % (mx->row);
          } while(ValInMatrix(gid, n) == 1 && k < (mx->row));
          if(k < mx->row){
            gid->data[i][j] = n;
            k++;
          }
          else
            continue;
        }
      }

      /*
      puts("Gid Matrix");
      PrintMatrix(gid);
      */

      /*step 2*/
      for(g = 0; g < gid->row; g++){ /*For aeach group */
        /* Estimate how many objects are inside the sub model without the group "g" */
        n = 0;
        for(i = 0; i < gid->row; i++){
          if(i != g){
            for(j = 0; j < gid->col; j++){
              if((int)getMatrixValue(gid, i, j) != -1)
                n++;
              else
                continue;
            }
          }
          else
            continue;
        }

        /*Allocate the submodel*/
        NewMatrix(&subX, n, mx->col);

        /* Estimate how many objects are inside the group "g" to predict*/
        n = 0;
        for(j = 0; j < gid->col; j++){
          if((int)getMatrixValue(gid, g, j) != -1)
            n++;
          else
            continue;
        }


        /*Allocate the */
        NewMatrix(&predictX, n, mx->col);


        /* copy the submodel values */

        for(i = 0, k = 0; i < gid->row; i++){
          if(i != g){
            for(j = 0; j < gid->col; j++){
              size_t a =  (size_t)gid->data[i][j]; /* get the row index */
              if(a != -1){
                for(n = 0; n < mx->col; n++){
                  subX->data[k][n] = mx->data[a][n];
                }
                k++;
              }
              else
                continue;
            }

          }
          else
            continue;
        }

        /* copy the objects to predict into predictmx*/
        for(j = 0, k = 0; j < gid->col; j++){
          size_t a = (size_t)gid->data[g][j];
          if(a != -1){
            for(n = 0; n < mx->col; n++){
              predictX->data[k][n] = mx->data[a][n];
            }
            k++;
          }
          else
            continue;
        }

        /*
        printf("Excuded the group number %u\n", (unsigned int)g);
        puts("Sub Model\nX:");
        PrintArray(subX);
        puts("Y:");
        PrintArray(subY);
        */

        NewPCAModel(&subm);

        PCA(subX, scaling, npc, subm, s);


        initDVector(&predictr2x);

        PCARSquared(predictX, subm, npc, &predictr2x);

        if(iterations_ == 0 && g == 0){
          for(i = 0; i < predictr2x->size; i++){
            DVectorAppend(r2, getDVectorValue(predictr2x, i));
          }
        }
        else{
          for(i = 0; i < predictr2x->size; i++){
            (*r2)->data[i] += predictr2x->data[i];
          }
        }

        DelDVector(&predictr2x);
        DelPCAModel(&subm);
        DelMatrix(&predictX);
        DelMatrix(&subX);
      }
      iterations_++;
    }
  }
  DelMatrix(&gid);

  /*Finalize the output by dividing for the number of iterations*/
  for(i = 0; i < (*r2)->size; i++){
    (*r2)->data[i] /= group*iterations;
  }
}

void GetResidualMatrix(matrix* mx, PCAMODEL* model, size_t pc, matrix** rmx)
{
  size_t i, j, k;
  ResizeMatrix(rmx, mx->row, mx->col);

  /*CENTERING */
  for(j = 0; j < mx->col; j++){
    for(i = 0; i < mx->row; i++){
      (*rmx)->data[i][j] = mx->data[i][j] - model->colaverage->data[j];
      /*setMatrixValue((*rmx), i, j, getMatrixValue(mx, i, j) - getDVectorValue(model->colaverage, j));*/
    }
  }

  if(model->colscaling->size > 0){
    for(j = 0; j < (*rmx)->col; j++){
      if(FLOAT_EQ(getDVectorValue(model->colscaling, j), 0, EPSILON)){
        for(i = 0; i< (*rmx)->row; i++){
          /*setMatrixValue((*rmx), i, j, 0);*/
          (*rmx)->data[i][j] = 0.f;
        }
      }
      else{
        for(i = 0; i < (*rmx)->row; i++){
          (*rmx)->data[i][j] /= model->colscaling->data[j];
          /*setMatrixValue((*rmx), i, j, getMatrixValue((*rmx), i, j)/getDVectorValue(model->colscaling, j));*/
        }
      }
    }
  }

   /* if the number of principal component selected is major of the permitted */
  if(pc > (*rmx)->col)
    pc = (*rmx)->col;    /* set the value to the max value */

  for(k = 0; k < pc; k++){
    for(i = 0; i < (*rmx)->row; i++){
      for(j = 0; j < (*rmx)->col; j++){
        (*rmx)->data[i][j] -= model->scores->data[i][k]*model->loadings->data[j][k];
        /*setMatrixValue((*rmx), i, j, getMatrixValue((*rmx), i, j)-(getMatrixValue(model->scores, i, k)*getMatrixValue(model->loadings, j, k)));*/
      }
    }
  }
}

void PrintPCA(PCAMODEL *m)
{
  size_t i, j;
  printf("Variance Explained\n");
  for(i = 0; i < m->varexp->size; i++){
    printf("PC%d: %.4f\n", (int)i+1, m->varexp->data[i]);
  }

  puts("Scores");
  for(j = 0; j < m->scores->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->scores->row; i++){
    for(j = 0; j < m->scores->col; j++){
      printf("%8.4f\t", m->scores->data[i][j]);
    }
    printf("\n");
  }

  puts("\nLoadings");
  for(j = 0; j < m->loadings->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->loadings->row; i++){
    for(j = 0; j < m->loadings->col; j++){
      printf("%8.4f\t", m->loadings->data[i][j]);
    }
    printf("\n");
  }

  puts("\nDmodX");
  for(j = 0; j < m->loadings->col; j++){
    printf("   PC %d\t", (int)j);
  }
  printf("\n");

  for(i = 0; i < m->dmodx->row; i++){
    for(j = 0; j < m->dmodx->col; j++){
      printf("%8.4f\t", m->dmodx->data[i][j]);
    }
    printf("\n");
  }
}
