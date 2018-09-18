//516030910259 Xinpeng Liu
/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]){
    if (N==32){
      int i, j, p, q, tmp, k;
      for (i = 0; i < N/8; i++)
        for (j = 0; j < M/8; j++)
          for (p=0;p<8;p++){
            for (q=0;q<8;q++){
              if (j*8+q==i*8+p){
                tmp=A[i*8+p][j*8+q];
                k=j*8+q;
                continue;
                // for elements A[k][k], postpone to avoid miss
              }
              B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
            }
            if (i==j) B[k][k]=tmp;
            // for elements A[k][k], postpone to avoid miss
          }
    }// divide matrix by 8x8
    else if (N==64){
      int i, j, p, q;
      int t[4];
      for (i=0;i<N;i+=8){
        for (p=0;p<4;p++){
          for (q=0;q<8;q++)
            B[p][q+8]=A[i+p][i+q];
          for (q=0;q<8;q++)
            B[p][q+16]=A[i+p+4][i+q];
        }
        for (p=0;p<8;p++){
          for (q=0;q<4;q++)
            B[i+p][i+q]=B[q][p+8];
          for (q=4;q<8;q++)
            B[i+p][i+q]=B[q-4][p+16];
        }
      }
      //for diagram matrixes, use other matrixes as temporary space
      for (i = 0; i < N/8; i++)
        for (j = 0; j < M/8; j++){
          if (i==j) continue;
          for (p=0;p<4;p++){
            for (q=0;q<4;q++)
              B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
            for (q=4;q<8;q++)
              B[j*8+q-4][i*8+p+4]=A[i*8+p][j*8+q];
          }//use the top right 4x4 matrix of B as temporary space
          for (p=0;p<4;p++){
            for (q=4;q<8;q++)
              t[q-4]=B[j*8+p][i*8+q];
              //store the pth line of the top right matrix in t
            for (q=4;q<8;q++)
              B[j*8+p][i*8+q]=A[i*8+q][j*8+p];
              //fill the pth line of the top right matrix with right values
            for (q=4;q<8;q++)
              B[j*8+p+4][i*8+q-4]=t[q-4];
              //fill the bottom left matrix with t
          }
          for (p=4;p<8;p++)
            for (q=4;q<8;q++)
              B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
              //fill the bottom right matrix with right values
        }
    }//divide matrix by 8x8, and change the order of handling 4x4 matrixes
    else {
      int i, j, p, q, tmp, k;
      for (i = 0; i <= N/16; i++)
        for (j = 0; j <= M/16; j++)
          for (p = 0; p < 16 && i*16+p<N; p++){
            for (q = 0; q < 16 && j*16+q<M; q++){
              if (j*16+q == i*16+p){
                tmp=A[i*16+p][j*16+q];
                k=j*16+q;
                continue;
                // for elements A[k][k], postpone to avoid miss
              }
              B[j*16+q][i*16+p] = A[i*16+p][j*16+q];
            }
            if (i==j) B[k][k]=tmp;
            // for elements A[k][k], postpone to avoid miss
          }
    }//divide matrix by 16x16
}

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

 /*
  * trans - A simple baseline transpose function, not optimized for the cache.
  */
 char trans_desc[] = "Simple row-wise scan transpose";
 void trans(int M, int N, int A[N][M], int B[M][N])
 {
     if (N==32){
       int i, j, p, q, tmp, k;
       for (i = 0; i < N/8; i++)
         for (j = 0; j < M/8; j++)
           for (p=0;p<8;p++){
             for (q=0;q<8;q++){
               if (j*8+q==i*8+p){
                 tmp=A[i*8+p][j*8+q];
                 k=j*8+q;
                 continue;
               }
               B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
             }
             if (i==j) B[k][k]=tmp;
           }
     }
     else if (N==64){
       int i, j, p, q;
       int t[8];
       for (i = 0; i < N/8; i++)
         for (j = 0; j < M/8; j++){
           for (p=0;p<4;p++){
             for (q=0;q<4;q++)
               B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
             for (q=4;q<8;q++){
               t[q-4]=A[i*8][j*8+q];
               t[q]=A[i*8+1][j*8+q];
             }
           }
           for (p=4;p<8;p++)
             for (q=0;q<4;q++)
               B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
           for (p=4;p<8;p++)
             for (q=4;q<8;q++)
               B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
           for (q=4;q<8;q++){
             B[j*8+q][i*8]=t[q-4];
             B[j*8+q][i*8+1]=t[q];
           }
           for (p=2;p<4;p++)
             for (q=4;q<8;q++)
               B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
         }
     }
     else {
       int i, j, p, q, tmp, k;
       for (i = 0; i <= N/16; i++)
         for (j = 0; j <= M/16; j++)
           for (p = 0; p < 16 && i*16+p<N; p++){
             for (q = 0; q < 16 && j*16+q<M; q++){
               if (j*16+q == i*16+p){
                 tmp=A[i*16+p][j*16+q];
                 k=j*16+q;
                 continue;
               }
               B[j*16+q][i*16+p] = A[i*16+p][j*16+q];
             }
             if (i==j) B[k][k]=tmp;
           }
     }
 }


  /*
   * Slightly changed.
   */
  char trans_sc[] = "Slightly changed";
  void trans_scf(int M, int N, int A[N][M], int B[M][N])
  {
      if (N==32){
        int i, j, p, q, tmp, k;
        for (i = 0; i < N/8; i++)
          for (j = 0; j < M/8; j++)
            for (p=0;p<8;p++){
              for (q=0;q<8;q++){
                if (j*8+q==i*8+p){
                  tmp=A[i*8+p][j*8+q];
                  k=j*8+q;
                  continue;
                }
                B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
              }
              if (i==j) B[k][k]=tmp;
            }
      }
      else if (N==64){
        int i, j, p, q;
        int t[4];
        for (i=0;i<N;i+=8){
          for (p=0;p<4;p++){
            for (q=0;q<8;q++)
              B[p][q+8]=A[i+p][i+q];
            for (q=0;q<8;q++)
              B[p][q+16]=A[i+p+4][i+q];
          }
          for (p=0;p<8;p++){
            for (q=0;q<4;q++)
              B[i+p][i+q]=B[q][p+8];
            for (q=4;q<8;q++)
              B[i+p][i+q]=B[q-4][p+16];
          }
        }
        for (i = 0; i < N/8; i++)
          for (j = 0; j < M/8; j++){
            if (i==j) continue;
            for (p=0;p<4;p++){
              for (q=0;q<4;q++)
                B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
              for (q=4;q<8;q++)
                B[j*8+q-4][i*8+p+4]=A[i*8+p][j*8+q];
            }
            for (p=0;p<4;p++){
              for (q=4;q<8;q++)
                t[q-4]=B[j*8+p][i*8+q];
              for (q=4;q<8;q++)
                B[j*8+p][i*8+q]=A[i*8+q][j*8+p];
              for (q=4;q<8;q++)
                B[j*8+p+4][i*8+q-4]=t[q-4];
            }
            for (p=4;p<8;p++)
              for (q=4;q<8;q++)
                B[j*8+q][i*8+p]=A[i*8+p][j*8+q];
          }
      }
      else {
        int i, j, p, q, tmp, k;
        for (i = 0; i <= N/16; i++)
          for (j = 0; j <= M/16; j++)
            for (p = 0; p < 16 && i*16+p<N; p++){
              for (q = 0; q < 16 && j*16+q<M; q++){
                if (j*16+q == i*16+p){
                  tmp=A[i*16+p][j*16+q];
                  k=j*16+q;
                  continue;
                }
                B[j*16+q][i*16+p] = A[i*16+p][j*16+q];
              }
              if (i==j) B[k][k]=tmp;
            }
      }
  }

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);
    registerTransFunction(trans_scf, trans_sc);

}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}
