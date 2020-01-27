/*
**      Author: yutato.gz@gmail.com
**      Date:   20 Apire 2019
**      File:   hmm.c
**      Purpose: functions used for HMM. 
*/
#include <stdio.h>
#include <stdlib.h>
#include "hmm.h"

static double *dvector(int nl, int nh)
{
    double *v;

    v = (double *)calloc((unsigned)(nh - nl + 1), sizeof(double));
    if (!v)
        printf("allocation failure in dvector()");
    return v - nl;
}

static void free_dvector(double *v, int nl, int nh)
{
    free((char*) (v+nl));
}

static double **dmatrix(int nrl, int nrh, int ncl, int nch) 
{
    int i;
    double **m;

    while(!(m = (double **)calloc((unsigned)(nrh - nrl + 1), sizeof(double *))));
    if (!m)
        printf("allocation failure 1 in dmatrix()");
    m -= nrl;

    for (i = 0; i < nrh; i++)
    {
        m[i] = (double *)calloc((unsigned)(nch - ncl + 1), sizeof(double));
        if (!m[i])
            printf("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}

static void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
    int i;
    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

/**@brief   初始化函数 */
void hmm_init(HMM *p_hmm,char* f_hmm){
    printf("Init HMM.\n");
    // 读取hmm文件
    FILE *fp = fopen(f_hmm, "r");
    if(!fp)
        printf("file open failed");
    // 获取观测序列长度
    fscanf(fp, "n_obser= %d\n", &p_hmm->n_obser);
    printf("n_obser = %d\n", p_hmm->n_obser);
    // 获取隐含状态数目
    fscanf(fp, "n_states= %d\n", &p_hmm->n_states);
    printf("n_states = %d\n", p_hmm->n_states);
    // 获取状态转移矩阵
    fscanf(fp, "A:\n");
    printf("trans_mat = \n");
    p_hmm->A = (double **)dmatrix(1, p_hmm->n_states, 1, p_hmm->n_states);
    for (int i = 0; i < p_hmm->n_states; i++)
    {
        for (int j = 0; j < p_hmm->n_states; j++)
        {
            fscanf(fp, "%lf", &(p_hmm->A[i][j]));
            printf("  A[%d][%d] = %lf",i,j,p_hmm->A[i][j]);
        }
        printf("\n");
        fscanf(fp, "\n");
    }
    // 获取混淆矩阵
    fscanf(fp, "B:\n");
    printf("confusion_mat =\n");
    p_hmm->B = (double **)dmatrix(1, p_hmm->n_states, 1, p_hmm->n_obser);
    for (int j = 0; j < p_hmm->n_states; j++)
    {
        for (int k = 0; k < p_hmm->n_obser; k++)
        {
            fscanf(fp, "%lf", &(p_hmm->B[j][k]));
            printf("  B[%d][%d] = %lf",j,k,p_hmm->B[j][k]);
        }
        printf("\n");
        fscanf(fp, "\n");
    }
    // 获取初始概率向量
    fscanf(fp, "pi:\n");
    printf("start_pprob =\n");
    p_hmm->PI = (double *)dvector(1, p_hmm->n_states);
    for (int i = 0; i < p_hmm->n_states; i++)
    {
        fscanf(fp, "%lf", &(p_hmm->PI[i]));
        printf("  PI[%d] = %lf",i,p_hmm->PI[i]);
    }
    printf("\n");
    fclose(fp);
}

/**@brief   前向算法 */
void hmm_forward_algorithm(HMM *p_hmm,int *O,int T)
{
    printf("Start Forward algorithm.\n");
    float pprob = 0.0;
    int n_states = p_hmm->n_states;
    int n_obser  = p_hmm->n_obser;

    // 初始化前向概率矩阵
    float a[T][n_states];
    for(int t=0;t<T;t++)
        for(int i=0;i<n_states;i++)
            a[t][i] = 0;
    // 计算初值
    for(int i=0;i<n_states;i++){
        a[0][i] = p_hmm->PI[i]*p_hmm->B[i][O[0]];
        printf("  a[1][%d] = %.5f",i+1,a[0][i]);
    }
    printf("\n");
    // 递推
    for(int t=1;t<T;t++){
        for(int i=0;i<n_states;i++){
            float sum_i = 0;
            for(int j=0;j<n_states;j++){
                sum_i += a[t-1][j]*p_hmm->A[j][i];
            }
            a[t][i] = sum_i*p_hmm->B[i][O[t]];
            printf("  a[%d][%d] = %.5f",t+1,i+1, a[t][i]);
        }
        printf("\n");
    }
    // 终止
    for(int i=0;i<n_states;i++)
        pprob += a[T-1][i];

    printf("P(O|HMM) = %.5f\n",pprob);
    printf("End Forward algorithm.\n");
}

/**@brief   维特比算法 */
void hmm_viterbi_algorithm(HMM *p_hmm, int *O, int T)
{
    printf("Start Viterbi algorithm.\n");
    // 储存路径的概率
    double **delta = (double **)dmatrix(1, p_hmm->n_obser, 1, p_hmm->n_states);
    // 储存最优路径的节点
    
    
    
    
}
