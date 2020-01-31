/*
**      Author: yutato.gz@gmail.com
**      Date:   20 Apire 2019
**      File:   hmm.h
**      Purpose: datastructures used for HMM. 
**      Version: v0.1
*/

#ifndef hmm_h
#define hmm_h

typedef struct
{
    int n_states;   // 隐含状态数目
    int n_obser;    // 观测值数目
    double **A;     // 状态转移矩阵
    double **B;     // 混淆矩阵
    double *PI;     // 初始概率向量
} HMM;

void hmm_init(HMM *p_hmm,char *f_hmm);
void hmm_forward_algorithm(HMM *p_hmm, int *O, int T);
void hmm_viterbi_algorithm(HMM *p_hmm, int *O, int T);

#endif
