#include <stdio.h>
#include "hmm.h"

int main(){
	printf("Start Main.\n");
    
    static HMM hmm_s;
    char *hmm_model = "t2.hmm";
    hmm_init(&hmm_s, hmm_model);

    //int obser_seq[3] = {0,1,0};
    int obser_seq[10] = {0,0,0,0,1,0,1,1,1,1};

    hmm_forward_algorithm(&hmm_s, obser_seq, 10);

    return 0;
}