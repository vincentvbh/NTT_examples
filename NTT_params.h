#ifndef NTT_PARAMS_H
#define NTT_PARAMS_H

#define ARRAY_N 256

#define NTT_N 256u
#define LOGNTT_N 8

// Q1
#define Q1 7681
#define omegaQ1 7146
// invomegaQ1 = omegaQ^{-1} mod Q1
#define invomegaQ1 7480
// RmodQ1 = 2^16 mod^{+-} Q1
#define RmodQ1 (-3593)
// Q1prime = -Q1^{-1} mod^{+-} 2^16
#define Q1prime 7679
// invNQ1 = NTT_N^{-1} mod Q1
#define invNQ1 7651


#endif





