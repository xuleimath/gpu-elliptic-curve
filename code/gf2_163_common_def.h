#ifndef COMMON_DEF_H
#define COMMON_DEF_H

/*********************************************************************************
   Consider only ECDLP over GF(2^163)
   each field element is saved in integer array, specifically, unsigned int[6]
   
   irreducible polynomial is fixed: z^163+z^8+z^2+z+1
                                    {0x8,0x0,0x0,0x0,0x0,0x107}
**********************************************************************************/



typedef unsigned int poly163;
typedef unsigned int poly326;

typedef unsigned int elem163;

typedef char op_state;

const char op_success = 1;
const char op_failure = 0;

#endif


