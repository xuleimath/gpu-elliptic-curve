#ifndef COMMON_DEF_H
#define COMMON_DEF_H

/*********************************************************************************
   Consider only ECDLP over GF(2^163)
   each field element is saved in integer array, specifically, unsigned int[6]
   
   irreducible polynomial is fixed: z^163+z^8+z^2+z+1
                                    {0x8,0x0,0x0,0x0,0x0,0x107}
**********************************************************************************/


extern "C"
{


typedef unsigned int int163;

typedef unsigned int poly163;
typedef unsigned int poly326;

typedef unsigned int elem163;

typedef unsigned int pt163;

typedef unsigned int op_state;


#define SIZE 6

// each finite field element and big integer occupy 6 32-bits words
#define WORDNUM 6

// the order of base point is N =0x 04 00000000 0000000 00020108 A2E0CC0D 99F8A5EF
// use our 28-bits limb, N = 0400000 0000000 0000000 20108A2 E0CC0D9 9F8A5EF
#define N0 0x9F8A5EF
#define N1 0xE0CC0D9
#define N2 0x20108A2
#define N3 0x0000000
#define N4 0x0000000
#define N5 0x0400000

// the irreducible polynomial is p(x) = z^163+z^8+z^2+z+1
// use 32-bits limb, P = 0x 0008 0000 0000 0000 0000 0107
#define P0 0x0107
#define P1 0x0000
#define P2 0x0000
#define P3 0x0000
#define P4 0x0000
#define P5 0x0008

// mask, get the highest 4 bits of a 32-bits word
#define HIGH4BITS 0xF0000000

// mask, the 29th bit
#define _29BIT 0x10000000

#define op_success  1
#define op_failure  0

};

#endif




