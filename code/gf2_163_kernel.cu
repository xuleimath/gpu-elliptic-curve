#ifndef GF2_163_CU
#define GF2_163_CU

/*
   irreducible polynomial is fixed: z^163+z^8+z^2+z+1
   {0x8,0x0,0x0,0x0,0x0,0x107}
*/
#include <device_functions.h>
#include <cutil_inline.h>


#include "gf2_163_common_def.h"

/*
	multiply lhs by rhs, where rhs=0,1.
	the result is lhs or 0
*/
__device__ void
poly_mulbit(poly163 *dst, const poly163 *lhs, unsigned int rhs)
{
	unsigned int t = rhs * 0xFFFFFFFF;
	
	dst[0] = lhs[0] & t;
	dst[1] = lhs[1] & t;
	dst[2] = lhs[2] & t;
	dst[3] = lhs[3] & t;
	dst[4] = lhs[4] & t;
	dst[5] = lhs[5] & t;
}

__device__ unsigned int 
brev(unsigned int a)
{
	a = ((a >>  1) & 0x55555555) + ((a & 0x55555555) <<  1);
	a = ((a >>  2) & 0x33333333) + ((a & 0x33333333) <<  2);
	a = ((a >>  4) & 0x0F0F0F0F) + ((a & 0x0F0F0F0F) <<  4);
	a = ((a >>  8) & 0x00FF00FF) + ((a & 0x00FF00FF) <<  8);
	a = ( a >> 16              ) + ( a               << 16);
	return a;

}

__device__ bool 
poly_isone(const poly163 *src)
{
	unsigned int t = src[0] | src[1] | src[2] | src[3] | src[4];

	return (t == 0) && (src[5] == 1);
}

__device__ unsigned int
poly_degree(const poly163 *src)
{
	int i = 5;
	while(i >= 0)
		if(src[i] == 0)
			--i;

	if(i < 0)
		return 0;
	
 	//return 8*i + 31 - __brev(src[i]);
}

__device__ void 
poly_add(poly163 *dst, const poly163 *lhs, const poly163 *rhs)
{
	dst[0] = lhs[0] ^ rhs[0];
	dst[1] = lhs[1] ^ rhs[1];
	dst[2] = lhs[2] ^ rhs[2];
	dst[3] = lhs[3] ^ rhs[3];
	dst[4] = lhs[4] ^ rhs[4];
	dst[5] = lhs[5] ^ rhs[5];
};

__device__ void 
poly_selfadd(poly163 *dst, const poly163 *rhs)
{
	dst[0] ^= rhs[0];
	dst[1] ^= rhs[1];
	dst[2] ^= rhs[2];
	dst[3] ^= rhs[3];
	dst[4] ^= rhs[4];
	dst[5] ^= rhs[5];
}

// add the irreducible polynomial to dst
__device__ void
poly_selfaddire(poly163 *dst)
{
	dst[0] ^= 0x107;
	dst[5] ^= 0x8;

}


__device__ void
poly_mul1(unsigned int *c, unsigned int a, unsigned int b)
{
	unsigned int hi, lo, t;
	unsigned int A[8];
	A[0] = 0;
	A[1] = a;
	A[2] = A[1] << 1;
	A[3] = A[2] ^ A[1];
	A[4] = A[2] << 1;
	A[5] = A[4] ^ A[1];
	A[6] = A[3] << 1;
	A[7] = A[6] ^ A[1];
	lo = A[b & 7]; t = A[(b >> 3) & 7]; hi = t >> 29; lo ^= t << 3;
	t = A[(b >> 6) & 7]; hi ^= t >> 26; lo ^= t << 6;
	t = A[(b >> 9) & 7]; hi ^= t >> 23; lo ^= t << 9;
	t = A[(b >> 12) & 7]; hi ^= t >> 20; lo ^= t << 12;
	t = A[(b >> 15) & 7]; hi ^= t >> 17; lo ^= t << 15;
	t = A[(b >> 18) & 7]; hi ^= t >> 14; lo ^= t << 18;
	t = A[(b >> 21) & 7]; hi ^= t >> 11; lo ^= t << 21;
	t = A[(b >> 24) & 7]; hi ^= t >> 8; lo ^= t << 24;
	t = A[(b >> 27) & 7]; hi ^= t >> 5; lo ^= t << 27;
	t = A[b >> 30]; hi ^= t >> 2; lo ^= t << 30;
	if (a >> 31) hi ^= ((b & 0xb6db6db6UL) >> 1);
	if ((a >> 30) & 1) hi ^= ((b & 0x24924924UL) >> 2);
	c[0] = lo;    c[1] = hi;

}

__device__ void 
poly_mul3 (unsigned int *c, const unsigned int *a, const unsigned int *b)
{
	unsigned int d0[2], d1[2], d2[2], d01[2], d02[2], d12[2];

	poly_mul1(d0, a[0], b[0]);
	poly_mul1(d1, a[1], b[1]);
	poly_mul1(d2, a[2], b[2]);
	poly_mul1(d01, a[0]^a[1], b[0]^b[1]);
	poly_mul1(d02, a[0]^a[2], b[0]^b[2]);
	poly_mul1(d12, a[1]^a[2], b[1]^b[2]);


	c[0] = d0[0];
	c[1] = d0[1] ^ d01[0] ^ d1[0] ^ d0[0];
	c[2] = d01[1] ^ d1[1] ^ d0[1] ^ d02[0] ^ d2[0] ^ d0[0] ^ d1[0];
	c[3] = d02[1] ^ d2[1] ^ d0[1] ^ d1[1] ^ d12[0] ^ d1[0] ^ d2[0];
	c[4] = d12[1] ^ d1[1] ^ d2[1] ^ d2[0];
	c[5] = d2[1];

}

__device__ void 
poly_mul6(unsigned int *c, const unsigned int *a, const unsigned int *b)
{
	unsigned int hs0[3], hs1[3];
	unsigned int hl2[6];

	hs0[0] = a[0] ^ a[3];   
	hs0[1] = a[1] ^ a[4];
	hs0[2] = a[2] ^ a[5];
	hs1[0] = b[0] ^ b[3];   
	hs1[1] = b[1] ^ b[4];
	hs1[2] = b[2] ^ b[5];

	poly_mul3(c, a, b);   
	poly_mul3(c+6, a+3, b+3); 
	poly_mul3(hl2, hs0, hs1);  

	hl2[0] = hl2[0] ^ c[0] ^ c[6];
	hl2[1] = hl2[1] ^ c[1] ^ c[7];
	hl2[2] = hl2[2] ^ c[2] ^ c[8];
	hl2[3] = hl2[3] ^ c[3] ^ c[9];
	hl2[4] = hl2[4] ^ c[4] ^ c[10];
	hl2[5] = hl2[5] ^ c[5] ^ c[11];

	c[3] ^= hl2[0];
	c[4] ^= hl2[1];
	c[5] ^= hl2[2];
	c[6] ^= hl2[3];
	c[7] ^= hl2[4];
	c[8] ^= hl2[5];
}

__device__ void
poly_mul(poly326 *dst, const poly163 *lhs, const poly163 *rhs)
{
	poly_mul6(dst, lhs, rhs);
}

// mod  z^163+z^8+z^2+z+1
__device__ void
poly_mod(poly163 *dst, const poly326 *src)
{
	unsigned int tp[11];
	
	// save polynomial src to be reduced in tp
	#pragma unroll
	for (int i=0; i<12; ++i)
	{
		tp[i] = src[i];
	}

	// reduce tp
	#pragma unroll
	for (int i=10; i>=6; ++i)
	{
		unsigned int t = tp[i];
		tp[i-6] = tp[i-6] ^ (t << 29) ^ (t << 30) ^ (t << 31);
		tp[i-5] = tp[i-5] ^ (t >> 3)  ^ (t >> 2)  ^ ( t >> 1) ^(t << 5);
		tp[i-4] = tp[i-4] ^ (t >> 27);
	}
	unsigned int t = (tp[5] >> 3);
	tp[0] = tp[0] ^ (t << 8) ^ (t << 2) ^ (t << 1) ^ t;
	tp[1] = tp[1] ^ (t >> 24) ^ (t >> 30) ^ (t >> 31);
	tp[5] &= 0x07;

	// save tp to result polynomial dst
	#pragma  unroll
	for (int i=0; i<6; ++i)
	{
		dst[i] =tp[i];
	}
}


__device__ void
poly_divx(poly163 *dst, const poly163 *src)
{
	dst[0] = (src[0] >> 1) ^ (src[1] & 0x01);
	dst[1] = (src[1] >> 1) ^ (src[2] & 0x01);
	dst[2] = (src[2] >> 1) ^ (src[3] & 0x01);
	dst[3] = (src[3] >> 1) ^ (src[4] & 0x01);
	dst[4] = (src[4] >> 1) ^ (src[5] & 0x01);
	dst[5] = (src[5] >> 1);
}

__device__ void
poly_selfdivx(poly163 *src)
{
	src[0] >>= 1;	src[0] ^= (src[1] & 0x01);
	src[1] >>= 1;	src[1] ^= (src[2] & 0x01);
	src[2] >>= 1;	src[2] ^= (src[3] & 0x01);
	src[3] >>= 1;	src[3] ^= (src[4] & 0x01);
	src[4] >>= 1;	src[4] ^= (src[5] & 0x01);
	src[5] >>= 1;
}

// try divide src by x
/* 
   if src is divisible by x, change src and return op_success
   else keep src unchanged and return op_failure
*/
__device__ op_state
poly_tryselfdivx(poly163 *src)
{
	if ((src[0] & 0x1) == 0)	// divisible by x
	{
		src[0] >>= 1;	src[0] ^= (src[1] & 0x01);
		src[1] >>= 1;	src[1] ^= (src[2] & 0x01);
		src[2] >>= 1;	src[2] ^= (src[3] & 0x01);
		src[3] >>= 1;	src[3] ^= (src[4] & 0x01);
		src[4] >>= 1;	src[4] ^= (src[5] & 0x01);
		src[5] >>= 1;

		return op_success;
	}
	else						// not divisible by x
		return op_failure;
}

// dst * src = 1 mod z^163+z^8+z^2+z+1
__device__ void
poly_modinv(poly163 *dst, const poly163 *src)
{
	// initialization
	poly163 u[6];		// src is saved in u
	poly163 v[6] = {0x8, 0x0, 0x0, 0x0, 0x0, 0x107};	// irreducible polynomial z^163+z^8+z^2+z+1

	poly163 g1[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x1};
	poly163 g2[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
	
	#pragma unroll
	for (int i=0; i<6; ++i)
	{
		u[i] = src[i];
	}

	// main loop
	while (!poly_isone(u) && !poly_isone(v))
	{
		while(poly_tryselfdivx(u))
		{
			if((g1[0] & 0x1) == 0)
				poly_selfdivx(g1);
			else
			{
				poly_selfaddire(g1);
				poly_selfdivx(g1);
			}
		}
		while (poly_tryselfdivx(v))
		{
			if((g2[0] & 0x1) == 0)
				poly_selfdivx(g2);
			else
			{
				poly_selfaddire(g2);
				poly_selfdivx(g2);
			}
		}
	}
}

#endif

