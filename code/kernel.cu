/*
for each unsigned int, low 28 bits are used to represent an limb of big integer
*/
#include <cutil.h>
#include <cuda_runtime.h>
#include <device_functions.h>
#include "common_def.h"

/************************************************************************/
/*					Integer Operations                                  */
/************************************************************************/

__device__ void 
int_add(int163 *dst, const int163 *lhs, const int163 *rhs)
{
	dst[0] = lhs[0] + rhs[0];
	dst[1] = lhs[1] + rhs[1];
	dst[2] = lhs[2] + rhs[2];
	dst[3] = lhs[3] + rhs[3];
	dst[4] = lhs[4] + rhs[4];
	dst[5] = lhs[5] + rhs[5];

	// treat possible carry
	unsigned int t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;
}

__device__ void
int_addmod(int163 *dst, const int163 *lhs, const int163 *rhs)
{
	dst[0] = lhs[0] + rhs[0];
	dst[1] = lhs[1] + rhs[1];
	dst[2] = lhs[2] + rhs[2];
	dst[3] = lhs[3] + rhs[3];
	dst[4] = lhs[4] + rhs[4];
	dst[5] = lhs[5] + rhs[5];

	// treat possible carry
	unsigned int t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;

	// if dst too bigger than N, minus N. Maybe not fully reduced
	t = ((dst[5] & HIGH4BITS) == 0) * 0xFFFFFFFF;
	dst[5] -= (t & N5); dst[5] --;
	dst[4] |= _29BIT;	/*dst[4] -= (t & N4);*/ dst[4] --;	// note N3 == N4 == 0
	dst[3] |= _29BIT;	/*dst[3] -= (t & N3);*/ dst[3] --;
	dst[2] |= _29BIT;	dst[2] -= (t & N2);     dst[2] --;
	dst[1] |= _29BIT;	dst[1] -= (t & N1);		dst[1] --;
	dst[0] |= _29BIT;	dst[4] -= (t & N0);

	// treat possible carry on
	t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;
}

__device__ void
int_selfadd1(int163 *dst)
{
	dst[0] ++;

	// treat possible carry
	unsigned int t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;

	// if dst too bigger than N, minus N. Maybe not fully reduced
	t = ((dst[5] & HIGH4BITS) == 0) * 0xFFFFFFFF;
	dst[5] -= (t & N5); dst[5] --;
	dst[4] |= _29BIT;	/*dst[4] -= (t & N4);*/ dst[4] --;	// note N3 == N4 == 0
	dst[3] |= _29BIT;	/*dst[3] -= (t & N3);*/ dst[3] --;
	dst[2] |= _29BIT;	dst[2] -= (t & N2);     dst[2] --;
	dst[1] |= _29BIT;	dst[1] -= (t & N1);		dst[1] --;
	dst[0] |= _29BIT;	dst[4] -= (t & N0);

	// treat possible carry on
	t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;
}
__device__ void
int_double(int163 *dst, const int163 *rhs)
{
	dst[0] = (rhs[0] << 1);
	dst[1] = (rhs[1] << 1);
	dst[2] = (rhs[2] << 1);
	dst[3] = (rhs[3] << 1);
	dst[4] = (rhs[4] << 1);
	dst[5] = (rhs[5] << 1);

	// if dst too bigger than N, minus N. Maybe not fully reduced
	unsigned int t = ((dst[5] & HIGH4BITS) == 0) * 0xFFFFFFFF;
	dst[5] -= (t & N5); dst[5] --;
	dst[4] |= _29BIT;	/*dst[4] -= (t & N4);*/ dst[4] --;	// note N3 == N4 == 0
	dst[3] |= _29BIT;	/*dst[3] -= (t & N3);*/ dst[3] --;
	dst[2] |= _29BIT;	dst[2] -= (t & N2);     dst[2] --;
	dst[1] |= _29BIT;	dst[1] -= (t & N1);		dst[1] --;
	dst[0] |= _29BIT;	dst[4] -= (t & N0);

	// treat possible carry on
	t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;

}

__device__ void
int_sefldouble(int163 *dst)
{
	dst[0] <<= 1;
	dst[1] <<= 1;
	dst[2] <<= 1;
	dst[3] <<= 1;
	dst[4] <<= 1;
	dst[5] <<= 1;

	// treat possible carry
	unsigned int t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;
}

__device__ void
int_selfdoublemod(int163 *dst)
{
	dst[0] <<= 1;
	dst[1] <<= 1;
	dst[2] <<= 1;
	dst[3] <<= 1;
	dst[4] <<= 1;
	dst[5] <<= 1;

	// treat possible carry
	unsigned int t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;

	// if dst too bigger than N, minus N. Maybe not fully reduced
	t = ((dst[5] & HIGH4BITS) == 0) * 0xFFFFFFFF;
	dst[5] -= (t & N5); dst[5] --;
	dst[4] |= _29BIT;	/*dst[4] -= (t & N4);*/ dst[4] --;	// note N3 == N4 == 0
	dst[3] |= _29BIT;	/*dst[3] -= (t & N3);*/ dst[3] --;
	dst[2] |= _29BIT;	dst[2] -= (t & N2);     dst[2] --;
	dst[1] |= _29BIT;	dst[1] -= (t & N1);		dst[1] --;
	dst[0] |= _29BIT;	dst[4] -= (t & N0);

	// treat possible carry on
	t = (dst[0] & HIGH4BITS) >> 28;
	dst[1] += t;
	t = (dst[1] & HIGH4BITS) >> 28;
	dst[2] += t;
	t = (dst[2] & HIGH4BITS) >> 28;
	dst[3] += t;
	t = (dst[3] & HIGH4BITS) >> 28;
	dst[4] += t;
	t = (dst[4] & HIGH4BITS) >> 28;
	dst[5] += t;

}



/************************************************************************/
/*       Polynomial(Finite Field Elements) Operation                    */
/************************************************************************/


/*
hanmming weight of a polynomial
*/
__device__ unsigned int
poly_hw(const poly163 *src)
{
	return __popc(src[0]) + __popc(src[1]) + __popc(src[2]) + __popc(src[3]) + __popc(src[4]) + __popc(src[5]);
}

/*
multiply lhs by 0,1. rhs should be 0 or 1
the result is saved to dst
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

/*
multiply dst by 0,1. rhs should be 0 or 1
the result is saved to dst
*/
__device__ void
poly_selfmulbit(poly163 *dst, unsigned int rhs)
{
	unsigned int t = rhs * 0xFFFFFFFF;

	dst[0] &= t;
	dst[1] &= t;
	dst[2] &= t;
	dst[3] &= t;
	dst[4] &= t;
	dst[5] &= t;
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

/*
if poly src is ONE, return 1, else return 0
*/
__device__ unsigned int 
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

	// 	return 8*i + 31 - __brev(src[i]);
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

__device__ void
poly_addword(poly163 *dst, const poly163 *lhs, unsigned int rhs)
{
	dst[0] = lhs[0] ^ rhs;
}

__device__ void
poly_selfaddword(poly163 *dst, unsigned int rhs)
{
	dst[0] ^= rhs;
}


// add the irreducible polynomial to dst
__device__ void
poly_selfaddire(poly163 *dst)
{
	dst[0] ^= N0;
	dst[5] ^= N5;

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
	unsigned int tp[12];

	// save polynomial src to be reduced in tp
#pragma unroll
	for (int i=0; i<12; ++i)
	{
		tp[i] = src[i];
	}

	// reduce tp
#pragma unroll
	for (int i=10; i>=6; --i)
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
poly_mulmod(poly163 *dst, const poly163 *lhs, const poly163 *rhs)
{
	poly326 r[12];
	poly_mul6(r, lhs, rhs);
	poly_mod(dst, r);
}

// dst *= lhs
__device__ void 
poly_selfmulmod(poly163 *dst, const poly163 *lhs)
{
	poly163 t[6];
	t[0] = dst[0]; t[1] = dst[1]; t[2] = dst[2]; t[3] = dst[3]; t[4] = dst[4]; t[5] = dst[5]; 

	poly_mulmod(dst, t, lhs);
}

// t = t*t;
__device__ void 
poly_selfsqrmod(poly163 *dst)
{
	poly163 t[6];
	t[0] = dst[0]; t[1] = dst[1]; t[2] = dst[2]; t[3] = dst[3]; t[4] = dst[4]; t[5] = dst[5]; 

	poly_mulmod(dst, t, t);
}

// __device__ void
// tool_sqrword(unsigned int *dst, unsigned int src)
// {
// 	unsigned int sqrtbl[256];
// 	
// 	unsigned int hi, lo;
// 	
// 	dst[0] = sqrtbl[src&255];
// 	dst[0] = dst[0] | (sqrtbl[(src>>8)&255]<<16);
// 	dst[1] = sqrtbl[(src>>16)&255];
// 	dst[1] = dst[1] | (sqrtbl[(src>>24)&255]<<16);
// }
// 
// __device__ void
// poly_sqr(poly326 *dst, const poly163 *src)
// {
// 	unsigned int sqrtbl[256];
// 
// }
// 
// __device__ void
// poly_sqrmod(poly163 *dst, const poly163 *src)
// {
// }


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
// __device__ void
// poly_modinv(poly163 *dst, const poly163 *src)
// {
// 	// initialization
// 	poly163 u[6];		// src is saved in u
// 	poly163 v[6] = {0x8, 0x0, 0x0, 0x0, 0x0, 0x107};	// irreducible polynomial z^163+z^8+z^2+z+1
// 
// 	poly163 g1[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x1};
// 	poly163 g2[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
// 	
// 	#pragma unroll
// 	for (int i=0; i<6; ++i)
// 	{
// 		u[i] = src[i];
// 	}
// 
// 	// main loop
// 	while (!poly_isone(u) && !poly_isone(v))
// 	{
// 		while(poly_tryselfdivx(u))
// 		{
// 			if((g1[0] & 0x1) == 0)
// 				poly_selfdivx(g1);
// 			else
// 			{
// 				poly_selfaddire(g1);
// 				poly_selfdivx(g1);
// 			}
// 		}
// 		while (poly_tryselfdivx(v))
// 		{
// 			if((g2[0] & 0x1) == 0)
// 				poly_selfdivx(g2);
// 			else
// 			{
// 				poly_selfaddire(g2);
// 				poly_selfdivx(g2);
// 			}
// 		}
// 	}
// }

// dst * src = 1 mod z^163+z^8+z^2+z+1
/*
use fermat little theorem
*/
__device__ void
poly_modinv(poly163 *dst, const poly163 *src)
{

}









/************************************************************************/
/*          Elliptic Curve Points Operation                             */
/************************************************************************/
/*
judge whether pt1 == -pt2.
if pt1 == -pt2, return 1, else return 0
*/
#define ISNEGATIVE(pt1, pt2)  \
	(pt1[0]==pt2[0])&(pt1[1]==pt2[1])&(pt1[2]==pt2[2])&(pt1[3]==pt2[3])&(pt1[4]==pt2[4])&(pt1[5]==pt2[5]) & \
	((pt1[6]^pt1[0])==pt2[6])&((pt1[7]^pt1[1])==pt2[7])&((pt1[8]^pt1[2])==pt2[8])&((pt1[9]^pt1[3])==pt2[9])&((pt1[10]^pt1[4])==pt2[10])&((pt1[11]^pt1[5])==pt2[11])

/*
judge whether pt is identity.
if pt is identity, return 1, else return 0
*/
#define ISIDENTITY(pt)\
	(pt[0] == 0) & (pt[1] == 0) & (pt[2] == 0) & (pt[3] == 0) & (pt[4] == 0) & (pt[5] == 0) &\
	(pt[6] == 0) & (pt[7] == 0) & (pt[8] == 0) & (pt[9] == 0) & (pt[10] == 0) & (pt[11] == 0)

/*
judge whether pt1 == pt2.
if pt1 == pt2, return 1, else return 0
*/
#define ISEQUAL(pt1, pt2)\
	(pt1[0] == pt2[0]) & (pt1[1] == pt2[1]) & (pt1[2] == pt2[2]) & (pt1[3] == pt2[3]) & (pt1[4] == pt2[4]) & (pt1[5] == pt2[5]) & \
	(pt1[6] == pt2[6]) & (pt1[7] == pt2[7]) & (pt1[8] == pt2[8]) & (pt1[9] == pt2[9]) & (pt1[10] == pt2[10]) & (pt1[11] == pt2[11])

/*
scalar mulplication, k is restricted to 1 or 0.
the result is whether src or identity point
*/	
__device__ void
pt163_mulbit(pt163 *dst, const pt163 *src, unsigned int k)
{
	unsigned int t = (k == 0) * 0xFFFFFFFF;
	dst[0] = src[0] & t;
	dst[1] = src[1] & t;
	dst[2] = src[2] & t;
	dst[3] = src[3] & t;
	dst[4] = src[4] & t;
	dst[5] = src[5] & t;
	dst[6] = src[6] & t;
	dst[7] = src[7] & t;
	dst[8] = src[8] & t;
	dst[9] = src[9] & t;
	dst[10] = src[10] & t;
	dst[11] = src[11] & t;
}

__device__ void
pt163_add(pt163 *pt3, const pt163 *pt1, const pt163 *pt2)
{
	const unsigned int c1 = (ISIDENTITY(pt1)) | (ISIDENTITY(pt2));	// c1=1 iff pt1 or pt2 is identity
	const unsigned int c2 = (ISNEGATIVE(pt1, pt2));					// c2=1 iff pt1=-pt2
	const unsigned int c3 = (ISEQUAL(pt1, pt2));					// c3=1 iff pt1=pt2
	const unsigned int c4 = 1 - c3;									// c4=1 iff pt1~=pt2

	poly163 t1[6];
	poly163 t2[6];
	poly163 t3[6];

	/*
	lambda is calculated
	*/
	poly_mulbit(t1, pt2, c1 | c2 | c3);
	poly_selfadd(t1, pt1);
	poly_modinv(t2, t1);			// t2 = 1/(x_1+k_3 x_2). t1, t3 are free

	poly_mulmod(t1, pt1, pt1);		// t1 = x1^2
	poly_selfmulbit(t1, c3);		// t1 = k2 x1^2
	poly_mulbit(t3, pt2+6, 1-c3);	// t3 = k1 y2
	poly_selfadd(t1, t3);			// t1 = k1 y2 + k2 x1^2
	poly_selfadd(t1, pt1+6);		// t1 = y1 + k1 y2 + k2 x1^2. t3 is free

	poly_mulmod(t3, t1, t2);		// t3 = lambda. t1, t2 are free

	/*
	x3 is calculated
	*/
	poly_add(pt3, pt1, pt2);					// x3 = x1 + x2
	poly_selfadd(pt3, t3);						// x3 = lambda + x1 + x2
	poly_selfaddword(pt3, 1 - (c1 | c2));		// x3 = lambda + x1 + x2 + l
	poly_mulmod(t1, t3, t3);					// t1 = lambda^2. t2 is free
	poly_selfadd(pt3, t1);						// x3 = lambda^2 + lambda + x1 + x2 + l. t1, t2 are free

	/*
	y3 is calculated
	*/
	poly_mulmod(pt3+6, t3, pt3);	// y3 = lambda x3. t1, t2 are free

	poly_mulmod(t1, t3, pt1);		// t1 = lambda x1. t2 is free.
	poly_selfmulbit(t1, c3);		// t1 = m1 lambda x1. t2 is free.
	poly_selfadd(pt3+6, t1);		// y3 = m1 lambda x1 + lambda x3. t1, t2, t3 are free.

	poly_mulmod(t1, pt1, pt1);			// t1 = x1^2. t2, t3 are free.
	poly_selfmulbit(t1, c3 & (1-c1));			// t1 = m2 x1^2. t2 is free
	poly_selfadd(pt3+6, t1);		// y3 = m1 lambda x1 + m2 x1^2 + lambda x3. t1, t2, t3 are free

	poly_mulbit(t1, pt1+6, c1|c4);	// t1 = m3 y1. m3 = 1 iff pt1 or pt2 is identity or pt1!=pt2
	poly_selfadd(pt3+6, t1);		// y3 = m1 lambda x1 + m2 x1^2 + lambda x3 + m3 y1. t1, t2, t3 are free

	poly_mulbit(t1, pt2+6, c1);		// t1 = m4 y2. m4 = 1 iff pt1 or pt2 is identity
	poly_selfadd(pt3+6, t1);		// y3 = m1 lambda x1 + m2 x1^2 + lambda x3 + m3 y1 + m4 y2. t1, t2, t3 are free

	poly_mulbit(t1, pt3, 1-(c1|c2));// t1 = m5 x3. m5 = 0 iff  pt1 or pt2 is identity or pt1=-pt2
	poly_selfadd(pt3+6, t1);		// y3 = m1 lambda x1 + m2 x1^2 + lambda x3 + m3 y1 + m4 y2 + m5 x3
}



/************************************************************************/
/*               Attack Function                                        */
/************************************************************************/
#define MAX_ITERATION_NUM 0xfffff
#define CACHE_SIZE 2400				// cache size that are used to hold distinguished points and coefficients
#define THREADS_PER_BLOCK 512		// 16 warps in one block
#define BLOCKS_PER_GRID 60			// 60 blocks, and there are 30 SMs in GTX 280

/*
PT saved the initial points 
A, B saved the coefficient, and will be updated during iteration
dp_set is used to saved distinguished points and their coefficients, each unit occupy 24 words
*/
__global__ void iteration(pt163 *PT, int163 *A, int163 *B, unsigned int *dp_set)
{
	__shared__ unsigned int counter;
	counter = 0;
	__threadfence_block();

	// initialize start point
	pt163 pt[12];
	pt[0] = PT[0];

	// initialze a_i,b_i
	int163 a[6];
	a[0] = A[0];
	int163 b[6];
	b[0] = B[0];


	while(counter > 0)
	{
		if (0)	// if a point is distinguished point, save it to dp_set
		{
			unsigned int idx = atomicAdd(&counter, 1);

			if (idx > 0)
			{
				
			}
		}
	}
}

__global__ void testAdd(poly163 *des, const poly163 *lh, const poly163 *rh)
{
	poly_add(des, lh, rh);
}

int main()
{
	printf("Hello, this is the testbed for xulei's cuda ecdlp solver\n\n");

	poly163 *lhs = (poly163 *)malloc(SIZE * 200 *sizeof(unsigned int));
	poly163 *rhs = (poly163 *)malloc(SIZE * 200 *sizeof(unsigned int));
	poly163 *des = (poly163 *)malloc(SIZE * 200 *sizeof(unsigned int));

	poly163 *lhs_d;
	cudaMalloc((void**)&lhs_d, SIZE * 200 *sizeof(unsigned int));
	cudaMemcpy(lhs_d, lhs, SIZE * 200 *sizeof(unsigned int), cudaMemcpyHostToDevice);

	poly163 *rhs_d;
	cudaMalloc((void**)&rhs_d, SIZE * 200 *sizeof(unsigned int));
	cudaMemcpy(rhs_d, rhs, SIZE * 200 *sizeof(unsigned int), cudaMemcpyHostToDevice);

	poly163 *des_d;
	cudaMalloc((void**)&des_d, SIZE * 200 *sizeof(unsigned int));


	return 1;
}
