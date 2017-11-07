// C++ port of Nupic HTM with the aim of being lite and fast
//
// Copyright (c) 2017 Henk-Jan Lebbink
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero Public License version 3 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero Public License for more details.
//
// You should have received a copy of the GNU Affero Public License
// along with this program.  If not, see http://www.gnu.org/licenses.

#pragma once
#include <iostream>		// std::cout
#include <intrin.h>

#include "assert.ipp"

namespace tools
{
	namespace random
	{
		using namespace ::tools::assert;

		namespace priv
		{
			static __m128i current_random_number4 = _mm_setr_epi32(0xACE1ACE1, 0xACE1ACE2, 0xACE1ACE3, 0xACE1ACE4);
			static unsigned int current_random_number = 0xACE1ACE1;

			// not thread safe!
			static const bool USE_LFSR_RANDOM_GENERATOR = true; 

			// Galois LFSR 32; see http://en.wikipedia.org/wiki/Linear_feedback_shift_register
			inline unsigned int lfsr32_galois_ref(const unsigned int i)
			{
				//return (i >> 1) ^ (-(signed int)(i & 1u) & 0xD0000001u); // The cast to signed int is needed to prevent the warning that unairy minus on unsigned ints yields unsigned ints.
				unsigned int j = i;
				unsigned int lsb = static_cast<unsigned int>(-static_cast<signed int>(j & 1));   /* Get LSB (i.e., the output bit). */
				j >>= 1;                /* Shift register */
				j ^= lsb & 0xD0000001u;

				#if _DEBUG
				if (false) std::cout << "log_INFO:tools:lfsr32_galois: in " << std::hex << i << "; out " << j << std::endl;
				#endif
				return j;
			}

			inline __m128i lfsr32_galois_sse(const __m128i i)
			{
				const __m128i i1 = _mm_and_si128(i, _mm_set1_epi32(1u));
				const __m128i i2 = _mm_sign_epi32(i1, _mm_set1_epi32(-1));
				const __m128i i3 = _mm_and_si128(i2, _mm_set1_epi32(0xD0000001u));
				const __m128i i4 = _mm_srli_epi32(i, 1);
				return _mm_xor_si128(i4, i3);

				/*
				movdqa      xmm2,xmmword ptr [__xmm@00000001000000010000000100000001 (013F7C3540h)]
				movdqa      xmm3,xmmword ptr [__xmm@ffffffffffffffffffffffffffffffff (013F7C3560h)]
				movdqa      xmm4,xmmword ptr [__xmm@d0000001d0000001d0000001d0000001 (013F7C3550h)]

				//  (V)PSIGNB/(V)PSIGNW/(V)PSIGND negates each data element of the destination operand (the first operand)
				//  if the signed integer value of the corresponding data element in the source operand (the second operand)
				//  is less than zero. If the signed integer value of a data element in the source operand is positive, the
				//  corresponding data element in the destination operand is unchanged. If a data element in the source
				//  operand is zero, the corresponding data element in the destination operand is set to zero.

				movdqa      xmm1,xmm2		// xmm1 := xmm2
				pand        xmm1,xmm6		// xmm1 := xmm1&xmm6
				psrld       xmm6,1		// xmm6 := xmm6>>1
				psignd      xmm1,xmm3		// xmm1 := -xmm1
				movdqa      xmm0,xmm6		// xmm0 := xmm6
				pand        xmm1,xmm4		// xmm1 := xmm1 & xmm4
				movdqa      xmm6,xmm1		// xmm6 := xmm1
				pxor        xmm6,xmm0		// xmm6 = xmm6 ^xmm0
				*/
			}
		
			// rescale between MIN and MAX (inclusive)
			template<unsigned int MIN, unsigned int MAX>
			constexpr inline unsigned int rescale_incl(const unsigned int number)
			{
				return (number % (unsigned int)(MAX + 1 - MIN)) + MIN;
			}

			// rescale between MIN and max (inclusive)
			template<unsigned int MIN>
			constexpr inline unsigned int rescale_incl(const unsigned int number, const unsigned int max)
			{
				return (number % (unsigned int)(max + 1 - MIN)) + MIN;
			}

			// rescale between min and max (inclusive)
			inline unsigned int rescale_incl(const unsigned int number, const unsigned int min, const unsigned int max)
			{
				assert_msg((max + 1 - min) > 0, "rescale ERROR: (max + 1 - min) > 0: min = ", min, "; max = ", max);

				const unsigned int result = (number % (unsigned int)(max + 1 - min)) + min;
				assert_msg(result >= min, "rescale ERROR : number = ", number, "; min = ", min, "; result = ", result);
				assert_msg(result <= max, "rescale ERROR : number = ", number, "; max = ", max, "; result = ", result);
				return result;
				//45 8B 03             mov         r8d,dword ptr [number]  
				//33 D2                xor         edx, edx ; clear edx prior to div
				//41 8B C0             mov         eax, r8d
				//41 F7 F2             div         eax, r10d ;  r10d = max + 1 - min; afterwards: EDX=remainder, EAX=quotient;
			}
		}

		inline unsigned int next_rand(const unsigned int i)
		{
			return random::priv::lfsr32_galois_ref(i);
		}

		inline int rand_int32(const int min, const int max, unsigned int& random_number)
		{
			if (random::priv::USE_LFSR_RANDOM_GENERATOR)
			{
				const bool fortran_compatible = false;
				if (fortran_compatible)
				{
					const int b = max - min + 1;
					int m = random_number % b;
					m = (m < 0) ? m + b : m;
					const int result = min + m;

					random_number = next_rand(random_number);
					return result;
				}
				else
				{
					const int result = min + (random_number % (max - min + 1));
					random_number = next_rand(random_number);
					return result;
				}
			}
			else
			{
				return min + (rand() % (max - min + 1));
			}
		}

		inline int rand_int32(const int min, const int max)
		{
			return rand_int32(min, max, random::priv::current_random_number);
		}

		inline unsigned int rand_int32(const unsigned int min, const unsigned int max, unsigned int& random_number)
		{
			const unsigned int result = min + (random_number % (max - min + 1));
			random_number = next_rand(random_number);
			return result;
		}

		inline unsigned int rand_int32(const unsigned int min, const unsigned int max)
		{
			return rand_int32(min, max, random::priv::current_random_number);
		}

		inline unsigned int rand_int32(const unsigned int max)
		{
			return rand_int32(0U, max, priv::current_random_number);
		}

		// return a random float between 0 and max (inclusive); uses and updates the provided random number
		inline float rand_float(const float max, unsigned int& random_number)
		{
			const float result = static_cast<float>((random_number * max) / 0xFFFFFFFFUL);
			//if (true) std::cout << "log_INFO:tools:rand_float: current_random_number = " << static_cast<int>(current_random_number) << "; result = " << std::fixed << std::setw(13) << std::setprecision(10) << result << std::endl;
			random_number = next_rand(random_number);
			return result;
		}

		// return a random float between 0 and max
		inline float rand_float(const float max)
		{
			if (priv::USE_LFSR_RANDOM_GENERATOR)
			{
				return rand_float(max, priv::current_random_number);
			}
			else
			{
				const double r = static_cast<double>(rand());
				return static_cast<float>((r * max) / static_cast<double>(RAND_MAX));
			}
		}

		// return a random float between 0 and 1
		inline float rand_float()
		{
			return rand_float(1.0);
		}

		inline unsigned int rdrand32()
		{
			unsigned int i;
			if (false)
			{
				//_rdrand32_step(&i);
			}
			else
			{
				i = (static_cast<unsigned int>(rand()) << 0) |
					(static_cast<unsigned int>(rand()) << 15) |
					(static_cast<unsigned int>(0x4 & rand()) << 30);
			}
			return i;
		}
	}
}