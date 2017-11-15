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
#include <string>
#include <bitset>
#include <intrin.h>

#include "log.ipp"

namespace tools
{
	namespace allocator
	{
		//Optionally increment the provided i such that it becomes a multiple of 64.
		constexpr int multiple_64(const int i)
		{
			return ((i & 0b111111) == 0) ? i : i + (0b1000000 - (i & 0b111111));
		}
		constexpr int multiple_32(const int i)
		{
			return ((i & 0b11111) == 0) ? i : i + (0b100000 - (i & 0b11111));
		}
		constexpr int multiple_16(const int i)
		{
			return ((i & 0b1111) == 0) ? i : i + (0b10000 - (i & 0b1111));
		}
		constexpr int multiple_N(const int i, const int N)
		{
			if (N == 16) return multiple_16(i);
			if (N == 32) return multiple_32(i);
			if (N == 64) return multiple_64(i);
			return multiple_64(i);
		}
		template <typename T, int ALIGN = 64>
		class Allocator_AVX512
		{
		public:

			typedef T value_type;
			typedef value_type* pointer;
			typedef const value_type* const_pointer;
			typedef value_type& reference;
			typedef const value_type& const_reference;
			typedef std::size_t size_type;
			typedef std::ptrdiff_t difference_type;

			template<typename U>
			struct rebind
			{
				typedef Allocator_AVX512<U> other;
			};

			inline explicit Allocator_AVX512() {}
			inline ~Allocator_AVX512() {}
			inline explicit Allocator_AVX512(Allocator_AVX512 const&) {}
			template<typename U>
			inline Allocator_AVX512(Allocator_AVX512<U> const&) {}

			inline pointer address(reference r) { return &r; }
			inline const_pointer address(const_reference r) { return &r; }

			pointer allocate(size_type n, const void *hint = 0)
			{
				const auto n_bytes = multiple_N(static_cast<int>(n) * sizeof(T), ALIGN);
				void * ptr = _mm_malloc(n_bytes, ALIGN);
				
				#if _DEBUG
				#pragma warning( push )
				#pragma warning( disable : 810)
				unsigned long ptr_value = reinterpret_cast<unsigned long>(ptr);
				if ((ptr_value & 0b111111) != 0) log_ERROR("allocator:allocate: ptr ", std::bitset<64>(ptr_value), " is not properly aligned.\n");
				#pragma warning( pop ) 
				#endif

				return reinterpret_cast<pointer>(ptr);
			}
			void deallocate(pointer p, size_type n)
			{
				_mm_free(p);
			}
			inline size_type max_size() const
			{
				return std::numeric_limits<size_type>::max() / sizeof(T);
			}
			inline void construct(pointer p, const T& t) {
				new(p) T(t); 
			}
			inline void destroy(pointer p) { p->~T(); }

			inline bool operator==(Allocator_AVX512 const&) { return true; }
			inline bool operator!=(Allocator_AVX512 const& a) { return !operator==(a); }
		};
	}
}