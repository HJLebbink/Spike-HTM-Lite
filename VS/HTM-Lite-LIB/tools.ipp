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
#include <array>
#include <vector>
#include <tuple>
#include <iostream>		// std::cout
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setw

#include "..\Spike-Tools-Lib\assert.ipp"

#include "parameters.ipp"

namespace htm
{
	namespace tools
	{
		//Return the number of blocks needed with 16 elements in each block; tail block included.
		constexpr int n_blocks_16(const int i)
		{
			return (i >> 4) + ((i & 0b1111) > 0); // add an optional extra block for the tail
		}
		//Return the number of blocks needed with 16 elements in each block; tail block included.
		constexpr int n_blocks_32(const int i)
		{
			return (i >> 5) + ((i & 0b11111) > 0); // add an optional extra block for the tail
		}
		//Return the number of blocks needed with 64 elements in each block; tail block included.
		constexpr int n_blocks_64(const int i)
		{
			return (i >> 6) + ((i & 0b111111) > 0); // add an optional extra block for the tail
		}

		//Optionally increment the provided i such that it becomes a multiple of 16.
		constexpr int multiple_16(const int i)
		{
			return ((i & 0b1111) == 0) ? i : i + (0b10000 - (i & 0b1111));
		}
		//Optionally increment the provided i such that it becomes a multiple of 64.
		constexpr int multiple_64(const int i)
		{
			return ((i & 0b111111) == 0) ? i : i + (0b1000000 - (i & 0b111111));
		}

		constexpr int get_global_cell_id(const int delay_and_cell_id)
		{
			return delay_and_cell_id & 0x1FFFFFFF;
		}
		constexpr int get_delay(const int delay_and_cell_id)
		{
			return static_cast<int>(static_cast<unsigned int>(delay_and_cell_id) >> (32 - 3));
		}
		constexpr int create_delay_and_cell_id(const int global_cell_id, const int delay)
		{
			return (static_cast<unsigned int>(delay) << (32 - 3)) | global_cell_id;
		}



		template <int N_COLUMNS>
		int get_random_column()
		{
			return rand_int32(0, N_COLUMNS - 1);
		}

		template <int N_CELLS_PC>
		int8_t get_random_cell()
		{
			return rand_int32(0, N_CELLS_PC - 1);
		}

		template <int N_CELLS>
		int get_random_global_cell()
		{
			return rand_int32(0, N_CELLS - 1);
		}


		template <int N_CELLS>
		int get_random_global_cell_include(
			const std::vector<int>& global_cell_ids)
		{
			int n_global_cell_ids = static_cast<int>(global_cell_ids.size());
			if (n_global_cell_ids == 0)
			{
				if (false) ::tools::log::log_WARNING("tools::get_random_global_cell_include: could not find random cell. include_cells is empty");
				return get_random_global_cell<N_CELLS>();
			}
			return global_cell_ids[tools::rand_int32(0, n_global_cell_ids-1)];
		}

		template <int N_COLUMNS, int N_BITS_CELL>
		constexpr std::tuple<int8_t, int> global_2_local_cell(
			const int global_cell_idx)
		{
			const int8_t cell = global_cell_idx & ((1 << N_BITS_CELL) - 1);
			const int column = global_cell_idx >> N_BITS_CELL;

			#if _DEBUG
			const auto N_CELLS_PC = 1 << N_BITS_CELL;
			if ((cell < 0) || (cell >= N_CELLS_PC)) ::tools::log::log_ERROR("tools:global_2_local_cell: invalid cell = ", static_cast<int>(cell));
			if ((column < 0) || (column >= N_COLUMNS)) ::tools::log::log_ERROR("tools:global_2_local_cell: invalid column = ", column);
			#endif
			return std::tuple<int8_t, int>(cell, column);
		}


		template <int N_BITS_CELL>
		constexpr int local_2_global_cell(
			const int8_t cell,
			const int column)
		{
			#if _DEBUG
			if (cell >= (1 << N_BITS_CELL)) ::tools::log::log_ERROR("tools::local_2_global_cell: cell ", static_cast<int>(cell), " is too large.");
			#endif
			return (column << N_BITS_CELL) | cell;
		}

		template <int N_CELLS_PC>
		void test_global_cell()
		{
			for (auto column_i = 0; column_i < 10; ++column_i)
			{
				for (auto cell_i = 0; cell_i < N_CELLS_PC; ++cell_i)
				{
					auto global_cell_id = local_2_global_cell(cell_i, column_i);
					std::cout << "log_INFO: column " << column_i << "; cell " << cell_i << "; global cell id " << global_cell_id << "." << std::endl;
				}
			}
		}

		template <class T>
		void select_random(const int size, std::vector<T>& inout)
		{
			const auto size_inout = inout.size();
			for (auto i = 0; i < size; ++i)
			{
				const int r = ::tools::rand_int32(i, size_inout);
				const T tmp = inout[i];
				inout[i] = inout[r];
				inout[r] = tmp;
			}
		}


		void add(std::vector<int>& a, const std::vector<int>& b)
		{
			for (int i = 0; i < static_cast<int>(b.size()); ++i) a[i] += b[i];
		}
		void clear(std::vector<int>& a)
		{
			for (int i = 0; i < static_cast<int>(a.size()); ++i) a[i] = 0;
		}

	}
}