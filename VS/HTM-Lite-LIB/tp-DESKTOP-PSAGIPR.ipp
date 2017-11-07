#pragma once

#include "tbb/tbb.h"

#include <algorithm>	// std::min
#include <limits>		// std::numeric_limits
#include <tuple>
#include <array>
#include <vector>
#include <map>
#include <set>

//#include "mmintrin.h"  // mmx
#include "emmintrin.h"  // sse
#include "pmmintrin.h"  // sse3
#include "tmmintrin.h"  // ssse3
#include "smmintrin.h"  // sse4.1
#include "nmmintrin.h"  // sse4.2
#include "immintrin.h"  // avx, avx2, avx512, FP16C, KNCNI, FMA
#include "zmmintrin.h"
#include "intrin.h"

#include "constants.ipp"
#include "tools.ipp"
#include "print.ipp"

#include "../../Tools/Tools-Lib/log.ipp"
#include "../../Tools/Tools-Lib/assert.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Temporal Memory/Pooler
	namespace tp
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;

		//HTM Temporal Pooler private methods
		namespace priv
		{
			//Get the number of active synapses of the provided segment in the provided column
			namespace get_number_active_DD_synapses
			{
				inline int get_number_active_DD_synapses_ref(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold)
				{
					const auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[segment];
					const auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];

					int n_active_synapses = 0;
					for (auto synapse_i = 0; synapse_i < column.dd_synapse_count[segment]; ++synapse_i)
					{
						if (dd_synapse_permanence_segment[synapse_i] >= permanence_threshold)
						{
							const auto global_cell_id = dd_synapse_origin_segment[synapse_i];
							const bool b = active_cells_all[global_cell_id];
							if (b) n_active_synapses++;
						}
					}
					if (false) log_INFO_DEBUG("TP:get_number_active_DD_synapses: column ", column.id, "; segment ", segment, " has ", column.dd_synapse_count[segment], " synapses and ", n_active_synapses, " active synapses.");
					return n_active_synapses;
				}

				inline int get_number_active_DD_synapses_avx(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold)
				{
					const __m256 * permanence = reinterpret_cast<const __m256 *>(column.dd_synapse_permanence[segment].data());
					const __m256 permanence_threshold_avx = _mm256_set1_ps(permanence_threshold);
					const __m256i * origin_ptr = reinterpret_cast<const __m256i *>(column.dd_synapse_origin[segment].data());

					int n_active_synapses = 0;
					const int n_synapses = column.dd_synapse_count[segment];
					const int n_blocks = (n_synapses >> 3) + (((n_synapses & 0b111) > 0) ? 1 : 0);
					for (int block = 0; block < n_blocks; ++block)
					{
						const int cmp_result1 = _mm256_movemask_ps(_mm256_cmp_ps(permanence[block], permanence_threshold_avx, _CMP_GE_OS));
						const __m256i origin1 = origin_ptr[block];

						if (((cmp_result1 & 0b00000001) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 0)]) n_active_synapses++;
						if (((cmp_result1 & 0b00000010) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 1)]) n_active_synapses++;
						if (((cmp_result1 & 0b00000100) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 2)]) n_active_synapses++;
						if (((cmp_result1 & 0b00001000) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 3)]) n_active_synapses++;
						if (((cmp_result1 & 0b00010000) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 4)]) n_active_synapses++;
						if (((cmp_result1 & 0b00100000) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 5)]) n_active_synapses++;
						if (((cmp_result1 & 0b01000000) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 6)]) n_active_synapses++;
						if (((cmp_result1 & 0b10000000) != 0) && active_cells_all[_mm256_extract_epi32(origin1, 7)]) n_active_synapses++;

					}
					if (false) log_INFO_DEBUG("TP:get_number_active_DD_synapses_avx: column ", column.id, "; segment ", segment, " has ", column.dd_synapse_count[segment], " synapses and ", n_active_synapses, " active synapses.");

					#if _DEBUG
					const int n_active_synapses_ref = get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold);
					if (n_active_synapses != n_active_synapses_ref) log_ERROR("TP:get_number_active_DD_synapses_avx:: UNEQUAL number active synapses ", n_active_synapses, " != ", n_active_synapses_ref);
					#endif

					return n_active_synapses;
				}

				inline int get_number_active_DD_synapses_avx512(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold)
				{
					const __m512 * permanence = reinterpret_cast<const __m512 *>(column.dd_synapse_permanence[segment].data());
					const __m512i * origin = reinterpret_cast<const __m512i *>(column.dd_synapse_origin[segment].data());
					const char * active_cells_ptr = active_cells_all.data();

					const __m512 permanence_threshold_avx = _mm512_set1_ps(permanence_threshold);
					const __m512i zeros = _mm512_setzero_epi32();
					const __m512i ones = _mm512_set1_epi32(1);

					__m512i n_active_synapses = zeros;

					const int n_synapses = column.dd_synapse_count[segment];
					const int n_blocks = n_synapses >> 4;

					if (n_blocks > 1) log_INFO("TP:get_number_active_DD_synapses_avx512: n_blocks = ", n_blocks);

					// approx. 40% of time is spend in this loop
					#pragma ivdep // ICC thinks there is a waw dependency but there isn't one.
					for (int block = 0; block < n_blocks; ++block) 
					{
						const __mmask16 cmp_result = _mm512_cmp_ps_mask(permanence[block], permanence_threshold_avx, _CMP_GE_OS);
						// not sure whether the branch misprediction is worse than calling a superflous _mm512_add_epi32 and _mm512_and_epi32
						// which are hidden behind the latency of the gather.

						//if (cmp_result != 0) 
						{
							// subtle bug: may read upto 3 bytes beyond the end of sensor_activity
							const __m512i sensors = _mm512_mask_i32gather_epi32(zeros, cmp_result, origin[block], active_cells_ptr, 1);
							// we read 4 bytes but we are only interested in the lowest position bit
							n_active_synapses = _mm512_add_epi32(n_active_synapses, _mm512_and_epi32(sensors, ones));
						}
					}
					int n_active_synapses_int = _mm512_reduce_add_epi32(n_active_synapses);

					const int tail = n_synapses & 0b1111;
					if (tail > 0)
					{
						const auto permanence_tail = reinterpret_cast<const float *>(column.dd_synapse_permanence[segment].data());
						const auto origin_tail = reinterpret_cast<const int *>(column.dd_synapse_origin[segment].data());
						for (auto synapse_i = (n_blocks << 4) + 1; synapse_i < n_synapses; ++synapse_i)
						{
							if (permanence_tail[synapse_i] >= permanence_threshold)
							{
								if (active_cells_all[origin_tail[synapse_i]]) n_active_synapses_int++;
							}
						}
					}

					if (false) log_INFO_DEBUG("TP:get_number_active_DD_synapses_avx: column ", column.id, "; segment ", segment, " has ", column.dd_synapse_count[segment], " synapses and ", n_active_synapses_int, " active synapses.");

					#if _DEBUG
					const int n_active_synapses_ref = get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold);
					if (n_active_synapses_int != n_active_synapses_ref) log_ERROR("TP:get_number_active_DD_synapses_avx:: UNEQUAL number active synapses ", n_active_synapses_int, " != ", n_active_synapses_ref);
					#endif

					return n_active_synapses_int;
				}

				inline std::tuple<int, int> get_number_active_DD_synapses_2x_avx512(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold1,
					const float permanence_threshold2)
				{
					const __m512 * permanence = reinterpret_cast<const __m512 *>(column.dd_synapse_permanence[segment].data());
					const __m512i * origin = reinterpret_cast<const __m512i *>(column.dd_synapse_origin[segment].data());
					const char * active_cells_ptr = active_cells_all.data();

					const __m512 permanence_threshold1_avx = _mm512_set1_ps(permanence_threshold1);
					const __m512 permanence_threshold2_avx = _mm512_set1_ps(permanence_threshold2);
					const __m512i zeros = _mm512_setzero_epi32();
					const __m512i ones = _mm512_set1_epi32(1);

					__m512i n_active_synapses1 = zeros;
					__m512i n_active_synapses2 = zeros;

					const int n_synapses = column.dd_synapse_count[segment];
					const int n_blocks = n_synapses >> 4;

					for (int block = 0; block < n_blocks; ++block)
					{
						const __m512 p = permanence[block];
						const __mmask16 cmp_result1 = _mm512_cmp_ps_mask(p, permanence_threshold1_avx, _CMP_GE_OS);
						const __mmask16 cmp_result2 = _mm512_cmp_ps_mask(p, permanence_threshold2_avx, _CMP_GE_OS);
						const __mmask16 cmp_result = cmp_result1 | cmp_result2;
						
						// not sure whether the branch misprediction is worse than calling a superflous _mm512_add_epi32 and _mm512_and_epi32
						// which are hidden behind the latency of the gather.
						//if (cmp_result != 0) 
						{
							// subtle bug: may read upto 3 bytes beyond the end of sensor_activity
							const __m512i sensors = _mm512_mask_i32gather_epi32(zeros, cmp_result, origin[block], active_cells_ptr, 1);
							// we read 4 bytes but we are only interested in the lowest position bit

							n_active_synapses1 = _mm512_add_epi32(n_active_synapses1, _mm512_mask_and_epi32(zeros, cmp_result1, sensors, ones));
							n_active_synapses2 = _mm512_add_epi32(n_active_synapses2, _mm512_mask_and_epi32(zeros, cmp_result2, sensors, ones));
						}
					}
					int n_active_synapses1_int = _mm512_reduce_add_epi32(n_active_synapses1);
					int n_active_synapses2_int = _mm512_reduce_add_epi32(n_active_synapses2);

					const int tail = n_synapses & 0b1111;
					if (tail > 0)
					{
						const auto permanence_tail = reinterpret_cast<const float *>(column.dd_synapse_permanence[segment].data());
						const auto origin_tail = reinterpret_cast<const int *>(column.dd_synapse_origin[segment].data());
						for (auto synapse_i = (n_blocks << 4) + 1; synapse_i < n_synapses; ++synapse_i)
						{
							if (permanence_tail[synapse_i] >= permanence_threshold1)
							{
								if (active_cells_all[origin_tail[synapse_i]]) n_active_synapses1_int++;
							}
							if (permanence_tail[synapse_i] >= permanence_threshold2)
							{
								if (active_cells_all[origin_tail[synapse_i]]) n_active_synapses2_int++;
							}
						}
					}

					if (false) log_INFO_DEBUG("TP:get_number_active_DD_synapses_avx: column ", column.id, "; segment ", segment, " has ", column.dd_synapse_count[segment], " synapses and ", n_active_synapses1_int, " active synapses.");

					#if _DEBUG
					const int n_active_synapses1_ref = get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold1);
					if (n_active_synapses1_int != n_active_synapses1_ref) log_ERROR("TP:get_number_active_DD_synapses_avx:: UNEQUAL number active synapses ", n_active_synapses1_int, " != ", n_active_synapses1_ref);
					const int n_active_synapses2_ref = get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold2);
					if (n_active_synapses2_int != n_active_synapses2_ref) log_ERROR("TP:get_number_active_DD_synapses_avx:: UNEQUAL number active synapses ", n_active_synapses2_int, " != ", n_active_synapses2_ref);
					#endif

					return std::make_tuple(n_active_synapses1_int, n_active_synapses2_int);
				}

				// assumes BitsetCell is a bit-array instead of a byte-array: consider using bit-array instead of byte-array for BitsetCell, but it still wont fit L3 cache (8MB)
				/*
					L1 Data cache = 32 KB, 64 B / line, 8 - WAY.
					L1 Instruction cache = 32 KB, 64 B / line, 8 - WAY.
					L2 cache = 256 KB, 64 B / line, 4 - WAY
					L3 cache = 8 MB, 64 B / line, 16 - WAY

					L1 Data Cache Latency = 4 cycles for simple access via pointer
					L1 Data Cache Latency = 5 cycles for access with complex address calculation(size_t n, *p; n = p[n]).
					L2 Cache Latency = 12 cycles
					L3 Cache Latency = 42 cycles(core 0) (i7 - 6700 Skylake 4.0 GHz)
					RAM Latency = 42 cycles + 51 ns(i7 - 6700 Skylake)
				*/
				inline int get_number_active_DD_synapses_avx512_bitarray(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold)
				{
					const __m512 * permanence = reinterpret_cast<const __m512 *>(column.dd_synapse_permanence[segment].data());
					const __m512i * origin = reinterpret_cast<const __m512i *>(column.dd_synapse_origin[segment].data());
					const char * active_cells_ptr = active_cells_all.data();

					const __m512 permanence_threshold_avx = _mm512_set1_ps(permanence_threshold);
					const __m512i zeros = _mm512_setzero_epi32();
					const __m512i ones = _mm512_set1_epi32(1);

					__m512i n_active_synapses = zeros;

					const int n_synapses = column.dd_synapse_count[segment];
					const int n_blocks = n_synapses >> 4;
					for (int block = 0; block < n_blocks; ++block)
					{
						const __mmask16 cmp_result = _mm512_cmp_ps_mask(permanence[block], permanence_threshold_avx, _CMP_GE_OS);
						if (cmp_result != 0)
						{
							const __m512i bit_addr = origin[block];
							const __m512i byte_addr = _mm512_srai_epi32(bit_addr, 3);
							const __m512i pos_in_byte = _mm512_and_epi32(bit_addr, _mm512_set1_epi32(0b111));
							
							// subtle bug: may read upto 3 bytes beyond the end of sensor_activity
							const __m512i sensors = _mm512_srav_epi32(_mm512_mask_i32gather_epi32(zeros, cmp_result, byte_addr, active_cells_ptr, 1), pos_in_byte);
							n_active_synapses = _mm512_add_epi32(n_active_synapses, _mm512_and_epi32(sensors, ones));
						}
					}
					int n_active_synapses_int = _mm512_reduce_add_epi32(n_active_synapses);

					const int tail = n_synapses & 0b1111;
					if (tail > 0)
					{
						const auto permanence_tail = reinterpret_cast<const float *>(column.dd_synapse_permanence[segment].data());
						const auto origin_tail = reinterpret_cast<const int *>(column.dd_synapse_origin[segment].data());
						for (auto synapse_i = (n_blocks << 4) + 1; synapse_i < n_synapses; ++synapse_i)
						{
							if (permanence_tail[synapse_i] >= permanence_threshold)
							{
								const int bit_addr = origin_tail[synapse_i];
								const int byte_addr = bit_addr >> 3;
								const int pos_in_byte = bit_addr & 0b111;

								const int sensor = active_cells_all[byte_addr] >> pos_in_byte;
								if ((sensor & 1) == 1) n_active_synapses_int++;
							}
						}
					}

					if (false) log_INFO_DEBUG("TP:get_number_active_DD_synapses_avx: column ", column.id, "; segment ", segment, " has ", column.dd_synapse_count[segment], " synapses and ", n_active_synapses_int, " active synapses.");

					#if _DEBUG
					const int n_active_synapses_ref = get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold);
					if (n_active_synapses_int != n_active_synapses_ref) log_ERROR("TP:get_number_active_DD_synapses_avx:: UNEQUAL number active synapses ", n_active_synapses_int, " != ", n_active_synapses_ref);
					#endif

					return n_active_synapses_int;
				}

				// Get the number of active synapses of the provided segment in the provided column
				inline int d(
					const Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS> & active_cells_all,
					const float permanence_threshold)
				{
					#if _DEBUG
					if (true)
					{
						const auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[segment];
						const auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];

						// number of synapses in use
						const auto n_synapses = column.dd_synapse_count[segment];

						if (false) log_INFO("TP:get_number_active_DD_synapses: A: column ", column.id, "; segment ", segment, "; n_synapses = ", n_synapses, ".");

						if ((n_synapses < 0) || (n_synapses > P::TP_N_DD_SYNAPSES)) // invalid number of synapses
						{
							log_ERROR("TP:get_number_active_DD_synapses: B: column ", column.id, "; segment ", segment, "; n_synapses = ", n_synapses, " which is invalid.");
						}

						for (auto synapse_i = 0; synapse_i < n_synapses; ++synapse_i)
						{
							const int global_cell_id = dd_synapse_origin_segment[synapse_i];
							const float permanence = dd_synapse_permanence_segment[synapse_i];

							if ((global_cell_id < 0) || (global_cell_id >= P::N_CELLS)) // invalid global_cell_id
							{
								log_ERROR("TP:get_number_active_DD_synapses: C: column ", column.id, "; segment ", segment, "; synapse ", synapse_i, "; has invalid origin cell = ", global_cell_id, "; permanence = ", permanence, ".");
							}
							if (permanence == P::TP_DD_PERMANENCE_INVALID)
							{
								log_ERROR("TP:get_number_active_DD_synapses: D: column ", column.id, "; segment ", segment, "; synapse ", synapse_i, "; origin = ", global_cell_id, "; has invalid permanence = ", permanence, ".");
							}
						}
					}
					#endif

					if (ARCH == arch_t::X64) return get_number_active_DD_synapses_ref(column, segment, active_cells_all, permanence_threshold);
					if (ARCH == arch_t::AVX) return get_number_active_DD_synapses_avx(column, segment, active_cells_all, permanence_threshold);
					if (ARCH == arch_t::AVX512) return get_number_active_DD_synapses_avx512(column, segment, active_cells_all, permanence_threshold);
				}
			}

			// return true when the provided column is predicted
			inline bool activate_predicted_column(
				//in
				const BitsetTiny<P::N_CELLS_PC>& prev_predictive_cells,
				const BitsetTiny<P::N_CELLS_PC>& prev_matching_cells,
				const bool active_column,
				//out
				BitsetTiny<P::N_CELLS_PC>& active_cells,
				BitsetTiny<P::N_CELLS_PC>& winner_cells)
			{
				#if _DEBUG
				if (false)
				{
					if (active_column)
					{
						log_INFO<false, false>("TP:activate_correctly_predictive_cells: active_column:", active_column, ", prev_predictive_cells:");
						print::print_bitset(prev_predictive_cells);
					}
				}
				if (false)
				{
					if (active_column)
					{
						log_INFO<false, false>("TP:activate_correctly_predictive_cells: active_column:", active_column, ", prev_matching_cells:");
						print::print_bitset(prev_matching_cells);
					}
				}
				#endif

				bool predicted_column;

				if (active_column)
				{
					predicted_column = prev_predictive_cells.any();
					if (false) log_INFO_DEBUG("TP:activate_correctly_predictive_cells: active column; predicted ", predicted_column);

					tools::copy(active_cells, prev_predictive_cells);
					tools::copy(winner_cells, prev_predictive_cells);
				}
				else
				{
					predicted_column = false;
					active_cells.reset();
					winner_cells.reset();
				}
				return predicted_column;
			}

			inline int8_t least_used_cell(const Column& column)
			{
				std::array<int, P::N_CELLS_PC> counter = { 0 };

				for (auto segment_i = 0; segment_i < column.dd_segment_count; ++segment_i)
				{
					const auto cell_i = column.dd_segment_destination[segment_i];
					counter[cell_i]++;
				}

				int best_cell = 0;
				int selected_count = counter[0];

				for (auto cell_i = 1; cell_i < P::N_CELLS_PC; ++cell_i)
				{
					if (counter[cell_i] < selected_count)
					{
						selected_count = counter[cell_i];
						best_cell = cell_i;
					}
				}
				return static_cast<int8_t>(best_cell);
			}

			inline std::tuple<bool, int, bool, int> best_matching_cell(
				const Column& column,
				const BitsetCell<P::N_CELLS>& active_cells_all)
			{
				//get the best matching cell from the provided column given the active_cells
				//if this code does not find a segment, then it is likely that a new segment
				//will be created
				#if _DEBUG
				if (false)
				{
					log_INFO<false, false>("TP:best_matching_cell: column ", column.id, "; active_cells_all:");
					print::print_bitset(active_cells_all);
				}
				#endif

				int best_cell = -1;
				int best_segment = -1;

				int highest_num_active_synapses = 0;
				for (auto segment_i = 0; segment_i < column.dd_segment_count; ++segment_i)
				{
					const auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment_i];
					int num_active_synapses = 0;

					for (auto synapse_i = 0; synapse_i < column.dd_synapse_count[segment_i]; ++synapse_i)
					{
						const auto global_cell_id = dd_synapse_origin_segment[synapse_i];
						const bool b = active_cells_all[global_cell_id];
						if (b) num_active_synapses++;

						if (false) log_INFO_DEBUG("TP:best_matching_cell: column ", column.id, "; segment_i ", segment_i, "; synapse_i = ", synapse_i, "; is active = ", b, "; global_cell_id = ", global_cell_id);
					}
					if (num_active_synapses > highest_num_active_synapses)
					{
						highest_num_active_synapses = num_active_synapses;
						best_segment = segment_i;
					}
					if (false) log_INFO_DEBUG("TP:best_matching_cell: column ", column.id, "; segment_i ", segment_i, "; num_active_synapses = ", num_active_synapses);
				}

				if (false) if (best_segment != -1) log_INFO_DEBUG("TP:best_matching_cell: column ", column.id, "; best_segment = ", best_segment, "; highest_num_active_synapses = ", highest_num_active_synapses);

				if (best_segment != -1)
				{
					best_cell = column.dd_segment_destination[best_segment];
					if (false) log_INFO_DEBUG("TP:best_matching_cell: column ", column.id, "; found a best segment: best_segment = ", best_segment, "; best_cell = ", best_cell);
				}
				if (best_cell == -1)
				{
					best_cell = least_used_cell(column);
					if (false) log_INFO_DEBUG("TP:best_matching_cell: column ", column.id, "; still no cell found: get the least used cell: best_cell = ", best_cell);
				}
				return std::tuple<bool, int, bool, int>(best_cell != -1, best_cell, best_segment != -1, best_segment);
			}

			inline void select_cell_to_learn_on(
				//in
				Column& column, // column cannot be readonly due to the random number generator
				const int segment,
				const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
				const int select_size,
				//out
				std::vector<int>& selected_cells) // assumes that selected_cells has sufficient capacity (size >= select_size)
			{
				//pick global cell ids from winner_cells that do not already
				//have a pathway to segment in column, if not enough winner cells
				//exists, take non winner cells

				const std::vector<int>& winner_cells = prev_winner_cells_all._data;
				const int n_winner_cells = static_cast<int>(winner_cells.size());

				#if _DEBUG
				if (false)
				{
					log_INFO<false, false>("TP:select_cell_to_learn_on: column ", column.id, "; prev_winner_cells_all: ");
					print::print_int_array(winner_cells, n_winner_cells);
				}
				#endif

				const auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];
				const auto n_synapses = column.dd_synapse_count[segment];

				int selected_cells_count = 0;
				int winner_cell_i = 0;

				while (selected_cells_count < select_size)
				{
					int global_cell_id;
					if (winner_cell_i < n_winner_cells) // get a random cell from winner cells
					{
						const int random_i = htm::tools::random::rand_int32(0, n_winner_cells - 1);
						global_cell_id = winner_cells[random_i];

						if (false) log_INFO_DEBUG("TP:select_cell_to_learn_on: column ", column.id, "; segment ", segment, "; random_i = ", random_i, "; global_cell_id = ", global_cell_id, "; current_random_number = ", static_cast<int>(tools::random::priv::current_random_number), "; min = ", winner_cell_i, "; max = ", n_winner_cells - 1);
						winner_cell_i++;
					}
					else // not enough winner cells, get a random cell
					{
						global_cell_id = htm::tools::random::rand_int32(0, P::N_CELLS - 1, column.random_number);
					}

					bool already_present = false;

					for (auto synapse_i = 0; synapse_i < selected_cells_count; ++synapse_i)
					{// search the already selected cells if it is already present
						if (selected_cells[synapse_i] == global_cell_id)
						{
							already_present = true; // found it
							break;
						}
					}
					if (!already_present)
					{
						for (auto synapse_i = 0; synapse_i < n_synapses; ++synapse_i)
						{// search the existing synapses whether the cell already has a pathway to this segment
							if (dd_synapse_origin_segment[synapse_i] == global_cell_id)
							{
								already_present = true; // found it
								break;
							}
						}
					}
					if (!already_present) // add the cell if not already present
					{
						selected_cells[selected_cells_count] = global_cell_id;
						selected_cells_count++;
					}
				}
			}
			
			//Grow new synapses on the provided segment in the provided column
			namespace grow_DD_synapses
			{
				inline void grow_DD_synapses_v1(
					Column& column,
					const int segment,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all)
				{
					const int new_number_synapses = P::TP_DD_MAX_NEW_SYNAPSE_COUNT;
					assert_msg(segment < column.dd_segment_count, "TP:replace_DD_synapses: segment", segment, " is too large");

					auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[segment];
					auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];
					assert_msg(dd_synapse_permanence_segment.size() == dd_synapse_origin_segment.size(), "TP:replace_DD_synapses: bug A.");

					const int old_capacity = dd_synapse_permanence_segment.size();
					const int old_size = column.dd_synapse_count[segment];
					assert_msg(old_size <= old_capacity, "TP:replace_DD_synapses: bug B.");

					//1] assume we are not going to replace any synapses, thus we will be adding 
					// new_number_synapses synapses. Allocate space for these synapses if needed

					if (old_size + new_number_synapses >= old_capacity)
					{	// not enough capacity
						dd_synapse_permanence_segment.resize(old_size + new_number_synapses);
						dd_synapse_origin_segment.resize(old_size + new_number_synapses);
					}

					//2] find indices that will be overwritten with new values
					std::vector<int> indices_to_update = std::vector<int>(new_number_synapses, -1);

					int n_synapses_kept = 0;

					//3] count the number of synapses that can be replaced
					int n_synapses_replaced = 0;
					for (auto synapse_i = 0; synapse_i < old_size; ++synapse_i)
					{
						//if (true) // strangely enough it works better to replace existing synapses only
						//if (dd_synapse_permanence_segment[synapse_i] < P::TP_DD_PERMANENCE_THRESHOLD) // replace synapses that are not connected: works badly
						//if (dd_synapse_permanence_segment[synapse_i] < 0.1f) // replace synapses that have been punished somewhat: works reasonable well
						if (dd_synapse_permanence_segment[synapse_i] < 0.15f) // replace synapses that have been punished somewhat: works reasonable well
						//if (dd_synapse_permanence_segment[synapse_i] = 0.0f) // replace synapses that have been punished much: does not work
						{
							indices_to_update[n_synapses_replaced] = synapse_i;
							n_synapses_replaced++;
							if (n_synapses_replaced >= new_number_synapses) break;
						}
						else
						{
							n_synapses_kept++;
						}
					}

					const int n_new_synapses = new_number_synapses - n_synapses_replaced;
					if (n_new_synapses > 0)
					{
						//log_INFO("TP:replace_DD_synapses: n_new_synapses = ", n_new_synapses);

						for (int i = n_synapses_replaced; i < new_number_synapses; ++i)
						{
							indices_to_update[i] = old_size + i;
						}
					}

					const int n_synapses_update = n_synapses_replaced + n_new_synapses;

					if (false) log_INFO_DEBUG("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; n_new_synapses = ", n_new_synapses, ".");

					//2] get new cells from which we will get input
					std::vector<int> selected_cells(n_synapses_update);
					select_cell_to_learn_on(column, segment, prev_winner_cells_all, n_synapses_update, selected_cells);

					#if _DEBUG
					if (false)
					{
						log_INFO<false, false>("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, ": cells selected for new synapses:");
						print::print_int_array(selected_cells, n_new_synapses);
					}
					#endif

					//3] create space for the new synapses
					if (dd_synapse_permanence_segment.size() < n_synapses_update)
					{
						for (int i = 0; i < htm::tools::multiple_16(n_new_synapses); ++i)
						{
							dd_synapse_permanence_segment.push_back(P::TP_DD_PERMANENCE_INVALID);
						}
					}

					//4] update the synapses with n_new_synapses new cells,
					for (auto counter = 0; counter < n_synapses_update; ++counter)
					{
						const int synapse_i = indices_to_update[counter];
						dd_synapse_permanence_segment[synapse_i] = P::TP_DD_PERMANENCE_INIT;
						dd_synapse_origin_segment[synapse_i] = selected_cells[counter];
					}
					column.dd_synapse_count[segment] = n_synapses_update;

					//if (n_synapses_replaced < new_number_synapses) log_INFO("TP:replace_DD_synapses: n_synapses_kept = ", n_synapses_kept, "; n_synapses_replaced = ", n_synapses_replaced, "; new n_synapses_update = ", n_synapses_update);



					#if _DEBUG
					if (false) log_INFO("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; done replacing ", n_new_synapses, " synapses.");
					if (false)
					{
						log_INFO("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; dd_synapses:");
						print::print_dd_synapses(column);
					}
					#endif
				}

				inline void grow_DD_synapses_v2(
					Column& column,
					const int segment,
					const BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& prev_active_segments,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all)
				{
					//if (prev_active_segments.count() > 0) log_INFO("TP:grow_DD_synapses_v2: prev_active_segments = ", prev_active_segments.count());

					const int new_number_synapses = P::TP_DD_MAX_NEW_SYNAPSE_COUNT - prev_active_segments.count();
					if ((new_number_synapses) < 1) return;

					assert_msg(segment < column.dd_segment_count, "TP:replace_DD_synapses: segment", segment, " is too large");

					auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[segment];
					auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];
					assert_msg(dd_synapse_permanence_segment.size() == dd_synapse_origin_segment.size(), "TP:replace_DD_synapses: bug A.");

					const int old_capacity = dd_synapse_permanence_segment.size();
					const int old_size = column.dd_synapse_count[segment];
					assert_msg(old_size <= old_capacity, "TP:replace_DD_synapses: bug B.");

					//1] assume we are not going to replace any synapses, thus we will be adding 
					// new_number_synapses synapses. Allocate space for these synapses

					if (old_size + new_number_synapses >= old_capacity)
					{	// not enough capacity
						dd_synapse_permanence_segment.resize(old_size + new_number_synapses);
						dd_synapse_origin_segment.resize(old_size + new_number_synapses);
					}

					//2] find indices that will be overwritten with new values
					std::vector<int> indices_to_update = std::vector<int>(new_number_synapses, -1);

					int n_synapses_kept = 0;

					//3] count the number of synapses that can be replaced
					int n_synapses_replaced = 0;
					for (auto synapse_i = 0; synapse_i < old_size; ++synapse_i)
					{
						//if (true) // strangely enough it works better to replace existing synapses only
						//if (dd_synapse_permanence_segment[synapse_i] < P::TP_DD_PERMANENCE_THRESHOLD) // replace synapses that are not connected: works badly
						//if (dd_synapse_permanence_segment[synapse_i] < 0.1f) // replace synapses that have been punished somewhat: works reasonable well
						if (dd_synapse_permanence_segment[synapse_i] < 0.15f) // replace synapses that have been punished somewhat: works reasonable well
						//if (dd_synapse_permanence_segment[synapse_i] = 0.0f) // replace synapses that have been punished much: does not work
						{
							indices_to_update[n_synapses_replaced] = synapse_i;
							n_synapses_replaced++;
							if (n_synapses_replaced >= new_number_synapses) break;
						}
						else
						{
							n_synapses_kept++;
						}
					}
					//if (n_synapses_replaced > 0) log_INFO("TP:replace_DD_synapses: n_synapses_kept = ", n_synapses_kept, "; n_synapses_replaced = ", n_synapses_replaced);

					const int n_new_synapses = new_number_synapses - n_synapses_replaced;
					if (n_new_synapses > 0)
					{
						for (int i = n_synapses_replaced; i < new_number_synapses; ++i)
						{
							indices_to_update[i] = old_size + i;
						}
					}

					const int n_synapses_update = n_synapses_replaced + n_new_synapses;

					if (false) log_INFO_DEBUG("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; n_new_synapses = ", n_new_synapses, ".");

					//2] get new cells from which we will get input
					std::vector<int> selected_cells(n_synapses_update);
					select_cell_to_learn_on(column, segment, prev_winner_cells_all, n_synapses_update, selected_cells);

					#if _DEBUG
					if (false)
					{
						log_INFO<false, false>("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, ": cells selected for new synapses:");
						print::print_int_array(selected_cells, n_new_synapses);
					}
					#endif

					//3] create space for the new synapses
					if (dd_synapse_permanence_segment.size() < n_synapses_update)
					{
						for (int i = 0; i < htm::tools::multiple_16(n_new_synapses); ++i)
						{
							dd_synapse_permanence_segment.push_back(P::TP_DD_PERMANENCE_INVALID);
						}
					}

					//4] update the synapses with n_new_synapses new cells,
					for (auto counter = 0; counter < n_synapses_update; ++counter)
					{
						const int synapse_i = indices_to_update[counter];
						dd_synapse_permanence_segment[synapse_i] = P::TP_DD_PERMANENCE_INIT;
						dd_synapse_origin_segment[synapse_i] = selected_cells[counter];
					}
					column.dd_synapse_count[segment] = n_synapses_update;

					#if _DEBUG
					if (false) log_INFO("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; done replacing ", n_new_synapses, " synapses.");
					if (false)
					{
						log_INFO("TP:replace_DD_synapses: column ", column.id, "; segment ", segment, "; dd_synapses:");
						print::print_dd_synapses(column);
					}
					#endif
				}
				
				inline void grow_DD_synapses_org(
					Column& column,
					const int segment,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all) 
				{
					/*
						@classmethod
						def _growSynapses(cls, connections, random, segment, nDesiredNewSynapes, prevWinnerCells, initialPermanence, maxSynapsesPerSegment) :
						"""
						Creates nDesiredNewSynapes synapses on the segment passed in if
						possible, choosing random cells from the previous winner cells that are
						not already on the segment.

						: param connections : (Object)Connections instance for the tm
						: param random : (Object)TM object used to generate random numbers
						: param segment : (int)Segment to grow synapses on.
						: param nDesiredNewSynapes : (int)Desired number of synapses to grow
						: param prevWinnerCells : (list)Winner cells in `t-1`
						:param initialPermanence : (float)Initial permanence of a new synapse.
						"""
						candidates = list(prevWinnerCells)

						for synapse in connections.synapsesForSegment(segment) :
							i = binSearch(candidates, synapse.presynapticCell)
							if i != -1 :
								del candidates[i]

						nActual = min(nDesiredNewSynapes, len(candidates))

						# Check if we're going to surpass the maximum number of synapses.
						overrun = connections.numSynapses(segment) + nActual - maxSynapsesPerSegment
						if overrun > 0:
							cls._destroyMinPermanenceSynapses(connections, random, segment, overrun, prevWinnerCells)

						# Recalculate in case we weren't able to destroy as many synapses as needed.
						nActual = min(nActual, maxSynapsesPerSegment - connections.numSynapses(segment))

						for _ in range(nActual):
							i = random.getUInt32(len(candidates))
							connections.createSynapse(segment, candidates[i], initialPermanence)
							del candidates[i]
					*/


				}

				inline void d(
					Column& column,
					const int segment,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all) 
				{
					//grow_DD_synapses_org(column, segment, prev_winner_cells_all);
					grow_DD_synapses_v1(column, segment, prev_winner_cells_all);
					//grow_DD_synapses_v2(column, segment, prev_winner_cells_all);
				}
			}
			
			inline void create_DD_segment(
				Column& column,
				const int8_t cell,
				const BitsetSparse<P::N_CELLS>& prev_winner_cells_all)
			{
				const int new_segment_i = column.dd_segment_count;
				if (new_segment_i >= P::TP_N_DD_SEGMENTS_MAX)
				{
					log_WARNING("TP:create_new_DD_segment: too many segments for column ", column.id, ", ignoring new segment for cell ", static_cast<int>(cell), ".");
					//TODO: remove the segment that is least recently used.
				}
				else
				{
					if (false) log_INFO_DEBUG("TP:create_new_DD_segment: column ", column.id, "; adding a new segment (", new_segment_i, ") to cell ", static_cast<int>(cell), "; new_segment_idx = ", new_segment_i, ".");

					const int n_new_synapses = P::TP_DD_MAX_NEW_SYNAPSE_COUNT;

					#pragma region Resize Synapses
					// this is the only place where synapse space is created
					assert_msg(column.dd_synapse_permanence.size() == column.dd_synapse_origin.size(), "TP:create_new_DD_segment: Bug A.");
					if (column.dd_synapse_permanence.size() <= new_segment_i)
					{
						const int n_new_synapsed_capacity = htm::tools::multiple_16(n_new_synapses); // multiple of 16 needed for vectorization;
						const auto empty_permanence = std::vector<float>(n_new_synapsed_capacity, P::TP_DD_PERMANENCE_INVALID); // init with invalid number for debugging purposes
						column.dd_synapse_permanence.push_back(empty_permanence);
						const auto empty_origin = std::vector<int>(n_new_synapsed_capacity, P::TP_DD_SYNAPSE_ORIGIN_INVALID); // init with invalid number for debugging purposes
						column.dd_synapse_origin.push_back(empty_origin);

						column.dd_synapse_count.resize(new_segment_i + 1);
						column.dd_segment_destination.resize(new_segment_i + 1);
					}
					#pragma endregion

					auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[new_segment_i];
					auto& dd_synapse_origin_segment = column.dd_synapse_origin[new_segment_i];

					column.dd_segment_count++;
					column.dd_segment_destination[new_segment_i] = cell;

					std::vector<int> selected_cells_local(n_new_synapses);
					select_cell_to_learn_on(column, new_segment_i, prev_winner_cells_all, n_new_synapses, selected_cells_local);

					for (auto synapse_i = 0; synapse_i < n_new_synapses; ++synapse_i)
					{
						dd_synapse_permanence_segment[synapse_i] = P::TP_DD_PERMANENCE_INIT;
						dd_synapse_origin_segment[synapse_i] = selected_cells_local[synapse_i];
					}
					column.dd_synapse_count[new_segment_i] = n_new_synapses;

					if (false) print::print_dd_synapses(column);
				}
			}

			template <bool LEARN>
			inline void burst_column(
				Column& column,
				const BitsetCell<P::N_CELLS>& prev_active_cells_all,
				const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
				//inout
				BitsetTiny<P::N_CELLS_PC>& active_cells,
				BitsetTiny<P::N_CELLS_PC>& winner_cells)
			{
				if (false) log_INFO_DEBUG("TP:burst_column: column ", column.id, " bursts.");

				active_cells.set(); //burst. mark all cells as active

				const std::tuple<bool, int, bool, int> tup = priv::best_matching_cell(column, prev_active_cells_all);
				const bool found_cell = std::get<0>(tup);
				const int best_cell = std::get<1>(tup);
				const bool found_segment = std::get<2>(tup);
				const int best_segment = std::get<3>(tup);

				if (found_cell)
				{
					if (false) log_INFO_DEBUG("TP:burst_column: column ", column.id, "; cell ", best_cell, " is added to winner cells.");
					winner_cells.set(best_cell);
				}
				const bool has_prev_winner_cells = prev_winner_cells_all.any();

				if (false) log_INFO_DEBUG("TP:burst_column: column ", column.id, "; found_segment = ", found_segment, "; has_prev_winner_cells = ", has_prev_winner_cells, "; n_prev_winner_cells = ", prev_winner_cells_all.count(), ".");

				if (!found_segment && has_prev_winner_cells)
				{
					#if _DEBUG
					if (false)
					{
						log_INFO<false, false>("TP:burst_column:prev_winner_cells: ");
						print::print_bitset(prev_winner_cells_all);
					}
					#endif
					create_DD_segment(column, best_cell, prev_winner_cells_all);
				}
			}

			//Change the permanence values of DD synapses on the provided segment in the provided column
			namespace adapt_segment
			{
				//#pragma omp declare simd notinbranch
				inline float min_simd(float a, float b) { return a < b ? a : b; }

				//#pragma omp declare simd notinbranch
				inline float max_simd(float a, float b) { return a > b ? a : b; }

				inline void adapt_segment_ref(
					Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS>& prev_active_cells_all,
					const float permanence_inc,
					const float permanence_dec)
				{
					const auto& dd_synapse_origin_segment = column.dd_synapse_origin[segment];
					auto& dd_synapse_permanence_segment = column.dd_synapse_permanence[segment];

					for (auto synapse_i = 0; synapse_i < column.dd_synapse_count[segment]; ++synapse_i)
					{
						const auto global_cell_id = dd_synapse_origin_segment[synapse_i];
						const float old_permanence = dd_synapse_permanence_segment[synapse_i];
						const bool b = prev_active_cells_all[global_cell_id];

						dd_synapse_permanence_segment[synapse_i] = (b)
							? min_simd(1.0f, old_permanence + permanence_inc)
							: max_simd(0.0f, old_permanence - permanence_dec);
					}
				}

				inline void adapt_segment_avx512(
					Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS>& prev_active_cells_all,
					const float permanence_inc,
					const float permanence_dec)
				{
					__m512 * permanence = reinterpret_cast<__m512 *>(column.dd_synapse_permanence[segment].data());
					const __m512i * origin = reinterpret_cast<const __m512i *>(column.dd_synapse_origin[segment].data());
					const char * cells_all_ptr = prev_active_cells_all.data();

					const __m512 inc = _mm512_set1_ps(permanence_inc);
					const __m512 dec = _mm512_set1_ps(permanence_dec);
					const __m512i one_i = _mm512_set1_epi32(1);
					const __m512i zero_i = _mm512_setzero_epi32();
					const __m512 one_f = _mm512_set1_ps(1.0);
					const __m512 zero_f = _mm512_setzero();

					const int n_synapses = column.dd_synapse_count[segment];
					const int n_blocks = n_synapses >> 4;
					const int tail = n_synapses & 0b1111;

					//bug here

					for (int block = 0; block < n_blocks; ++block)
					{
						const __m512 p_old = permanence[block];
						const __mmask16 k = _mm512_cmpeq_epi32_mask(_mm512_and_epi32(_mm512_i32gather_epi32(origin[block], cells_all_ptr, 1), one_i), one_i);
						permanence[block] = _mm512_mask_blend_ps(k, _mm512_min_ps(one_f, _mm512_add_ps(p_old, inc)), _mm512_max_ps(zero_f, _mm512_sub_ps(p_old, dec)));
					}
					if (tail > 0) {
						const __mmask16 m = (1 << tail) - 1;
						const __m512 p_old = permanence[n_blocks];
						const __mmask16 k = _mm512_cmpeq_epi32_mask(_mm512_and_epi32(_mm512_mask_i32gather_epi32(zero_i, m, origin[n_blocks], cells_all_ptr, 1), one_i), one_i);
						permanence[n_blocks] = _mm512_mask_blend_ps(k, _mm512_min_ps(one_f, _mm512_add_ps(p_old, inc)), _mm512_max_ps(zero_f, _mm512_sub_ps(p_old, dec)));
					}
				}

				inline void d(
					Column& column,
					const int segment,
					const BitsetCell<P::N_CELLS>& prev_active_cells_all,
					const float permanence_inc,
					const float permanence_dec)
				{
					if (false) log_INFO_DEBUG("TP:adapt_segment: column ", column.id, "; segment ", segment);

					if (ARCH == arch_t::X64) return adapt_segment_ref(column, segment, prev_active_cells_all, permanence_inc, permanence_dec);
					if (ARCH == arch_t::AVX) return adapt_segment_ref(column, segment, prev_active_cells_all, permanence_inc, permanence_dec);
					if (ARCH == arch_t::AVX512) return adapt_segment_ref(column, segment, prev_active_cells_all, permanence_inc, permanence_dec);
				}
			}

			//Punishes the segments that incorrectly predicted a column to be active.
			inline void punish_predicted_column(
				Column& column,
				const BitsetCell<P::N_CELLS>& prev_active_cells_all,
				const BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& prev_matching_segments)
			{
				if (P::TP_DD_PREDICTED_SEGMENT_DEC > 0.0)
				{
					//appreciate that this branch is taken for non active columns, thus we cannot parallelize on active columns
					//punish the segments that predict the activity of a cell active
					for (auto i = 0; i < prev_matching_segments.count(); ++i)
					{
						auto segment_i = prev_matching_segments.get(i);
						//const auto cell = column.dd_segment_destination[segment_i];
						//if (predicted_inactive_cells.get(cell))
						{
							adapt_segment::d(column, segment_i, prev_active_cells_all, -P::TP_DD_PREDICTED_SEGMENT_DEC, 0.0);
							//if (false) log_INFO("TP:learn_on_segments: column ", column.id, "; segment ", segment_i, " (from cell ", static_cast<int>(cell), ") is punished.");
						}
					}
				}
			}

			inline void learn_on_segments(
				Column& column,
				const BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& prev_active_segments,
				const BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& learning_segments,
				const BitsetCell<P::N_CELLS>& prev_active_cells_all,
				const BitsetTiny<P::N_CELLS_PC>& winner_cells,
				const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
				const BitsetTiny<P::N_CELLS_PC>& predicted_inactive_cells,
				const BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& prev_matching_segments)
			{
				#if _DEBUG
				if (false)
				{
					if (prev_active_segments.any())
					{
						log_INFO<false, false>("TP:learn_on_segments: column ", column.id, ": prev_active_segments = ");
						print::print_bitset(prev_active_segments);
					}
				}
				if (false)
				{
					if (learning_segments.any())
					{
						log_INFO<false, false>("TP:learn_on_segments: column ", column.id, ": learning_segments = ");
						print::print_bitset(learning_segments);
					}
				}
				#endif
				if (!(prev_active_segments.none() && learning_segments.none()))
				{
					for (auto segment_i = 0; segment_i < column.dd_segment_count; ++segment_i)
					{
						const bool is_learning_segment = learning_segments.contains(segment_i);
						const bool prev_active_segment = prev_active_segments.contains(segment_i);

						if (prev_active_segment || is_learning_segment)
						{
							const auto cell = column.dd_segment_destination[segment_i];
							const bool is_from_winner_cell = winner_cells.get(cell);

							if (false) log_INFO_DEBUG("TP:learn_on_segments: column ", column.id, "; segment ", segment_i, ": learning segment OR segment is previously active.");

							if (is_learning_segment || is_from_winner_cell)
							{
								adapt_segment::d(column, segment_i, prev_active_cells_all, P::TP_DD_PERMANENCE_INC, P::TP_DD_PERMANENCE_DEC);
								if (false) log_INFO_DEBUG("TP:learn_on_segments: column ", column.id, "; cell ", static_cast<int>(cell), "; segment ", segment_i, ": learning segment OR from winner cell,");
							}
							if (is_learning_segment)
							{
								grow_DD_synapses::d(column, segment_i, prev_winner_cells_all);
								//grow_DD_synapses::grow_DD_synapses_v2(column, segment_i, prev_active_segments, prev_winner_cells_all);
								if (false) log_INFO_DEBUG("TP:learn_on_segments: column", column.id, "; cell ", static_cast<int>(cell), "; segment", segment_i, " is a learning segment.");
							}
						}
					}
				}
				punish_predicted_column(column, prev_active_cells_all, prev_matching_segments);
			}

			namespace activate_cells
			{
				inline void activate_cells_per_column(
					Column& column,
					//in
					const bool learn,
					const bool active_column,
					const BitsetCell<P::N_CELLS>& prev_active_cells_all,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
					const Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& prev_predictive_cells_all_2D,
					//out
					BitsetTiny<P::N_CELLS_PC>& active_cells,
					BitsetTiny<P::N_CELLS_PC>& winner_cells)
				{
					const int column_i = column.id;

					tools::swap(column.prev_active_segments, column.active_segments);
					//tools::swap(column.prev_matching_segments, column.matching_segments);
					tools::swap(column.prev_matching_cells, column.matching_cells);

					//for each column
					//	if column is active and has active distal dendrite segments
					//		call activatePredictedColumn
					//	if column is active and doesn't have active distal dendrite segments
					//		call burstColumn
					//	if column is inactive and has matching distal dendrite segments
					//		call punishPredictedColumn

					const bool predicted_column = priv::activate_predicted_column(
						//in
						prev_predictive_cells_all_2D[column_i],
						column.prev_matching_cells,
						active_column,
						//out
						active_cells,
						winner_cells);

					if (active_column && !predicted_column)
					{
						priv::burst_column<LEARN>(
							column,
							//in 
							prev_active_cells_all,
							prev_winner_cells_all,
							//inout
							active_cells,
							winner_cells);
					}
					if (!active_column && learn)
					{
						priv::punish_predicted_column(
							column, 
							column.active_segments, 
							column.matching_segments, 
							prev_active_cells_all, 
							prev_winner_cells_all);
					}
				}

				using namespace std::placeholders;  // For `_1`

				class ApplyFoo
				{
					Network& network;
					const bool learn;
					const Bitset<P::N_COLUMNS>& active_columns;
					const BitsetCell<P::N_CELLS>& prev_active_cells_all;
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all;
					const Bitset2<P::N_COLUMNS, P::N_CELLS_PC> prev_predictive_cells_all_2D;
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& active_cells_all_2D;
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& winner_cells_all_2D;
				public:
					void operator()(const tbb::blocked_range<int>& r) const
					{
						for (int column_i = r.begin(); column_i < r.end(); ++column_i)
						{
							activate_cells_per_column(network[column_i], learn, active_columns[column_i], prev_active_cells_all, prev_winner_cells_all, prev_predictive_cells_all_2D, active_cells_all_2D[column_i], winner_cells_all_2D[column_i]);
						}
					}
					ApplyFoo(
						Network& network,
						const bool learn,
						const Bitset<P::N_COLUMNS>& active_columns,
						const BitsetCell<P::N_CELLS>& prev_active_cells_all,
						const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
						const Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& prev_predictive_cells_all_2D,
						Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& active_cells_all_2D,
						Bitset2<P::N_COLUMNS, P::N_CELLS_PC>& winner_cells_all_2D)
					:
						network(network),
						learn(learn),
						active_columns(active_columns),
						prev_active_cells_all(prev_active_cells_all),
						prev_winner_cells_all(prev_winner_cells_all),
						prev_predictive_cells_all_2D(prev_predictive_cells_all_2D),
						active_cells_all_2D(active_cells_all_2D),
						winner_cells_all_2D(winner_cells_all_2D)
					{}
				};

				template <bool LEARN>
				inline void d(
					Network& network,
					//in
					const Bitset<P::N_COLUMNS>& active_columns,
					const BitsetCell<P::N_CELLS>& prev_active_cells_all,
					const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
					const BitsetCell<P::N_CELLS>& prev_predictive_cells_all,
					//out
					BitsetCell<P::N_CELLS>& active_cells_all,
					BitsetCell<P::N_CELLS>& winner_cells_all)
				{
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> prev_predictive_cells_all_2D;
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> active_cells_all_2D;
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> winner_cells_all_2D;

					tools::copy(prev_predictive_cells_all_2D, prev_predictive_cells_all);

					if (PARALLEL)
					{
						auto c = ApplyFoo(network, LEARN, active_columns, prev_active_cells_all, prev_winner_cells_all, prev_predictive_cells_all_2D, active_cells_all_2D, winner_cells_all_2D);
						tbb::parallel_for(tbb::blocked_range<int>(0, P::N_COLUMNS, P::N_COLUMNS / 64), c);
					}
					else
					{
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							activate_cells_per_column(
								network[column_i],
								LEARN,
								active_columns[column_i],
								prev_active_cells_all,
								prev_winner_cells_all,
								prev_predictive_cells_all_2D,
								//out
								active_cells_all_2D[column_i],
								winner_cells_all_2D[column_i]);
						}
					}

					tools::copy(active_cells_all, active_cells_all_2D);
					tools::copy(winner_cells_all, winner_cells_all_2D);
				}
			}

			namespace activate_dendrites
			{
				inline void compute_activity(
					//in
					const Column& column,
					const BitsetCell<P::N_CELLS>& active_cells_all,
					const float permanence_threshold,
					const int synapse_threshold,
					//out
					BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& active_segments,
					BitsetTiny<P::N_CELLS_PC>& predictive_cells)
				{
					if (false) log_INFO_DEBUG("TP:compute_activity: column ", column.id, "; number of active cells = ", active_cells_all.count());

					predictive_cells.reset();
					active_segments.reset();

					for (auto segment_i = 0; segment_i < column.dd_segment_count; ++segment_i)
					{
						const int n_active_synapses = get_number_active_DD_synapses::d(column, segment_i, active_cells_all, permanence_threshold);

						if (false) log_INFO_DEBUG("TP:compute_activity: column ", column.id, "; segment ", segment_i, ": n_active_synapes = ", n_active_synapses, ", synapse_threshold = ", synapse_threshold, ", permanence_threshold = ", permanence_threshold);
						if (n_active_synapses >= synapse_threshold)
						{
							active_segments.add(segment_i);
							const auto cell = column.dd_segment_destination[segment_i];
							predictive_cells.set(cell);
						}
					}
				}

				inline void compute_activity2(
					//in
					const Column& column,
					const BitsetCell<P::N_CELLS>& active_cells_all,
					const float permanence_threshold1,
					const int synapse_threshold1,
					const float permanence_threshold2,
					const int synapse_threshold2,
					//out
					BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& active_segments1,
					BitsetTiny<P::N_CELLS_PC>& predictive_cells1,
					BitsetSparse<P::TP_N_DD_SEGMENTS_MAX>& active_segments2,
					BitsetTiny<P::N_CELLS_PC>& predictive_cells2)
				{
					predictive_cells1.reset();
					active_segments1.reset();
					predictive_cells2.reset();
					active_segments2.reset();

					for (auto segment_i = 0; segment_i < column.dd_segment_count; ++segment_i)
					{
						const auto tup = get_number_active_DD_synapses::get_number_active_DD_synapses_2x_avx512(column, segment_i, active_cells_all, permanence_threshold1, permanence_threshold2);
						const int n_active_synapses1 = std::get<0>(tup);
						const int n_active_synapses2 = std::get<1>(tup);

						if (false) log_INFO_DEBUG("TP:compute_activity: column ", column.id, "; segment ", segment_i, ": n_active_synapses1 = ", n_active_synapses1, ", synapse_threshold1 = ", synapse_threshold1, ", permanence_threshold1 = ", permanence_threshold1);
						if (n_active_synapses1 >= synapse_threshold1)
						{
							active_segments1.add(segment_i);
							const auto cell = column.dd_segment_destination[segment_i];
							predictive_cells1.set(cell);
						}
						if (n_active_synapses2 >= synapse_threshold2)
						{
							active_segments2.add(segment_i);
							const auto cell = column.dd_segment_destination[segment_i];
							predictive_cells2.set(cell);
						}
					}
				}

				template <bool LEARN>
				inline void d(
					Network& network,
					//in
					const BitsetCell<P::N_CELLS>& active_cells_all,
					//out
					BitsetCell<P::N_CELLS>& predictive_cells_all)
				{
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> predictive_cells_all_2D;

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						Column& column = network[column_i];

						if (false)
						{
							compute_activity(
								//in
								column,
								active_cells_all,
								P::TP_DD_PERMANENCE_THRESHOLD,
								P::TP_DD_ACTIVATION_THRESHOLD,
								//out
								column.active_segments,
								predictive_cells_all_2D[column_i]);

							compute_activity(
								//in
								column,
								active_cells_all,
								0.0,
								P::MIN_DD_ACTIVATION_THRESHOLD,
								//out
								column.matching_segments,
								column.matching_cells);
						}
						else
						{
							// compute_activity2 is slightly slower than compute_activity, thus this branch is almost twice as fast.
							compute_activity2(
								//in
								column,
								active_cells_all,
								P::TP_DD_PERMANENCE_THRESHOLD,
								P::TP_DD_ACTIVATION_THRESHOLD,
								0.0,
								P::MIN_DD_ACTIVATION_THRESHOLD,
								//out
								column.active_segments,
								predictive_cells_all_2D[column_i],
								column.matching_segments,
								column.matching_cells);
						}
					}
					tools::copy(predictive_cells_all, predictive_cells_all_2D);
				}
			}
		}

		template <bool LEARN>
		inline void compute_tp(
			const int t,
			Network& network,
			const Bitset<P::N_COLUMNS>& active_columns,
			const BitsetCell<P::N_CELLS>& prev_active_cells_all,
			const BitsetSparse<P::N_CELLS>& prev_winner_cells_all,
			const BitsetCell<P::N_CELLS>& prev_predictive_cells_all,
			// out
			BitsetCell<P::N_CELLS>& active_cells_all,
			BitsetCell<P::N_CELLS>& winner_cells_all,
			BitsetCell<P::N_CELLS>& predictive_cells_all)
		{
			#if _DEBUG
			if (false)
			{
				log_INFO<false, false>("TP:compute_tp: prev_predictive_cells: ");
				print::print_active_cells(prev_predictive_cells_all);
			}
			if (false)
			{
				log_INFO<false, false>("TP:compute_tp: prev_winner_cells: ");
				//print::print_active_cells(prev_winner_cells_all);
			}
			if (false)
			{
				log_INFO<false, false>("TP:compute_tp: prev_active_cells: ");
				print::print_active_cells(prev_active_cells_all);
			}
			#endif
			
			priv::activate_cells::d<LEARN>(
				network,
				active_columns,
				prev_active_cells_all,
				prev_winner_cells_all,
				prev_predictive_cells_all,
				//out
				active_cells_all,
				winner_cells_all);

			#if _DEBUG
			if (true)
			{
				log_INFO<false, false>("TP:compute_tp: active_cells_all: t = ", t, ":");
				print::print_active_cells(active_cells_all);
			}
			if (true)
			{
				log_INFO<false, false>("TP:compute_tp: winner_cells_all: t = ", t, ":");
				print::print_active_cells(winner_cells_all);
			}
			#endif

			//barrier: active cells of all columns have to be known

			priv::activate_dendrites::d<LEARN>(
				network, 
				active_cells_all,
				//out
				predictive_cells_all);
			
			#if _DEBUG
			if (true)
			{
				log_INFO<false, false>("TP:compute_tp: predictive_cells_all: t = ", t, ":");
				print::print_active_cells(predictive_cells_all);
			}
			if (false)
			{
				log_INFO("TP:compute_tp: n_active_cells = ", active_cells_all.count(), "; segment activation threshold = ", P::TP_DD_ACTIVATION_THRESHOLD, "; synapse permanence threshold = ", P::TP_DD_PERMANENCE_THRESHOLD, "; segment activity is: ");
				print::print_segment_activity(network, active_cells_all, P::TP_DD_PERMANENCE_THRESHOLD);
			}
			if ((false) && (t > 7))
			{
				log_INFO("TP:compute_tp: all synapses: t = ", t, ":");
				print::print_dd_synapses(network);
			}
			#endif
		}
	}
}