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
#include <algorithm>	// std::min
#include <limits>		// std::numeric_limits
#include <tuple>
#include <array>
#include <vector>
#include <algorithm>
#include <functional>

#include "../Spike-Tools-Lib/log.ipp"
#include "../Spike-Tools-Lib/assert.ipp"

#include "constants.ipp"
#include "print.ipp"
#include "types.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Spacial Pooler
	namespace sp 
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;
		using namespace htm::types;

		//HTM Spacial Pooler private methods
		namespace priv
		{
			namespace inhibit_columns
			{
				template <typename P>
				void active_columns_ref1(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const int inhibition_top,
					Bitset<P::N_COLUMNS>& active_columns)
				{
					for (int column_i1 = 0; column_i1 < P::N_COLUMNS; ++column_i1)
					{
						const float oa = boosted_overlap[column_i1];
						int sum = 0;

						for (int column_i2 = 0; column_i2 < P::N_COLUMNS; ++column_i2)
						{
							sum += (boosted_overlap[column_i2] > oa);
						}
						active_columns[column_i1] = sum < inhibition_top;
					}
				}

				// faster than ref1
				template <typename P>
				void active_columns_ref2(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const int inhibition_top,
					Bitset<P::N_COLUMNS>& active_columns)
				{
					active_columns.reset();

					for (int i = 0; i < inhibition_top; ++i)
					{
						float f = -10;
						int best_i = -1;

						for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							if ((boosted_overlap[column_i] > f) && !active_columns[column_i])
							{
								f = boosted_overlap[column_i];
								best_i = column_i;
							}
						}
						if (best_i != -1) active_columns[best_i] = true;
					}
					#if _DEBUG
					Bitset<P::N_COLUMNS> active_columns2;
					active_columns_ref1(boosted_overlap, inhibition_top, active_columns2);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2[column_i] != active_columns[column_i])
					{
						log_INFO("SP:active_columns_ref2: active_columns=");
						print::print_bitset(active_columns);
						log_INFO("SP:active_columns_ref2: active_columns2=");
						print::print_bitset(active_columns2);
						::tools::log::log_ERROR("SP:active_columns_ref2: column ", column_i, " is not unequal (inhibition_top=", inhibition_top, ")");
					}
					#endif
				}

				// faster than ref2
				template <typename P>
				void active_columns_ref3(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const int inhibition_top,
					Bitset<P::N_COLUMNS>& active_columns)
				{
					active_columns.reset();
					std::vector<float> tmp = std::vector<float>(boosted_overlap);

					for (int i = 0; i < inhibition_top; ++i)
					{
						float f = -10;
						int best_i = -1;

						for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							if (tmp[column_i] > f)
							{
								f = boosted_overlap[column_i];
								best_i = column_i;
							}
						}
						if (best_i != -1)
						{
							active_columns[best_i] = true;
							tmp[best_i] = -1;
						}
					}

					#if _DEBUG
					Bitset<P::N_COLUMNS> active_columns2;
					active_columns_ref1(boosted_overlap, inhibition_top, active_columns2);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2[column_i] != active_columns[column_i])
					{
						log_INFO("SP:active_columns_ref3: active_columns=");
						print::print_bitset(active_columns);
						log_INFO("SP:active_columns_ref3: active_columns2=");
						print::print_bitset(active_columns2);
						log_ERROR("SP:active_columns_ref3: column ", column_i, " is not unequal (inhibition_top=", inhibition_top, ")");
					}
					#endif
				}

				// faster than ref3 when n_columsn > 2000; for smaller number of columns (500) speed is comparable to ref1
				template <typename P>
				void active_columns_ref4(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const int inhibition_top,
					Bitset<P::N_COLUMNS>& active_columns)
				{
					active_columns.reset();
					std::vector<float> boosted_overlap_copy = std::vector<float>(boosted_overlap);

					std::nth_element(boosted_overlap_copy.begin(), boosted_overlap_copy.begin() + inhibition_top, boosted_overlap_copy.end(), std::greater<float>());
					const float nth = boosted_overlap_copy[inhibition_top];

					for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						if (boosted_overlap[column_i] > nth) active_columns[column_i] = true;
					}

					#if _DEBUG
					if (false) // it seems ref1 does not always yield the same results as ref4
					{
						Bitset<P::N_COLUMNS> active_columns2;
						active_columns_ref1(boosted_overlap, inhibition_top, active_columns2);
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2[column_i] != active_columns[column_i])
						{
							log_INFO("SP:active_columns_ref4: active_columns=");
							print::print_bitset(active_columns);
							log_INFO("SP:active_columns_ref4: active_columns2=");
							print::print_bitset(active_columns2);
							log_ERROR("SP:active_columns_ref4: column ", column_i, " is not unequal (inhibition_top=", inhibition_top, ";nth=", nth, ")");
						}
					}
					#endif
				}

				template <typename P>
				void d(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const Dynamic_Param& param,
					Bitset<P::N_COLUMNS>& active_columns)
				{
					const int inhibition_top = static_cast<int>(P::N_COLUMNS * param.SP_LOCAL_AREA_DENSITY);

					if (P::SP_GLOBAL_INHIBITION)
					{
						//active_columns_ref1<P>(boosted_overlap, inhibition_top, active_columns);
						active_columns_ref4<P>(boosted_overlap, inhibition_top, active_columns);
					}
					else
					{
						log_ERROR("SP::active_columns: not implemented yet.");
					}

					#if _DEBUG
					if (false) if (active_columns.count() != inhibition_top) log_WARNING("SP:compute_sp: number_of_active_columns = ", active_columns.count(), " is not equal to inhibition_top = ", inhibition_top);
					#endif
				}
			}

			namespace calc_overlap
			{
				// Why is this an inefficient algorithm: for example, assume 1024 columns,
				// each column has 256 synapses and there are 400 sensor cells. The 400 cells 
				// are stored in 40 bytes, which fits in one L1 cache line (which is 64 bytes = 512 bits),
				// every columns is going to read 256 bytes, that is, 256KB will be read in total, 
				// thus every byte (of the 40) will on average be read 6553 times! Since only 5 percent 
				// of the sensors are active, it may be 20 times faster to do a scatter instead of a gather.

				template <typename P>
				void calc_overlap_ref(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out
					std::vector<int>& overlaps) //size = P::N_COLUMNS
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const auto& column = layer[column_i];

						int overlap = 0;
						for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
						{
							if (column.pd_synapse_permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
							{
								const auto origin_sensor = column.pd_synapse_origin[synapse_i];
								overlap += sensor_activity.get(origin_sensor); // deadly gather here!
							}
						}
						if (false) log_INFO_DEBUG("SP:calc_overlap_ref: column ", column.id, " has overlap = ", overlap, ".\n");

						overlaps[column_i] = (overlap < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap;
					}
				}

				template <typename P>
				void calc_overlap_scatter(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out
					std::vector<int>& overlaps) //size = P::N_COLUMNS
				{
					for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
					{
						if (sensor_activity.get(sensor_i))
						{
							const auto& destination_columns = layer.sp_pd_destination_column[sensor_i];
							const auto& permanences = layer.sp_pd_synapse_permanence[sensor_i];

							for (auto i = 0; i < destination_columns.size(); ++i)
							{
								if (permanences[i] > param.SP_PD_CONNECTED_THRESHOLD)
								{
									const auto column_i = destination_columns[i];
									overlaps[column_i]++; //deadly scatter here!
								}
							}
						}
					}
				}

				__m512i get_sensors_epi32(
					const __mmask16 mask,
					const __m512i origin,
					const void * active_sensors_ptr)
				{
					const __m512i int_addr = _mm512_srli_epi32(origin, 5);
					const __m512i pos_in_int = _mm512_and_epi32(origin, _mm512_set1_epi32(0b11111));
					const __m512i sensor_int = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), mask, int_addr, active_sensors_ptr, 4);
					const __m512i sensors = _mm512_srlv_epi32(sensor_int, pos_in_int);
					return _mm512_and_epi32(sensors, _mm512_set1_epi32(1));
				}

				__m512i get_sensors_epi32(
					const __mmask16 mask,
					const __m512i origin_epi32,
					const __m512i active_sensors)
				{
					const __m512i int_addr = _mm512_srli_epi32(origin_epi32, 5);
					const __m512i pos_in_int = _mm512_and_epi32(origin_epi32, _mm512_set1_epi32(0b11111));
					const __m512i sensor_int = _mm512_mask_permutexvar_epi32(_mm512_setzero_epi32(), mask, int_addr, active_sensors);
					const __m512i sensors = _mm512_srlv_epi32(sensor_int, pos_in_int);
					return _mm512_and_epi32(sensors, _mm512_set1_epi32(1));
				}

				__m512i get_sensors_epi16(
					const __mmask32 mask,
					const __m512i origin_epi32,
					const __m512i active_sensors)
				{
					const __m512i short_addr = _mm512_srli_epi16(origin_epi32, 4);
					const __m512i pos_in_int = _mm512_and_si512(origin_epi32, _mm512_set1_epi16(0b1111));
					const __m512i sensor_int = _mm512_mask_permutexvar_epi16(_mm512_setzero_epi32(), mask, short_addr, active_sensors);
					const __m512i sensors = _mm512_srlv_epi16(sensor_int, pos_in_int);
					return _mm512_and_si512(sensors, _mm512_set1_epi16(1));
				}

				template <typename P>
				void calc_overlap_avx512(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out 
					std::vector<int>& overlaps) //size = P::N_COLUMNS
				{
					const __m512i connected_threshold_epi8 = _mm512_set1_epi8(param.SP_PD_CONNECTED_THRESHOLD);
					auto active_sensors_ptr = sensor_activity.data();
					const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const auto& column = layer[column_i];
						auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_permanence.data());
						auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_origin.data());

						__m512i overlap = _mm512_setzero_epi32();

						for (int block = 0; block < n_blocks; ++block)
						{
							const __mmask64 connected_mask_64 = _mm512_cmp_epi8_mask(permanence_epi8_ptr[block], connected_threshold_epi8, _MM_CMPINT_NLE);
							for (int i = 0; i < 4; ++i)
							{
								const __mmask16 mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
								if (mask_16 != 0) overlap = _mm512_add_epi32(overlap, get_sensors_epi32(mask_16, origin_epi32_ptr[(block * 4) + i], active_sensors_ptr));
							}
						}
						const int overlap_int = _mm512_reduce_add_epi32(overlap);
						overlaps[column_i] = (overlap_int < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap_int;

						if (false) log_INFO_DEBUG("SP:calc_overlap_avx512: column ", column.id, " has overlap = ", overlaps[column_i], ".\n");
					}
					#if _DEBUG
					std::vector<int> overlaps_ref = std::vector<int>(P::N_COLUMNS);
					priv::calc_overlap::calc_overlap_ref(layer, param, sensor_activity, overlaps_ref);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const int overlap_ref = overlaps_ref[column_i];
						const int overlap_avx512 = overlaps[column_i];
						if (overlap_ref != overlap_avx512) log_ERROR("SP:calc_overlap_avx512:: UNEQUAL: column ", column_i, "; overlap ref ", overlap_ref, " != avx512 ", overlap_avx512, ".\n");
					}
					#endif
				}

				//P::N_SENSORS < 512
				template <typename P>
				void calc_overlap_avx512_small_epi32(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out 
					std::vector<int>& overlaps)
				{
					assert_msg(P::N_SENSORS < 512, "ERROR: calc_overlap_avx512_small: N_SENSORS is larger than 512");

					const __m512i connected_threshold_epi8 = _mm512_set1_epi8(param.SP_PD_CONNECTED_THRESHOLD);
					const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

					std::array<int, 16> t = { 0 };
					for (int i = 0; i < sensor_activity.N_BLOCKS; ++i)
					{
						t[i] = sensor_activity._data[i];
					}
					const __m512i active_sensors_simd = _mm512_load_epi32(t.data());

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const auto& column = layer[column_i];
						auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_permanence.data());
						auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_origin.data());

						__m512i overlap_epi16 = _mm512_setzero_epi32();

						for (int block = 0; block < n_blocks; ++block)
						{
							const __mmask64 connected_mask_64 = _mm512_cmp_epi8_mask(permanence_epi8_ptr[block], connected_threshold_epi8, _MM_CMPINT_NLE);
							for (int i = 0; i < 4; ++i)
							{
								const __mmask16 mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
								//if (mask_16 != 0) // slower with this check
								{
									overlap_epi16 = _mm512_add_epi32(overlap_epi16, get_sensors_epi32(mask_16, origin_epi32_ptr[(block * 4) + i], active_sensors_simd));
								}
							}
						}
						const int overlap_int = _mm512_reduce_add_epi32(overlap_epi16);
						overlaps[column_i] = (overlap_int < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap_int;
					}

					#if _DEBUG
					std::vector<int> overlaps_ref = std::vector<int>(P::N_COLUMNS);
					priv::calc_overlap::calc_overlap_ref(layer, param, sensor_activity, overlaps_ref);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const int overlap_ref = overlaps_ref[column_i];
						const int overlap_avx512 = overlaps[column_i];
						if (overlap_ref != overlap_avx512) log_ERROR("SP:calc_overlap_avx512_small_epi32:: UNEQUAL: column ", column_i, "; overlap ref ", overlap_ref, " != avx512 ", overlap_avx512, ".\n");
					}
					#endif
				}
				// not much faster than calc_overlap_avx512_small_epi32
				template <typename P>
				void calc_overlap_avx512_small_epi16(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out 
					std::vector<int>& overlaps)
				{
					assert_msg(P::N_SENSORS < 512, "ERROR: calc_overlap_avx512_small: N_SENSORS is larger than 512");

					const __m512i connected_threshold_epi8 = _mm512_set1_epi8(param.SP_PD_CONNECTED_THRESHOLD);
					const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

					std::array<int, 16> t = { 0 };
					for (int i = 0; i < sensor_activity.N_BLOCKS; ++i)
					{
						t[i] = sensor_activity._data[i];
					}
					const __m512i active_sensors_simd = _mm512_load_epi32(t.data());

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const auto& column = layer[column_i];
						auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_permanence.data());
						auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(column.pd_synapse_origin.data());

						__m512i overlap_epu16_AB = _mm512_setzero_si512(); // contains 32 overlap values of 16bits
						__m512i overlap_epu16_CD = _mm512_setzero_si512(); // contains 32 overlap values of 16bits

						for (int block = 0; block < n_blocks; ++block)
						{
							const __mmask64 connected_mask_64 = _mm512_cmp_epi8_mask(permanence_epi8_ptr[block], connected_threshold_epi8, _MM_CMPINT_NLE);

							const __mmask32 mask_32_A = static_cast<__mmask32>(connected_mask_64 >> (0 * 32));
							const __m512i origin_epi32_A = origin_epi32_ptr[(block * 4) + 0];
							const __m512i origin_epi32_B = origin_epi32_ptr[(block * 4) + 1];
							const __m512i origin_epu16_A = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_A));
							const __m512i origin_epu16_B = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_B));
							const __m512i origin_epu16_AB = _mm512_shuffle_i64x2(origin_epu16_A, origin_epu16_B, 0b01000100);
							overlap_epu16_AB = _mm512_adds_epu16(overlap_epu16_AB, get_sensors_epi16(mask_32_A, origin_epu16_AB, active_sensors_simd));

							const __mmask32 mask_32_B = static_cast<__mmask32>(connected_mask_64 >> (1 * 32));
							const __m512i origin_epi32_C = origin_epi32_ptr[(block * 4) + 2];
							const __m512i origin_epi32_D = origin_epi32_ptr[(block * 4) + 3];
							const __m512i origin_epu16_C = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_C));
							const __m512i origin_epu16_D = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_D));
							const __m512i origin_epu16_CD = _mm512_shuffle_i64x2(origin_epu16_C, origin_epu16_D, 0b01000100);
							overlap_epu16_CD = _mm512_adds_epu16(overlap_epu16_CD, get_sensors_epi16(mask_32_B, origin_epu16_CD, active_sensors_simd));
						}
						const __m512i overlap_epu16 = _mm512_adds_epu16(overlap_epu16_AB, overlap_epu16_CD);
						const __m512i overlap_epi32 = _mm512_add_epi32(_mm512_and_epi32(overlap_epu16, _mm512_set1_epi32(0xFFFF)), _mm512_srli_epi32(overlap_epu16, 16));
						const int overlap_int = _mm512_reduce_add_epi32(overlap_epi32);

						overlaps[column_i] = (overlap_int < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap_int;
					}

					#if _DEBUG
					std::vector<int> overlaps_ref = std::vector<int>(P::N_COLUMNS);
					priv::calc_overlap::calc_overlap_ref(layer, param, sensor_activity, overlaps_ref);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const int overlap_ref = overlaps_ref[column_i];
						const int overlap_avx512 = overlaps[column_i];
						if (overlap_ref != overlap_avx512) log_ERROR("SP:calc_overlap_avx512_small_epi16:: UNEQUAL: column ", column_i, "; overlap ref ", overlap_ref, " != avx512 ", overlap_avx512, ".\n");
					}
					#endif
				}

				template <typename P>
				void d(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity,
					//out
					std::vector<int>& overlaps) //size = P::N_COLUMNS
				{
					if (SP_GATHER)
					{
						if (P::ARCH == arch_t::X64) calc_overlap_ref(layer, param, sensor_activity, overlaps);
						if (P::ARCH == arch_t::AVX512)
							if (P::N_SENSORS < 512)
							{
								calc_overlap_avx512_small_epi16(layer, param, sensor_activity, overlaps);
								//calc_overlap_avx512_small_epi32(layer, param, sensor_activity, overlaps);
							}
							else
								calc_overlap_avx512(layer, param, sensor_activity, overlaps);
					}
					else
					{
						if (P::ARCH == arch_t::X64) calc_overlap_scatter(layer, param, sensor_activity, overlaps);
						if (P::ARCH == arch_t::AVX512) calc_overlap_scatter(layer, param, sensor_activity, overlaps);
					}
				}
			}

			namespace update_synapses
			{
				// assume b is positive
				Permanence add_saturate(Permanence a, Permanence b)
				{
					const int result = static_cast<int>(a) + static_cast<int>(b);
					return (result > 127) ? 127 : result;
				}
				// assume b is positive
				Permanence sub_saturate(Permanence a, Permanence b)
				{
					const int result = static_cast<int>(a) - static_cast<int>(b);
					return (result < -128) ? -128 : result;
				}

				template <typename P>
				void update_synapses_ref(
					Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset<P::N_COLUMNS>& active_columns,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity)
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						if (active_columns[column_i])
						{
							auto& column = layer[column_i];

							for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
							{
								const auto sensor_i = column.pd_synapse_origin[synapse_i];
								const int old_permanence = column.pd_synapse_permanence[synapse_i];
								const int increment = (sensor_activity.get(sensor_i)) ? param.SP_PD_PERMANENCE_INC : -param.SP_PD_PERMANENCE_DEC;
								const int new_permanence = std::min(127, std::max(-128, old_permanence + increment));

								if (false) log_INFO_DEBUG("SP:update_synapses: inc: column ", column_i, "; synpase ", synapse_i, ": old permanence ", old_permanence, "; new permanence = ", new_permanence, ".");
								column.pd_synapse_permanence[synapse_i] = static_cast<Permanence>(new_permanence);
							}
						}
					}
				}

				template <typename P>
				void update_synapses_scatter(
					Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset<P::N_COLUMNS>& active_columns,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity)
				{
					for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
					{
						const int increment = (sensor_activity.get(sensor_i)) ? param.SP_PD_PERMANENCE_INC : -param.SP_PD_PERMANENCE_DEC;

						auto& permanence = layer.sp_pd_synapse_permanence[sensor_i];
						const auto& destination_columns = layer.sp_pd_destination_column[sensor_i];

						for (auto synapse_i = 0; synapse_i < destination_columns.size(); ++synapse_i)
						{
							const int column_i = destination_columns[synapse_i];
							if (active_columns.get(column_i))
							{
								const int old_permanence = permanence[synapse_i];
								const int new_permanence = std::min(127, std::max(-128, old_permanence + increment));

								if (false) log_INFO_DEBUG("SP:update_synapses_scatter: inc: column ", column_i, "; synpase ", synapse_i, ": old permanence ", old_permanence, "; new permanence = ", new_permanence, ".");
								permanence[synapse_i] = new_permanence;
							}
						}
					}
				}

				template <typename P>
				void d(
					Layer<P>& layer,
					const Dynamic_Param& param,
					const Bitset<P::N_COLUMNS>& active_columns,
					const Bitset_Compact<P::N_SENSORS>& sensor_activity)
				{
					if (SP_GATHER)
					{
						if (P::ARCH == arch_t::X64) update_synapses_ref(layer, param, active_columns, sensor_activity);
						if (P::ARCH == arch_t::AVX512) update_synapses_ref(layer, param, active_columns, sensor_activity);
					}
					else
					{
						if (P::ARCH == arch_t::X64) update_synapses_scatter(layer, param, active_columns, sensor_activity);
						if (P::ARCH == arch_t::AVX512) update_synapses_scatter(layer, param, active_columns, sensor_activity);
					}
				}
			}

		/* 
			Updates the duty cycles for each column. The OVERLAP duty cycle is a moving
			average of the number of inputs which overlapped with the each column. The
			ACTIVITY duty cycles is a moving average of the frequency of activation for
			each column.

			overlap: An int vector containing the overlap score for each column.
			The overlap score for a column is defined as the number
			of synapses in a "connected state" (connected synapses)
			that are connected to input bits which are turned on.

			active_columns. An bit array containing the the active columns,
			the set of columns which survived inhibition
		*/
			template <typename P>
			void update_duty_cycles(
				Layer<P>& layer,
				const int column_i,
				const int overlap,
				const bool active_column)
			{
				//Checked: 20-09-17
				const int period = (P::SP_DUTY_CYCLE_PERIOD > layer.iteration_num) ? layer.iteration_num : P::SP_DUTY_CYCLE_PERIOD;
				layer.overlap_duty_cycles[column_i] = ((layer.overlap_duty_cycles[column_i] * (period - 1) + ((overlap > 0) ? 1 : 0)) / period);
				layer.active_duty_cycles[column_i] = ((layer.active_duty_cycles[column_i] * (period - 1) + ((active_column) ? 1 : 0)) / period);
			}

			template <typename P>
			void update_boost_factors(
				const Layer<P>& layer,
				Column<P>& column,
				const Dynamic_Param& param,
				const float active_duty_cycle)
			{
				/*
				Update the boost factors for all columns.The boost factors are used to
				increase the overlap of inactive columns to improve their chances of
				becoming active. and hence encourage participation of more columns in the
				learning process.
				*/
				if (P::SP_GLOBAL_INHIBITION)
				{
					const float target_desity = param.SP_LOCAL_AREA_DENSITY;
					const float new_boost_factor = exp((target_desity - active_duty_cycle) * P::SP_BOOST_STRENGTH);
					if (false) log_INFO_DEBUG("SP:update_boost_factors: column ", column.id, "; old boost factor ", column.boost_factor, "; new boost factor ", new_boost_factor, ".\n");
					column.boost_factor = new_boost_factor;
				}
				else
				{
				/*
					UInt inhibitionArea = pow((Real) (2 * inhibitionRadius_ + 1), (Real) columnDimensions_.size());
					inhibitionArea = min(inhibitionArea, numColumns_);
					targetDensity = ((Real) numActiveColumnsPerInhArea_) / inhibitionArea;
					targetDensity = min(targetDensity, (Real) 0.5);
					boostFactors_[i] = exp((targetDensity - activeDutyCycles_[i]) * boostStrength_);
				*/
				}
			}

			///This method increases the permanence values of synapses of columns whose
			//activity level has been too low.Such columns are identified by having an
			//overlap duty cycle that drops too much below those of their peers.The
			//permanence values for such columns are increased.
			template <typename P>
			void bump_up_weak_columns(
				Column<P>& column,
				const Dynamic_Param& param,
				const float overlap_duty_cycle,
				const float min_overlap_duty_cycle)
			{
				if (false) log_INFO("SP:bump_up_weak_columns: column ", column.id, "; overlap_duty_cycle = ", overlap_duty_cycle, "; min_overlap_duty_cycle = ", min_overlap_duty_cycle, ".\n");

				if (overlap_duty_cycle < min_overlap_duty_cycle) // the provided column is a weak column
				{
					if (false) log_INFO("SP:bump_up_weak_columns: column ", column.id, " is bumped up; overlap_duty_cycle = ", overlap_duty_cycle, "; min_overlap_duty_cycle = ", min_overlap_duty_cycle, ".\n");

					for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
					{
						if (column.pd_synapse_permanence[synapse_i] < param.SP_PD_CONNECTED_THRESHOLD)
						{
							column.pd_synapse_permanence[synapse_i] += param.SP_PD_PERMANENCE_INC_WEAK;
						}
					}
				}
			}

			template <typename P>
			void update_inhibition_radius(
				const Layer<P>& layer,
				Column<P>& column)
			{
				if (P::SP_GLOBAL_INHIBITION)
				{
					//TODO
				}
				else
				{
					//TODO
				}
			}

			template <typename P>
			void update_min_duty_cycles(
				Layer<P>& layer,
				const int column_i)
			{
				if (P::SP_GLOBAL_INHIBITION)
				{	
					//Real maxOverlapDutyCycles = *max_element(overlapDutyCycles_.begin(), overlapDutyCycles_.end());
					//fill(minOverlapDutyCycles_.begin(), minOverlapDutyCycles_.end(), minPctOverlapDutyCycles_ * maxOverlapDutyCycles);

					if (column_i == 0)
					{
						const float max_overlap_duty_cycles = max(layer.overlap_duty_cycles);
						layer.min_overlap_duty_cycles[column_i] = P::SP_MIN_PCT_OVERLAP_DUTY_CYCLES * max_overlap_duty_cycles;
						if (false) log_INFO_DEBUG("SP:update_min_duty_cycles: min_overlap_duty_cycles of column ", column_i, " = ", layer.min_overlap_duty_cycles[column_i], ".");
					}
					else
					{
						layer.min_overlap_duty_cycles[column_i] = layer.min_overlap_duty_cycles[0];
					}
				}
				else
				{
					::tools::log::log_ERROR("SP:update_min_duty_cycles: not implemented yet.");
				}
			}

			template <typename P>
			bool is_update_round(const Layer<P>& layer)
			{
				return (layer.iteration_num % P::SP_DUTY_CYCLE_PERIOD) == 0;
			}
		}

		template <bool LEARN, typename P>
		void compute_sp(
			const Bitset_Compact<P::N_SENSORS>& sensor_activity,
			Layer<P>& layer,
			const Dynamic_Param& param,
			//out
			Bitset<P::N_COLUMNS>& active_columns)
		{
			//local variables
			auto overlap_local = std::vector<int>(P::N_COLUMNS);
			auto boosted_overlap_local = std::vector<float>(P::N_COLUMNS);

			#if _DEBUG
			if (false) log_INFO("SP:compute_sp: sensor activity IN:", print::print_sensor_activity(sensor_activity, P::N_SENSORS_DIM1));
			#endif

			layer.iteration_num++;
			if (LEARN) layer.iteration_learn_num++;

			priv::calc_overlap::d(layer, param, sensor_activity, overlap_local);

			// update the boost factors
			//#pragma ivdep // ignore write after write dependency in rand_float
			for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
			{
				auto& column = layer[column_i];
				const int overlap = overlap_local[column_i];

				//add a small random number seems to improve learning speed
				const float r = random::rand_float(0.1f, column.random_number);
				boosted_overlap_local[column_i] = ((LEARN) ? (overlap * column.boost_factor) : overlap) + r;
			}

			#if _DEBUG
			if (false) log_INFO("SP:compute_sp: overlap:", print::print_int_array(overlap_local, P::N_COLUMNS));
			if (false) log_INFO("SP:compute_sp: boosted_overlap:", print::print_float_array(boosted_overlap_local, P::N_COLUMNS));
			#endif

			priv::inhibit_columns::d<P>(boosted_overlap_local, param, active_columns);

			if (LEARN)
			{
				priv::update_synapses::d(layer, param, active_columns, sensor_activity);

				if (true)
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						auto& column = layer[column_i];
						priv::update_duty_cycles(layer, column_i, overlap_local[column_i], active_columns[column_i]);
						priv::update_boost_factors(layer, column, param, layer.active_duty_cycles[column_i]);

						if (false)
						{
							priv::bump_up_weak_columns(column, param, layer.overlap_duty_cycles[column_i], layer.min_overlap_duty_cycles[column_i]);
							if (priv::is_update_round(layer))
							{
								priv::update_inhibition_radius(layer, column);
								priv::update_min_duty_cycles(layer, column_i);
							}
						}
					}
					//TODO: 1] update layer_boost(column_i)
					//TODO: 2] increase permanence if input overlap is small
				}
			}
			
			#if _DEBUG
			if (false) log_INFO("SP:compute_sp: active columns OUT:", print::print_bitset(active_columns));
			#endif
		}
	}
}
