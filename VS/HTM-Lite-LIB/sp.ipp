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

#include "..\Spike-Tools-Lib\log.ipp"
#include "..\Spike-Tools-Lib\assert.ipp"

#include "parameters.ipp"
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
					typename Layer<P>::Active_Columns& active_columns)
				{
					for (int column_i1 = 0; column_i1 < P::N_COLUMNS; ++column_i1)
					{
						const float oa = boosted_overlap[column_i1];
						int sum = 0;

						for (int column_i2 = 0; column_i2 < P::N_COLUMNS; ++column_i2)
						{
							sum += (boosted_overlap[column_i2] > oa);
						}
						active_columns.set(column_i1, sum < inhibition_top);
					}
				}

				// faster than ref1
				template <typename P>
				void active_columns_ref2(
					const std::vector<float>& boosted_overlap, //N_COLUMNS
					const int inhibition_top,
					typename Layer<P>::Active_Columns& active_columns)
				{
					active_columns.reset();

					for (int i = 0; i < inhibition_top; ++i)
					{
						float f = -10;
						int best_i = -1;

						for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							if ((boosted_overlap[column_i] > f) && !active_columns.get(column_i))
							{
								f = boosted_overlap[column_i];
								best_i = column_i;
							}
						}
						if (best_i != -1) active_columns.set(best_i);
					}
					#if _DEBUG
					Layer<P>::Active_Columns active_columns2;
					active_columns_ref1(boosted_overlap, inhibition_top, active_columns2);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2.get(column_i) != active_columns.get(column_i))
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
					typename Layer<P>::Active_Columns& active_columns)
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
							active_columns.set(best_i);
							tmp[best_i] = -1;
						}
					}

					#if _DEBUG
					Layer<P>::Active_Columns active_columns2;
					active_columns_ref1(boosted_overlap, inhibition_top, active_columns2);
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2.get(column_i) != active_columns.get(column_i))
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
					typename Layer<P>::Active_Columns& active_columns)
				{
					active_columns.clear_all();
					std::vector<float> boosted_overlap_copy = std::vector<float>(boosted_overlap);

					std::nth_element(boosted_overlap_copy.begin(), boosted_overlap_copy.begin() + inhibition_top, boosted_overlap_copy.end(), std::greater<float>());
					const float nth = boosted_overlap_copy[inhibition_top];

					for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						if (boosted_overlap[column_i] > nth) active_columns.set(column_i, true);
					}

					#if _DEBUG
					if (false) // it seems ref1 does not always yield the same results as ref4
					{
						Layer<P>::Active_Columns active_columns2;
						active_columns_ref1<P>(boosted_overlap, inhibition_top, active_columns2);
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i) if (active_columns2.get(column_i) != active_columns.get(column_i))
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
					typename Layer<P>::Active_Columns& active_columns)
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
				namespace synapse_backward
				{
					//Calculate overlap iterator over sensors
					template <typename P>
					void calc_overlap_sb_ref(
						const Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Sensors& active_sensors,
						//out
						std::vector<int>& overlaps) //size = P::N_COLUMNS
					{
						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							if (active_sensors.get(sensor_i))
							{
								const auto& destination_columns = layer.sp_pd_destination_column_sb[sensor_i];
								const auto& permanences = layer.sp_pd_synapse_permanence_sb[sensor_i];

								for (auto i = 0; i < layer.sp_pd_synapse_count_sb[sensor_i]; ++i)
								{
									if (permanences[i] > P::SP_PD_PERMANENCE_THRESHOLD)
									{
										const auto column_i = destination_columns[i];
										overlaps[column_i]++; //deadly scatter here!
									}
								}
							}
						}
					}
				}
				namespace synapse_forward
				{
					// Why is this an inefficient algorithm: for example, assume 1024 columns,
					// each column has 256 synapses and there are 400 sensor cells. The 400 cells 
					// are stored in 40 bytes, which fits in one L1 cache line (which is 64 bytes = 512 bits),
					// every columns is going to read 256 bytes, that is, 256KB will be read in total, 
					// thus every byte (of the 40) will on average be read 6553 times! Since only 5 percent 
					// of the sensors are active, it may be 20 times faster to do a scatter instead of a gather.

					template <typename P>
					void calc_overlap_sf_ref(
						const Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Sensors& active_sensors,
						//out
						std::vector<int>& overlaps) //size = P::N_COLUMNS
					{
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							const auto& permanence = layer.sp_pd_synapse_permanence_sf[column_i];
							const auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];

							int overlap = 0;
							for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
							{
								if (permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
								{
									const auto origin_sensor = synapse_origin[synapse_i];
									overlap += active_sensors.get(origin_sensor); // deadly gather here!
								}
							}
							if (false) log_INFO_DEBUG("SP:calc_overlap_sf_ref: column ", column_i, " has overlap = ", overlap, ".\n");

							overlaps[column_i] = (overlap < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap;
						}
					}

					__m512i get_sensors_epi32(
						const __mmask16 mask,
						const __m512i origin_epi32,
						const void * active_sensors_ptr)
					{
						const __m512i int_addr = _mm512_srli_epi32(origin_epi32, 5);
						const __m512i pos_in_int = _mm512_and_epi32(origin_epi32, _mm512_set1_epi32(0b11111));
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
						const __m512i origin_epi16,
						const __m512i active_sensors)
					{
						const __m512i short_addr = _mm512_srli_epi16(origin_epi16, 4);
						const __m512i pos_in_int = _mm512_and_si512(origin_epi16, _mm512_set1_epi16(0b1111));
						const __m512i sensor_int = _mm512_mask_permutexvar_epi16(_mm512_setzero_epi32(), mask, short_addr, active_sensors);
						const __m512i sensors = _mm512_srlv_epi16(sensor_int, pos_in_int);
						return _mm512_and_si512(sensors, _mm512_set1_epi16(1));
					}

					template <typename P>
					void calc_overlap_sf_avx512(
						const Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Sensors& active_sensors,
						//out 
						std::vector<int>& overlaps) //size = P::N_COLUMNS
					{
						const __m512i connected_threshold_epi8 = _mm512_set1_epi8(P::SP_PD_PERMANENCE_THRESHOLD);
						auto active_sensors_ptr = active_sensors.data();
						const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_permanence_sf[column_i].data());
							auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_origin_sensor_sf[column_i].data());

							__m512i overlap_epi32 = _mm512_setzero_epi32();

							for (int block = 0; block < n_blocks; ++block)
							{
								//const __m512i permanence_epi8 = _mm512_stream_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values
								const __m512i permanence_epi8 = _mm512_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values

								const __mmask64 connected_mask_64 = _mm512_cmpgt_epi8_mask(permanence_epi8, connected_threshold_epi8);
								for (int i = 0; i < 4; ++i)
								{
									const __mmask16 mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
									if (mask_16 != 0)
									{
										//const __m512i origin_epi32 = _mm512_stream_load_si512(&origin_epi32_ptr[(block * 4) + i]);
										const __m512i origin_epi32 = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + i]);
										overlap_epi32 = _mm512_add_epi32(overlap_epi32, get_sensors_epi32(mask_16, origin_epi32, active_sensors_ptr));
									}
								}
							}
							const int overlap_int = _mm512_reduce_add_epi32(overlap_epi32);
							if (false) log_INFO_DEBUG("SP:calc_overlap_avx512: column ", column_i, " has overlap = ", overlaps[column_i], ".\n");
							overlaps[column_i] = (overlap_int < P::SP_STIMULUS_THRESHOLD) ? 0 : overlap_int;
						}
						#if _DEBUG
						std::vector<int> overlaps_ref = std::vector<int>(P::N_COLUMNS);
						calc_overlap_sf_ref(layer, param, active_sensors, overlaps_ref);
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
					void calc_overlap_avx512_sf_small_epi32(
						const Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Sensors& active_sensors,
						//out 
						std::vector<int>& overlaps)
					{
						assert_msg(P::N_SENSORS < 512, "ERROR: calc_overlap_avx512_small: N_SENSORS is larger than 512");

						const __m512i connected_threshold_epi8 = _mm512_set1_epi8(P::SP_PD_PERMANENCE_THRESHOLD);
						std::array<int, 16> t = { 0 };
						for (int i = 0; i < active_sensors.N_BLOCKS; ++i)
						{
							t[i] = active_sensors._data[i];
						}
						const __m512i active_sensors_simd = _mm512_load_epi32(t.data());
						const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_permanence[column_i].data());
							auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_origin[column_i].data());

							__m512i overlap_epi32 = _mm512_setzero_epi32();

							for (int block = 0; block < n_blocks; ++block)
							{
								//const __m512i permanence_epi8 = _mm512_stream_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values
								const __m512i permanence_epi8 = _mm512_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values
								const __mmask64 connected_mask_64 = _mm512_cmpgt_epi8_mask(permanence_epi8, connected_threshold_epi8);
								for (int i = 0; i < 4; ++i)
								{
									const __mmask16 mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
									//if (mask_16 != 0) // slower with this check
									{
										//const __m512i origin_epi32 = _mm512_stream_load_si512(&origin_epi32_ptr[(block * 4) + i]);
										const __m512i origin_epi32 = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + i]);
										overlap_epi32 = _mm512_add_epi32(overlap_epi32, get_sensors_epi32(mask_16, origin_epi32, active_sensors_simd));
									}
								}
							}
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
							if (overlap_ref != overlap_avx512) log_ERROR("SP:calc_overlap_avx512_small_epi32:: UNEQUAL: column ", column_i, "; overlap ref ", overlap_ref, " != avx512 ", overlap_avx512, ".\n");
						}
						#endif
					}
					// not much faster than calc_overlap_avx512_small_epi32
					template <typename P>
					void calc_overlap_avx512_sf_small_epi16(
						const Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Sensors& active_sensors,
						//out 
						std::vector<int>& overlaps)
					{
						assert_msg(P::N_VISIBLE_SENSORS < 512, "ERROR: calc_overlap_avx512_small: N_SENSORS is larger than 512");

						const __m512i connected_threshold_epi8 = _mm512_set1_epi8(P::SP_PD_PERMANENCE_THRESHOLD);

						std::array<int, 16> t = { 0 };
						for (int i = 0; i < active_sensors.N_BLOCKS; ++i)
						{
							t[i] = active_sensors._data[i];
						}
						const __m512i active_sensors_simd = _mm512_load_epi32(t.data());
						const int n_blocks = htm::tools::n_blocks_64(P::SP_N_PD_SYNAPSES);

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_permanence_sf[column_i].data());
							auto origin_epi32_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_synapse_origin_sensor_sf[column_i].data());

							__m512i overlap_epu16_AB = _mm512_setzero_si512(); // contains 32 overlap values of 16bits
							__m512i overlap_epu16_CD = _mm512_setzero_si512(); // contains 32 overlap values of 16bits

							for (int block = 0; block < n_blocks; ++block)
							{
								//const __m512i permanence_epi8 = _mm512_stream_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values
								const __m512i permanence_epi8 = _mm512_load_si512(&permanence_epi8_ptr[block]); //load 64 permanence values
								const __mmask64 connected_mask_64 = _mm512_cmpgt_epi8_mask(permanence_epi8, connected_threshold_epi8);

								const __mmask32 mask_32_A = static_cast<__mmask32>(connected_mask_64 >> (0 * 32));
								const __m512i origin_epi32_A = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + 0]);
								const __m512i origin_epi32_B = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + 1]);
								const __m512i origin_epu16_A = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_A));
								const __m512i origin_epu16_B = _mm512_castsi256_si512(_mm512_cvtepi32_epi16(origin_epi32_B));
								const __m512i origin_epu16_AB = _mm512_shuffle_i64x2(origin_epu16_A, origin_epu16_B, 0b01000100);
								overlap_epu16_AB = _mm512_adds_epu16(overlap_epu16_AB, get_sensors_epi16(mask_32_A, origin_epu16_AB, active_sensors_simd));

								const __mmask32 mask_32_B = static_cast<__mmask32>(connected_mask_64 >> (1 * 32));
								const __m512i origin_epi32_C = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + 2]);
								const __m512i origin_epi32_D = _mm512_load_si512(&origin_epi32_ptr[(block * 4) + 3]);
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
						calc_overlap_sf_ref(layer, param, active_sensors, overlaps_ref);
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							const int overlap_ref = overlaps_ref[column_i];
							const int overlap_avx512 = overlaps[column_i];
							if (overlap_ref != overlap_avx512) log_ERROR("SP:calc_overlap_avx512_small_epi16:: UNEQUAL: column ", column_i, "; overlap ref ", overlap_ref, " != avx512 ", overlap_avx512, ".\n");
						}
						#endif
					}
				}

				template <typename P>
				void d(
					const Layer<P>& layer,
					const Dynamic_Param& param,
					const typename Layer<P>::Active_Sensors& active_sensors,
					//out
					std::vector<int>& overlaps) //size = P::N_COLUMNS
				{
					if (P::SP_SYNAPSE_FORWARD)
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) synapse_forward::calc_overlap_sf_ref(layer, param, active_sensors, overlaps);
						if (architecture_switch(P::ARCH) == arch_t::AVX512)
						{
							if (P::N_SENSORS < 512)
							{
								synapse_forward::calc_overlap_avx512_sf_small_epi16(layer, param, active_sensors, overlaps);
								//indexed_by_column::calc_overlap_avx512_ic_small_epi32(layer, param, active_sensors, overlaps);
							}
							else
							{
								synapse_forward::calc_overlap_sf_avx512(layer, param, active_sensors, overlaps);
							}
						}
					}
					else
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) synapse_backward::calc_overlap_sb_ref(layer, param, active_sensors, overlaps);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) synapse_backward::calc_overlap_sb_ref(layer, param, active_sensors, overlaps);
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

				namespace synapse_backward
				{
					//Update synapses iterator over sensors.
					template <typename P>
					void update_synapses_is_ref(
						Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Columns& active_columns,
						const typename Layer<P>::Active_Sensors& active_sensors)
					{
						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							const auto& destination_columns = layer.sp_pd_destination_column_sb[sensor_i];
							auto& permanence = layer.sp_pd_synapse_permanence_sb[sensor_i];
							const bool sensor_is_active = active_sensors.get(sensor_i);

							for (auto synapse_i = 0; synapse_i < layer.sp_pd_synapse_count_sb[sensor_i]; ++synapse_i)
							{
								const auto column_i = destination_columns[synapse_i];
								if (active_columns.get(column_i))
								{
									const Permanence old_permanence = permanence[synapse_i];
									const Permanence new_permanence = (sensor_is_active)
										? add_saturate(old_permanence, param.SP_PD_PERMANENCE_INC)
										: sub_saturate(old_permanence, param.SP_PD_PERMANENCE_DEC);

									if (false) log_INFO_DEBUG("SP:update_synapses_is_ref: inc: column ", column_i, "; synpase ", synapse_i, ": old permanence ", old_permanence, "; new permanence = ", new_permanence, ".");
									permanence[synapse_i] = new_permanence;
								}
							}
						}
					}

					template <typename P>
					void update_synapses_is_avx512(
						Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Columns& active_columns,
						const typename Layer<P>::Active_Sensors& active_sensors)
					{
						//TODO: Code has a bug!

						#if _DEBUG
						using vector_type = std::vector<Permanence, types::priv::Allocator<Permanence>>;
						auto permanence_org = std::vector<vector_type>(P::N_SENSORS);
						auto permanence_ref = std::vector<vector_type>(P::N_SENSORS);

						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							permanence_org[sensor_i] = vector_type(layer.sp_pd_synapse_permanence[sensor_i]);
						}
						update_synapses_ref(layer, param, active_columns, active_sensors);

						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							permanence_ref[sensor_i] = layer.sp_pd_synapse_permanence[sensor_i];
							layer.sp_pd_synapse_permanence[sensor_i] = permanence_org[sensor_i];
						}
						#endif

						auto active_columns_ptr = reinterpret_cast<const __m512i *>(active_columns.data());

						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							auto destination_columns_epi32_ptr = reinterpret_cast<const __m512i *>(layer.sp_pd_destination_column_sb[sensor_i].data());
							auto permanence_epi8_ptr = reinterpret_cast<__m128i *>(layer.sp_pd_synapse_permanence_sb[sensor_i].data());

							const __m128i inc_epi8 = _mm_set1_epi8(active_sensors.get(sensor_i) ? param.SP_PD_PERMANENCE_INC : -param.SP_PD_PERMANENCE_DEC);

							const int n_blocks = tools::n_blocks_16(layer.sp_pd_synapse_count[sensor_i]);
							for (int block = 0; block < n_blocks; ++block)
							{
								const __m512i destination_columns_epi32 = destination_columns_epi32_ptr[block];

								const __m512i int_addr = _mm512_srli_epi32(destination_columns_epi32, 2);
								const __m512i sensor_int = _mm512_i32gather_epi32(int_addr, active_columns_ptr, 4);
								const __m512i byte_pos_in_int = _mm512_and_epi32(destination_columns_epi32, _mm512_set1_epi32(0b11));
								const __m512i bit_mask_epi32 = _mm512_sllv_epi32(_mm512_set1_epi32(1), byte_pos_in_int);
								const __mmask16 mask_16 = _mm512_cmpeq_epi32_mask(_mm512_and_epi32(sensor_int, bit_mask_epi32), bit_mask_epi32);

								const __m128i old_permanence = permanence_epi8_ptr[block];
								permanence_epi8_ptr[block] = _mm_mask_adds_epu8(old_permanence, mask_16, old_permanence, inc_epi8);
							}
						}

						#if _DEBUG
						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							const auto& ref = permanence_ref[sensor_i];
							const auto& org = permanence_org[sensor_i];
							for (auto synapse_i = 0; synapse_i < layer.sp_pd_synapse_count[sensor_i]; ++synapse_i)
							{
								if (ref[synapse_i] != org[synapse_i])
								{
									log_ERROR("TP:update_synapses_avx512:: UNEQUAL permanence for synapse_i ", synapse_i, ": ref ", static_cast<int>(ref[synapse_i]), "; avx512 ", static_cast<int>(org[synapse_i]));
								}
							}
						}
						#endif
					}
				}
				namespace indexed_by_column
				{
					//Update synapses iterator over column. When iterating over the columns, the sparsity of the culumns yields in less memory access.
					template <typename P>
					void update_synapses_ic_ref(
						Layer<P>& layer,
						const Dynamic_Param& param,
						const typename Layer<P>::Active_Columns& active_columns,
						const typename Layer<P>::Active_Sensors& active_sensors)
					{
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							if (active_columns.get(column_i))
							{
								const auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];
								auto& permanence = layer.sp_pd_synapse_permanence_sf[column_i];

								for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
								{
									const auto sensor_i = synapse_origin[synapse_i];
									const bool sensor_is_active = active_sensors.get(sensor_i);

									const Permanence old_permanence = permanence[synapse_i];
									const Permanence new_permanence = (sensor_is_active)
										? add_saturate(old_permanence, param.SP_PD_PERMANENCE_INC)
										: sub_saturate(old_permanence, param.SP_PD_PERMANENCE_DEC);

									if (false) log_INFO_DEBUG("SP:update_synapses_ic_ref: inc: column ", column_i, "; synpase ", synapse_i, ": old permanence ", old_permanence, "; new permanence = ", new_permanence, ".");
									permanence[synapse_i] = new_permanence;
								}
							}
						}
					}
				}

				template <typename P>
				void d(
					Layer<P>& layer,
					const Dynamic_Param& param,
					const typename Layer<P>::Active_Columns& active_columns,
					const typename Layer<P>::Active_Sensors& active_sensors)
				{
					if (P::SP_SYNAPSE_FORWARD)
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) indexed_by_column::update_synapses_ic_ref(layer, param, active_columns, active_sensors);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) indexed_by_column::update_synapses_ic_ref(layer, param, active_columns, active_sensors);
					}
					else
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) synapse_backward::update_synapses_is_ref(layer, param, active_columns, active_sensors);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) synapse_backward::update_synapses_is_ref(layer, param, active_columns, active_sensors);
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
				layer.sp_overlap_duty_cycles[column_i] = ((layer.sp_overlap_duty_cycles[column_i] * (period - 1) + ((overlap > 0) ? 1 : 0)) / period);
				layer.sp_active_duty_cycles[column_i] = ((layer.sp_active_duty_cycles[column_i] * (period - 1) + ((active_column) ? 1 : 0)) / period);
			}

			template <typename P>
			void update_boost_factors(
				Layer<P>& layer,
				const int column_i,
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
					if (false) log_INFO_DEBUG("SP:update_boost_factors: column ", column_i, "; old boost factor ", layer.boost_factor[column_i], "; new boost factor ", new_boost_factor, ".\n");
					layer.boost_factor[column_i] = new_boost_factor;
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
				Layer<P>& layer,
				const Dynamic_Param& param,
				const std::vector<float>& overlap_duty_cycle,
				const std::vector<float>& min_overlap_duty_cycle)
			{
				//if (false) log_INFO("SP:bump_up_weak_columns: column ", column_i, "; overlap_duty_cycle = ", overlap_duty_cycle, "; min_overlap_duty_cycle = ", min_overlap_duty_cycle, ".\n");

				for (auto sensor_i = 0; sensor_i < P::SP_N_PD_SYNAPSES; ++sensor_i)
				{
					const auto& destination = layer.sp_pd_destination_column_sb[sensor_i];
					auto& permanance = layer.sp_pd_synapse_permanence_sb[sensor_i];

					for (auto synapse_i = 0; synapse_i < layer.sp_pd_synapse_count_sb[sensor_i]; ++synapse_i)
					{
						const auto column_i = destination[synapse_i];
						if (overlap_duty_cycle[column_i] < min_overlap_duty_cycle[column_i]) // the provided column is a weak column
						{
							if (permanance[synapse_i] <= P::SP_PD_PERMANENCE_THRESHOLD)
							{
								permanance[synapse_i] += param.SP_PD_PERMANENCE_INC_WEAK;
							}
						}
					}
				}
			}

			template <typename P>
			void update_inhibition_radius(
				const Layer<P>& layer)
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
				Layer<P>& layer)
			{
				if (P::SP_GLOBAL_INHIBITION)
				{
					//Real maxOverlapDutyCycles = *max_element(overlapDutyCycles_.begin(), overlapDutyCycles_.end());
					//fill(minOverlapDutyCycles_.begin(), minOverlapDutyCycles_.end(), minPctOverlapDutyCycles_ * maxOverlapDutyCycles);

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						if (column_i == 0)
						{
							const float max_overlap_duty_cycles = max(layer.sp_overlap_duty_cycles);
							layer.sp_min_overlap_duty_cycles[column_i] = P::SP_MIN_PCT_OVERLAP_DUTY_CYCLES * max_overlap_duty_cycles;
							if (false) log_INFO_DEBUG("SP:update_min_duty_cycles: min_overlap_duty_cycles of column ", column_i, " = ", layer.sp_min_overlap_duty_cycles[column_i], ".");
						}
						else
						{
							layer.sp_min_overlap_duty_cycles[column_i] = layer.sp_min_overlap_duty_cycles[0];
						}
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
			const typename Layer<P>::Active_Sensors& active_sensors,
			Layer<P>& layer,
			const Dynamic_Param& param,
			//out
			typename Layer<P>::Active_Columns& active_columns)
		{
			//local variables
			auto overlap_local = std::vector<int>(P::N_COLUMNS);
			auto boosted_overlap_local = std::vector<float>(P::N_COLUMNS);

			layer.iteration_num++;
			if (LEARN) layer.iteration_learn_num++;

			priv::calc_overlap::d(layer, param, active_sensors, overlap_local);

			// update the boost factors
			//#pragma ivdep // ignore write after write dependency in rand_float
			for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
			{
				const int overlap = overlap_local[column_i];

				//add a small random number seems to improve learning speed
				const float r = random::rand_float(0.1f, layer.random_number[column_i]);
				boosted_overlap_local[column_i] = ((LEARN) ? (overlap * layer.boost_factor[column_i]) : overlap) + r;
			}

			#if _DEBUG
			if (false) log_INFO("SP:compute_sp: overlap:", print::print_int_array(overlap_local, P::N_COLUMNS));
			if (false) log_INFO("SP:compute_sp: boosted_overlap:", print::print_float_array(boosted_overlap_local, P::N_COLUMNS));
			#endif

			priv::inhibit_columns::d<P>(boosted_overlap_local, param, active_columns);

			if (LEARN)
			{
				priv::update_synapses::d(layer, param, active_columns, active_sensors);

				if (true)
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						priv::update_duty_cycles(layer, column_i, overlap_local[column_i], active_columns.get(column_i));
						priv::update_boost_factors(layer, column_i, param, layer.sp_active_duty_cycles[column_i]);
					}
					if (false)
					{
						priv::bump_up_weak_columns(layer, param, layer.sp_overlap_duty_cycles, layer.sp_min_overlap_duty_cycles);
						if (priv::is_update_round(layer))
						{
							priv::update_inhibition_radius(layer);
							priv::update_min_duty_cycles(layer);
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