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
#include <map>
#include <set>

#include "..\Spike-Tools-Lib\log.ipp"
#include "..\Spike-Tools-Lib\assert.ipp"
#include "..\Spike-Tools-Lib\random.ipp"

#include "parameters.ipp"
#include "print.ipp"
#include "types.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Temporal Memory/Pooler
	namespace tp
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;
		using namespace htm::types;
		using namespace htm::tools;

		//HTM Temporal Pooler private methods
		namespace priv
		{
			namespace activate_cells
			{
				template <typename P>
				void select_delay_and_cell_to_learn_on(
					//in
					Layer<P>& layer, // column cannot be readonly due to the random number generator
					const int column_i,
					const int segment_i,
					const Layer<P>::Winner_Cells& winner_cells,
					const int select_size,
					//out
					std::vector<int>& selected_delay_and_cells) // assumes that selected_cells has sufficient capacity (size >= select_size)
				{
					assert_msg(segment_i < layer.dd_segment_count[column_i], "TP:select_cell_to_learn_on: segment_i=", segment_i + " is too large. dd_segment_count=", layer.dd_segment_count[column_i]);

					//pick global cell ids from winner_cells that do not already
					//have a pathway to segment in column, if not enough winner cells
					//exists, take non winner cells

					const std::vector<int>& delay_and_winner_cells_vector = winner_cells.get_sparse_history();
					const int n_winner_cells = static_cast<int>(delay_and_winner_cells_vector.size());

					unsigned int random_number = layer.random_number[column_i];

					const auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];
					const auto n_synapses = layer.dd_synapse_count[column_i][segment_i];

					int selected_cells_count = 0;
					int winner_cell_i = 0;

					while (selected_cells_count < select_size)
					{
						int delay_and_cell_id;
						if (winner_cell_i < n_winner_cells) // get a random cell from winner cells
						{
							const int random_i = random::rand_int32(0, n_winner_cells - 1, random_number);
							delay_and_cell_id = delay_and_winner_cells_vector[random_i];
							if (false) log_INFO_DEBUG("TP:select_cell_to_learn_on: column ", column_i, "; segment_i ", segment_i, "; global_cell_id = ", get_global_cell_id(delay_and_cell_id), "; current_random_number = ", static_cast<int>(random::priv::current_random_number), "; min = ", winner_cell_i, "; max = ", n_winner_cells - 1, "\n");
							winner_cell_i++;
						}
						else // not enough winner cells, get a random cell
						{
							const int global_cell_id = random::rand_int32(0, P::N_CELLS - 1, random_number);
							const int delay = 1;
							delay_and_cell_id = create_delay_and_cell_id(global_cell_id, delay);
						}

						bool already_present = false;

						for (auto synapse_i = 0; synapse_i < selected_cells_count; ++synapse_i)
						{// search the already selected cells if it is already present
							if (selected_delay_and_cells[synapse_i] == delay_and_cell_id)
							{
								already_present = true; // found it
								break;
							}
						}
						if (!already_present)
						{
							for (auto synapse_i = 0; synapse_i < n_synapses; ++synapse_i)
							{// search the existing synapses whether the cell already has a pathway to this segment
								if (get_global_cell_id(dd_synapse_delay_origin_segment[synapse_i]) == delay_and_cell_id)
								{
									already_present = true; // found it
									break;
								}
							}
						}
						if (!already_present) // add the cell if not already present
						{
							selected_delay_and_cells[selected_cells_count] = delay_and_cell_id;
							selected_cells_count++;
						}
					}

					layer.random_number[column_i] = random_number;
				}

				//Change the permanence values of DD synapses on the provided segment in the provided column
				namespace adapt_segment
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
					void adapt_segment_ref(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Permanence permanence_dec)
					{
						const auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];
						auto& dd_synapse_permanence_segment = layer.dd_synapse_permanence[column_i][segment_i];

						for (auto synapse_i = 0; synapse_i < layer.dd_synapse_count[column_i][segment_i]; ++synapse_i)
						{
							const auto delay_cell_id = dd_synapse_delay_origin_segment[synapse_i];
							const int global_cell_id = get_global_cell_id(delay_cell_id);
							const int delay = get_delay(delay_cell_id);
							const bool b = active_cells.get(global_cell_id, delay);
							if (!b) dd_synapse_permanence_segment[synapse_i] = sub_saturate(dd_synapse_permanence_segment[synapse_i], permanence_dec);
						}
					}

					template <typename P>
					void adapt_segment_ref(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Permanence permanence_inc,
						const Permanence permanence_dec)
					{
						const auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];
						auto& dd_synapse_permanence_segment = layer.dd_synapse_permanence[column_i][segment_i];
						
						for (auto synapse_i = 0; synapse_i < layer.dd_synapse_count[column_i][segment_i]; ++synapse_i)
						{
							const Permanence old_permanence = dd_synapse_permanence_segment[synapse_i];
							if (old_permanence > P::TP_DD_CONNECTED_THRESHOLD)
							{
								const auto delay_and_cell_id = dd_synapse_delay_origin_segment[synapse_i];
								const int global_cell_id = get_global_cell_id(delay_and_cell_id);
								const int delay = get_delay(delay_and_cell_id);
								const bool b = active_cells.get(global_cell_id, delay); // deadly gather here!

								dd_synapse_permanence_segment[synapse_i] = (b)
									? add_saturate(old_permanence, permanence_inc)
									: sub_saturate(old_permanence, permanence_dec);
							}
						}
					}

					__mmask16 get_sensors_mask(
						const __mmask16 mask,
						const __m512i delay_and_origin,
						const void * active_sensors_ptr)
					{
						const __m512i global_cell_id = _mm512_and_epi32(delay_and_origin, _mm512_set1_epi32(0x1FFFFFFF));
						const __m512i byte_pos_in_int = _mm512_and_epi32(delay_and_origin, _mm512_set1_epi32(0b11));
						const __m512i delay = _mm512_srli_epi32(delay_and_origin, 32 - 3);

						const __m512i int_addr = _mm512_srli_epi32(global_cell_id, 2);
						const __m512i sensor_int = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), mask, int_addr, active_sensors_ptr, 4);
						const __m512i pos_in_int = _mm512_slli_epi32(byte_pos_in_int, 3);
						const __m512i delay_shift = _mm512_add_epi32(delay, pos_in_int);
						const __m512i delay_mask = _mm512_sllv_epi32(_mm512_set1_epi32(1), delay_shift);
						return _mm512_cmpeq_epi32_mask(_mm512_and_epi32(sensor_int, delay_mask), delay_mask);
					}

					template <typename P>
					void adapt_segment_avx512(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Permanence permanence_inc,
						const Permanence permanence_dec)
					{
						#if _DEBUG
						const bool compare_to_ref = true;
						std::vector<Permanence, types::priv::Allocator> permanence_org;
						std::vector<Permanence, types::priv::Allocator> permanence_ref;
						if (compare_to_ref)
						{
							permanence_org = layer.dd_synapse_permanence[column_i][segment_i];
							adapt_segment_ref(layer, column_i, segment_i, active_cells, permanence_inc, permanence_dec);
							permanence_ref = layer.dd_synapse_permanence[column_i][segment_i];
							layer.dd_synapse_permanence[column_i][segment_i] = permanence_org;
						}
						#endif

						auto permanence_epi8_ptr = reinterpret_cast<__m512i *>(layer.dd_synapse_permanence[column_i][segment_i].data());
						auto delay_origin_epi32_ptr = reinterpret_cast<const __m512i *>(layer.dd_synapse_delay_origin[column_i][segment_i].data());
						auto active_cells_ptr = active_cells.data();

						const __m512i connected_threshold_epi8 = _mm512_set1_epi8(P::TP_DD_CONNECTED_THRESHOLD);
						const __m512i inc_epi8 = _mm512_set1_epi8(permanence_inc);
						const __m512i dec_epi8 = _mm512_set1_epi8(-permanence_dec);

						const int n_blocks = tools::n_blocks_64(layer.dd_synapse_count[column_i][segment_i]);

						for (int block = 0; block < n_blocks; ++block)
						{
							const __m512i old_permanence_epi8 = permanence_epi8_ptr[block]; //load 64 permanence values
							const __mmask64 connected_mask_64 = _mm512_cmpgt_epi8_mask(old_permanence_epi8, connected_threshold_epi8);
							__mmask64 active_cells_mask_64 = 0;
							{
								for (int i = 0; i < 4; ++i)
								{
									const __mmask16 connected_mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
									if (connected_mask_16 != 0)
									{
										const __mmask64 tmp_mask_64 = get_sensors_mask(connected_mask_16, delay_origin_epi32_ptr[(block * 4) + i], active_cells_ptr);
										active_cells_mask_64 |= tmp_mask_64 << (i * 16);
									}
								}
							}
							const __m512i inc_mask = _mm512_mask_blend_epi8(active_cells_mask_64, dec_epi8, inc_epi8);
							permanence_epi8_ptr[block] = _mm512_mask_adds_epi8(old_permanence_epi8, connected_mask_64, old_permanence_epi8, inc_mask);
						}

						#if _DEBUG
						if (compare_to_ref)
						{
							const auto& permanence_avx512 = layer.dd_synapse_permanence[column_i][segment_i];

							for (auto synapse_i = 0; synapse_i < layer.dd_synapse_count[column_i][segment_i]; ++synapse_i)
							{
								if (permanence_ref[synapse_i] != permanence_avx512[synapse_i])
								{
									log_ERROR("TP:adapt_segment_avx512:: UNEQUAL permanence for synapse_i ", synapse_i, ": ref ", static_cast<int>(permanence_ref[synapse_i]), "; avx512 ", static_cast<int>(permanence_avx512[synapse_i]));
								}
							}
						}
						#endif
					}

					template <typename P>
					void d(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Permanence permanence_inc,
						const Permanence permanence_dec)
					{
						assert_msg(segment_i < layer.dd_segment_count[column_i], "TP:grow_DD_synapses: segment_i=", segment_i + " is too large. dd_segment_count=", layer.dd_segment_count[column_i]);
						if (false) log_INFO_DEBUG("TP:adapt_segment: column ", column_i, "; segment_i ", segment_i);

						if (architecture_switch(P::ARCH) == arch_t::X64) return adapt_segment_ref(layer, column_i, segment_i, active_cells, permanence_inc, permanence_dec);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) return adapt_segment_avx512(layer, column_i, segment_i, active_cells, permanence_inc, permanence_dec);
					}

					template <typename P>
					void d(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Permanence permanence_dec)
					{
						assert_msg(segment_i < layer.dd_segment_count[column_i], "TP:grow_DD_synapses: segment_i=", segment_i, " is too large. dd_segment_count=", layer.dd_segment_count[column_i]);
						if (false) log_INFO_DEBUG("TP:adapt_segment: column ", column_i, "; segment_i ", segment_i);
						return adapt_segment_ref(layer, column_i, segment_i, active_cells, permanence_dec);
					}
				}

				//Grow new synapses on the provided segment in the provided column
				namespace grow_DD_synapses
				{
					template <typename P>
					void grow_DD_synapses_ref(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const int n_desired_new_synapses,
						const Layer<P>::Winner_Cells& winner_cells)
					{
						#if _DEBUG
						if (false) log_INFO("TP:grow_DD_synapses_ref: before: column ", column_i, "; segment_i ", segment_i, "; dd_synapses:", print::print_dd_synapses(layer, column_i));
						#endif
						assert_msg(segment_i < layer.dd_segment_count[column_i], "TP:grow_DD_synapses_ref: segment=", segment_i, " is too large; dd_segment_count=", layer.dd_segment_count[column_i]);

						if (n_desired_new_synapses <= 0) return; // nothing to do

						auto& dd_synapse_permanence_segment = layer.dd_synapse_permanence[column_i][segment_i];
						auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];
						assert_msg(dd_synapse_permanence_segment.size() == dd_synapse_delay_origin_segment.size(), "TP:grow_DD_synapses_ref: bug A.");

						const int old_size = layer.dd_synapse_count[column_i][segment_i];

						//find indices that will be overwritten with new values
						std::vector<int> indices_to_update = std::vector<int>(n_desired_new_synapses, -1);

						//log_INFO_DEBUG("TP:grow_DD_synapses_ref: Start: old_size=", old_size, "; n_desired_new_synapses=", n_desired_new_synapses);

						// a new synapse can either be on the place of an existing synapse that has a low permanence, recycled synapse, or
						// a new synapse can be a newly created synapse, we may need to allocate space for such a synapse and test if the 
						// total number of synapses is not above TP_N_DD_SEGMENTS_MAX 

						int k = 0; // number of newly created synapses
						int j = 0; // index of synapse that will be a new synapse
						int i = 0; // counter of new synapses (added + recycled)
						while (i < n_desired_new_synapses)
						{
							//log_INFO_DEBUG("TP:grow_DD_synapses_ref: i=", i, "; j=", j, "; k=", k);

							if (j < old_size)
							{
								if (dd_synapse_permanence_segment[j] <= 0)
								{ // recycle synapse with index j
									indices_to_update[i] = j;
									i++;
								}
							}
							else 
							{ // create a new synapse with index j
								indices_to_update[i] = j;
								i++;
								k++;
							}
							j++;
						}

						#if _DEBUG
						for (int i = 0; i < n_desired_new_synapses; ++i) assert_msg(indices_to_update[i] != -1, "TP:grow_DD_synapses_ref: bug D.");
						#endif

						int overflow = 0;
						if (k > 0)
						{ // test if we need to allocate space for new synapses
							if ((old_size + k) > P::TP_N_DD_SYNAPSES_MAX) 
							{// ajust k and j for the situation in which the segment is overflowing
								overflow = (old_size + k) - P::TP_N_DD_SYNAPSES_MAX;
								assert_msg(overflow >= 0, "TP:grow_DD_synapses_ref: overflow cannot be negative; overflow=", overflow);
								k -= overflow;
								j -= overflow;
							}
							if (k > 0)
							{
								const int old_capacity = static_cast<int>(dd_synapse_permanence_segment.size());
								assert_msg(old_size <= old_capacity, "TP:grow_DD_synapses_ref: bug B.");
								const int new_size = old_size + k;
								if (new_size > old_capacity)
								{
									dd_synapse_permanence_segment.resize(htm::tools::multiple_64(new_size), P::TP_DD_CONNECTED_THRESHOLD);
									dd_synapse_delay_origin_segment.resize(htm::tools::multiple_64(new_size), P::TP_DD_SYNAPSE_ORIGIN_INVALID);
								}
								layer.dd_synapse_count[column_i][segment_i] = new_size;

								assert_msg(new_size <= dd_synapse_permanence_segment.size(), "TP:grow_DD_synapses_ref: invalid new_size=", new_size, "; while dd_synapse_permanence_segment.size()=", dd_synapse_permanence_segment.size());
								assert_msg(new_size <= dd_synapse_delay_origin_segment.size(), "TP:grow_DD_synapses_ref: invalid new_size=", new_size, "; while dd_synapse_origin_segment.size()=", dd_synapse_delay_origin_segment.size());
								assert_msg(new_size <= P::TP_N_DD_SYNAPSES_MAX, "TP:grow_DD_synapses_ref: invalid new_size=", new_size, "; TP_N_DD_SYNAPSES_MAX=", P::TP_N_DD_SYNAPSES_MAX);
							}
						}

						const int n_exact_new_synapses = n_desired_new_synapses - overflow;
						assert_msg(n_exact_new_synapses <= n_desired_new_synapses, "TP:grow_DD_synapses_ref: n_exact_new_synapses=", n_exact_new_synapses, "; n_desired_new_synapses=", n_desired_new_synapses);

						if (n_exact_new_synapses <= 0)
						{
							if (false) log_INFO_DEBUG("TP:grow_DD_synapses_ref: ", column_i, "; segment_i ", segment_i, "; segment is full. n_desired_new_synapses=", n_desired_new_synapses, "; overflow=", overflow);
						}
						else
						{
							std::vector<int> selected_delay_and_cells(n_exact_new_synapses);
							select_delay_and_cell_to_learn_on(layer, column_i, segment_i, winner_cells, n_exact_new_synapses, selected_delay_and_cells);
							
							for (int i = 0; i < n_exact_new_synapses; ++i)
							{
								const int synapse_i = indices_to_update[i];
								dd_synapse_permanence_segment[synapse_i] = P::TP_DD_PERMANENCE_INIT;
								dd_synapse_delay_origin_segment[synapse_i] = selected_delay_and_cells[i];
							}
						}

						#if _DEBUG
						for (auto synapse_i = 0; synapse_i < layer.dd_synapse_count[column_i][segment_i]; ++synapse_i)
						{
							const auto origin = get_global_cell_id(dd_synapse_delay_origin_segment[synapse_i]);
							assert_msg(origin >= 0, "TP:grow_DD_synapses_ref: invalid origin ", origin);
							assert_msg(origin < P::N_CELLS, "TP:grow_DD_synapses_ref: invalid origin ", origin);
						}
						if (false) log_INFO("TP:grow_DD_synapses_ref: column ", column_i, "; segment_i ", segment_i, "; done replacing ", n_desired_new_synapses, " synapses.");
						if (false) log_INFO("TP:grow_DD_synapses_ref: after: column ", column_i, "; segment_i ", segment_i, "; dd_synapses:\n", print::print_dd_synapses(layer, column_i));
						#endif
					}

					template <typename P>
					void d(
						Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const int n_desired_new_synapses,
						const Layer<P>::Winner_Cells& winner_cells)
					{
						assert_msg(segment_i < layer.dd_segment_count[column_i], "TP:grow_DD_synapses: segment_i=", segment_i + " is too large. dd_segment_count=", layer.dd_segment_count[column_i]);
						grow_DD_synapses_ref(layer, column_i, segment_i, n_desired_new_synapses, winner_cells);
					}
				}

				template <typename P>
				int get_least_used_cell(
					Layer<P>& layer,
					const int column_i)
				{
					std::array<int, P::N_CELLS_PC> counter = { 0 };
					const auto& segment_destination = layer.dd_segment_destination[column_i];

					for (auto segment_i = 0; segment_i < layer.dd_segment_count[column_i]; ++segment_i)
					{
						counter[segment_destination[segment_i]]++;
					}

					unsigned int random_number = layer.random_number[column_i];
					int best_cell = 0;
					int selected_count = counter[0];
					int num_tied_cells = 1;

					for (auto cell_i = 1; cell_i < P::N_CELLS_PC; ++cell_i)
					{
						const int count = counter[cell_i];
						if (count < selected_count)
						{
							selected_count = counter[cell_i];
							best_cell = cell_i;
							num_tied_cells = 1;
						}
						else if (count == selected_count)
						{
							num_tied_cells++;
						}
					}

					if (num_tied_cells == 1)
					{
						return best_cell;
					}
					else
					{
						const int rand_index = ::tools::random::rand_int32(0, num_tied_cells, random_number);
						int i = 0;
						for (auto cell_i = 1; cell_i < P::N_CELLS_PC; ++cell_i)
						{
							if (counter[cell_i] == selected_count)
							{
								if (i == rand_index) return cell_i;
								i++;
							}
						}
						return best_cell;
					}
					layer.random_number[column_i] = random_number;
				}

				template <typename P>
				void create_DD_segment(
					Layer<P>& layer,
					const int column_i,
					const int time,
					const int n_desired_new_synapses,
					const Dynamic_Param& param,
					const int8_t cell,
					const Layer<P>::Winner_Cells& winner_cells)
				{
					if (n_desired_new_synapses <= 0) return; // nothing to do

					auto& synapse_permanence = layer.dd_synapse_permanence[column_i];
					auto& synapse_delay_origin = layer.dd_synapse_delay_origin[column_i];
					auto& synapse_count = layer.dd_synapse_count[column_i];


					int new_segment_i = -1;

					if (layer.dd_segment_count[column_i] >= P::TP_N_DD_SEGMENTS_MAX) // no empty segments available, use the least recent used one
					{
						const auto& active_time = layer.dd_synapse_active_time[column_i];

						int least_recent_used_time = active_time[0];
						new_segment_i = 0;

						for (int segment_i = 1; segment_i < layer.dd_segment_count[column_i]; ++segment_i)
						{
							const int time = active_time[segment_i];
							if (time < least_recent_used_time)
							{
								least_recent_used_time = time;
								new_segment_i = segment_i;
							}
						}
						if (false) log_INFO("TP:create_DD_segment: column ", column_i, ", has no segment slots left, recycling segment ", new_segment_i, " with least recent used time ", least_recent_used_time);
						// cleanup the old synapses
						synapse_permanence[new_segment_i].clear();
						synapse_delay_origin[new_segment_i].clear();
						synapse_count[new_segment_i] = 0;
					}
					else
					{
						new_segment_i = layer.dd_segment_count[column_i];
						layer.dd_segment_count[column_i]++;
					}

					assert_msg(new_segment_i != -1, "TP:create_DD_segment: error A; new_segment_i = ", new_segment_i);
					if (false) log_INFO_DEBUG("TP:create_new_DD_segment: column ", column_i, "; adding a new segment (", new_segment_i, ") to cell ", static_cast<int>(cell), "; new_segment_idx = ", new_segment_i, ".");

					const int n_new_synapses = n_desired_new_synapses;

					#pragma region Resize Synapses
					// this is the only place where synapse space is created
					assert_msg(synapse_permanence.size() == synapse_delay_origin.size(), "TP:create_new_DD_segment: Bug A.");
					if (synapse_permanence.size() <= new_segment_i)
					{
						const int new_capacity = new_segment_i + 1;
						synapse_permanence.resize(new_capacity);
						synapse_delay_origin.resize(new_capacity);
						synapse_count.resize(new_capacity, 0);

						layer.dd_segment_destination[column_i].resize(new_capacity);
						layer.dd_synapse_active_time[column_i].resize(new_capacity);
					}
					#pragma endregion

					auto& dd_synapse_permanence_segment = synapse_permanence[new_segment_i];
					auto& dd_synapse_delay_origin_segment = synapse_delay_origin[new_segment_i];

					const int n_new_synapsed_capacity = htm::tools::multiple_64(n_new_synapses); // multiple of 16 needed for vectorization;
					dd_synapse_permanence_segment.resize(n_new_synapsed_capacity, P::TP_DD_CONNECTED_THRESHOLD);
					dd_synapse_delay_origin_segment.resize(n_new_synapsed_capacity, P::TP_DD_SYNAPSE_ORIGIN_INVALID); // init with invalid number for debugging purposes

					layer.dd_segment_destination[column_i][new_segment_i] = cell;

					std::vector<int> selected_delay_and_cells(n_new_synapses, layer.random_number[column_i]);
					select_delay_and_cell_to_learn_on(layer, column_i, new_segment_i, winner_cells, n_new_synapses, selected_delay_and_cells);

					for (auto synapse_i = 0; synapse_i < n_new_synapses; ++synapse_i)
					{
						dd_synapse_permanence_segment[synapse_i] = P::TP_DD_PERMANENCE_INIT;
						dd_synapse_delay_origin_segment[synapse_i] = selected_delay_and_cells[synapse_i];
					}

					synapse_count[new_segment_i] = n_new_synapses;
					assert_msg(n_new_synapses <= dd_synapse_permanence_segment.size(), "TP:create_DD_segment: n_new_synapses=", n_new_synapses, "; while dd_synapse_permanence_segment.size()=", dd_synapse_permanence_segment.size());
					assert_msg(n_new_synapses <= dd_synapse_delay_origin_segment.size(), "TP:create_DD_segment: n_new_synapses=", n_new_synapses, "; while dd_synapse_origin_segment.size()=", dd_synapse_delay_origin_segment.size());

					layer.dd_synapse_active_time[column_i][new_segment_i] = time;
				}

				template <bool LEARN, typename P>
				void burst_column(
					Layer<P>& layer,
					const int column_i,
					const int time,
					const Dynamic_Param& param,
					//in
					const Layer<P>::Active_Cells& active_cells,
					const Layer<P>::Winner_Cells& winner_cells,
					//out
					Bitset_Tiny<P::N_CELLS_PC>& current_active_cells,
					Bitset_Tiny<P::N_CELLS_PC>& current_winner_cells)
				{
					if (false) log_INFO_DEBUG("TP:burst_column: column ", column_i, " bursts.");

					const auto& prev_matching_segments = layer.matching_dd_segments[column_i].prev();

					current_active_cells.set_all(); //burst!

					// find the best matching segment
					const auto tup = prev_matching_segments.highest_activity();
					const int best_matching_segment = std::get<0>(tup);
					const int best_segment_activity = std::get<1>(tup);

					const int winner_cell = (best_matching_segment != -1)
						? layer.dd_segment_destination[column_i][best_matching_segment]
						: get_least_used_cell(layer, column_i);

					current_winner_cells.clear_all();
					current_winner_cells.set(winner_cell, true);

					if (LEARN)
					{
						if (best_matching_segment != -1) // found a best matching segment
						{
							// Learn on the best matching segment.
							adapt_segment::d(layer, column_i, best_matching_segment, active_cells, param.TP_DD_PERMANENCE_INC, param.TP_DD_PERMANENCE_DEC);

							const int n_grow_desired = param.TP_DD_MAX_NEW_SYNAPSE_COUNT - best_segment_activity;
							if (false) log_INFO("TP:burst_column: column ", column_i, " bursts. Found best segment ", best_matching_segment, "; n_grow_desired = ", n_grow_desired);
							grow_DD_synapses::d(layer, column_i, best_matching_segment, n_grow_desired, winner_cells);
						}
						else // No matching segments found. Grow a new segment and learn on it.
						{
							const int n_grow_exact = std::min(param.TP_DD_MAX_NEW_SYNAPSE_COUNT, winner_cells.prev().count());
							if (false) log_INFO("TP:burst_column: column ", column_i, " bursts. Not found best segment; n_grow_exact = ", n_grow_exact);
							create_DD_segment(layer, column_i, time, n_grow_exact, param, winner_cell, winner_cells);
						}
					}
				}

				template <bool LEARN, typename P>
				void activate_predicted_column(
					Layer<P>& layer,
					const int column_i,
					const Dynamic_Param& param,
					//in
					const Layer<P>::Active_Cells& active_cells,
					const Layer<P>::Winner_Cells& winner_cells,
					//out
					Bitset_Tiny<P::N_CELLS_PC>& current_active_cells,
					Bitset_Tiny<P::N_CELLS_PC>& current_winner_cells)
				{
					current_active_cells.clear_all();
					current_winner_cells.clear_all();

					const auto& prev_active_segments = layer.active_dd_segments[column_i].prev();
					const auto& segment_destination = layer.dd_segment_destination[column_i];

					for (int i = 0; i < prev_active_segments.count(); ++i)
					{
						const int segment_i = prev_active_segments.get_id(i);
						const auto cell = segment_destination[segment_i];
						current_active_cells.set(cell, true);
						current_winner_cells.set(cell, true);

						if (LEARN)
						{
							adapt_segment::d(layer, column_i, segment_i, active_cells, param.TP_DD_PERMANENCE_INC, param.TP_DD_PERMANENCE_DEC);

							const int segment_activity = prev_active_segments.get_activity(i);
							const int n_grow_desired = param.TP_DD_MAX_NEW_SYNAPSE_COUNT - segment_activity;
							if (n_grow_desired > 0) grow_DD_synapses::d(layer, column_i, segment_i, n_grow_desired, winner_cells);
						}
					}
				}

				//Punishes the segments that incorrectly predicted a column to be active.
				template <typename P>
				void punish_predicted_column(
					Layer<P>& layer,
					const int column_i,
					const Dynamic_Param& param,
					const Layer<P>::Active_Cells& active_cells)
				{
					//if (param.TP_DD_PREDICTED_SEGMENT_DEC > 0.0)
					{
						const auto& prev_matching_segments = layer.matching_dd_segments[column_i].prev();
						for (auto i = 0; i < prev_matching_segments.count(); ++i)
						{
							const auto segment_i = prev_matching_segments.get_id(i);
							adapt_segment::d(layer, column_i, segment_i, active_cells, param.TP_DD_PREDICTED_SEGMENT_DEC);
						}
					}
				}

				template <bool LEARN, typename P>
				void activate_cells_per_column(
					Layer<P>& layer,
					const int column_i,
					const int time,
					const Dynamic_Param& param,
					//in
					const bool is_active_column,
					const Layer<P>::Active_Cells& active_cells,
					const Layer<P>::Winner_Cells& winner_cells,
					//out
					Bitset_Tiny<P::N_CELLS_PC>& current_active_cells,
					Bitset_Tiny<P::N_CELLS_PC>& current_winner_cells)
				{
					if (is_active_column)
					{
						//const bool is_predicted_column = column.active_segments.any_past();
						const bool is_predicted_column = layer.active_dd_segments[column_i].prev().any();
						if (is_predicted_column)
						{
							activate_predicted_column<LEARN>(
								layer,
								column_i,
								param,
								//in
								active_cells,
								winner_cells,
								//out
								current_active_cells,
								current_winner_cells);
						}
						else
						{
							burst_column<LEARN>(
								layer,
								column_i,
								time,
								param,
								//in 
								active_cells,
								winner_cells,
								//out
								current_active_cells,
								current_winner_cells);
						}
					}
					else
					{
						current_active_cells.clear_all();
						current_winner_cells.clear_all();
						if (LEARN) punish_predicted_column(
							layer,
							column_i,
							param,
							active_cells);
					}
				}
				
				template <bool LEARN, typename P>
				void d(
					Layer<P>& layer,
					const int time,
					const Dynamic_Param& param,
					//in
					const Layer<P>::Active_Columns& active_columns,
					//inout
					Layer<P>::Active_Cells& active_cells,
					Layer<P>::Winner_Cells& winner_cells)
				{
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> active_cells_all_2D;
					Bitset2<P::N_COLUMNS, P::N_CELLS_PC> winner_cells_all_2D;

					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						layer.active_dd_segments[column_i].advance_time();
						layer.matching_dd_segments[column_i].advance_time();

						activate_cells_per_column<LEARN>(
							layer,
							column_i,
							time,
							param,
							//in
							active_columns.get(column_i),
							active_cells,
							winner_cells,
							//out
							active_cells_all_2D[column_i],
							winner_cells_all_2D[column_i]);
					}

					active_cells.set_current(active_cells_all_2D);
					copy(winner_cells.current(), winner_cells_all_2D);
				}
			}

			namespace activate_dendrites
			{
				//Get the number of active and matching synapses of the provided segment in the provided column
				namespace count_active_potential_DD_synapses
				{
					template <typename P>
					std::tuple<int, int> count_active_potential_DD_synapses_ref(
						const Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Dynamic_Param& param)
					{
						const auto& dd_synapse_permanence_segment = layer.dd_synapse_permanence[column_i][segment_i];
						const auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];

						int n_potential_synapses = 0;
						int n_active_synapses = 0;

						for (auto synapse_i = 0; synapse_i < layer.dd_synapse_count[column_i][segment_i]; ++synapse_i)
						{
							const Permanence permanence = dd_synapse_permanence_segment[synapse_i];
							if (permanence > P::TP_DD_CONNECTED_THRESHOLD)
							{
								const auto delay_and_cell_id = dd_synapse_delay_origin_segment[synapse_i];
								const int global_cell_id = get_global_cell_id(delay_and_cell_id);
								const int delay = get_delay(delay_and_cell_id) - 1; // can we remove the minus one here: very confusing
								if (active_cells.get(global_cell_id, delay)) // deadly gather here!
								{
									n_potential_synapses++;
									n_active_synapses += (permanence > param.TP_DD_PERMANENCE_THRESHOLD);
								}
							}
						}
						return std::make_tuple(n_active_synapses, n_potential_synapses);
					}

					__m512i get_sensors_epi32(
						const __mmask16 mask,
						const __m512i delay_and_origin_epi32,
						const void * active_sensors_ptr)
					{
						const __m512i global_cell_id = _mm512_and_epi32(delay_and_origin_epi32, _mm512_set1_epi32(0x1FFFFFFF));
						const __m512i byte_pos_in_int = _mm512_and_epi32(delay_and_origin_epi32, _mm512_set1_epi32(0b11));
						const __m512i delay_epi32 = _mm512_sub_epi32(_mm512_srli_epi32(delay_and_origin_epi32, 32 - 3), _mm512_set1_epi32(1));

						const __m512i int_addr = _mm512_srli_epi32(global_cell_id, 2);
						const __m512i sensor_int = _mm512_mask_i32gather_epi32(_mm512_setzero_epi32(), mask, int_addr, active_sensors_ptr, 4);
						const __m512i pos_in_int = _mm512_slli_epi32(byte_pos_in_int, 3);
						const __m512i delay_shift = _mm512_add_epi32(delay_epi32, pos_in_int);
						return _mm512_and_epi32(_mm512_srlv_epi32(sensor_int, delay_shift), _mm512_set1_epi32(1));
					}

					/*
					L1 Data cache = 32 KB, 64 B / line, 8-WAY.
					L1 Instruction cache = 32 KB, 64 B / line, 8-WAY.
					L2 cache = 256 KB, 64 B / line, 4-WAY
					L3 cache = 8 MB, 64 B / line, 16-WAY

					L1 Data Cache Latency = 4 cycles for simple access via pointer
					L1 Data Cache Latency = 5 cycles for access with complex address calculation
					L2 Cache Latency = 12 cycles
					L3 Cache Latency = 42 cycles (i7 - 6700 Skylake 4.0 GHz)
					RAM Latency = 42 cycles + 51 ns (i7 - 6700 Skylake)
					*/
					template <typename P>
					std::tuple<int, int> count_active_potential_DD_synapses_avx512(
						const Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Dynamic_Param& param)
					{
						auto permanence_epi8_ptr = reinterpret_cast<const __m512i *>(layer.dd_synapse_permanence[column_i][segment_i].data());
						auto delay_origin_epi32_ptr = reinterpret_cast<const __m512i *>(layer.dd_synapse_delay_origin[column_i][segment_i].data());
						auto active_cells_ptr = active_cells.data();
						
						const __m512i connected_threshold_epi8 = _mm512_set1_epi8(P::TP_DD_CONNECTED_THRESHOLD);
						const __m512i active_threshold_simd = _mm512_set1_epi8(param.TP_DD_PERMANENCE_THRESHOLD);

						__m512i n_potential_synapses = _mm512_setzero_si512();
						__m512i n_active_synapses = _mm512_setzero_si512();

						const int n_synapses = layer.dd_synapse_count[column_i][segment_i];
						const int n_blocks = htm::tools::n_blocks_64(n_synapses);

						for (int block = 0; block < n_blocks; ++block)
						{
							const __m512i permanence_epi8 = permanence_epi8_ptr[block]; //load 64 permanence values
							const __mmask64 connected_mask_64 = _mm512_cmpgt_epi8_mask(permanence_epi8, connected_threshold_epi8);
							const __mmask64 active_mask_64 = _mm512_cmpgt_epi8_mask(permanence_epi8, active_threshold_simd);

							for (int i = 0; i < 4; ++i)
							{
								const __mmask16 connected_mask_16 = static_cast<__mmask16>(connected_mask_64 >> (i * 16));
								if (connected_mask_16 != 0)
								{
									const __m512i sensors_epi32 = get_sensors_epi32(connected_mask_16, delay_origin_epi32_ptr[(block * 4) + i], active_cells_ptr);
									n_potential_synapses = _mm512_add_epi32(n_potential_synapses, sensors_epi32);
									const __mmask16 active_mask_16 = static_cast<__mmask16>(active_mask_64 >> (i * 16));
									n_active_synapses = _mm512_mask_add_epi32(n_active_synapses, active_mask_16, n_active_synapses, sensors_epi32);
								}
							}
						}

						const int n_potential_synapses_int = _mm512_reduce_add_epi32(n_potential_synapses);
						const int n_active_synapses_int = _mm512_reduce_add_epi32(n_active_synapses);

						if (false) log_INFO_DEBUG("TP:count_active_potential_DD_synapses_avx512: column ", column_i, "; segment_i ", segment_i, " has ", n_synapses, " synapses and ", n_active_synapses_int, " active synapses.\n");

						#if _DEBUG
						const auto tup = count_active_potential_DD_synapses_ref(layer, column_i, segment_i, active_cells, param);
						const int n_active_synapses_ref = std::get<0>(tup);
						const int n_potential_synapses_ref = std::get<1>(tup);

						if (n_active_synapses_int != n_active_synapses_ref)
							log_ERROR("TP:count_active_potential_DD_synapses_avx512:: UNEQUAL number active synapses: AVX512 ", n_active_synapses_int, " != ref ", n_active_synapses_ref, ".\n");
						else
							//log_INFO("TP:count_active_potential_DD_synapses_avx512:: EQUAL number active synapses.\n");

						if (n_potential_synapses_int != n_potential_synapses_ref) log_ERROR("TP:count_active_potential_DD_synapses_avx512:: UNEQUAL number potential synapses: AVX512 ", n_potential_synapses_int, " != ref ", n_potential_synapses_ref, ".\n");
						
						#endif

						return std::make_tuple(n_active_synapses_int, n_potential_synapses_int);
					}

					// Get the number of active and connected synapses of the provided segment in the provided column
					template <typename P>
					std::tuple<int, int> d(
						const Layer<P>& layer,
						const int column_i,
						const int segment_i,
						const Layer<P>::Active_Cells& active_cells,
						const Dynamic_Param& param)
					{
						#if _DEBUG
						if (true)
						{
							const auto& dd_synapse_permanence_segment = layer.dd_synapse_permanence[column_i][segment_i];
							const auto& dd_synapse_delay_origin_segment = layer.dd_synapse_delay_origin[column_i][segment_i];

							// number of synapses in use
							const auto n_synapses = layer.dd_synapse_count[column_i][segment_i];

							if (false) log_INFO("TP:count_active_potential_DD_synapses: A: column ", column_i, "; segment_i ", segment_i, "; n_synapses = ", n_synapses, ".\n");

							if ((n_synapses < 0) || (n_synapses > P::TP_N_DD_SYNAPSES_MAX)) // invalid number of synapses
							{
								log_ERROR("TP:count_active_potential_DD_synapses: B: column ", column_i, "; segment_i ", segment_i, "; n_synapses = ", n_synapses, " which is invalid. P::TP_N_DD_SYNAPSES_MAX = ", P::TP_N_DD_SYNAPSES_MAX, ".\n");
							}

							for (auto synapse_i = 0; synapse_i < n_synapses; ++synapse_i)
							{
								const int delay_and_cell_id = dd_synapse_delay_origin_segment[synapse_i];
								const int global_cell_id = get_global_cell_id(delay_and_cell_id);
								const int delay = get_delay(delay_and_cell_id);
								const Permanence permanence = dd_synapse_permanence_segment[synapse_i];

								if ((global_cell_id < 0) || (global_cell_id >= P::N_CELLS)) // invalid global_cell_id
								{
									log_ERROR("TP:count_active_potential_DD_synapses: C: column ", column_i, "; segment_i ", segment_i, "; synapse ", synapse_i, "; has invalid origin cell = ", global_cell_id, "; permanence = ", static_cast<int>(permanence), ".\n");
								}
							}
						}
						#endif

						switch (architecture_switch(P::ARCH))
						{
							case arch_t::X64: return count_active_potential_DD_synapses_ref(layer, column_i, segment_i, active_cells, param);
							case arch_t::AVX512: return count_active_potential_DD_synapses_avx512(layer, column_i, segment_i, active_cells, param);
							default: return count_active_potential_DD_synapses_ref(layer, column_i, segment_i, active_cells, param);
						}
					}
				}
			
				template <bool LEARN, typename P>
				void d(
					Layer<P>& layer,
					const int time,
					//in
					const Layer<P>::Active_Cells& active_cells,
					const Dynamic_Param& param)
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						auto& active_segments_current = layer.active_dd_segments[column_i].current();
						auto& matching_segments_current = layer.matching_dd_segments[column_i].current();

						active_segments_current.reset();
						matching_segments_current.reset();

						for (auto segment_i = 0; segment_i < layer.dd_segment_count[column_i]; ++segment_i)
						{
							const auto tup = count_active_potential_DD_synapses::d(layer, column_i, segment_i, active_cells, param);
							const int n_active_synapses = std::get<0>(tup);
							const int n_potential_synapses = std::get<1>(tup);

							if (n_potential_synapses > param.TP_MIN_DD_ACTIVATION_THRESHOLD)
							{
								matching_segments_current.add(segment_i, n_potential_synapses);
							}
							if (n_active_synapses > param.TP_DD_SEGMENT_ACTIVE_THRESHOLD)
							{
								active_segments_current.add(segment_i, n_active_synapses);
								if (LEARN) layer.dd_synapse_active_time[column_i][segment_i] = time;
							}
						}
					}
				}
			}
		}
		
		template <bool LEARN, typename P>
		void compute_tp(
			Layer<P>& layer,
			const int time,
			const Dynamic_Param& param,
			//in
			const Layer<P>::Active_Columns& active_columns,
			// inout
			Layer<P>::Active_Cells& active_cells,
			Layer<P>::Winner_Cells& winner_cells)
		{
			#if _DEBUG
			//if (false) log_INFO("TP:compute_tp: prev_winner_cells: ", print::print_active_cells(winner_cells.prev()));
			//if (false) log_INFO("TP:compute_tp: prev_active_cells: ", print::print_active_cells(active_cells.prev()));
			#endif

			priv::activate_cells::d<LEARN>(
				layer,
				time,
				param,
				//in
				active_columns,
				//inout
				active_cells,
				winner_cells);

			#if _DEBUG
			//if (false) log_INFO("TP:compute_tp: active_cells current: time = ", time, ":", print::print_active_cells(active_cells.current()));
			//if (false) log_INFO("TP:compute_tp: winner_cells current: time = ", time, ":", print::print_active_cells(winner_cells.current()));
			#endif

			priv::activate_dendrites::d<LEARN>(
				layer,
				time,
				active_cells,
				param);

			#if _DEBUG
			if (false) log_INFO("TP:compute_tp: all synapses: time = ", time, ":", print::print_dd_synapses(layer));
			#endif
		}
	}
}