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
#include <iostream>		// std::cout
#include <iomanip>		// std::setw
#include <tuple>
#include <array>
#include <vector>

#include "constants.ipp"
#include "tools.ipp"
#include "types.ipp"
#include "print.ipp"
#include "encoder.ipp"
#include "sp.ipp"
#include "tp.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Layer
	namespace layer
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;
		using namespace htm::types;

		//HTM layer private methods
		namespace priv
		{
			//return the number of times a sensor is predicteed.
			// if the predicted sensor influx is ABOVE (not equal) this threshold, the sensor is said to be active.
			void get_predicted_sensor_activity(
				const Layer& layer,
				const int sensor_threshold,
				const Dynamic_Param& param,
				//out
				Bitset_Compact<P::N_SENSORS>& predicted_sensor)
			{
				std::vector<int> predicted_sensor_activity = std::vector<int>(P::N_SENSORS, 0);

				if (SP_GATHER)
				{
					for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
					{
						const Column& column = layer[column_i];
						const bool column_is_predicted = column.active_dd_segments.any_current();
						if (column_is_predicted)
						{
							for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
							{
								if (column.pd_synapse_permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
								{
									const auto sensor_i = column.pd_synapse_origin[synapse_i]; // gather
									predicted_sensor_activity[sensor_i]++;
								}
							}
						}
					}
				}
				else
				{
					if (true)
					{	// somewhat slow method
						Bitset<P::N_COLUMNS> predicted_columns;

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							const Column& column = layer[column_i];
							predicted_columns[column_i] = column.active_dd_segments.any_current();
						}

						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							int sensor_activity = 0;

							const auto& permanence = layer.sp_pd_synapse_permanence[sensor_i];
							const auto& destination_columns = layer.sp_pd_destination_column[sensor_i];

							for (auto synapse_i = 0; synapse_i < permanence.size(); ++synapse_i)
							{
								if (permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
								{
									const int column_i = destination_columns[synapse_i];
									if (predicted_columns.get(column_i))
									{
										sensor_activity++;
										if (sensor_activity > sensor_threshold) break;
									}
								}
							}
							predicted_sensor_activity[sensor_i] = sensor_activity;
						}
					}
					else
					{
						/*
						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							const Column& column = layer[column_i];
							if (column.active_dd_segments.any_current())
							{
								for (auto sensor_i : layer.sp_pd_origin_sensor[column_i])
								{
									const auto& permanence = layer.sp_pd_synapse_permanence[sensor_i][column_i];
									if (permanence > param.SP_PD_CONNECTED_THRESHOLD)
									{

									}
								}
							}
						}

						for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
						{
							int sensor_activity = 0;

							const auto& permanence = layer.sp_pd_synapse_permanence[sensor_i];
							const auto& destination_columns = layer.sp_pd_destination_column[sensor_i];


							for (auto synapse_i = 0; synapse_i < permanence.size(); ++synapse_i)
							{
								if (permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
								{
									const int column_i = destination_columns[synapse_i];
									if (predicted_columns.get(column_i))
									{
										sensor_activity++;
										if (sensor_activity > sensor_threshold) break;
									}
								}
							}
							predicted_sensor_activity[sensor_i] = sensor_activity;
						}
						*/
					}
				}

				predicted_sensor.reset();
				for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
				{
					if (predicted_sensor_activity[sensor_i] > sensor_threshold) predicted_sensor.set(sensor_i);
				}
			}

			int calc_mismatch(
				const int t,
				const Dynamic_Param& param,
				const std::vector<Bitset_Compact<P::N_SENSORS>>& data,
				const Layer& layer)
			{
				Bitset_Compact<P::N_SENSORS> actual_sensors;
				Bitset_Compact<P::N_SENSORS> predicted_sensors;

				const int sensor_threshold = 0; // if the predicted sensor influx is ABOVE (not equal) this threshold, the sensor is said to be active.
				get_predicted_sensor_activity(layer, sensor_threshold, param, predicted_sensors);
				encoder::get_sensor_activity(t + 1, data, actual_sensors);

				int mismatch = 0;
				for (auto i = 0; i < P::N_SENSORS; ++i)
				{
					if (predicted_sensors.get(i) != actual_sensors.get(i)) mismatch++;
				}
				return mismatch;
			}

			void add_sensor_noise(
				Bitset_Compact<P::N_SENSORS>& sensor_activity)
			{
				const int RAND_PERCENT = 1;
				int i = random::rand_int32(1, 200 * RAND_PERCENT);
				while (i < P::N_SENSORS)
				{
					sensor_activity.negate(i);
					i += random::rand_int32(1, 200 * RAND_PERCENT);;
				}
			}

			void load_inferred_sensor_activity(
				const Layer& layer,
				const Dynamic_Param& param,
				const Bitset<P::N_COLUMNS>& active_columns,
				// out
				std::vector<int>& inferred_sensor_activity)
			{
				for (auto i = 0; i < P::N_SENSORS; ++i) inferred_sensor_activity[i] = 0;

				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					if (active_columns[column_i])
					{
						const Column& column = layer[column_i];

						for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
						{
							if (column.pd_synapse_permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
							{
								const auto sensor_idx = column.pd_synapse_origin[synapse_i];
								inferred_sensor_activity[sensor_idx]++;
							}
						}
					}
				}
			}

			void get_projected_boost_factors(
				const Layer& layer,
				const Dynamic_Param& param,
				std::vector<float>& projected_boost_factor) // N_SENSORS
			{
				std::vector<int> counter = std::vector<int>(P::N_SENSORS, 0);

				for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
				{
					projected_boost_factor[sensor_i] = 0.001;
				}

				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					const Column& column = layer[column_i];
					for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
					{
						if (column.pd_synapse_permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
						{
							const auto sensor_idx = column.pd_synapse_origin[synapse_i];
							projected_boost_factor[sensor_idx] += column.boost_factor;
							counter[sensor_idx]++;
						}
					}
				}
				for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
				{
					if (counter[sensor_i] > 0)
					{
						projected_boost_factor[sensor_i] = projected_boost_factor[sensor_i] / counter[sensor_i];
					}
				}
			}

			void show_progress(
				const int t,
				const Layer& layer,
				const Dynamic_Param& param,
				const std::vector<Bitset_Compact<P::N_SENSORS>>& data,
				const Bitset<P::N_COLUMNS>& active_columns)
			{
				const int progress_frequency = 1;
				const int sensor_threshold = 1;
				Bitset_Compact<P::N_SENSORS> sensor_activity_local1;
				Bitset_Compact<P::N_SENSORS> sensor_activity_local2;

				if ((t % progress_frequency) == 0)
				{
					if (true) //print predicted sensor activity
					{
						std::cout << "=====" << std::endl;

						get_predicted_sensor_activity(layer, sensor_threshold, param, sensor_activity_local1);
						encoder::get_sensor_activity(t + 1, data, sensor_activity_local2);

						std::cout << "at t = " << t << ": predicted sensor activity at (future) t = " << (t + 1) << ":" << std::endl;
						std::cout << std::setw(P::N_SENSORS_DIM1) << "predicted";
						std::cout << " | ";
						std::cout << std::setw(P::N_SENSORS_DIM1) << "correct";
						std::cout << " | ";
						std::cout << std::setw(P::N_SENSORS_DIM1) << "mismatch";
						std::cout << std::endl;
						std::cout << print::print_sensor_activity2(sensor_activity_local1, sensor_activity_local2, P::N_SENSORS_DIM1);
					}

					// print the boost values
					if (false) log_INFO("boost factors = \n", print::print_boost_factors(layer, P::N_COLUMNS/20));

					if (false)
					{
						auto projected_boost_factors = std::vector<float>(P::N_SENSORS);
						get_projected_boost_factors(layer, param, projected_boost_factors);
						log_INFO("projected boost factors = \n", print::print_float_array(projected_boost_factors, P::N_SENSORS_DIM1));
					}
					std::cout << "=============================================================" << std::endl;
				}
			}
		}

		template <bool LEARN>
		int run(
			const std::vector<Bitset_Compact<P::N_SENSORS>>& data,
			Layer& layer,
			const Dynamic_Param& param)
		{
			Bitset_Compact<P::N_SENSORS> sensor_activity_local;
			Bitset<P::N_COLUMNS> active_columns_local;

			const int display_interval = 10;
			int total_mismatch = 0;

			for (auto times = 0; times < param.n_times; times++)
			{
				// reset the layer
				layer.init(param);

				int mismatch = 0;

				for (int time = 0; time < param.n_time_steps; ++time)
				{
					//swap(layer.prev_active_cells_all, layer.active_cells_all);
					layer.active_cells.advance_time();
					layer.winner_cells.advance_time();

					encoder::get_sensor_activity(time, data, sensor_activity_local);
					//priv::add_sensor_noise(sensor_activity_local);

					#if _DEBUG
					if (false) log_INFO("layer:run: sensor activity IN:\n", print::print_sensor_activity(sensor_activity_local, P::N_SENSORS_DIM1));
					#endif

					sp::compute_sp<LEARN>(
						sensor_activity_local,
						layer,
						param,
						active_columns_local);

					tp::compute_tp<LEARN>(
						layer,
						time,
						param,
						//in
						active_columns_local,
						//inout
						layer.active_cells,
						layer.winner_cells);

					#if _DEBUG
					if (false) log_INFO("layer:run: active columns at t = ", time, ":\n", print::print_active_columns(active_columns_local, static_cast<int>(std::sqrt(P::N_COLUMNS))));
					//if (false) log_INFO("layer:run: active cells at t = ", time, ": ", print::print_active_cells(layer.active_cells.current()));
					if (false) log_INFO("layer:run: dd_synapes at t = ", time, ": ", print::print_dd_synapses(layer));
					if (false) log_INFO("layer:run: pd_synapes at t = ", time, ": ", print::print_pd_synapses(layer));
					#endif

					const int current_mismatch = priv::calc_mismatch(time, param, data, layer);
					mismatch += current_mismatch;
					total_mismatch += current_mismatch;
					
					if (!param.quiet)
					{
						if (time == 0) std::cout << "layer:run: total mismatch: ";
						if (((time % display_interval) == 0) && (time > 0))
						{
							const float average_mismatch = static_cast<float>(mismatch) / display_interval;
							std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
							mismatch = 0;
						}
					}
					if (param.progress) priv::show_progress(time, layer, param, data, active_columns_local);
				}
				if (!param.quiet) std::cout << std::endl;
			}
			return total_mismatch;
		}
	}
}