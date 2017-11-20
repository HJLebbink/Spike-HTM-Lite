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

#include "parameters.ipp"
#include "tools.ipp"
#include "types.ipp"
#include "print.ipp"
#include "datastream.ipp"
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
		using namespace htm::datastream;

		//HTM layer private methods
		namespace priv
		{
			namespace get_predicted_sensors
			{
				namespace synapse_backward
				{
					//Return the number of times a sensor is predicteed. If the predicted sensor 
					//influx is ABOVE (not equal) this threshold, the sensor is said to be active.
					//Synapse backwards
					template <typename P>
					void get_predicted_sensors_sb_ref(
						const Layer_Fluent<P>& layer_fluent,
						const Layer_Persisted<P>& layer,
						const Dynamic_Param& param,
						//out
						typename Layer_Fluent<P>::Active_Visible_Sensors& predicted_sensor)
					{
						std::vector<int> predicted_sensor_activity = std::vector<int>(P::N_VISIBLE_SENSORS, 0);

						Layer_Fluent<P>::Active_Columns predicted_columns;

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							predicted_columns.set(column_i, layer_fluent.active_dd_segments[column_i].any_current());
						}

						for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
						{
							int sensor_activity = 0;

							const auto& permanence = layer.sp_pd_synapse_permanence_sb[sensor_i];
							const auto& destination_columns = layer.sp_pd_destination_column_sb[sensor_i];

							for (auto synapse_i = 0; synapse_i < layer.sp_pd_synapse_count_sb[sensor_i]; ++synapse_i)
							{
								if (permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
								{
									const int column_i = destination_columns[synapse_i];
									if (predicted_columns.get(column_i))
									{
										sensor_activity++;
										if (sensor_activity > param.sensor_threshold) break;
									}
								}
							}
							predicted_sensor_activity[sensor_i] = sensor_activity;
						}

						predicted_sensor.clear_all();
						for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
						{
							if (predicted_sensor_activity[sensor_i] > param.sensor_threshold) predicted_sensor.set(sensor_i, true);
						}
					}
				}

				namespace synapse_forward
				{
					// predicts sensors one steps into the future
					template <typename P>
					void get_predicted_sensors_sf_ref(
						const Layer_Fluent<P>& layer_fluent,
						Layer_Persisted<P>& layer,//TODO: should have been const but ICC does not allow it
						const Dynamic_Param& param,
						//out
						typename Layer_Fluent<P>::Active_Visible_Sensors& predicted_visible_sensor)
					{
						//log_INFO("get_predicted_sensors_sf_ref: active_sensors:\n", print::print_active_sensors<P>(active_sensors, param.n_visible_sensors_dim1));

						auto predicted_visible_sensor_activity = std::vector<int>(P::N_VISIBLE_SENSORS, 0);

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							const bool column_is_predicted = layer_fluent.active_dd_segments[column_i].any_current();
							if (column_is_predicted)
							{
								const auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];
								const auto& synapse_permanence = layer.sp_pd_synapse_permanence_sf[column_i];

								for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
								{
									if (synapse_permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
									{
										const auto sensor_i = synapse_origin[synapse_i];
										if (sensor_i < P::N_VISIBLE_SENSORS)
										{
											predicted_visible_sensor_activity[sensor_i]++;
										}
									}
								}
							}
						}
						for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
						{
							predicted_visible_sensor.set(sensor_i, (predicted_visible_sensor_activity[sensor_i] > param.sensor_threshold));
						}
					}
					
					// predicts sensor multi time steps into the future, 
					template <typename P>
					void get_predicted_sensors_sf_multifuture_ref(
						const int time,
						const typename Layer_Fluent<P>::Active_Sensors& active_sensors,
						Layer_Fluent<P>& layer_fluent,
						Layer_Persisted<P>& layer, // should have been const but ICC does not allow it
						const Dynamic_Param& param,
						//out
						std::vector<typename Layer_Fluent<P>::Active_Visible_Sensors>& predicted_visible_sensor)
					{
						const int n_futures = static_cast<int>(predicted_visible_sensor.size());
						auto predicted_sensor_activity = std::vector<int>(P::N_SENSORS, 0);
						typename Layer_Fluent<P>::Active_Sensors predicted_sensors = typename Layer_Fluent<P>::Active_Sensors(active_sensors);
						typename Layer_Fluent<P>::Active_Columns predicted_columns;

						//log_INFO("get_predicted_sensors_sf_multifuture_ref: active_sensors:\n", print::print_active_sensors<P>(layer.active_sensors, param.n_visible_sensors_dim1));

						for (auto future_i = 0; future_i < n_futures; ++future_i)
						{
							assert_msg(false, "NOT implemented yet");
							//TODO: calling one_step destroys the fluent state of the layer: only the (const) persisted state is needed
							//priv::one_step<false>(predicted_sensors, layer_fluent, layer, time, param);

							for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
							{
								const bool column_is_predicted = layer_fluent.active_dd_segments[column_i].any_current();
								//const bool column_is_predicted = predicted_columns.get(column_i);
								if (column_is_predicted)
								{
									const auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];
									const auto& synapse_permanence = layer.sp_pd_synapse_permanence_sf[column_i];

									for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
									{
										if (synapse_permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
										{
											const auto sensor_i = synapse_origin[synapse_i];
											predicted_sensor_activity[sensor_i]++;
										}
									}
								}
							}
							if (future_i < (n_futures - 1))
							{
								for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
								{
									predicted_sensors.set(sensor_i, (predicted_sensor_activity[sensor_i] > param.sensor_threshold));
								}
							}
							for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
							{
								predicted_visible_sensor[future_i].set(sensor_i, (predicted_sensor_activity[sensor_i] > param.sensor_threshold));
							}
						}
					}

					// predicts sensors one steps into the future
					template <typename P>
					void get_predicted_sensors_sf_ref(
						const int time,
						Layer_Fluent<P>& layer_fluent,
						Layer_Persisted<P>& layer, //TODO: should have been const but ICC does not allow it
						const Dynamic_Param& param,
						//out
						std::vector<typename Layer_Fluent<P>::Active_Visible_Sensors>& predicted_visible_sensor)
					{
						if (predicted_visible_sensor.size() == 1)
							get_predicted_sensors_sf_ref(layer_fluent, layer, param, predicted_visible_sensor[0]);
							//get_predicted_sensors_sf_multifuture_ref(time, layer_fluent.active_sensors, layer_fluent, layer, param, predicted_visible_sensor);
						else
							get_predicted_sensors_sf_multifuture_ref(time, layer_fluent.active_sensors, layer_fluent, layer, param, predicted_visible_sensor);
					}

					template <typename P>
					void get_predicted_sensors_sf_future_X(
						Layer_Fluent<P>& layer_fluent,
						const Layer_Persisted<P>& layer,
						const int future,
						const Dynamic_Param& param,
						//out
						typename Layer_Fluent<P>::Active_Visible_Sensors& predicted_visible_sensor)
					{
						std::vector<int> predicted_sensor_activity = std::vector<int>(P::N_VISIBLE_SENSORS, 0);
						auto& active_cells = layer.active_cells;

						for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
						{
							bool column_is_predicted = false;
							{
								const int n_segments = layer.dd_segment_count[column_i];
								const auto& permanence_segment = layer.dd_synapse_permanence_sf[column_i];
								const auto& delay_origin_segment = layer.dd_synapse_delay_origin_sf[column_i];
								const auto& synapse_count_segment = layer.dd_synapse_count_sf[column_i];

								for (auto segment_i = 0; segment_i < n_segments; ++segment_i)
								{
									const auto& dd_synapse_permanence_segment = permanence_segment[segment_i];
									const auto& dd_synapse_delay_origin_segment = delay_origin_segment[segment_i];
									const int n_synpases = synapse_count_segment[segment_i];

									int n_active_synapses = 0;

									for (auto synapse_i = 0; synapse_i < n_synpases; ++synapse_i)
									{
										const Permanence permanence = dd_synapse_permanence_segment[synapse_i];
										if (permanence > P::TP_DD_PERMANENCE_THRESHOLD)
										{
											const auto delay_and_cell_id = dd_synapse_delay_origin_segment[synapse_i];
											const int delay = tools::get_delay(delay_and_cell_id) - future;
											if (delay >= 0)
											{
												const int global_cell_id = tools::get_global_cell_id(delay_and_cell_id);
												if (active_cells.get(global_cell_id, delay)) // deadly gather here!
												{
													n_active_synapses++;
												}
											}
										}
									}
									if (n_active_synapses > param.TP_DD_SEGMENT_ACTIVE_THRESHOLD)
									{
										column_is_predicted = true;
										break; // no need to check any other segments in this column
									}
								}
							}
							if (column_is_predicted)
							{
								const auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];
								const auto& synapse_permanence = layer.sp_pd_synapse_permanence_sf[column_i];

								for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
								{
									if (synapse_permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
									{
										const auto sensor_i = synapse_origin[synapse_i];
										if (sensor_i < P::N_VISIBLE_SENSORS)
										{
											predicted_sensor_activity[sensor_i]++;
										}
									}
								}
							}

							for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
							{
								predicted_visible_sensor.set(sensor_i, (predicted_sensor_activity[sensor_i] > param.sensor_threshold));
							}
						}
					}
				}
				
				template <typename P>
				void d(
					const int time,
					Layer_Fluent<P>& layer_fluent,
					Layer_Persisted<P>& layer,//TODO: should have been const but ICC does not allow it
					const Dynamic_Param& param,
					//out
					std::vector<typename Layer_Fluent<P>::Active_Visible_Sensors>& predicted_sensors)
				{
					if (P::SP_SYNAPSE_FORWARD)
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) synapse_forward::get_predicted_sensors_sf_ref(time, layer_fluent, layer, param, predicted_sensors);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) synapse_forward::get_predicted_sensors_sf_ref(time, layer_fluent, layer, param, predicted_sensors);
					}
					else
					{
						if (architecture_switch(P::ARCH) == arch_t::X64) synapse_backward::get_predicted_sensors_sb_ref(layer_fluent, layer, param, predicted_sensors[0]);
						if (architecture_switch(P::ARCH) == arch_t::AVX512) synapse_backward::get_predicted_sensors_sb_ref(layer_fluent, layer, param, predicted_sensors[0]);
					}
				}
			}

			template <typename P>
			void calc_mismatch(
				const int time,
				Layer_Fluent<P>& layer_fluent,
				Layer_Persisted<P>& layer,//TODO: should have been const but ICC does not allow it
				const Dynamic_Param& param,
				const DataStream<P>& datastream,
				//out
				std::vector<int>& mismatch)
			{
				const int n_futures = static_cast<int>(mismatch.size());

				typename std::vector<Layer_Fluent<P>::Active_Sensors> actual_sensors(n_futures);
				typename std::vector<Layer_Fluent<P>::Active_Visible_Sensors> predicted_sensors(n_futures);

				datastream.future_sensors(actual_sensors);
				get_predicted_sensors::d(time, layer_fluent, layer, param, predicted_sensors);

				for (auto future = 0; future < n_futures; ++future)
				{
					int mismatch_counter = 0;
					if (datastream.sensors_predictable(future))
					{
						const auto& predicted = predicted_sensors[future];
						const auto& actual = actual_sensors[future];

						//TODO: the folling loop can be done by xoring the Active_Visible_Sensors
						for (auto i = 0; i < P::N_VISIBLE_SENSORS; ++i)
						{
							if (predicted.get(i) != actual.get(i)) mismatch_counter++;
						}
					}
					mismatch[future] = mismatch_counter;
				}
			}

			template <typename P>
			void load_inferred_sensor_activity(
				const Layer_Persisted<P>& layer,
				const Dynamic_Param& param,
				const typename Layer_Fluent<P>::Active_Columns& active_columns,
				// out
				std::vector<int>& inferred_sensor_activity)
			{
				for (auto i = 0; i < P::N_VISIBLE_SENSORS; ++i) inferred_sensor_activity[i] = 0;

				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					if (active_columns[column_i])
					{
						const Column& column = layer[column_i];

						for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
						{
							if (column.pd_synapse_permanence[synapse_i] > P::SP_PD_PERMANENCE_THRESHOLD)
							{
								const auto sensor_idx = column.pd_synapse_origin[synapse_i];
								inferred_sensor_activity[sensor_idx]++;
							}
						}
					}
				}
			}

			template <typename P>
			void show_input_and_prediction(
				const int time,
				Layer_Fluent<P>& layer_fluent,
				Layer_Persisted<P>& layer,//TODO: should have been const but ICC does not allow it
				const int n_futures,
				const Dynamic_Param& param,
				const DataStream<P>& datastream,
				const typename Layer_Fluent<P>::Active_Columns& active_columns)
			{
				if ((time % param.show_input_and_prediction_interval) == 0)
				{
					if (true) //print predicted visible sensor activity
					{
						typename std::vector<Layer_Fluent<P>::Active_Visible_Sensors> active_visible_sensors(n_futures);
						typename std::vector<Layer_Fluent<P>::Active_Sensors> active_sensors(n_futures);

						std::cout << "=====" << std::endl;

						get_predicted_sensors::d(time, layer_fluent, layer, param, active_visible_sensors);
						datastream.future_sensors(active_sensors);

						for (auto future = 0; future < n_futures; ++future)
						{
							std::cout << "at t = " << time << ": predicted sensor activity at (future " << future << ") t = " << (time + future) << ":" << std::endl;
							std::cout << std::setw(param.n_visible_sensors_dim1) << "predicted";
							std::cout << " | ";
							std::cout << std::setw(param.n_visible_sensors_dim1) << "correct";
							std::cout << " | ";
							std::cout << std::setw(param.n_visible_sensors_dim1) << "mismatch";
							std::cout << std::endl;
							std::cout << print::print_visible_sensor_activity2<P>(active_visible_sensors[0], active_sensors[0], param.n_visible_sensors_dim1);
						}
					}

					// print the boost values
					if (false) log_INFO("boost factors = \n", print::print_boost_factors(layer, P::N_COLUMNS / 20));

					std::cout << "=============================================================" << std::endl;
				}
			}

			template <typename P>
			void show_mismatch(
				const Layer_Persisted<P>& layer,
				const int time,
				const Dynamic_Param& param,
				std::vector<int>& mismatch)
			{
				const int n_futures = static_cast<int>(mismatch.size());
				if (n_futures == 1)
				{
					if (time == 0) std::cout << "layer:show_mismatch: (future=" << n_futures << "):";
					if (((time % param.show_mismatch_interval) == 0) && (time > 0))
					{
						const float average_mismatch = static_cast<float>(mismatch[0]) / param.show_mismatch_interval;
						std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
						tools::clear(mismatch);
					}
				}
				else
				{
					if (time == 0) std::cout << "layer:show_mismatch: (future=" << n_futures << "):" << std::endl;
					if (((time % param.show_mismatch_interval) == 0) && (time > 0))
					{
						std::cout << "time " << time << ":";
						for (auto future = 0; future < n_futures; ++future)
						{
							const float average_mismatch = static_cast<float>(mismatch[future]) / param.show_mismatch_interval;
							std::cout <<"\t"<< std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch ;
						}
						std::cout << std::endl;
						tools::clear(mismatch);
					}
				}
			}

			template <bool LEARN, typename P>
			void one_step(
				const typename Layer_Fluent<P>::Active_Sensors& active_sensors,
				Layer_Fluent<P>& layer_fluent, 
				Layer_Persisted_C<LEARN, P>& layer,
				const int time,
				const Dynamic_Param& param)
			{
				layer_fluent.active_cells.advance_time();
				layer_fluent.winner_cells.advance_time();

				sp::compute_sp<LEARN>(
					active_sensors,
					layer_fluent,
					layer,
					param,
					//out
					layer_fluent.active_columns);

				tp::compute_tp<LEARN>(
					layer_fluent,
					layer,
					time,
					param,
					//in
					layer_fluent.active_columns,
					//inout
					layer_fluent.active_cells,
					layer_fluent.winner_cells);

				#if _DEBUG
				if (false) log_INFO("layer:run: active columns at t = ", time, ":\n", print::print_active_columns<P>(layer.fluent.active_columns, static_cast<int>(std::sqrt(P::N_COLUMNS))), "\n");
				if (false) log_INFO("layer:run: dd_synapes at t = ", time, ": ", print::print_dd_synapses(layer), "\n");
				#endif
			}
		}

		//Reset the provided layer with the properties as provided in param
		template <typename P>
		void init(Layer_Fluent<P>& layer_fluent, Layer_Persisted<P>& layer, const Dynamic_Param& param)
		{
			// reset pd synapses
			if (P::SP_SYNAPSE_FORWARD)
			{
				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					auto& synapse_origin = layer.sp_pd_synapse_origin_sensor_sf[column_i];
					auto& synapse_permanence = layer.sp_pd_synapse_permanence_sf[column_i];
					unsigned int random_number = layer_fluent.random_number[column_i];

					for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
					{
						const int random_sensor = random::rand_int32(0, P::N_SENSORS - 1, random_number);
						synapse_origin[synapse_i] = random_sensor;
						synapse_permanence[synapse_i] = param.SP_PD_PERMANENCE_INIT;
					}
					layer_fluent.random_number[column_i] = random_number;
				}
			}
			else
			{
				for (auto sensor_i = 0; sensor_i < P::N_SENSORS; ++sensor_i)
				{
					layer.sp_pd_destination_column_sb[sensor_i].clear();
					layer.sp_pd_synapse_permanence_sb[sensor_i].clear();
					layer.sp_pd_synapse_count_sb[sensor_i] = 0;
				}
				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
					{
						const int random_sensor = random::rand_int32(0, P::N_SENSORS - 1, layer_fluent.random_number[column_i]);

						auto& destination = layer.sp_pd_destination_column_sb[random_sensor];
						auto& permanence = layer.sp_pd_synapse_permanence_sb[random_sensor];
						const int old_size = layer.sp_pd_synapse_count_sb[random_sensor];
						const int new_size = old_size + 1;

						if (destination.size() <= new_size)
						{
							const int new_capacity = tools::multiple_16(new_size);
							destination.resize(new_capacity);
							permanence.resize(new_capacity);
						}

						destination[old_size] = column_i;
						permanence[old_size] = param.SP_PD_PERMANENCE_INIT;
						layer.sp_pd_synapse_count_sb[random_sensor] = new_size;
					}
				}
			}

			//init permanence values
			for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
			{
				// reset dd synapses
				layer.dd_segment_count[column_i] = 0;
				layer.dd_synapse_count_sf[column_i].clear();
				layer.dd_synapse_permanence_sf[column_i].clear();
				layer.dd_synapse_delay_origin_sf[column_i].clear();
				layer_fluent.dd_synapse_active_time[column_i].clear();

				// reset activity
				layer_fluent.active_dd_segments[column_i].reset();
				layer_fluent.matching_dd_segments[column_i].reset();
			}

			// reset global state
			layer_fluent.active_cells.reset();
			layer_fluent.winner_cells.reset();
		}

		template <typename P>
		void one_step(
			const typename Layer_Fluent<P>::Active_Sensors& active_sensors,
			Layer_Fluent<P>& layer_fluent, Layer_Persisted<P>& layer,
			const int time,
			const Dynamic_Param& param)
		{
			if (param.learn)
				priv::one_step<true>(active_sensors, layer_fluent, layer, time, param);
			else
				priv::one_step<false>(active_sensors, layer_fluent, layer, time, param);
		}

		template <typename P>
		void display_info(
			const DataStream<P>& datastream,
			Layer_Fluent<P>& layer_fluent,
			Layer_Persisted<P>& layer,//TODO: should have been const but ICC does not allow it
			const int time,
			const Dynamic_Param& param,
			const std::vector<int>& current_mismatch,
			std::vector<int>& mismatch)
		{
			tools::add(mismatch, current_mismatch);

			if (!param.quiet)
			{
				if (param.show_mismatch_interval > 0)
				{
					priv::show_mismatch(layer, time, param, mismatch);
				}
				if (param.show_input_and_prediction_interval > 0)
				{
					priv::show_input_and_prediction(time, layer_fluent, layer, 1, param, datastream, layer_fluent.active_columns);
				}
			}
		}

		//Run the provided layer once, update steps as provided in param
		template <typename P>
		void run(
			const DataStream<P>& datastream,
			const Dynamic_Param& param,
			Layer_Fluent<P>& layer_fluent, Layer_Persisted<P>& layer,
			//out
			std::vector<int>& prediction_mismatch)
		{
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			tools::clear(prediction_mismatch);

			auto mismatch = std::vector<int>(n_futures, 0);
			auto current_mismatch = std::vector<int>(n_futures, 0);

			for (auto time = 0; time < param.n_time_steps; ++time)
			{
				datastream.current_sensors(layer_fluent.active_sensors);
				encoder::add_sensor_noise<P>(layer_fluent.active_sensors);
				one_step(layer_fluent.active_sensors, layer_fluent, layer, time, param);

				if (n_futures > 0)
				{
					layer::priv::calc_mismatch(time, layer_fluent, layer, param, datastream, current_mismatch);
					tools::add(prediction_mismatch, current_mismatch);
					layer::display_info(datastream, layer_fluent, layer, time, param, current_mismatch, mismatch);
				}
				datastream.advance_time();
			}
			if (!param.quiet) std::cout << std::endl;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P>
		void run_multiple_times(
			const DataStream<P>& datastream,
			Layer_Fluent<P>& layer_fluent, Layer_Persisted<P>& layer,
			const Dynamic_Param& param,
			//out
			std::vector<int>& prediction_mismatch)
		{
			tools::clear(prediction_mismatch);
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			auto mismatch = std::vector<int>(n_futures);
			for (auto i = 0; i < param.n_times; ++i)
			{
				init(layer_fluent, layer, param);
				datastream.reset_time();
				run(datastream, param, layer_fluent, layer, mismatch);
				tools::add(prediction_mismatch, mismatch);
			}
		}
	}
}