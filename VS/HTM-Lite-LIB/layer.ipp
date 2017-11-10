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
			template <typename P>
			void get_predicted_sensors(
				const Layer<P>& layer,
				const int sensor_threshold,
				const Dynamic_Param& param,
				//out
				Layer<P>::Active_Visible_Sensors& predicted_sensor)
			{
				std::vector<int> predicted_sensor_activity = std::vector<int>(P::N_VISIBLE_SENSORS, 0);

				Layer<P>::Active_Columns predicted_columns;

				for (auto column_i = 0; column_i < P::N_COLUMNS; ++column_i)
				{
					const auto& column = layer[column_i];
					predicted_columns.set(column_i, column.active_dd_segments.any_current());
				}

				for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
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

				predicted_sensor.clear_all();
				for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
				{
					if (predicted_sensor_activity[sensor_i] > sensor_threshold) predicted_sensor.set(sensor_i, true);
				}
			}

			template <typename P>
			int calc_mismatch(
				const int t,
				const Dynamic_Param& param,
				const std::vector<Layer<P>::Active_Visible_Sensors>& data,
				const Layer<P>& layer)
			{
				Layer<P>::Active_Sensors active_sensors;
				Layer<P>::Active_Visible_Sensors predicted_sensors;

				const int sensor_threshold = 0; // if the predicted sensor influx is ABOVE (not equal) this threshold, the sensor is said to be active.
				get_predicted_sensors(layer, sensor_threshold, param, predicted_sensors);
				encoder::get_active_sensors<P>(t + 1, data, active_sensors);

				int mismatch = 0;
				//TODO: the folling loop can be done by xoring the Active_Visible_Sensors
				for (auto i = 0; i < P::N_VISIBLE_SENSORS; ++i)
				{
					if (predicted_sensors.get(i) != active_sensors.get(i)) mismatch++;
				}
				return mismatch;
			}

			template <typename P>
			void add_sensor_noise(
				Layer<P>::Active_Sensors& active_sensors)
			{
				if (P::SENSOR_NOISE_PERCENT > 0)
				{
					const int RAND_PERCENT = P::SENSOR_NOISE_PERCENT;
					int i = random::rand_int32(1, 200 * RAND_PERCENT);
					while (i < P::N_VISIBLE_SENSORS)
					{
						active_sensors.negate(i);
						i += random::rand_int32(1, 200 * RAND_PERCENT);;
					}
				}
			}

			template <typename P>
			void load_inferred_sensor_activity(
				const Layer<P>& layer,
				const Dynamic_Param& param,
				const Layer<P>::Active_Columns& active_columns,
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
							if (column.pd_synapse_permanence[synapse_i] > param.SP_PD_CONNECTED_THRESHOLD)
							{
								const auto sensor_idx = column.pd_synapse_origin[synapse_i];
								inferred_sensor_activity[sensor_idx]++;
							}
						}
					}
				}
			}

			template <typename P>
			void show_progress(
				const int t,
				const Layer<P>& layer,
				const Dynamic_Param& param,
				const std::vector<Layer<P>::Active_Visible_Sensors>& data,
				const Layer<P>::Active_Columns& active_columns)
			{
				const int progress_frequency = 1;
				const int sensor_threshold = 1;
				Layer<P>::Active_Visible_Sensors active_visible_sensors;
				Layer<P>::Active_Sensors active_sensors;

				if ((t % progress_frequency) == 0)
				{
					if (true) //print predicted visible sensor activity
					{
						std::cout << "=====" << std::endl;

						get_predicted_sensors(layer, sensor_threshold, param, active_visible_sensors);
						encoder::get_active_sensors<P>(t + 1, data, active_sensors);

						std::cout << "at t = " << t << ": predicted sensor activity at (future) t = " << (t + 1) << ":" << std::endl;
						std::cout << std::setw(param.n_visible_sensors_dim1) << "predicted";
						std::cout << " | ";
						std::cout << std::setw(param.n_visible_sensors_dim1) << "correct";
						std::cout << " | ";
						std::cout << std::setw(param.n_visible_sensors_dim1) << "mismatch";
						std::cout << std::endl;
						std::cout << print::print_visible_sensor_activity2<P>(active_visible_sensors, active_sensors, param.n_visible_sensors_dim1);
					}

					// print the boost values
					if (false) log_INFO("boost factors = \n", print::print_boost_factors(layer, P::N_COLUMNS / 20));

					std::cout << "=============================================================" << std::endl;
				}
			}

			template <bool LEARN, typename P>
			void one_step(
				const Layer<P>::Active_Sensors& active_sensors,
				Layer<P>& layer,
				const int time,
				const Dynamic_Param& param)
			{
				layer.active_cells.advance_time();
				layer.winner_cells.advance_time();

				sp::compute_sp<LEARN>(
					active_sensors,
					layer,
					param,
					//out
					layer.active_columns);

				tp::compute_tp<LEARN>(
					layer,
					time,
					param,
					//in
					layer.active_columns,
					//inout
					layer.active_cells,
					layer.winner_cells);

				#if _DEBUG
				if (false) log_INFO("layer:run: active columns at t = ", time, ":\n", print::print_active_columns<P>(layer.active_columns, static_cast<int>(std::sqrt(P::N_COLUMNS))), "\n");
				if (false) log_INFO("layer:run: dd_synapes at t = ", time, ": ", print::print_dd_synapses(layer), "\n");
				#endif
			}
		}

		//Reset the provided layer with the properties as provided in param
		template <typename P>
		void init(Layer<P>& layer, const Dynamic_Param& param)
		{
			for (auto sensor_i = 0; sensor_i < P::N_VISIBLE_SENSORS; ++sensor_i)
			{
				layer.sp_pd_destination_column[sensor_i].clear();
				layer.sp_pd_synapse_permanence[sensor_i].clear();
			}

			//init permanence values
			for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i)
			{
				auto& column = layer.data[column_i];
				column.id = column_i;

				// reset pd synapses
				for (auto synapse_i = 0; synapse_i < P::SP_N_PD_SYNAPSES; ++synapse_i)
				{
					const int random_sensor = random::rand_int32(0, P::N_VISIBLE_SENSORS - 1, column.random_number);
					layer.sp_pd_destination_column[random_sensor].push_back(column_i);
					layer.sp_pd_synapse_permanence[random_sensor].push_back(param.SP_PD_PERMANENCE_INIT);

					if (false) log_INFO("Layer::init: column ", column_i, "; synapse ", synapse_i, "; pd_synapse_origin ", random_sensor, "; current_random_number ", static_cast<int>(random::priv::current_random_number));
				}

				// reset dd synapses
				column.dd_segment_count = 0;
				column.dd_synapse_count.clear();
				column.dd_synapse_permanence.clear();
				column.dd_synapse_delay_origin.clear();
				column.dd_synapse_active_time.clear();

				// reset activity
				column.active_dd_segments.reset();
				column.matching_dd_segments.reset();
			}

			// reset global state
			layer.active_cells.reset();
			layer.winner_cells.reset();
		}


		template <typename P>
		void one_step(
			const Layer<P>::Active_Sensors& active_sensors,
			Layer<P>& layer,
			const int time,
			const Dynamic_Param& param)
		{
			if (param.learn)
				priv::one_step<true>(active_sensors, layer, time, param);
			else
				priv::one_step<false>(active_sensors, layer, time, param);
		}


		//Run the provided layer once, update steps as provided in param
		template <typename P>
		int run(
			const std::vector<Layer<P>::Active_Visible_Sensors>& data,
			Layer<P>& layer,
			const Dynamic_Param& param)
		{
			int total_mismatch = 0;
			int mismatch = 0;
			int current_mismatch = 0;

			for (int time = 0; time < param.n_time_steps; ++time)
			{
				encoder::get_active_sensors<P>(time, data, layer.active_sensors);
				priv::add_sensor_noise<P>(layer.active_sensors);
				one_step(layer.active_sensors, layer, time, param);

				if (param.progress_display_interval > 0)
				{
					current_mismatch = priv::calc_mismatch(time, param, data, layer);
				}
				total_mismatch += current_mismatch;

				if (!param.quiet)
				{
					mismatch += current_mismatch;

					if (time == 0) std::cout << "layer:run: total mismatch: ";
					if (((time % param.progress_display_interval) == 0) && (time > 0))
					{
						const float average_mismatch = static_cast<float>(mismatch) / param.progress_display_interval;
						std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
						mismatch = 0;
					}
				}
				if (param.progress) priv::show_progress(time, layer, param, data, layer.active_columns);
			}
			if (!param.quiet) std::cout << std::endl;
			return total_mismatch;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P>
		int run_multiple_times(
			const std::vector<Layer<P>::Active_Visible_Sensors>& data,
			Layer<P>& layer,
			const Dynamic_Param& param)
		{
			int mismatch = 0;
			for (auto i = 0; i < param.n_times; ++i)
			{
				init(layer, param);
				mismatch += run(data, layer, param);
			}
			return mismatch;
		}
	}
}