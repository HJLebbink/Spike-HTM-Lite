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
#include <vector>

#include "parameters.ipp"
#include "tools.ipp"
#include "types.ipp"
#include "print.ipp"
#include "layer.ipp"
#include "datastream.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Layer
	namespace network
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;
		using namespace htm::types;
		using namespace htm::datastream;

		//HTM layer private methods
		namespace priv
		{
			template <typename P1, typename P2>
			void one_step(
				const std::array<Dynamic_Param, 2>& param,
				const int time,
				Layer<P1>& layer1,
				Layer<P2>& layer2)
			{
				static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");
				{
					//2] concat previous layer 2 columns to become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_sensor_activity<P1>(layer1.active_sensors, param[0].n_visible_sensors_dim1), "\n");
					layer::one_step(layer1.active_sensors, layer1, time, param[0]);
				}
				{
					for (int i = 0; i < P2::N_HIDDEN_SENSORS; ++i)
					{
						layer2.active_sensors.set(P2::N_VISIBLE_SENSORS + i, layer1.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_sensor_activity<P2>(layer2.active_sensors, param[1].n_visible_sensors_dim1), "\n");
					layer::one_step(layer2.active_sensors, layer2, time, param[1]);
				}
			}

			template <typename P1, typename P2, typename P3>
			void one_step(
				const std::array<Dynamic_Param, 3>& param, 
				const int time,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				Layer<P3>& layer3)
			{
				{
					//Previous layer 2 columns become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_sensor_activity<P1>(layer1.active_sensors, param[0].n_visible_sensors_dim1), "\n");
					layer::one_step(layer1.active_sensors, layer1, time, param[0]);
				}
				{
					//Previous layer 1 columns become visible sensors in layer 2
					for (int i = 0; i < P2::N_VISIBLE_SENSORS; ++i)
					{
						layer2.active_sensors.set(i, layer1.active_columns.get(i));
					}
					//Previous layer 3 columns become hidden sensors in layer 2
					for (int i = 0; i < P2::N_HIDDEN_SENSORS; ++i)
					{
						layer2.active_sensors.set(P2::N_VISIBLE_SENSORS + i, layer3.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_sensor_activity<P2>(layer2.active_sensors, param[1].n_visible_sensors_dim1), "\n");
					layer::one_step(layer2.active_sensors, layer2, time, param[1]);
				}
				{
					//Previous layer 2 columns become hidden sensors in layer 3
					for (int i = 0; i < P3::N_VISIBLE_SENSORS; ++i)
					{
						layer3.active_sensors.set(i, layer2.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer3:\n", print::print_sensor_activity<P3>(layer3.active_sensors, param[2].n_visible_sensors_dim1), "\n");
					layer::one_step(layer3.active_sensors, layer3, time, param[2]);
				}
			}
		}

		template <
			int N_VISIBLE_SENSORS,
			int N_COLUMNS_L1, int N_BITS_CELL_L1, int HISTORY_L1,
			int N_COLUMNS_L2, int N_BITS_CELL_L2, int HISTORY_L2,
			arch_t ARCH>
		struct network_2Layer
		{
			static const int N_VISIBLE_SENSORS_L1 = N_VISIBLE_SENSORS;
			static const int N_HIDDEN_SENSORS_L1 = N_COLUMNS_L2;

			static const int N_VISIBLE_SENSORS_L2 = N_COLUMNS_L1;
			static const int N_HIDDEN_SENSORS_L2 = 0;

			using P_L1 = Static_Param<N_COLUMNS_L1, N_BITS_CELL_L1, N_VISIBLE_SENSORS_L1, N_HIDDEN_SENSORS_L1, HISTORY_L1, ARCH>;
			using P_L2 = Static_Param<N_COLUMNS_L2, N_BITS_CELL_L2, N_VISIBLE_SENSORS_L2, N_HIDDEN_SENSORS_L2, HISTORY_L2, ARCH>;
		};

		template <
			int N_VISIBLE_SENSORS,
			int N_COLUMNS_L1, int N_BITS_CELL_L1, int HISTORY_L1,
			int N_COLUMNS_L2, int N_BITS_CELL_L2, int HISTORY_L2,
			int N_COLUMNS_L3, int N_BITS_CELL_L3, int HISTORY_L3,
			arch_t ARCH>
		struct network_3Layer
		{
			static const int N_VISIBLE_SENSORS_L1 = N_VISIBLE_SENSORS;
			static const int N_HIDDEN_SENSORS_L1 = N_COLUMNS_L2;

			static const int N_VISIBLE_SENSORS_L2 = N_COLUMNS_L1;
			static const int N_HIDDEN_SENSORS_L2 = N_COLUMNS_L3;

			static const int N_VISIBLE_SENSORS_L3 = N_COLUMNS_L2;
			static const int N_HIDDEN_SENSORS_L3 = 0;

			using P_L1 = Static_Param<N_COLUMNS_L1, N_BITS_CELL_L1, N_VISIBLE_SENSORS_L1, N_HIDDEN_SENSORS_L1, HISTORY_L1, ARCH>;
			using P_L2 = Static_Param<N_COLUMNS_L2, N_BITS_CELL_L2, N_VISIBLE_SENSORS_L2, N_HIDDEN_SENSORS_L2, HISTORY_L2, ARCH>;
			using P_L3 = Static_Param<N_COLUMNS_L3, N_BITS_CELL_L3, N_VISIBLE_SENSORS_L3, N_HIDDEN_SENSORS_L3, HISTORY_L3, ARCH>;
		};

		template <typename P1, typename P2>
		int run(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 2>& param,
			Layer<P1>& layer1,
			Layer<P2>& layer2)
		{
			int total_mismatch = 0;
			int mismatch = 0;
			int current_mismatch = 0;

			for (int time = 0; time < param[0].n_time_steps; ++time)
			{
				datastream.current_sensors(layer1.active_sensors);
				encoder::add_sensor_noise<P1>(layer1.active_sensors);

				priv::one_step(param, time, layer1, layer2);

				if (param[0].progress_display_interval > 0)
				{
					current_mismatch = layer::priv::calc_mismatch(time, param[0], datastream, layer1);
				}
				total_mismatch += current_mismatch;

				if (!param[0].quiet)
				{
					mismatch += current_mismatch;

					if (time == 0) std::cout << "layer:run: total mismatch: ";
					if (((time % param[0].progress_display_interval) == 0) && (time > 0))
					{
						const float average_mismatch = static_cast<float>(mismatch) / param[0].progress_display_interval;
						std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
						mismatch = 0;
					}
				}
				if (param[0].progress) layer::priv::show_progress(time, layer1, param[0], datastream, layer1.active_columns);

				datastream.advance_time();
			}
			if (!param[0].quiet) std::cout << std::endl;
			return total_mismatch;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2>
		int run_multiple_times(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 2>& param,
			Layer<P1>& layer1,
			Layer<P2>& layer2)
		{
			int mismatch = 0;
			for (auto i = 0; i < param[0].n_times; ++i)
			{
				layer::init(layer1, param[0]);
				layer::init(layer2, param[1]);
				mismatch += run(datastream, param, layer1, layer2);
			}
			return mismatch;
		}

		template <typename P1, typename P2, typename P3>
		int run(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 3>& param,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			Layer<P3>& layer3)
		{
			int total_mismatch = 0;
			int mismatch = 0;
			int current_mismatch = 0;

			for (int time = 0; time < param[0].n_time_steps; ++time)
			{
				datastream.current_sensors(layer1.active_sensors);
				encoder::add_sensor_noise<P1>(layer1.active_sensors);

				priv::one_step(param, time, layer1, layer2, layer3);

				if (param[0].progress_display_interval > 0)
				{
					current_mismatch = layer::priv::calc_mismatch(time, param[0], datastream, layer1);
				}
				total_mismatch += current_mismatch;

				if (!param[0].quiet)
				{
					mismatch += current_mismatch;

					if (time == 0) std::cout << "layer:run: total mismatch: ";
					if (((time % param[0].progress_display_interval) == 0) && (time > 0))
					{
						const float average_mismatch = static_cast<float>(mismatch) / param[0].progress_display_interval;
						std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
						mismatch = 0;
					}
				}
				if (param[0].progress) layer::priv::show_progress(time, layer1, param[0], datastream, layer1.active_columns);

				datastream.advance_time();
			}
			if (!param[0].quiet) std::cout << std::endl;
			return total_mismatch;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2, typename P3>
		int run_multiple_times(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 3>& param,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			Layer<P3>& layer3)
		{
			int mismatch = 0;
			for (auto i = 0; i < param[0].n_times; ++i)
			{
				layer::init(layer1, param[0]);
				layer::init(layer2, param[1]);
				layer::init(layer3, param[2]);
				mismatch += run(datastream, param, layer1, layer2, layer3);
			}
			return mismatch;
		}
	}
}