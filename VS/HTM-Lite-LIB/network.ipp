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

#include "constants.ipp"
#include "tools.ipp"
#include "types.ipp"
#include "print.ipp"
#include "layer.ipp"
#include "encoder.ipp"

//Hierarchical Temporal Memory (HTM)
namespace htm
{
	//HTM Layer
	namespace network
	{
		using namespace ::tools::log;
		using namespace ::tools::assert;
		using namespace htm::types;

		//HTM layer private methods
		namespace priv
		{
			template <typename P1, typename P2>
			void one_step(
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				const int time,
				const Dynamic_Param& param1,
				const Dynamic_Param& param2)
			{
				static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");
				{
					//2] concat previous layer 2 columns to become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_sensor_activity<P1>(layer1.active_sensors, param1.n_visible_sensors_dim1), "\n");
					layer::one_step(layer1.active_sensors, layer1, time, param1);
				}
				{
					for (int i = 0; i < P2::N_HIDDEN_SENSORS; ++i)
					{
						layer2.active_sensors.set(P2::N_VISIBLE_SENSORS + i, layer1.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_sensor_activity<P2>(layer2.active_sensors, param2.n_visible_sensors_dim1), "\n");
					layer::one_step(layer2.active_sensors, layer2, time, param2);
				}
			}

			template <typename P1, typename P2, typename P3>
			void one_step(
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				Layer<P3>& layer3,
				const int time,
				const Dynamic_Param& param1,
				const Dynamic_Param& param2,
				const Dynamic_Param& param3)
			{
				//static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");
				{
					//Previous layer 2 columns become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2.active_columns.get(i));
					}
					if (true) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_sensor_activity<P1>(layer1.active_sensors, param1.n_visible_sensors_dim1), "\n");
					layer::one_step(layer1.active_sensors, layer1, time, param1);
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
					if (true) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_sensor_activity<P2>(layer2.active_sensors, param2.n_visible_sensors_dim1), "\n");
					layer::one_step(layer2.active_sensors, layer2, time, param2);
				}
				{
					//Previous layer 2 columns become hidden sensors in layer 3
					for (int i = 0; i < P3::N_VISIBLE_SENSORS; ++i)
					{
						layer3.active_sensors.set(i, layer2.active_columns.get(i));
					}
					if (true) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer3:\n", print::print_sensor_activity<P3>(layer3.active_sensors, param3.n_visible_sensors_dim1), "\n");
					layer::one_step(layer3.active_sensors, layer3, time, param3);
				}
			}

			template <typename P1, typename P2>
			int run(
				const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				const Dynamic_Param& param1,
				const Dynamic_Param& param2)
			{
				int total_mismatch = 0;
				int mismatch = 0;
				int current_mismatch = 0;

				for (int time = 0; time < param1.n_time_steps; ++time)
				{
					encoder::get_active_sensors<P1>(time, data, layer1.active_sensors);
					layer::priv::add_sensor_noise<P1>(layer1.active_sensors);

					one_step(layer1, layer2, time, param1, param2);

					if (param1.progress_display_interval > 0)
					{
						current_mismatch = layer::priv::calc_mismatch(time, param1, data, layer1);
					}
					total_mismatch += current_mismatch;

					if (!param1.quiet)
					{
						mismatch += current_mismatch;

						if (time == 0) std::cout << "layer:run: total mismatch: ";
						if (((time % param1.progress_display_interval) == 0) && (time > 0))
						{
							const float average_mismatch = static_cast<float>(mismatch) / param1.progress_display_interval;
							std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
							mismatch = 0;
						}
					}
					if (param1.progress) layer::priv::show_progress(time, layer1, param1, data, layer1.active_columns);
				}
				if (!param1.quiet) std::cout << std::endl;
				return total_mismatch;
			}

			template <typename P1, typename P2, typename P3>
			int run(
				const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				Layer<P3>& layer3,
				const Dynamic_Param& param1,
				const Dynamic_Param& param2,
				const Dynamic_Param& param3)
			{
				//static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");

				int total_mismatch = 0;
				int mismatch = 0;
				int current_mismatch = 0;

				for (int time = 0; time < param1.n_time_steps; ++time)
				{
					encoder::get_active_sensors<P1>(time, data, layer1.active_sensors);
					layer::priv::add_sensor_noise<P1>(layer1.active_sensors);

					one_step(layer1, layer2, layer3, time, param1, param2, param3);

					if (param1.progress_display_interval > 0)
					{
						current_mismatch = layer::priv::calc_mismatch(time, param1, data, layer1);
					}
					total_mismatch += current_mismatch;

					if (!param1.quiet)
					{
						mismatch += current_mismatch;

						if (time == 0) std::cout << "layer:run: total mismatch: ";
						if (((time % param1.progress_display_interval) == 0) && (time > 0))
						{
							const float average_mismatch = static_cast<float>(mismatch) / param1.progress_display_interval;
							std::cout << " " << std::setw(5) << std::setfill(' ') << std::setprecision(2) << average_mismatch;
							mismatch = 0;
						}
					}
					if (param1.progress) layer::priv::show_progress(time, layer1, param1, data, layer1.active_columns);
				}
				if (!param1.quiet) std::cout << std::endl;
				return total_mismatch;
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
			const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2)
		{
			return priv::run(data, layer1, layer2, param1, param2);
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2>
		int run_multiple_times(
			const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2)
		{
			int mismatch = 0;
			for (auto i = 0; i < param1.n_times; ++i)
			{
				layer::init(layer1, param1);
				layer::init(layer2, param2);
				mismatch += run(data, layer1, layer2, param1, param2);
			}
			return mismatch;
		}

		template <typename P1, typename P2, typename P3>
		int run(
			const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			Layer<P3>& layer3,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2,
			const Dynamic_Param& param3)
		{
			return priv::run(data, layer1, layer2, layer3, param1, param2, param3);
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2, typename P3>
		int run_multiple_times(
			const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			Layer<P3>& layer3,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2,
			const Dynamic_Param& param3)
		{
			int mismatch = 0;
			for (auto i = 0; i < param1.n_times; ++i)
			{
				layer::init(layer1, param1);
				layer::init(layer2, param2);
				layer::init(layer3, param3);
				mismatch += run(data, layer1, layer2, layer3, param1, param2, param3);
			}
			return mismatch;
		}
	}
}