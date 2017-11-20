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
				Layer_Fluent<P1>& layer1_fluent, Layer_Persisted<P1>& layer1,
				Layer_Fluent<P2>& layer2_fluent, Layer_Persisted<P2>& layer2)
			{
				static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");
				{
					//2] concat previous layer 2 columns to become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1_fluent.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2_fluent.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_active_sensors<P1>(layer1_fluent.active_sensors, param[0].n_visible_sensors_dim1), "\n");
					layer::one_step(layer1_fluent.active_sensors, layer1_fluent, layer1, time, param[0]);
				}
				{
					for (int i = 0; i < P2::N_HIDDEN_SENSORS; ++i)
					{
						layer2_fluent.active_sensors.set(P2::N_VISIBLE_SENSORS + i, layer1_fluent.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_active_sensors<P2>(layer2_fluent.active_sensors, param[1].n_visible_sensors_dim1), "\n");
					layer::one_step(layer2_fluent.active_sensors, layer2_fluent, layer2, time, param[1]);
				}
			}

			template <typename P1, typename P2, typename P3>
			void one_step(
				const std::array<Dynamic_Param, 3>& param, 
				const int time,
				Layer_Fluent<P1>& layer1_fluent, Layer_Persisted<P1>& layer1,
				Layer_Fluent<P2>& layer2_fluent, Layer_Persisted<P2>& layer2,
				Layer_Fluent<P3>& layer3_fluent, Layer_Persisted<P3>& layer3)
			{
				{
					//Previous layer 2 columns become hidden sensors in layer 1
					for (int i = 0; i < P1::N_HIDDEN_SENSORS; ++i)
					{
						layer1_fluent.active_sensors.set(P1::N_VISIBLE_SENSORS + i, layer2_fluent.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer1:\n", print::print_active_sensors<P1>(layer1_fluent.active_sensors, param[0].n_visible_sensors_dim1), "\n");
					layer::one_step(layer1_fluent.active_sensors, layer1_fluent, layer1, time, param[0]);
				}
				{
					//Previous layer 1 columns become visible sensors in layer 2
					for (int i = 0; i < P2::N_VISIBLE_SENSORS; ++i)
					{
						layer2_fluent.active_sensors.set(i, layer1_fluent.active_columns.get(i));
					}
					//Previous layer 3 columns become hidden sensors in layer 2
					for (int i = 0; i < P2::N_HIDDEN_SENSORS; ++i)
					{
						layer2_fluent.active_sensors.set(P2::N_VISIBLE_SENSORS + i, layer3_fluent.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer2:\n", print::print_active_sensors<P2>(layer2_fluent.active_sensors, param[1].n_visible_sensors_dim1), "\n");
					layer::one_step(layer2_fluent.active_sensors, layer2_fluent, layer2, time, param[1]);
				}
				{
					//Previous layer 2 columns become hidden sensors in layer 3
					for (int i = 0; i < P3::N_VISIBLE_SENSORS; ++i)
					{
						layer3_fluent.active_sensors.set(i, layer2_fluent.active_columns.get(i));
					}
					if (false) log_INFO_DEBUG("network:run: active sensors at t = ", time, ": Layer3:\n", print::print_active_sensors<P3>(layer3_fluent.active_sensors, param[2].n_visible_sensors_dim1), "\n");
					layer::one_step(layer3_fluent.active_sensors, layer3_fluent, layer3, time, param[2]);
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
		void run(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 2>& param,
			Layer_Fluent<P1>& layer1_fluent, Layer_Persisted<P1>& layer1, 
			Layer_Fluent<P2>& layer2_fluent, Layer_Persisted<P2>& layer2,
			//out
			std::vector<int>& prediction_mismatch)
		{
			const auto param0 = param[0];
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			tools::clear(prediction_mismatch);

			auto mismatch = std::vector<int>(n_futures, 0);
			auto current_mismatch = std::vector<int>(n_futures, 0);

			for (auto time = 0; time < param0.n_time_steps; ++time)
			{
				datastream.current_sensors(layer1_fluent.active_sensors);
				encoder::add_sensor_noise<P1>(layer1_fluent.active_sensors);
				priv::one_step(param, time, layer1_fluent, layer1, layer2_fluent, layer2);

				if (n_futures > 0)
				{
					layer::priv::calc_mismatch(time, layer1_fluent, layer1, param0, datastream, current_mismatch);
					tools::add(prediction_mismatch, current_mismatch);
					layer::display_info(datastream, layer1_fluent, layer1, time, param0, current_mismatch, mismatch);
				}
				datastream.advance_time();
			}
			if (!param0.quiet) std::cout << std::endl;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2>
		void run_multiple_times(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 2>& param,
			Layer_Fluent<P1> layer1_fluent, Layer_Persisted<P1>& layer1,
			Layer_Fluent<P2> layer2_fluent, Layer_Persisted<P2>& layer2,
			//out
			std::vector<int>& prediction_mismatch)
		{
			tools::clear(prediction_mismatch);
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			auto mismatch = std::vector<int>(n_futures);
			for (auto i = 0; i < param[0].n_times; ++i)
			{
				layer::init(layer1_fluent, layer1, param[0]);
				layer::init(layer2_fluent, layer2, param[1]);
				run(datastream, param, layer1_fluent, layer1, layer2_fluent, layer2, mismatch);
				tools::add(prediction_mismatch, mismatch);
			}
		}

		template <typename P1, typename P2, typename P3>
		void run(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 3>& param,
			Layer_Fluent<P1>& layer1_fluent, Layer_Persisted<P1>& layer1,
			Layer_Fluent<P2>& layer2_fluent, Layer_Persisted<P2>& layer2,
			Layer_Fluent<P3>& layer3_fluent, Layer_Persisted<P3>& layer3,
			std::vector<int>& prediction_mismatch)
		{
			const auto param0 = param[0];
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			tools::clear(prediction_mismatch);

			auto mismatch = std::vector<int>(n_futures, 0);
			auto current_mismatch = std::vector<int>(n_futures, 0);

			for (auto time = 0; time < param0.n_time_steps; ++time)
			{
				datastream.current_sensors(layer1_fluent.active_sensors);
				encoder::add_sensor_noise<P1>(layer1_fluent.active_sensors);
				priv::one_step(param, time, layer1_fluent, layer1, layer2_fluent, layer2, layer3_fluent, layer3);

				if (n_futures > 0)
				{
					layer::priv::calc_mismatch(time, layer1_fluent, layer1, param0, datastream, current_mismatch);
					tools::add(prediction_mismatch, current_mismatch);
					layer::display_info(datastream, layer1_fluent, layer1, time, param0, current_mismatch, mismatch);
				}
				datastream.advance_time();
			}
			if (!param0.quiet) std::cout << std::endl;
		}

		//Run the provided layer a number of times, update steps as provided in param
		template <typename P1, typename P2, typename P3>
		void run_multiple_times(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 3>& param,
			Layer_Fluent<P1>& layer1_fluent, Layer_Persisted<P1>& layer1,
			Layer_Fluent<P2>& layer2_fluent, Layer_Persisted<P2>& layer2,
			Layer_Fluent<P3>& layer3_fluent, Layer_Persisted<P3>& layer3,
			std::vector<int>& prediction_mismatch)
		{
			tools::clear(prediction_mismatch);
			const int n_futures = static_cast<int>(prediction_mismatch.size());
			auto mismatch = std::vector<int>(n_futures);
			for (auto i = 0; i < param[0].n_times; ++i)
			{
				layer::init(layer1_fluent, layer1, param[0]);
				layer::init(layer2_fluent, layer2, param[1]);
				layer::init(layer3_fluent, layer3, param[2]);
				run(datastream, param, layer1_fluent, layer1, layer2_fluent, layer2, layer3_fluent, layer3, mismatch);
				tools::add(prediction_mismatch, mismatch);
			}
		}
	}
}