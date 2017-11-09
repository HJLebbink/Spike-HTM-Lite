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

		}

		template <
			int N_COLUMNS_L1, int N_BITS_CELL_L1, int N_VISIBLE_SENSORS_L1, int HISTORY_L1,
			int N_COLUMNS_L2, int N_BITS_CELL_L2, int HISTORY_L2,
			arch_t ARCH = arch_t::X64>
		struct network_2Layer
		{
			static const int N_VISIBLE_SENSORS_L2 = N_COLUMNS_L1;
			static const int N_HIDDEN_SENSORS_L1 = N_COLUMNS_L2;
			static const int N_HIDDEN_SENSORS_L2 = 0;

			using P_L1 = Static_Param<N_COLUMNS_L1, N_BITS_CELL_L1, N_VISIBLE_SENSORS_L1, N_HIDDEN_SENSORS_L1, HISTORY_L1, ARCH>;
			using P_L2 = Static_Param<N_COLUMNS_L2, N_BITS_CELL_L2, N_VISIBLE_SENSORS_L2, N_HIDDEN_SENSORS_L2, HISTORY_L2, ARCH>;
		};

		template <typename P1, typename P2>
		int run(
			const std::vector<Layer<P1>::Active_Sensors>& data,
			Layer<P1>& layer1,
			Layer<P2>& layer2,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2)
		{
			static_assert(P1::N_COLUMNS == P2::N_SENSORS, "ERROR: layer1 and layer2 are not matched.");
			const bool LEARN = true;

			layer::init(layer1, param1);
			layer::init(layer2, param2);

			for (int time = 0; time < param1.n_time_steps; ++time)
			{
				//1] get active sensors
				encoder::get_active_sensors<P1>(time, data, layer1.active_sensors);
				layer::priv::add_sensor_noise<P1>(layer1.active_sensors);

				//2] concat previous l2 columns to sensors


				// learn L1
				layer::priv::one_step<LEARN>(layer1.active_sensors, layer1, time, param1);
				// Learn L2
				layer::priv::one_step<LEARN>(layer1.active_columns, layer2, time, param2);
			}

			return 0;
		}
	}
}