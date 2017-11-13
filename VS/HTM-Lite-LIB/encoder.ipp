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
#include <string>
#include <vector>
#include <array>
#include <iostream>		// std::cout
#include <fstream>

#include "..\Spike-Tools-Lib\log.ipp"

#include "tools.ipp"
#include "parameters.ipp"
#include "types.ipp"

namespace htm
{
	namespace encoder
	{
		using namespace ::tools::log;
		using namespace htm::types;

		template <typename P>
		std::vector<Layer<P>::Active_Visible_Sensors> encode_pass_through(
			const std::string& filename,
			const Dynamic_Param& param)
		{
			auto data = std::vector<Layer<P>::Active_Visible_Sensors>();

			std::ifstream input(filename);
			if (!input.good())
			{
				log_WARNING("encoder:encode_pass_through: could not open file ", filename, ".\n");
				return data;
			}

			std::string str;
			int line = 0;

			bool endFile = false;
			while (!endFile)
			{
				Layer<P>::Active_Visible_Sensors item;
				auto pos = 0;

				for (auto i1 = 0; i1 < param.n_visible_sensors_dim1; ++i1)
				{
					std::getline(input, str);
					if (input.bad())
					{
						log_ERROR("encoder::encode_pass_through: i1 = ", i1, "; D2 = ", param.n_visible_sensors_dim2, "; line ", line, ".\n");
					}
					else if (input.eof())
					{
						endFile = true;
					}
					else
					{
						auto strLength = str.length();

						if (strLength == 0) // read optional empty space
						{
							i1--;
						}
						else
						{
							if (false) log_INFO("encoder::encode_pass_through: i1 = ", i1, "; D2 = ", param.n_visible_sensors_dim2, "; line ", line, "; content = ", str, ".\n");
							line++;

							if (str.length() < param.n_visible_sensors_dim1) log_WARNING("encoder:encode_pass_through: str ", str, " is smaller than D1 = ", param.n_visible_sensors_dim1, ".\n");

							for (auto i2 = 0; i2 < std::min(static_cast<int>(str.length()), param.n_visible_sensors_dim1); ++i2)
							{
								switch (str[i2])
								{
									case '0': item.set(pos, false); pos++; break;
									case '1': item.set(pos, true); pos++; break;
									default:
										log_WARNING("encoder:encode_pass_through: found ", str[i2], ".\n");
										break;
								}
							}
						}
					}
				}
				if (pos > 0) // something is present in item: the item is not empty.
				{
					data.push_back(item);
				}
			}
			if (true) log_INFO_DEBUG("encoder::encode_pass_through: done reading ", data.size(), " items from file ", filename, ".\n");
			input.close();
			return data;
		}
	
		template <typename P>
		void add_sensor_noise(
			Layer<P>::Active_Sensors& active_sensors)
		{
			if (P::SP_SENSOR_NOISE_PERCENT > 0)
			{
				const int range = std::floorf(200.0f / P::SP_SENSOR_NOISE_PERCENT) - 1;
				int i = random::rand_int32(0, range);
				while (i < P::N_VISIBLE_SENSORS)
				{
					active_sensors.negate(i);
					i += random::rand_int32(0, range);
				}
			}
		}
	}
}