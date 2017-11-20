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
#include <random>

#include "tools.ipp"
#include "parameters.ipp"
#include "types.ipp"
#include "encoder.ipp"


namespace htm
{
	namespace datastream
	{
		using namespace ::tools::log;
		using namespace htm::types;

		template <typename P>
		class DataStream
		{
		private:

			using data_type = typename Layer<P>::Active_Visible_Sensors;

			bool use_file_data = false;

			mutable int time = 0;
			std::vector<data_type> file_data;

			mutable int sequence_i = 0;
			mutable int pos_in_sequence = 0;
			mutable int sequence_i_next = 0;
			mutable int pos_in_sequence_next = 0;
			std::vector<std::vector<data_type>> sequences;

			mutable std::default_random_engine random_engine;
			mutable std::uniform_int_distribution<unsigned int> random_number_dist;
			mutable std::uniform_int_distribution<int> random_sequence_dist;

			data_type create_random(float sparcity, unsigned int random_number) const
			{
				data_type result;
				result.clear_all();

				const int range = static_cast<int>(std::floorf(2.0f / sparcity) - 1);

				int i = random::rand_int32(0, range, random_number);
				while (i < P::N_VISIBLE_SENSORS)
				{
					result.set(i, true);
					i += random::rand_int32(0, range, random_number);
				}
				return result;
			}

			std::tuple<int, int> next_time(int pos_in_sequence, int sequence_i) const
			{
				int pos_in_sequence_next = pos_in_sequence + 1;
				int sequence_i_next = sequence_i;

				if (pos_in_sequence_next >= static_cast<int>(this->sequences[sequence_i_next].size()))
				{ // start a new random sequence
					pos_in_sequence_next = 0;
					sequence_i_next = this->random_sequence_dist(this->random_engine);
				}
				return std::make_tuple(pos_in_sequence_next, sequence_i_next);
			}

			std::tuple<int, int> future_time(int pos_in_sequence, int sequence_i, int future) const
			{
				int pos_in_sequence_future = this->pos_in_sequence;
				int sequence_i_future = this->sequence_i;

				for (auto i = 0; i < future; ++i)
				{
					const auto tup = next_time(pos_in_sequence_future, sequence_i_future);
					pos_in_sequence_future = std::get<0>(tup);
					sequence_i_future = std::get<1>(tup);
				}
				return std::make_tuple(pos_in_sequence_future, sequence_i_future);
			}

			void advance_time_current() const
			{
				if (this->use_file_data)
				{
					this->time++;
				}
				else
				{
					this->sequence_i = this->sequence_i_next;
					this->pos_in_sequence = this->pos_in_sequence_next;
				}
			}
			void update_time_next() const 
			{
				if (this->use_file_data)
				{
					// nothing to do;
				}
				else
				{
					const auto tup = next_time(this->pos_in_sequence, this->sequence_i);
					this->pos_in_sequence_next = std::get<0>(tup);
					this->sequence_i_next = std::get<1>(tup);
				}
			}

		public:

			void advance_time() const
			{
				this->advance_time_current();
				this->update_time_next();
			}
			void reset_time() const
			{
				if (this->use_file_data)
				{
					this->time = 0;
				}
				else
				{
					this->sequence_i = this->random_sequence_dist(this->random_engine);
					this->pos_in_sequence = 0;
					this->update_time_next();
				}
			}
			void load_from_file(const std::string& filename, const Dynamic_Param& param)
			{
				this->use_file_data = true;
				this->file_data = encoder::encode_pass_through<P>(filename, param);
				this->reset_time();
			}
			void generate_random_NxR(float sparsity, int n_sequences, int sequence_length)
			{
				std::random_device r;
				this->random_engine = std::default_random_engine(r());
				this->random_sequence_dist = std::uniform_int_distribution<int>(0, n_sequences-1);
				this->random_number_dist = std::uniform_int_distribution<unsigned int>(0, 0xFFFFFFFF);

				this->use_file_data = false;
				for (int sequence_i = 0; sequence_i < n_sequences; ++sequence_i)
				{
					std::vector<data_type> sequence;
					for (int data_i = 0; data_i < sequence_length; ++data_i)
					{
						sequence.push_back(create_random(sparsity, this->random_number_dist(this->random_engine)));
					}
					this->sequences.push_back(sequence);
				}
				this->reset_time();
			}
			bool sensors_predictable(int future) const 
			{
				if (this->use_file_data)
				{
					return true;
				}
				else
				{
					const auto tup = future_time(this->pos_in_sequence_next, this->sequence_i_next, future);
					const int pos_in_sequence_future = std::get<1>(tup);
					return (pos_in_sequence_future != 0);
				}
			}
			bool next_sensors_predictable() const
			{
				return sensors_predictable(1);
			}
			void current_sensors(typename Layer<P>::Active_Sensors& sensor_activity) const
			{
				if (this->use_file_data)
				{
					const int time_step_max = static_cast<int>(this->file_data.size());
					const int i = this->time % time_step_max;
					copy_partial(sensor_activity, this->file_data[i]);
				}
				else
				{
					copy_partial(sensor_activity, this->sequences[this->sequence_i][this->pos_in_sequence]);
				}
			}

			// load the future (next) active sensors
			void future_sensors(typename Layer<P>::Active_Sensors& sensor_activity, int future) const
			{
				if (this->use_file_data)
				{
					const int time_step_max = static_cast<int>(this->file_data.size());
					const int i = (this->time + future) % time_step_max;
					copy_partial(sensor_activity, this->file_data[i]);
				}
				else
				{
					if (future == 0)
					{
						copy_partial(sensor_activity, this->sequences[this->sequence_i][this->pos_in_sequence]);
					}
					else if (future == 1)
					{
						copy_partial(sensor_activity, this->sequences[this->sequence_i_next][this->pos_in_sequence_next]);
					}
					else
					{
						const auto tup = future_time(this->pos_in_sequence_next, this->sequence_i_next, future);
						const int sequence_i_future = std::get<0>(tup);
						const int pos_in_sequence_future = std::get<1>(tup);
						copy_partial(sensor_activity, this->sequences[sequence_i_future][pos_in_sequence_future]);
					}
				}
			}

			// load the future (next) active sensors
			void future_sensors(typename Layer<P>::Active_Sensors& sensor_activity) const
			{
				future_sensors(sensor_activity, 1);
			}

			void future_sensors(typename std::vector<Layer<P>::Active_Sensors>& sensor_activity) const
			{
				const int n_futures = static_cast<int>(sensor_activity.size());
				for (auto future = 0; future < n_futures; ++future)
				{
					future_sensors(sensor_activity[future], future + 1);
				}
			}
		};
	}
}