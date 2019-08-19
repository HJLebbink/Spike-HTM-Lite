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

#include "..\Spike-Tools-LIB\assert.ipp"
#include "..\Spike-Tools-LIB\allocator.ipp"

#include "parameters.ipp"
#include "tools.ipp"

namespace htm
{
	namespace types
	{
		using namespace ::tools::assert;

		namespace priv
		{
			template <typename T>
			//using Allocator = std::allocator<T>;
			using Allocator = ::tools::allocator::Allocator_AVX512<T>;

			template <int SIZE> struct basetype {};

			template <> struct basetype<2> { using type = int8_t; };
			template <> struct basetype<3> { using type = int8_t; };
			template <> struct basetype<4> { using type = int8_t; };
			template <> struct basetype<5> { using type = int8_t; };
			template <> struct basetype<6> { using type = int8_t; };
			template <> struct basetype<7> { using type = int8_t; };
			template <> struct basetype<8> { using type = int8_t; };

			template <> struct basetype<9> { using type = int16_t; };
			template <> struct basetype<10> { using type = int16_t; };
			template <> struct basetype<11> { using type = int16_t; };
			template <> struct basetype<12> { using type = int16_t; };
			template <> struct basetype<13> { using type = int16_t; };
			template <> struct basetype<14> { using type = int16_t; };
			template <> struct basetype<15> { using type = int16_t; };
			template <> struct basetype<16> { using type = int16_t; };

			template <> struct basetype<17> { using type = int32_t; };
			template <> struct basetype<18> { using type = int32_t; };
			template <> struct basetype<19> { using type = int32_t; };
			template <> struct basetype<20> { using type = int32_t; };
			template <> struct basetype<21> { using type = int32_t; };
			template <> struct basetype<22> { using type = int32_t; };
			template <> struct basetype<23> { using type = int32_t; };
			template <> struct basetype<24> { using type = int32_t; };

			template <> struct basetype<25> { using type = int32_t; };
			template <> struct basetype<26> { using type = int32_t; };
			template <> struct basetype<27> { using type = int32_t; };
			template <> struct basetype<28> { using type = int32_t; };
			template <> struct basetype<29> { using type = int32_t; };
			template <> struct basetype<30> { using type = int32_t; };
			template <> struct basetype<31> { using type = int32_t; };
			template <> struct basetype<32> { using type = int32_t; };

			template <> struct basetype<33> { using type = int64_t; };
			template <> struct basetype<34> { using type = int64_t; };
			template <> struct basetype<35> { using type = int64_t; };
			template <> struct basetype<36> { using type = int64_t; };
			template <> struct basetype<37> { using type = int64_t; };
			template <> struct basetype<38> { using type = int64_t; };
			template <> struct basetype<39> { using type = int64_t; };
			template <> struct basetype<40> { using type = int64_t; };
			template <> struct basetype<41> { using type = int64_t; };
			template <> struct basetype<42> { using type = int64_t; };
			template <> struct basetype<43> { using type = int64_t; };
			template <> struct basetype<44> { using type = int64_t; };
			template <> struct basetype<45> { using type = int64_t; };
			template <> struct basetype<46> { using type = int64_t; };
			template <> struct basetype<47> { using type = int64_t; };
			template <> struct basetype<48> { using type = int64_t; };
			template <> struct basetype<49> { using type = int64_t; };
			template <> struct basetype<50> { using type = int64_t; };
			template <> struct basetype<51> { using type = int64_t; };
			template <> struct basetype<52> { using type = int64_t; };
			template <> struct basetype<53> { using type = int64_t; };
			template <> struct basetype<54> { using type = int64_t; };
			template <> struct basetype<55> { using type = int64_t; };
			template <> struct basetype<56> { using type = int64_t; };
			template <> struct basetype<57> { using type = int64_t; };
			template <> struct basetype<58> { using type = int64_t; };
			template <> struct basetype<59> { using type = int64_t; };
			template <> struct basetype<60> { using type = int64_t; };
			template <> struct basetype<61> { using type = int64_t; };
			template <> struct basetype<62> { using type = int64_t; };
			template <> struct basetype<63> { using type = int64_t; };
			template <> struct basetype<64> { using type = int64_t; };
		}

		//========================================================================
		//Bitset: a bitset with an internal vector that stores only the set bit positions.
		template <int SIZE_IN>
		struct Bitset_Sparse
		{
			static constexpr int SIZE = SIZE_IN;
			std::vector<int, priv::Allocator<int>> _data;

			Bitset_Sparse()
			{
				this->_data;
			}
			void set(const std::vector<char>& bitset)
			{
				this->_data.clear();
				for (int i = 0; i < SIZE; ++i) if (bitset[i]) this->_data.push_back(i);
			}
			void set(const Bitset_Sparse<SIZE>& other)
			{
				this->_data.clear();
				if (!other.none()) this->_data = other._data;
			}
			//Test if no bit is set.
			bool none() const
			{
				return this->_data.size() == 0;
			}
			//Count bits set: Returns the number of bits in the bitset that are set (i.e., that have a value of one).
			int count() const
			{
				return static_cast<int>(this->_data.size());
			}
			//Returns the number of bits in the bitset.
			constexpr int size() const
			{
				return SIZE;
			}
			bool contains(const int i) const
			{
				return std::find(this->_data.begin(), this->_data.end(), i) != this->_data.end();
			}
			void set(const int i)
			{
				::tools::assert::assert_msg(((i >= 0) && (i < SIZE)), "types:BitsetSparse:set invalid i ", i);
				this->_data.push_back(i);
			}
			int get(const int i) const
			{
				return this->_data[i];
			}
			//Test if any bit is set.
			bool any() const
			{
				return !this->_data.empty();
			}
			void reset()
			{
				this->_data.clear();
			}
		};

		//========================================================================
		//Bitset Tiny: a bitset with a max size of 64.
		template <int SIZE_IN>
		struct Bitset_Tiny
		{
			static_assert(SIZE_IN <= 64, "ERROR:BitsetTiny: provided SIZE is too large; max SIZE is 64.");
			static constexpr int SIZE = SIZE_IN;

			using internal_type = typename types::priv::basetype<SIZE_IN>::type;
			static constexpr internal_type MASK = (SIZE_IN >= 64) ? 0xFFFFFFFFFFFFFFFF : (1ull << SIZE_IN) - 1;
			internal_type _data;

			// default constructor
			Bitset_Tiny()
			{
				//tools::log::log_INFO("tools::BitsetTiny: SIZE=", SIZE, "; type size=", sizeof(internal_type), "byte(s).");
				this->_data = 0;
			}
			bool get(int i) const
			{
				return this->_data & (1 << i);
			}
			// set all
			void set_all()
			{
				this->_data = MASK;
			}
			void set(const int i, const bool value)
			{
				::tools::assert::assert_msg(((i >= 0) && (i < SIZE)), "types:BitsetSparse:set invalid i ", i);
				if (value)
					this->_data |= (1 << i);
				else
					this->_data &= ~(1 << i);
			}
			void clear_all()
			{
				this->_data = 0;
			}
			bool any() const
			{
				return this->_data > 0;
			}
			bool empty() const
			{
				return this->_data == 0
			}
			int count() const
			{
				int counter = 0;
				for (int i = 0; i < SIZE; ++i) if (this->get(i)) counter++;
				return counter;
			}
		};

		//========================================================================
		//Generic Bitset: Uses an array of char (as boolean) needed for vectorization.
		template <int SIZE_IN>
		struct Bitset
		{
			static constexpr int SIZE = SIZE_IN;
			std::vector<char> _data;

			// default constructor
			Bitset()
			{
				this->_data = std::vector<char>(SIZE);
			}
			bool get(int i) const
			{
				return this->_data[i];
			}
			void set(int i)
			{
				this->_data[i] = true;
			}
			char * data()
			{
				return this->_data.data();
			}
			char const * data() const
			{
				return this->_data.data();
			}
			void reset()
			{
				std::fill(this->_data.begin(), this->_data.end(), 0);
			}
			bool any() const
			{
				for (int i = 0; i < SIZE; ++i) if (this->_data[i]) return true;
				return false;
			}
			//count the number of set bits.
			int count() const
			{
				int counter = 0;
				for (int i = 0; i < SIZE; ++i) if (this->_data[i]) counter++;
				return counter;
			}
		};

		//========================================================================
		//Bitset with 2 dimensions: max of second dimension is 64
		template <int SIZE1_IN, int SIZE2_IN>
		struct Bitset2
		{
			static constexpr int SIZE1 = SIZE1_IN;
			static constexpr int SIZE2 = SIZE2_IN;
			using content_t = Bitset_Tiny<SIZE2>;

			std::vector<content_t> _data;

			// default constructor
			Bitset2()
			{
				this->_data = std::vector<content_t>(SIZE1);
			}
			content_t& operator[] (int i)
			{
				return this->_data[i];
			}
			const content_t& operator[] (int i) const
			{
				return this->_data[i];
			}
		};

		//========================================================================
		//Bitset with a compact representation. Uses an array of int (with 32 booleans) needed for vectorization. 
		template <int SIZE_IN>
		struct Bitset_Compact
		{
			static constexpr int SIZE = SIZE_IN;
			static constexpr int N_BLOCKS = tools::n_blocks_32(SIZE);
			using t1 = priv::Allocator<int>;
			std::vector<int, t1> _data;

			// default constructor
			Bitset_Compact()
			{
				this->_data = std::vector<int, t1>(N_BLOCKS, 0);
			}
			bool get(int i) const
			{
				const int pos_in_array = i >> 5;
				const int pos_in_word = i & 0b11111;
				return (this->_data[pos_in_array] & (1 << pos_in_word)) != 0;
			}
			void set(const int i, const bool value)
			{
				const int pos_in_array = i >> 5;
				const int pos_in_word = i & 0b11111;

				if (value)
					this->_data[pos_in_array] |= (1 << pos_in_word);
				else
					this->_data[pos_in_array] &= ~(1 << pos_in_word);
			}
			void negate(int i)
			{
				const int pos_in_array = i >> 5;
				const int pos_in_word = i & 0b11111;
				this->_data[pos_in_array] ^= (1 << pos_in_word);
			}
			int * data()
			{
				return this->_data.data();
			}
			int const * data() const
			{
				return this->_data.data();
			}
			void clear_all()
			{
				std::fill(this->_data.begin(), this->_data.end(), 0);
			}
			bool any() const
			{
				for (int i = 0; i < N_BLOCKS; ++i) if (this->_data[i] > 0) return true;
				return false;
			}
			//count the number of set bits.
			int count() const
			{
				int result = 0;
				for (int i = 0; i < N_BLOCKS; ++i) result += _mm_popcnt_u32(static_cast<unsigned int>(this->_data[i]));
				return result;
			}
		};

		//========================================================================
		//Set with segments and their activity
		struct Segments_Set
		{
			std::vector<uint64_t> _data;

			//default constructor
			Segments_Set() = default;
			//copy constructor
			Segments_Set(const Segments_Set& other)
			{
				this->_data = other._data;
			}
			int count() const
			{
				return static_cast<int>(this->_data.size());
			}
			int get_id(const int i) const
			{
				return static_cast<int>(this->_data[i] >> 32);
			}
			int get_activity(const int i) const
			{
				return static_cast<int>(this->_data[i] & 0xFFFFFFFF);
			}
			void add(const int segment_i, const int activity)
			{
				this->_data.push_back(static_cast<uint64_t>(segment_i) << 32 | static_cast<uint64_t>(activity));
			}
			bool any() const
			{
				return !this->_data.empty();
			}
			std::tuple<int, int> highest_activity() const
			{
				int best_matching_segment = -1;
				int best_segment_activity = 0;
				for (int i = 0; i < this->count(); ++i)
				{
					const int activity = this->get_activity(i);
					if (activity > best_segment_activity)
					{
						best_segment_activity = activity;
						best_matching_segment = this->get_id(i);
					}
				}
				return std::make_tuple(best_matching_segment, best_segment_activity);
			}
			void reset()
			{
				this->_data.clear();
			}
		};

		//========================================================================
		//Add history to the provided type T
		template <typename T, int SIZE_IN>
		class History
		{
		public:
			static constexpr int SIZE = SIZE_IN;
			static_assert(SIZE > 1, "ERROR:History: SIZE of history has to be larger than 1.");
			using base_type = typename T;

		private:
			std::array<base_type, SIZE> _data;
			int _current_index;

			mutable bool _uptodate;
			mutable std::vector<int> _sparse_history;

			int get_prev_index(const int i) const
			{
				const int prev_index = i - 1;
				return (prev_index < 0) ? SIZE - 1 : prev_index;
			}
			int get_historic_index(const int h) const
			{
				const int hist_index = this->_current_index - h;
				return (hist_index < 0) ? hist_index + SIZE : hist_index;
			}

		public:
			// default constructor
			History()
			{
				this->_uptodate = false;
				this->_current_index = 0;
			}
			void advance_time()
			{
				this->_current_index++;
				if (this->_current_index >= SIZE) this->_current_index = 0;
				this->_uptodate = false;
			}
			const base_type& current() const
			{
				return this->_data[this->_current_index];
			}
			base_type& current()
			{
				return this->_data[this->_current_index];
			}
			// shorthand for prev(1)
			const base_type& prev() const
			{
				return this->_data[this->get_prev_index(this->_current_index)];
			}
			const base_type& prev(const int delay) const
			{
				if (delay == 1)
				{
					return this->prev();
				}
				else
				{
					return this->_data[this->get_historic_index(delay)];
				}
			}
			bool get(const int i, const int delay) const
			{
				return this->prev(delay).get(i);
			}
			void reset()
			{
				this->_uptodate = false;
				for (int i = 0; i < SIZE; ++i) this->_data[i].reset();
			}
			template <int SIZE1, int SIZE2>
			void set_current(const Bitset2<SIZE1, SIZE2>& current)
			{
				copy(this.current(), current);
			}

			//Any on all histories.
			bool any_history() const
			{
				int prev_index = this->_current_index;
				for (int i = 0; i < SIZE - 1; ++i)
				{
					prev_index = this->get_prev_index(prev_index);
					if (this->_data[prev_index].any()) return true;
				}
				return false;
			}
			//Any on the current.
			bool any_current() const
			{
				return this->current().any();
			}
			const std::vector<int>& get_sparse_history() const
			{
				if (!this->_uptodate)
				{
					this->_sparse_history.clear();

					int prev_index = this->_current_index;
					for (int i = 0; i < SIZE - 1; ++i)
					{
						prev_index = this->get_prev_index(prev_index);
						const auto& d = this->_data[prev_index]._data;
						for (int j = 0; j < d.size(); ++j)
						{
							this->_sparse_history.push_back(htm::tools::create_delay_and_cell_id(d[j], i + 1));
						}
					}
					this->_uptodate = true;
				}

				return this->_sparse_history;
			}
		};

		template <int SIZE_IN, int HISTORY_SIZE_IN>
		class Bitset_Hist8
		{
		public:
			static constexpr int SIZE = SIZE_IN;
			static constexpr int HISTORY_SIZE = HISTORY_SIZE_IN;
			static_assert(HISTORY_SIZE >= 1, "invalid HISTORY_SIZE_IN");
			static_assert(HISTORY_SIZE <= 8, "invalid HISTORY_SIZE_IN");

			std::vector<char> _data;

			// default constructor
			Bitset_Hist8()
			{
				this->_data = std::vector<char>(SIZE);
			}
			char& operator[] (int i)
			{
				return this->_data[i];
			}
			const char& operator[] (int i) const
			{
				return this->_data[i];
			}
			char * data()
			{
				return this->_data.data();
			}
			char const * data() const
			{
				return this->_data.data();
			}
			void reset()
			{
				std::fill(this->_data.begin(), this->_data.end(), 0);
			}
			void advance_time()
			{
				for (int i = 0; i < SIZE; ++i) this->_data[i] <<= 1;
			}

			bool get(const int i, const int delay) const
			{
				assert_msg(delay < HISTORY_SIZE, "ERROR:Bitset_Hist: invalid delay ", delay, "; HISTORY_SIZE = ", HISTORY_SIZE);
				return ((this->_data[i] & (1 << delay)) != 0);
			}

			template <int SIZE1, int SIZE2>
			void set_current(const Bitset2<SIZE1, SIZE2>& current)
			{
				auto counter = 0;
				for (auto column_i = 0; column_i < SIZE1; ++column_i)
				{
					const auto & a2 = current[column_i];
					for (auto cell_i = 0; cell_i < SIZE2; ++cell_i)
					{
						if (a2.get(cell_i)) this->_data[counter] |= 1;
						counter++;
					}
				}
			}

			bool any() const
			{
				for (int i = 0; i < SIZE; ++i) if (this->_data[i]) return true;
				return false;
			}
			//count the number of set bits.
			int count() const
			{
				int counter = 0;
				for (int i = 0; i < SIZE; ++i) if (this->_data[i]) counter++;
				return counter;
			}
		};

		//========================================================================
		// Fluent properties of a layer, the non learned states
		template <typename P>
		struct Layer_Fluent
		{
			using Active_Cells = Bitset_Hist8<P::N_CELLS, P::HISTORY_SIZE + 1>;
			using Winner_Cells = History<Bitset_Sparse<P::N_CELLS>, P::HISTORY_SIZE + 1>;
			using Active_Columns = Bitset_Compact<P::N_COLUMNS>;
			using Active_Sensors = Bitset_Compact<P::N_SENSORS>;
			using Active_Visible_Sensors = Bitset_Compact<P::N_VISIBLE_SENSORS>;
			using Active_Segments = History<Segments_Set, 2>;
			using Matching_Segments = History<Segments_Set, 2>;

			//Pseudo random number used only by the column.
			std::vector<unsigned int> random_number;

			Active_Cells active_cells; //32MB for 1M columns times History
			Winner_Cells winner_cells;

			//Segments that currently and previously are active.
			std::vector<Active_Segments> active_dd_segments = std::vector<Active_Segments>(P::N_COLUMNS);

			//Segments that currently and previously are matching.
			std::vector<Matching_Segments> matching_dd_segments = std::vector<Matching_Segments>(P::N_COLUMNS);

			//Last time step synapse was active.
			std::vector<std::vector<int>> dd_synapse_active_time = std::vector<std::vector<int>>(P::N_COLUMNS);


			#pragma region Used by SP only
			//Number of iterations.
			int iteration_num = 0;

			//Number of iterations while learning.
			int iteration_learn_num = 0;

			#pragma region Used by TP only
			Active_Sensors active_sensors;
			Active_Columns active_columns;
			#pragma endregion 

			//moving average denoting the frequency of column activation: used by boosting
			std::vector<float> sp_active_duty_cycles = std::vector<float>(P::N_COLUMNS, 0.0f);

			//moving average denoting the frequency of the column's overlap value being at least equal to the proximal segment activation threshold
			std::vector<float> sp_overlap_duty_cycles = std::vector<float>(P::N_COLUMNS, 0.0f);
			std::vector<float> sp_min_overlap_duty_cycles = std::vector<float>(P::N_COLUMNS, 0.0f);

			// default constructor
			Layer_Fluent()
			{
				this->random_number = std::vector<unsigned int>(P::N_COLUMNS);
				for (int column_i = 0; column_i < P::N_COLUMNS; ++column_i) this->random_number[column_i] = random::rdrand32();
			}
		};

		//========================================================================
		template <typename P>
		struct Layer_Persisted
		{
			//Number of segments this column currently has in use.
			std::vector<int> dd_segment_count = std::vector<int>(P::N_COLUMNS);

			//Cell to which this segment belongs to.
			std::vector<std::vector<int8_t>> dd_segment_destination = std::vector<std::vector<int8_t>>(P::N_COLUMNS);

			using t5 = std::vector<std::vector<Permanence, priv::Allocator<Permanence>>>;
			//Permanence of the synapses of the the provided segment index.
			std::vector<t5> dd_synapse_permanence_sf = std::vector<t5>(P::N_COLUMNS);

			using t6 = std::vector<std::vector<int, priv::Allocator<int>>>;
			//Originating cell id of the synapses of the provided segment index.
			std::vector<t6> dd_synapse_delay_origin_sf = std::vector<t6>(P::N_COLUMNS);

			//Number of synapses this column currently has in use; this int is smaller than TP_N_DD_SYNAPSES.
			std::vector<std::vector<int>> dd_synapse_count_sf = std::vector<std::vector<int>>(P::N_COLUMNS);

			//Current boost factor.
			std::vector<float> boost_factor = std::vector<float>(P::N_COLUMNS, 1.0f);

			#pragma region synapse forward
			using t1 = std::vector<Permanence, priv::Allocator<Permanence>>;
			//Proximal dendrite synapse permanence: iterate over columns: column pushes
			std::vector<t1> sp_pd_synapse_permanence_sf;

			using t2 = std::vector<int, priv::Allocator<int>>;
			//Proximal dendrite synapse origin cell ID: iterate over columns: column pushes
			std::vector<t2> sp_pd_synapse_origin_sensor_sf;
			#pragma endregion

			#pragma region synapse backward
			using t3 = std::vector<int, priv::Allocator<int>>;
			//Proximal dendrite synapse origin cell ID: iterate over sensors: sensor pulls
			std::vector<t3> sp_pd_destination_column_sb;

			using t4 = std::vector<Permanence, priv::Allocator<Permanence>>;
			//Proximal dendrite synapse permanence: iterate over sensors: sensor pulls
			std::vector<t4> sp_pd_synapse_permanence_sb;

			//Number of dendrite synapses: iterate over sensors: sensor pulls
			std::vector<int> sp_pd_synapse_count_sb;
			#pragma endregion

			// default constructor
			Layer_Persisted()
			{
				if (P::SP_SYNAPSE_FORWARD)
				{
					this->sp_pd_synapse_permanence_sf = std::vector<t1>(P::N_COLUMNS, t1(P::SP_N_PD_SYNAPSES, P::SP_PD_CONNECTED_THRESHOLD));
					this->sp_pd_synapse_origin_sensor_sf = std::vector<t2>(P::N_COLUMNS, t2(P::SP_N_PD_SYNAPSES, P::SP_PD_SYNAPSE_ORIGIN_INVALID));
				}
				else
				{
					this->sp_pd_destination_column_sb = std::vector<t3>(P::N_SENSORS);
					this->sp_pd_synapse_permanence_sb = std::vector<t4>(P::N_SENSORS);
					this->sp_pd_synapse_count_sb = std::vector<int>(P::N_SENSORS);
				}
			}
		};

		template <bool LEARN, typename P>
		//using Layer_Persisted_C = typename std::conditional<LEARN, Layer_Persisted<P>, const Layer_Persisted<P>>::type;
		using Layer_Persisted_C = typename std::conditional<LEARN, Layer_Persisted<P>, Layer_Persisted<P>>::type;


		#pragma region Copy

		template <int SIZE1, int SIZE2>
		void copy(Bitset_Compact<SIZE1 * SIZE2>& out, const Bitset2<SIZE1, SIZE2>& in)
		{
			out.clear_all();
			auto counter = 0;

			//TODO: this can be done by just copying in to the correct position in out

			for (auto column_i = 0; column_i < SIZE1; ++column_i)
			{
				const auto & a2 = in[column_i];
				for (auto cell_i = 0; cell_i < SIZE2; ++cell_i)
				{
					if (a2.get(cell_i)) out.set(counter);
					counter++;
				}
			}

			#if _DEBUG
			if (count(in) != out.count()) log_ERROR("tools::copy: BUG: count in ", count(in), "; count out ", out.count());
			#endif
		}

		template <int SIZE1, int SIZE2>
		void copy(Bitset_Sparse<SIZE1 * SIZE2>& out, const Bitset2<SIZE1, SIZE2>& in)
		{
			out.reset();
			auto counter = 0;

			for (auto column_i = 0; column_i < SIZE1; ++column_i)
			{
				const auto & a2 = in[column_i];
				for (auto cell_i = 0; cell_i < SIZE2; ++cell_i)
				{
					if (a2.get(cell_i)) out.set(counter);
					counter++;
				}
			}

			#if _DEBUG
			if (count(in) != out.count()) log_ERROR("tools::copy: BUG: count in ", count(in), "; count out ", out.count());
			#endif
		}

		template <int SIZE>
		void copy(Bitset<SIZE>& out, const Bitset<SIZE>& in)
		{
			out._data = in._data;
		}

		template <int SIZE>
		void copy(Bitset_Compact<SIZE>& out, const Bitset_Compact<SIZE>& in)
		{
			out._data = in._data;
		}

		template <int SIZE1, int SIZE2>
		void copy_partial(Bitset_Compact<SIZE1>& out, const Bitset_Compact<SIZE2>& in)
		{
			static_assert(SIZE2 <= SIZE1, "ERROR: copy_partial: Bitset in is larger than Bitset out.");
			for (int i = 0; i < Bitset_Compact<SIZE2>::N_BLOCKS; ++i)
			{
				out._data[i] = in._data[i];
			}
		}

		template <int SIZE>
		void copy(Bitset_Tiny<SIZE>& out, const Bitset_Tiny<SIZE>& in)
		{
			out._data = in._data;
		}

		template <int SIZE>
		void copy(Bitset_Sparse<SIZE>& out, const Bitset_Sparse<SIZE>& in)
		{
			if (in.none())
			{
				out.reset();
			}
			else
			{
				//tools::log::log_INFO("tools::copy: BitsetSparse in.size=", in.size());
				out.set(in);
			}
		}

		template <class T, int SIZE>
		void copy(std::array<T, SIZE>& __restrict out, const std::array<T, SIZE>& __restrict in)
		{
			for (auto i = 0; i < SIZE; ++i) out[i] = in[i];
		}

		#pragma endregion

		#pragma region Swap
		template <int SIZE>
		void swap(Bitset_Tiny<SIZE>& a, Bitset_Tiny<SIZE>& b)
		{
			std::swap(a._data, b._data);
		}

		template <int SIZE>
		void swap(Bitset_Sparse<SIZE>& a, Bitset_Sparse<SIZE>& b)
		{
			std::swap(a._data, b._data);
		}

		template <int SIZE>
		void swap(Bitset_Compact<SIZE>& a, Bitset_Compact<SIZE>& b)
		{
			std::swap(a._data, b._data);
		}

		void swap(Segments_Set& a, Segments_Set& b)
		{
			std::swap(a._data, b._data);
		}

		#pragma endregion

		#pragma region Misc
		template <class T>
		T max(const std::vector<T>& v)
		{
			return *std::max_element(v.begin(), v.end());
		}
		template <class T, int SIZE>
		void set(std::array<T, SIZE>& a, T value)
		{
			for (auto i = 0; i < SIZE; ++i) a[i] = value;
		}

		void clear(std::vector<int>& a)
		{
			for (auto i = 0; i < static_cast<int>(a.size()); ++i) a[i] = 0;
		}

		template <int SIZE1, int SIZE2>
		void reset(Bitset2<SIZE1, SIZE2>& a)
		{
			for (auto i1 = 0; i1 < SIZE1; ++i1) a[i1].reset();
		}

		template <int SIZE1, int SIZE2>
		bool is_empty(const Bitset2<SIZE1, SIZE2>& a)
		{
			for (auto i1 = 0; i1 < SIZE1; ++i1) if (a[i1].any()) return false;
			return true;
		}

		template <class T, int SIZE>
		T max(const std::array<T, SIZE>& a)
		{
			static_assert(SIZE > 0);
			T max = a[0];
			for (auto i = 1; i < SIZE; ++i) if (a[i] > max) max = a[i];
			return max;
		}

		template <int SIZE1, int SIZE2>
		int count(const Bitset2<SIZE1, SIZE2>& a)
		{
			int count = 0;
			for (auto i = 0; i < SIZE1; ++i) count += a[i].count();
			return count;
		}

		template <int SIZE1, int SIZE2>
		std::vector<int> convert_to_int_vector(const Bitset2<SIZE1, SIZE2>& cells)
		{
			std::vector<int> global_cell_ids;
			for (auto column_i = 0; column_i < SIZE1; ++column_i)
			{
				for (auto cell_i = 0; cell_i < SIZE2; ++cell_i)
				{
					if (cells[column_i][cell_i])
					{
						global_cell_ids.push_back(local_2_global_cell<P::N_BITS_CELL>(cell_i, column_i));
					}
				}
			}

			#if _DEBUG
			if (tools::count(cells) != global_cell_ids.size()) ::tools::log::log_ERROR("tools:convert_to_int_vector: BUG");
			#endif

			return global_cell_ids;
		}
		#pragma endregion
	}
}