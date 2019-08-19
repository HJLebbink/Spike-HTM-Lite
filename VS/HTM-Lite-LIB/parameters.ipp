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
#include <algorithm>
#include <iomanip>      // std::setw
#include <ratio>

//#include <immintrin.h> // _may_i_use_cpu_feature
#include <intrin.h>


#include "..\Spike-Tools-Lib\log.ipp"
#include "..\Spike-Tools-Lib\random.ipp"

#include "tools.ipp"

namespace htm
{
	using namespace ::tools::log;
	using namespace ::tools;

	//========================================================================
	enum arch_t
	{
		//Reference implementation, regular c++ code
		X64, 
		//Explicit use of AVX512 instructions (Skylake X)
		AVX512, 
		//Determine instruction set a runtime
		RUNTIME
	};

	arch_t architecture_switch(const arch_t arch)
	{
		if (arch == arch_t::X64) return arch_t::X64;
		if (arch == arch_t::AVX512) return arch_t::AVX512;
		if (arch == arch_t::RUNTIME)
		{
			#if __INTEL_COMPILER
			return (_may_i_use_cpu_feature(_FEATURE_AVX512F) == 1) ? arch_t::AVX512 : arch_t::X64;
			#else
			//log_WARNING("architecture_switch: could not determine desired architecture at runtime: assuming AVX512.\n");
			return arch_t::AVX512;
			#endif
		}
		return arch_t::X64;
	}

	using Permanence = int8_t; // = signed char

	template <int N, int D>
	float constexpr calc_noise_percentage() noexcept
	{
		if constexpr (N == 0) return 0;
		if constexpr (D == 0) return 0;
		return (static_cast<float>(N) / D) * 100;
	}


	//========================================================================
	template <
		int N_COLUMNS_IN,
		int N_BITS_CELL_IN,
		int N_VISIBLE_SENSORS_IN,
		int N_HIDDEN_SENSORS_IN = 0,
		int HISTORY_SIZE_IN = 1,
		arch_t ARCH_IN = arch_t::X64,
		typename SENSOR_NOISE_PERCENT_IN = std::ratio<0,100>
	>
	struct Static_Param
	{
		static_assert(N_COLUMNS_IN > 0, "ERROR: Parameters: provided N_COLUMNS is too small; min N_COLUMNS is 64.");
		static_assert((N_COLUMNS_IN & 0b111111) == 0, "ERROR: Parameters: provided N_COLUMNS is not a multiple of 64.");
		static_assert(N_BITS_CELL_IN > 0, "ERROR: Parameters: provided N_BITS_CELL is too small; min N_BITS_CELL is 1.");
		static_assert(N_BITS_CELL_IN < 7, "ERROR: Parameters: provided N_BITS_CELL is too large; max N_BITS_CELL is 6.");
		static_assert(N_VISIBLE_SENSORS_IN > 0, "ERROR: Parameters: provided N_VISIBLE_SENSORS_IN is too small; min N_VISIBLE_SENSORS_IN is 1.");
//		static_assert(HISTORY_SIZE_IN <= 1, "ERROR: Parameters: provided HISTORY_SIZE_IN is too small; min HISTORY_SIZE_IN is 1.");
//		static_assert(HISTORY_SIZE_IN >= 7, "ERROR: Parameters: provided HISTORY_SIZE_IN is too large; max HISTORY_SIZE_IN is 7.");

		//========================================================================
		//Number of columns in layer: Multiple of 64.
		static constexpr int N_COLUMNS = N_COLUMNS_IN;

		static constexpr int N_VISIBLE_SENSORS = N_VISIBLE_SENSORS_IN;

		static constexpr int N_HIDDEN_SENSORS = N_HIDDEN_SENSORS_IN;

		static constexpr int N_SENSORS = N_VISIBLE_SENSORS + N_HIDDEN_SENSORS;

		// History of 1 (minimum) means that there is one history, history of 7 means there are zeven histories.
		static constexpr int HISTORY_SIZE = HISTORY_SIZE_IN;

		//Whether AVX512 optimized code is used or a reference implementation
		static constexpr arch_t ARCH = ARCH_IN;

		//========================================================================
		#pragma region Spacial Pooler constants

		// number of potential synapses per proximal dendrite. Must be multiple of 64
		static constexpr int SP_N_PD_SYNAPSES = tools::multiple_64(4 * 64);
		static_assert((SP_N_PD_SYNAPSES & 0b111111) == 0, "ERROR: Parameters: provided SP_N_PD_SYNAPSES is not a multiple of 16.");

		//If the permanence value for a synapse is LARGER (NOT EQUAL) than this value,
		//the synapse is said to be active. (if the permanence is negative, then it is not connected)
		static constexpr Permanence SP_PD_PERMANENCE_THRESHOLD = -1;

		//This is a number specifying the minimum number of synapses that must be
		//on in order for a columns to turn ON. The purpose of this is to prevent
		//noise input from activating columns.
		static constexpr int SP_STIMULUS_THRESHOLD = 1;

		//If true, then during inhibition phase the winning columns are selected as
		//the most active columns from the region as a whole.Otherwise, the winning
		//columns are selected with respect to their local neighborhoods.
		static constexpr bool SP_GLOBAL_INHIBITION = true;

		//If a permanence is LARGER (NOT EQUAL) than this value, the synapse is said to be connected.
		static constexpr Permanence SP_PD_CONNECTED_THRESHOLD = -128;

		// Origin of ininstantiated (invalid) values, for debuggin purposes 
		static constexpr int SP_PD_SYNAPSE_ORIGIN_INVALID = -2;

		//A number between 0 and 1.0, used to set a floor on how often a column
		//should have at least stimulusThreshold active inputs. Periodically, each
		//column looks at the overlap duty cycle of all other columns within its
		//inhibition radius and sets its own internal minimal acceptable duty cycle
		//to : minPctDutyCycleBeforeInh * max(other columns' duty cycles).  On each
		//iteration, any column whose overlap duty cycle falls below this computed
		//value will get all of its permanence values boosted up by
		//synPermActiveInc. Raising all permanences in response to a sub-par duty
		//cycle before inhibition allows a cell to search for new inputs when
		//either its previously learned inputs are no longer ever active, or when
		//the vast majority of them have been "hijacked" by other columns.
		static constexpr float SP_MIN_PCT_OVERLAP_DUTY_CYCLES = 0.1f;

		//The period used to calculate duty cycles.Higher values make it take
		//longer to respond to changes in boost or synPerConnectedCell.Shorter
		//values make it more unstable and likely to oscillate.
		static constexpr int SP_DUTY_CYCLE_PERIOD = 10;

		//A number greater or equal than 0.0, used to control the strength of
		//boosting. No boosting is applied if it is set to 0. Boosting strength
		//increases as a function of boostStrength. Boosting encourages columns to
		//have similar activeDutyCycles as their neighbors, which will lead to more
		//efficient use of columns. However, too much boosting may also lead to
		//instability of SP outputs.
		static constexpr float SP_BOOST_STRENGTH = 0.25;

		//Percentage of noise added to sensors; use zero for no noise.
		static constexpr float SP_SENSOR_NOISE_PERCENT = calc_noise_percentage<SENSOR_NOISE_PERCENT_IN::num, SENSOR_NOISE_PERCENT_IN::den>();

		//Whether the spacial pooler is computed in a forward fashion 
		static constexpr bool SP_SYNAPSE_FORWARD = true;

		#pragma endregion
		//========================================================================
		#pragma region Temporal Pooler constants

		//Number of bits to encode the number of cells per column
		static constexpr int N_BITS_CELL = N_BITS_CELL_IN;

		//Number of cells per column.
		static constexpr int N_CELLS_PC = 1 << N_BITS_CELL;

		//Total number of cells in layer, maximum 0x1FFFFFFF
		static constexpr int N_CELLS = N_CELLS_PC * N_COLUMNS;

		static_assert((static_cast<int64_t>(N_COLUMNS) * N_CELLS_PC) <= 0x3FFFFFFF, "ERROR: Parameters: N_CELLS is too large; max N_CELLS is 0x1FFFFFFF.");

		//If the permanence value for a synapse is LARGER (NOT EQUAL) than this value, the synapse is said to be active. (negative values are not connected)
		static constexpr Permanence TP_DD_PERMANENCE_THRESHOLD = -1;

		//If a permanence is LARGER (NOT EQUAL) than this value, the synapse is said to be connected (but not yet active).
		static constexpr Permanence TP_DD_CONNECTED_THRESHOLD = -128;
		static_assert(TP_DD_CONNECTED_THRESHOLD < TP_DD_PERMANENCE_THRESHOLD, "ERROR");

		//The maximum number of segments per cell.
		static constexpr int TP_N_DD_SEGMENTS_MAX = 250;

		//The maximum number of synapses per segment.
		static constexpr int TP_N_DD_SYNAPSES_MAX = tools::multiple_64(250);

		static constexpr int TP_DD_SYNAPSE_ORIGIN_INVALID = -3;
		static constexpr int8_t TP_DD_SEGMENT_DESTINATION_INVALID = -4;

		static constexpr bool TP_SYNAPSE_FORWARD = true;

		#pragma endregion
		//========================================================================
		#pragma endregion
	};

	struct Dynamic_Param
	{
		//Whether learning is switched on or off
		bool learn = true;

		//Number of (time) steps the layer will be updated (run)
		int n_time_steps = 100;

		//How often the complete learning run will be repeated
		int n_times = 1;

		//Whether stuff will be written to stdout
		bool quiet = false;

		//Show the sensor input, the predicted sensor input and the mismatch; zero means off; one means every time step.
		int show_input_and_prediction_interval = 0;

		//Show the predicted mismatch as the average number of misspredicted sensors; zero means off; one means every time step.
		int show_mismatch_interval = 0;

		//Show the predicted mismatch with the number of futures; 1 means only one future prediction, 2 means to future predictions.
		int show_mismatch_n_futures = 1;

		int n_visible_sensors_dim1 = 20;
		int n_visible_sensors_dim2 = 20;

		#pragma region Spacial Pooler stuff

		// if the predicted sensor influx is ABOVE (not equal) this threshold, the sensor is said to be active.
		int sensor_threshold = 0;



		//Initial permanence of a new synapse
		Permanence SP_PD_PERMANENCE_INIT = 1;

		//Amount by which permanence of synapses are incremented during learning.
		Permanence SP_PD_PERMANENCE_INC = 15;

		//Amount by which permanences of synapses are decremented during learning.
		Permanence SP_PD_PERMANENCE_DEC = 10;

		//the permanence increment amount for columns that have not been recently active.
		Permanence SP_PD_PERMANENCE_INC_WEAK = 1;

		//The desired density of active columns within a local inhibition area
		//(the size of which is set by the internally calculated inhibitionRadius,
		//which is in turn determined from the average size of the connected
		//potential pools of all columns). The inhibition logic will insure that
		//at most N columns remain ON within a local inhibition area, where
		//N = localAreaDensity * (total number of columns in inhibition area).
		float SP_LOCAL_AREA_DENSITY = 0.05f;
		#pragma endregion

		#pragma region Temporal Pooler stuff
		//If the number of active connected synapses on a segment is at least
		//this threshold, the segment is said to be active.
		int TP_DD_SEGMENT_ACTIVE_THRESHOLD = 8;
		
		//If the number of potential synapses active on a
		//segment is at least this threshold, it is said to be "matching" and
		//is eligible for learning.
		//If the number of synapses active on a segment is at least this
		//threshold, it is selected as the best matching cell in a bursting column.
		int TP_MIN_DD_ACTIVATION_THRESHOLD = 4;

		//The maximum number of synapses added to a segment during learning.
		int TP_DD_MAX_NEW_SYNAPSE_COUNT = 16;

		//Initial permanence of a new synapse.
		Permanence TP_DD_PERMANENCE_INIT = 10;

		//Amount by which permanences of synapses are incremented during learning.
		Permanence TP_DD_PERMANENCE_INC = 20;

		//Amount by which permanences of synapses are decremented during learning.
		Permanence TP_DD_PERMANENCE_DEC = 19;

		//Amount by which segments are punished for incorrect predictions.
		//A good value is just a bit larger than (the column - level sparsity *
		//permanenceIncrement). So, if column - level sparsity is 2 % and 
		//permanenceIncrement is 0.01, this parameter should be something 
		//like 4 % * 0.01 = 0.0004).
		Permanence TP_DD_PREDICTED_SEGMENT_DEC = 10; // 10 * SP_LOCAL_AREA_DENSITY * TP_DD_PERMANENCE_INC;

		#pragma endregion

		static std::string header_str()
		{
			std::ostringstream result;

			result << "SP_LOCAL_AREA_DENSITY\t";

			result << "SP_PD_PERMANENCE_INIT\t";
			result << "SP_PD_PERMANENCE_INC\t";
			result << "SP_PD_PERMANENCE_DEC\t";
			result << "SP_PD_PERMANENCE_INC_WEAK\t";

			result << "TP_DD_PERMANENCE_INIT\t";
			result << "TP_DD_PERMANENCE_INC\t";
			result << "TP_DD_PERMANENCE_DEC\t";
			result << "TP_DD_PREDICTED_SEGMENT_DEC";

			result << "TP_DD_SEGMENT_ACTIVE_THRESHOLD\t";
			result << "TP_MIN_DD_ACTIVATION_THRESHOLD\t";
			result << "TP_DD_MAX_NEW_SYNAPSE_COUNT\t";

			return result.str();
		}

		std::string str() const
		{
			std::ostringstream result;
			result << std::fixed << std::setw(7) << std::setprecision(5) << this->SP_LOCAL_AREA_DENSITY << "\t";

			result << std::setw(2) << static_cast<int>(this->SP_PD_PERMANENCE_INIT) << "\t";
			result << std::setw(2) << static_cast<int>(this->SP_PD_PERMANENCE_INC) << "\t";
			result << std::setw(2) << static_cast<int>(this->SP_PD_PERMANENCE_DEC) << "\t";
			result << std::setw(2) << static_cast<int>(this->SP_PD_PERMANENCE_INC_WEAK) << "\t";

			result << std::setw(2) << static_cast<int>(this->TP_DD_PERMANENCE_INIT) << "\t";
			result << std::setw(2) << static_cast<int>(this->TP_DD_PERMANENCE_INC) << "\t";
			result << std::setw(2) << static_cast<int>(this->TP_DD_PERMANENCE_DEC) << "\t";
			result << std::setw(2) << static_cast<int>(this->TP_DD_PREDICTED_SEGMENT_DEC) << "\t";

			result << std::setw(2) << this->TP_DD_SEGMENT_ACTIVE_THRESHOLD << "\t";
			result << std::setw(2) << this->TP_MIN_DD_ACTIVATION_THRESHOLD << "\t";
			result << std::setw(2) << this->TP_DD_MAX_NEW_SYNAPSE_COUNT;

			return result.str();
		}
	};
}