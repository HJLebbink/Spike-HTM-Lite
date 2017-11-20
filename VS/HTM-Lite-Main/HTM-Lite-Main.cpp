// C++port of Nupic HTM with the aim of being lite and fast
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

#include <chrono>
#include <stdio.h>  // for printf

#include "..\Spike-Tools-Lib\log.ipp"
#include "..\Spike-Tools-Lib\timing.ipp"
#include "..\Spike-Tools-Lib\profiler.ipp"

#include "..\HTM-Lite-LIB\parameters.ipp"
#include "..\HTM-Lite-LIB\types.ipp"
#include "..\HTM-Lite-LIB\layer.ipp"
#include "..\HTM-Lite-LIB\datastream.ipp"
#include "..\HTM-Lite-LIB\swarm.ipp"
#include "..\HTM-Lite-LIB\network.ipp"

using namespace htm;
using namespace htm::types;
using namespace htm::datastream;

inline void test_1layer()
{
	// static properties: properties that need to be known at compile time:

	const int N_SENSORS_DIM1 = 20;
	const int N_SENSORS_DIM2 = 20;
	const int N_VISIBLE_SENSORS = N_SENSORS_DIM1 * N_SENSORS_DIM2;
	const int N_HIDDEN_SENSORS = 0;

	//const int N_BLOCKS = 8; // 512 columns: use sparsity 0.05 -> 25
	//const int N_BLOCKS = 4096; // 262144 columns: use sparsity of 0.005 -> 1310
	//const int N_BLOCKS = 16384; // 1048576 columns: use sparsity of 0.002 -> 2048
	const int N_BLOCKS = 1 * 8;
	const int N_COLUMNS = 64 * N_BLOCKS;
	const int N_BITS_CELL = 4;
	const int HISTORY_SIZE = 1;

	//const arch_t ARCH = arch_t::X64;
	const arch_t ARCH = arch_t::RUNTIME;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 1000;
	param1.n_times = 1;
	param1.n_visible_sensors_dim1 = N_SENSORS_DIM1;
	param1.n_visible_sensors_dim2 = N_SENSORS_DIM2;

	param1.show_mismatch_interval = 40;
	param1.show_mismatch_n_futures = 1;

	param1.show_input_and_prediction_interval = 0;

	using P = Static_Param<N_COLUMNS, N_BITS_CELL, N_VISIBLE_SENSORS, N_HIDDEN_SENSORS, HISTORY_SIZE, ARCH>;
	Layer<P> layer;
	DataStream<P> datastream;

	const bool load_from_file = true;
	if (load_from_file)
	{
		param1.SP_LOCAL_AREA_DENSITY = 0.10000;

		param1.SP_PD_PERMANENCE_INIT = 6;
		param1.SP_PD_PERMANENCE_INC = 7;
		param1.SP_PD_PERMANENCE_DEC = 23;
		param1.SP_PD_PERMANENCE_INC_WEAK = 18;

		param1.TP_DD_PERMANENCE_INIT = 12;
		param1.TP_DD_PERMANENCE_INC = 26;
		param1.TP_DD_PERMANENCE_DEC = 15;
		param1.TP_DD_PREDICTED_SEGMENT_DEC = 15;

		param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 22;
		param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 16;
		param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 27;

		const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
		//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
		datastream.load_from_file(input_filename, param1);
	}
	else
	{
		const float sparsity = 0.05f;
		const int n_sequences = 18;
		const int sequence_length = 3;

		param1.SP_LOCAL_AREA_DENSITY = 0.095;

		param1.SP_PD_PERMANENCE_INIT = 7;
		param1.SP_PD_PERMANENCE_INC = 14;
		param1.SP_PD_PERMANENCE_DEC = 28;
		param1.SP_PD_PERMANENCE_INC_WEAK = 16;

		param1.TP_DD_PERMANENCE_INIT = 18;
		param1.TP_DD_PERMANENCE_INC = 18;
		param1.TP_DD_PERMANENCE_DEC = 11;
		param1.TP_DD_PREDICTED_SEGMENT_DEC = 22;

		param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 23;
		param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 14;
		param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 24;

		datastream.generate_random_NxR(sparsity, n_sequences, sequence_length);
	}

	std::vector<int> prediction_mismatch(param1.show_mismatch_n_futures);
	htm::layer::run_multiple_times(datastream, layer, param1, prediction_mismatch);
}

inline void test_2layers()
{
	// static properties: properties that need to be known at compile time:

	const int N_SENSORS_DIM1 = 20;
	const int N_SENSORS_DIM2 = 20;
	const int N_VISIBLE_SENSORS_L1 = N_SENSORS_DIM1 * N_SENSORS_DIM2;

	const int N_BLOCKS_L1 = 2 * 8;
	const int N_COLUMNS_L1 = 64 * N_BLOCKS_L1;
	const int N_BITS_CELL_L1 = 4;
	const int HISTORY_SIZE_L1 = 7;

	const int N_BLOCKS_L2 = 1 * 8;
	const int N_COLUMNS_L2 = 64 * N_BLOCKS_L2;
	const int N_BITS_CELL_L2 = 4;
	const int HISTORY_SIZE_L2 = 7;

	const arch_t ARCH = arch_t::RUNTIME;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 3000;
	param1.n_times = 10;
	param1.n_visible_sensors_dim1 = N_SENSORS_DIM1;
	param1.n_visible_sensors_dim2 = N_SENSORS_DIM2;

	param1.show_input_and_prediction_interval = 1;
	param1.show_mismatch_interval = 200;

	param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 19;
	param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 14;
	param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 25;

	Dynamic_Param param2(param1);

	const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
	using Network_Config = network::network_2Layer<
		N_VISIBLE_SENSORS_L1,
		N_COLUMNS_L1, N_BITS_CELL_L1, HISTORY_SIZE_L1,
		N_COLUMNS_L2, N_BITS_CELL_L2, HISTORY_SIZE_L2,
		ARCH>;

	using P1 = Network_Config::P_L1;
	using P2 = Network_Config::P_L2;

	Layer<P1> layer1;
	Layer<P2> layer2;

	DataStream<P1> datastream;
	datastream.load_from_file(input_filename, param1);

	auto param = std::array<Dynamic_Param, 2>{param1, param2};
	std::vector<int> prediction_mismatch(1);
	htm::network::run_multiple_times(datastream, param, layer1, layer2, prediction_mismatch);
}

inline void test_3layers()
{
	// static properties: properties that need to be known at compile time:

	const int N_SENSORS_DIM1 = 20;
	const int N_SENSORS_DIM2 = 20;
	const int N_VISIBLE_SENSORS_L1 = N_SENSORS_DIM1 * N_SENSORS_DIM2;

	const int N_BLOCKS_L1 = 2 * 8;
	const int N_COLUMNS_L1 = 64 * N_BLOCKS_L1;
	const int N_BITS_CELL_L1 = 4;
	const int HISTORY_SIZE_L1 = 7;

	const int N_BLOCKS_L2 = 1 * 8;
	const int N_COLUMNS_L2 = 64 * N_BLOCKS_L2;
	const int N_BITS_CELL_L2 = 4;
	const int HISTORY_SIZE_L2 = 7;

	const int N_BLOCKS_L3 = 1 * 8;
	const int N_COLUMNS_L3 = 64 * N_BLOCKS_L3;
	const int N_BITS_CELL_L3 = 4;
	const int HISTORY_SIZE_L3 = 7;

	const arch_t ARCH = arch_t::RUNTIME;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 500;
	param1.n_times = 10;
	param1.n_visible_sensors_dim1 = N_SENSORS_DIM1;
	param1.n_visible_sensors_dim2 = N_SENSORS_DIM2;

	param1.show_input_and_prediction_interval = 0;
	param1.show_mismatch_interval = 20;

	Dynamic_Param param2(param1);
	Dynamic_Param param3(param1);

	const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";

	param1.SP_LOCAL_AREA_DENSITY = 0.05963;
	param1.SP_PD_PERMANENCE_INIT = 23;
	param1.SP_PD_PERMANENCE_INC = 8;
	param1.SP_PD_PERMANENCE_DEC = 30;
	param1.SP_PD_PERMANENCE_INC_WEAK = 26;
	param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 18;
	param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 21;
	param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 38;
	param1.TP_DD_PERMANENCE_INC = 26;
	param1.TP_DD_PERMANENCE_DEC = 21;
	param1.TP_DD_PREDICTED_SEGMENT_DEC = 7;

	param2.SP_LOCAL_AREA_DENSITY = 0.01196;
	param2.SP_PD_PERMANENCE_INIT = 4;
	param2.SP_PD_PERMANENCE_INC = 5;
	param2.SP_PD_PERMANENCE_DEC = 1;
	param2.SP_PD_PERMANENCE_INC_WEAK = 27;
	param2.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 21;
	param2.TP_MIN_DD_ACTIVATION_THRESHOLD = 5;
	param2.TP_DD_MAX_NEW_SYNAPSE_COUNT = 23;
	param2.TP_DD_PERMANENCE_INC = 14;
	param2.TP_DD_PERMANENCE_DEC = 30;
	param2.TP_DD_PREDICTED_SEGMENT_DEC = 23;

	param3.SP_LOCAL_AREA_DENSITY = 0.01781;
	param3.SP_PD_PERMANENCE_INIT = 25;
	param3.SP_PD_PERMANENCE_INC = 12;
	param3.SP_PD_PERMANENCE_DEC = 10;
	param3.SP_PD_PERMANENCE_INC_WEAK = 11;
	param3.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 7;
	param3.TP_MIN_DD_ACTIVATION_THRESHOLD = 18;
	param3.TP_DD_MAX_NEW_SYNAPSE_COUNT = 32;
	param3.TP_DD_PERMANENCE_INC = 9;
	param3.TP_DD_PERMANENCE_DEC = 10;
	param3.TP_DD_PREDICTED_SEGMENT_DEC = 8;

	using Network_Config = network::network_3Layer<
		N_VISIBLE_SENSORS_L1,
		N_COLUMNS_L1, N_BITS_CELL_L1, HISTORY_SIZE_L1,
		N_COLUMNS_L2, N_BITS_CELL_L2, HISTORY_SIZE_L2,
		N_COLUMNS_L3, N_BITS_CELL_L3, HISTORY_SIZE_L3,
		ARCH>;

	using P1 = Network_Config::P_L1;
	using P2 = Network_Config::P_L2;
	using P3 = Network_Config::P_L3;

	Layer<P1> layer1;
	Layer<P2> layer2;
	Layer<P3> layer3;

	DataStream<P1> datastream;
	datastream.load_from_file(input_filename, param1);

	auto param = std::array<Dynamic_Param, 3>{param1, param2, param3};
	std::vector<int> prediction_mismatch(1);
	htm::network::run_multiple_times(datastream, param, layer1, layer2, layer3, prediction_mismatch);
}

inline void test_swarm_1layer()
{
	// static properties: properties that need to be known at compile time:

	const int N_SENSOR_DIM1 = 20;
	const int N_SENSOR_DIM2 = 20;
	const int N_VISIBLE_SENSORS = N_SENSOR_DIM1 * N_SENSOR_DIM2;
	const int N_HIDDEN_SENSORS = 0;

	//const int N_BLOCKS = 9; // 1024 columns
	//const int N_BLOCKS = 4096; // 262144 columns
	//const int N_BLOCKS = 16384; // 1048576 columns
	const int N_BLOCKS = 4 * 8;
	const int N_COLUMNS = 64 * N_BLOCKS;
	const int N_BITS_CELL = 4;
	const int HISTORY_SIZE = 1;

	const arch_t ARCH = arch_t::RUNTIME;

	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 1000;
	param1.n_times = 1;
	param1.quiet = true;

	const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
	//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
	//htm::swarm::run_random(input_filename, param);

	// Swarm options
	const std::string outputfilename = "C:/Temp/VS/htm_ga_results_l1.csv";
	swarm::Swarm_Options options;
	options.population_size = 100;
	options.n_epochs = 100;
	options.mutation_rate = 0.1;
	options.quiet = false;

	std::ofstream file; // out file stream
	options.outputfile.open(outputfilename, std::ios::out | std::ios::trunc);
	if (options.outputfile.good())
		log_INFO("swarm:run_ga: writing results to file ", outputfilename, "\n");
	else
		log_WARNING("swarm:run_ga: could not write to file ", outputfilename, "\n");

	using P = Static_Param<N_COLUMNS, N_BITS_CELL, N_VISIBLE_SENSORS, N_HIDDEN_SENSORS, HISTORY_SIZE, ARCH>;
	auto param = std::array<Dynamic_Param, 1>{param1};
	DataStream<P> datastream;

	const bool load_from_file = true;
	if (load_from_file)
	{
		const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
		//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
		datastream.load_from_file(input_filename, param1);
	}
	else
	{
		const float sparsity = 0.05f;
		const int n_sequences = 16;
		const int sequence_length = 3;
		datastream.generate_random_NxR(sparsity, n_sequences, sequence_length);
	}
	htm::swarm::run_ga<P>(datastream, param, options);
}

inline void test_swarm_2layers()
{
	// static properties: properties that need to be known at compile time:

	//const int N_BLOCKS = 8; // 512 columns: use sparsity 0.05 -> 25
	//const int N_BLOCKS = 4096; // 262144 columns: use sparsity of 0.005 -> 1310
	//const int N_BLOCKS = 16384; // 1048576 columns: use sparsity of 0.002 -> 2048
	const int N_BLOCKS_L1 = 1 * 2 * 8;
	const int N_COLUMNS_L1 = 64 * N_BLOCKS_L1;
	const int N_BITS_CELL_L1 = 4;
	const int N_SENSORS_DIM1 = 20;
	const int N_SENSORS_DIM2 = 20;
	const int N_VISIBLE_SENSORS_L1 = N_SENSORS_DIM1 * N_SENSORS_DIM2;
	const int HISTORY_SIZE_L1 = 7;

	const int N_BLOCKS_L2 = 1 * 1 * 8;
	const int N_COLUMNS_L2 = 64 * N_BLOCKS_L2;
	const int N_BITS_CELL_L2 = 4;
	const int HISTORY_SIZE_L2 = 7;

	const arch_t ARCH = arch_t::RUNTIME;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 500;
	param1.n_times = 10;
	param1.n_visible_sensors_dim1 = N_SENSORS_DIM1;
	param1.n_visible_sensors_dim2 = N_SENSORS_DIM2;

	param1.show_input_and_prediction_interval = 0;
	param1.quiet = true;

	param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 19;
	param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 14;
	param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 25;

	Dynamic_Param param2(param1);

	// Swarm options
	const std::string outputfilename = "C:/Temp/VS/htm_ga_results_l2.csv";
	swarm::Swarm_Options options;
	options.population_size = 100;
	options.n_epochs = 100;
	options.mutation_rate = 0.1;
	options.quiet = false;

	std::ofstream file; // out file stream
	options.outputfile.open(outputfilename, std::ios::out | std::ios::trunc);
	if (options.outputfile.good())
		log_INFO("swarm:run_ga: writing results to file ", outputfilename, "\n");
	else
		log_WARNING("swarm:run_ga: could not write to file ", outputfilename, "\n");

	using Network_Config = network::network_2Layer<
		N_VISIBLE_SENSORS_L1,
		N_COLUMNS_L1, N_BITS_CELL_L1, HISTORY_SIZE_L1,
		N_COLUMNS_L2, N_BITS_CELL_L2, HISTORY_SIZE_L2,
		ARCH>;

	using P1 = Network_Config::P_L1;
	using P2 = Network_Config::P_L2;
	auto param = std::array<Dynamic_Param, 2>{param1, param2};
	DataStream<P1> datastream;

	const bool load_from_file = false;
	if (load_from_file)
	{
		const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
		//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
		datastream.load_from_file(input_filename, param1);
	}
	else
	{
		const float sparsity = 0.1f;
		const int n_sequences = 3;
		const int sequence_length = 3;
		datastream.generate_random_NxR(sparsity, n_sequences, sequence_length);
	}
	htm::swarm::run_ga<P1, P2>(datastream, param, options);
}

inline void test_swarm_3layers()
{
	// static properties: properties that need to be known at compile time:

	const int N_SENSORS_DIM1 = 20;
	const int N_SENSORS_DIM2 = 20;
	const int N_VISIBLE_SENSORS_L1 = N_SENSORS_DIM1 * N_SENSORS_DIM2;

	const int N_BLOCKS_L1 = 2 * 8;
	const int N_COLUMNS_L1 = 64 * N_BLOCKS_L1;
	const int N_BITS_CELL_L1 = 4;
	const int HISTORY_SIZE_L1 = 7;

	const int N_BLOCKS_L2 = 1 * 8;
	const int N_COLUMNS_L2 = 64 * N_BLOCKS_L2;
	const int N_BITS_CELL_L2 = 4;
	const int HISTORY_SIZE_L2 = 7;

	const int N_BLOCKS_L3 = 1 * 8;
	const int N_COLUMNS_L3 = 64 * N_BLOCKS_L3;
	const int N_BITS_CELL_L3 = 4;
	const int HISTORY_SIZE_L3 = 7;

	const arch_t ARCH = arch_t::RUNTIME;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param1;
	param1.learn = true;
	param1.n_time_steps = 500;
	param1.n_times = 10;
	param1.n_visible_sensors_dim1 = N_SENSORS_DIM1;
	param1.n_visible_sensors_dim2 = N_SENSORS_DIM2;

	param1.show_input_and_prediction_interval = 0;
	param1.quiet = true;

	param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 19;
	param1.TP_MIN_DD_ACTIVATION_THRESHOLD = 14;
	param1.TP_DD_MAX_NEW_SYNAPSE_COUNT = 25;

	Dynamic_Param param2(param1);
	Dynamic_Param param3(param1);

	// Swarm options
	const std::string outputfilename = "C:/Temp/VS/htm_ga_results_l3.csv";
	swarm::Swarm_Options options;
	options.population_size = 100;
	options.n_epochs = 100;
	options.mutation_rate = 0.1;
	options.quiet = false;

	std::ofstream file; // out file stream
	options.outputfile.open(outputfilename, std::ios::out | std::ios::trunc);
	if (options.outputfile.good())
		log_INFO("swarm:run_ga: writing results to file ", outputfilename, "\n");
	else
		log_WARNING("swarm:run_ga: could not write to file ", outputfilename, "\n");

	using Network_Config = network::network_3Layer<
		N_VISIBLE_SENSORS_L1,
		N_COLUMNS_L1, N_BITS_CELL_L1, HISTORY_SIZE_L1,
		N_COLUMNS_L2, N_BITS_CELL_L2, HISTORY_SIZE_L2,
		N_COLUMNS_L3, N_BITS_CELL_L3, HISTORY_SIZE_L3,
		ARCH>;

	using P1 = Network_Config::P_L1;
	using P2 = Network_Config::P_L2;
	using P3 = Network_Config::P_L3;
	auto param = std::array<Dynamic_Param, 3>{param1, param2, param3};
	DataStream<P1> datastream;

	const bool load_from_file = false;
	if (load_from_file)
	{
		const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
		//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
		datastream.load_from_file(input_filename, param1);
	}
	else
	{
		const float sparsity = 0.1f;
		const int n_sequences = 3;
		const int sequence_length = 3;
		datastream.generate_random_NxR(sparsity, n_sequences, sequence_length);
	}
	htm::swarm::run_ga<P1, P2, P3>(datastream, param, options);
}

int main()
{
	const auto start_time = std::chrono::system_clock::now();
	if (true) test_1layer();
	if (false) test_2layers();
	if (false) test_3layers();
	if (false) test_swarm_1layer();
	if (false) test_swarm_2layers();
	if (false) test_swarm_3layers();

	const auto end_time = std::chrono::system_clock::now();

	printf("Elapsed Time: %s\n", ::tools::timing::elapsed_time_str(start_time, end_time).c_str());
	printf("\n-------------------\n");
	printf("\nPress RETURN to finish:");
	getchar();
	return 0;
}