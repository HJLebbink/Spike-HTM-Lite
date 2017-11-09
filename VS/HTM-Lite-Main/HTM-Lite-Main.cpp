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

#include <chrono>
#include <stdio.h>  // for printf

#include "../Spike-Tools-Lib/log.ipp"
#include "../Spike-Tools-Lib/timing.ipp"
#include "../Spike-Tools-Lib/profiler.ipp"

#include "../HTM-Lite-LIB/constants.ipp"
#include "../HTM-Lite-LIB/types.ipp"
#include "../HTM-Lite-LIB/layer.ipp"
#include "../HTM-Lite-LIB/swarm.ipp"

using namespace htm;
using namespace htm::types;

inline void test_run()
{
	// static properties: properties that need to be known at compile time:

	//const int N_BLOCKS = 8; // 512 columns: use sparsity 0.05 -> 25
	//const int N_BLOCKS = 4096; // 262144 columns: use sparsity of 0.005 -> 1310
	//const int N_BLOCKS = 16384; // 1048576 columns: use sparsity of 0.002 -> 2048
	const int N_BLOCKS = 100 * 8;
	const int N_COLUMNS_LOCAL = 64 * N_BLOCKS;
	const int N_BITS_CELL_LOCAL = 4;
	const int SENSOR_DIM1_LOCAL = 20;
	const int SENSOR_DIM2_LOCAL = 20;
	const arch_t ARCH = arch_t::RUNTIME;
	const int HISTORY_SIZE = 8;

	// dynamic properties: properties that can be changed while the program is running.
	Dynamic_Param param;
	param.n_time_steps = 2000;
	param.n_times = 1;
	param.progress = false;
	param.TP_DD_SEGMENT_ACTIVE_THRESHOLD = 19;
	param.MIN_DD_ACTIVATION_THRESHOLD = 14;
	param.TP_DD_MAX_NEW_SYNAPSE_COUNT = 25;


	const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
	//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";

	using P = Static_Param<N_COLUMNS_LOCAL, N_BITS_CELL_LOCAL, SENSOR_DIM1_LOCAL, SENSOR_DIM2_LOCAL, ARCH, HISTORY_SIZE>;
	Layer<P> layer;
	auto data = encoder::encode_pass_through<P>(input_filename);

	const bool LEARN = true;
	htm::layer::run<LEARN>(data, layer, param);
}

inline void test_swarm()
{
	// static properties: properties that need to be known at compile time:

	//const int N_BLOCKS = 9; // 1024 columns
	//const int N_BLOCKS = 4096; // 262144 columns
	//const int N_BLOCKS = 16384; // 1048576 columns
	const int N_BLOCKS = 10 * 8;
	const int N_COLUMNS_LOCAL = 64 * N_BLOCKS;
	const int N_BITS_CELL_LOCAL = 4;
	const int SENSOR_DIM1_LOCAL = 20;
	const int SENSOR_DIM2_LOCAL = 20;
	const arch_t ARCH = arch_t::RUNTIME;
	const int HISTORY_SIZE = 8;

	Dynamic_Param param;
	param.n_time_steps = 1000;
	param.n_times = 1;

	const std::string input_filename = "../../Misc/data/ABBCBBA_20x20/input.txt";
	//const std::string input_filename = "../../Misc/data/AAAX_16x16/input.txt";
	//htm::swarm::run_random(input_filename, param);

	// Swarm options
	const std::string outputfilename = "C:/Temp/VS/htm_ga_results.csv";
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

	using P = Static_Param<N_COLUMNS_LOCAL, N_BITS_CELL_LOCAL, SENSOR_DIM1_LOCAL, SENSOR_DIM2_LOCAL, ARCH, HISTORY_SIZE>;
	htm::swarm::run_ga<P>(input_filename, param, options);
}

int main()
{
	const auto start_time = std::chrono::system_clock::now();
	test_run();
	//test_swarm();
	const auto end_time = std::chrono::system_clock::now();

	printf("Elapsed Time: %s\n", ::tools::timing::elapsed_time_str(start_time, end_time).c_str());
	printf("\n-------------------\n");
	printf("\nPress RETURN to finish:");
	getchar();
	return 0;
}