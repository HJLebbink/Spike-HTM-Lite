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
#include <iostream>
#include <fstream>
#include <random>
#include <mutex>
#include <iomanip>        // std::setw
#include <atomic>         // std::atomic_flag
#include <thread>         // std::thread
#include <algorithm>

#include "..\Spike-Tools-Lib\log.ipp"

#include "parameters.ipp"
#include "types.ipp"
#include "datastream.ipp"
#include "layer.ipp"

namespace htm
{
	namespace swarm
	{
		using namespace htm::types;
		using namespace htm::datastream;

		struct Swarm_Options
		{
			mutable std::mutex write_file;

			bool quiet;
			int population_size;
			int n_epochs;
			float mutation_rate;
			mutable std::ofstream outputfile;
			
			void writeline_outputfile(const std::string& str) const
			{
				if (outputfile.good())
				{
					write_file.lock();
					outputfile << str << std::endl;
					outputfile.flush();
					write_file.unlock();
				}
			}
		};

		template <int N_LAYERS>
		struct Configuration
		{
			bool flag;
			float mismatch;
			std::array<Dynamic_Param, N_LAYERS> param;

			static constexpr float MISMATCH_INVALID = 4533234.0f; // just a high random number

			Configuration(const std::array<Dynamic_Param, N_LAYERS>& param) :
				param(param)
			{
				this->mismatch = MISMATCH_INVALID;
				this->flag = false;
			}
			static std::string header_str()
			{
				std::ostringstream result;
				result << "mismatch";
				for (int i = 0; i < N_LAYERS; ++i)
				{
					result << "\t" << Dynamic_Param::header_str();
				}
				return result.str();
			}
			std::string str() const
			{
				std::ostringstream result;
				result << std::fixed << std::setw(7) << std::setprecision(3) << this->mismatch;
				for (int i = 0; i < N_LAYERS; ++i)
				{
					result << "\t" << param[i].str();
				}
				return result.str();
			}
			const Dynamic_Param& get_param(const int i) const
			{
				return this->param[i];
			}
		};

		namespace priv
		{
			namespace ga
			{
				template <int SIZE>
				struct Configuration_better
				{
					constexpr bool operator()(const Configuration<SIZE>& left, const Configuration<SIZE>& right) const
					{
						return (left.mismatch < right.mismatch);
					}
				};

				template <int SIZE> 
				std::vector<Configuration<SIZE>> get_best(const int selection_size, std::vector<Configuration<SIZE>>& pool)
				{
					if (pool.size() < selection_size) log_ERROR("swarm::get_best: pool is too small.\n");

					std::vector<Configuration<SIZE>> results;
					if (pool.size() == selection_size)
					{
						for (auto& p : pool) results.push_back(p);
					}
					else
					{
						std::nth_element(pool.begin(), pool.begin() + selection_size, pool.end(), Configuration_better<SIZE>());
						const float nth_element_mismatch = pool[selection_size].mismatch;
						if (false) log_INFO("swarm:get_best: nth_element_mismatch=", nth_element_mismatch,"\n");

						for (int i = 0; i < selection_size; ++i)
						{
							if (pool[i].mismatch <= nth_element_mismatch)
							{
								if (false) log_INFO("swarm:get_best:", i, ": ", pool[i].mismatch, "\n");
								results.push_back(pool[i]);
							}
						}
					}
					return results;
				}

				template <class T>
				T min_(T a, T b) { return a < b ? a : b; }

				template <class T>
				T max_(T a, T b) { return a > b ? a : b; }

				void mutate_float(
					float& value, const float min, const float max, const float mutation_rate,
					std::default_random_engine& engine)
				{
					std::uniform_real_distribution<float> mutation_dist(0.0f, 1.0f);
					if (mutation_dist(engine) < mutation_rate)
					{
						std::uniform_real_distribution<float> change_dist(-0.1f, 0.1f);
						const float change = change_dist(engine);
						const float new_value = min_(max, max_(min, value + (value * change)));
						if (false) log_INFO("swarm::mutate_float: change = ", change, "; old value = ", value, "; new_value = ", new_value, "\n");
						value = new_value;
					}
				}
				
				void mutate_int(
					int& value, const int min, const int max, const float mutation_rate,
					std::default_random_engine& engine)
				{
					std::uniform_real_distribution<float> mutation_dist(0.0f, 1.0f);
					if (mutation_dist(engine) < mutation_rate)
					{
						std::uniform_int_distribution<int> change_dist(-2, 2);
						const int change = change_dist(engine);
						const int new_value = min_(max, max_(min, value + change));
						if (false) log_INFO("swarm::mutate_int: change = ", change, "; old value = ", value, "; new_value = ", new_value, "\n");
						value = new_value;
					}
				}
				
				void mutate_byte(
					signed char& value, const int min, const int max, const float mutation_rate,
					std::default_random_engine& engine)
				{
					std::uniform_real_distribution<float> mutation_dist(0.0f, 1.0f);
					if (mutation_dist(engine) < mutation_rate)
					{
						std::uniform_int_distribution<int> change_dist(-2, 2);
						const int change = change_dist(engine);
						const signed char new_value = min_(max, max_(min, value + change));
						if (false) log_INFO("swarm::mutate_int: change = ", change, "; old value = ", value, "; new_value = ", new_value, "\n");
						value = new_value;
					}
				}

				Dynamic_Param cross_over(
					const Dynamic_Param& param1,
					const Dynamic_Param& param2,
					const float mutation_rate,
					std::default_random_engine& engine)
				{
					Dynamic_Param result(param1);
					std::uniform_int_distribution<int> bool_dist(0, 1);

					result.SP_PD_PERMANENCE_INIT = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INIT : param2.SP_PD_PERMANENCE_INIT;
					result.SP_PD_PERMANENCE_INC = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INC : param2.SP_PD_PERMANENCE_INC;
					result.SP_PD_PERMANENCE_DEC = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_DEC : param2.SP_PD_PERMANENCE_DEC;
					result.SP_PD_PERMANENCE_INC_WEAK = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INC_WEAK : param2.SP_PD_PERMANENCE_INC_WEAK;
					result.SP_LOCAL_AREA_DENSITY = (bool_dist(engine) == 1) ? param1.SP_LOCAL_AREA_DENSITY : param2.SP_LOCAL_AREA_DENSITY;

					result.TP_DD_SEGMENT_ACTIVE_THRESHOLD = (bool_dist(engine) == 1) ? param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD : param2.TP_DD_SEGMENT_ACTIVE_THRESHOLD;
					result.TP_MIN_DD_ACTIVATION_THRESHOLD = (bool_dist(engine) == 1) ? param1.TP_MIN_DD_ACTIVATION_THRESHOLD : param2.TP_MIN_DD_ACTIVATION_THRESHOLD;
					result.TP_DD_MAX_NEW_SYNAPSE_COUNT = (bool_dist(engine) == 1) ? param1.TP_DD_MAX_NEW_SYNAPSE_COUNT : param2.TP_DD_MAX_NEW_SYNAPSE_COUNT;

					result.TP_DD_PERMANENCE_INC = (bool_dist(engine) == 1) ? param1.TP_DD_PERMANENCE_INC : param2.TP_DD_PERMANENCE_INC;
					result.TP_DD_PERMANENCE_DEC = (bool_dist(engine) == 1) ? param1.TP_DD_PERMANENCE_DEC : param2.TP_DD_PERMANENCE_DEC;
					result.TP_DD_PREDICTED_SEGMENT_DEC = (bool_dist(engine) == 1) ? param1.TP_DD_PREDICTED_SEGMENT_DEC : param2.TP_DD_PREDICTED_SEGMENT_DEC;

					mutate_byte(result.SP_PD_PERMANENCE_INIT, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_INC, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_DEC, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_INC_WEAK, 0, 128, mutation_rate, engine);
					mutate_float(result.SP_LOCAL_AREA_DENSITY, 0, 0.1, mutation_rate, engine);

					mutate_int(result.TP_DD_SEGMENT_ACTIVE_THRESHOLD, 1, 100, mutation_rate, engine);
					mutate_int(result.TP_MIN_DD_ACTIVATION_THRESHOLD, 1, 100, mutation_rate, engine);
					mutate_int(result.TP_DD_MAX_NEW_SYNAPSE_COUNT, 1, 100, mutation_rate, engine);

					mutate_byte(result.TP_DD_PERMANENCE_INC, 0, 127, mutation_rate, engine);
					mutate_byte(result.TP_DD_PERMANENCE_DEC, 0, 127, mutation_rate, engine);
					mutate_byte(result.TP_DD_PREDICTED_SEGMENT_DEC, 0, 127, mutation_rate, engine);

					if (false)
					{
						log_INFO("swarm:cross_over: parent1=", param1.str() + "\n");
						log_INFO("swarm:cross_over: parent2=", param2.str() + "\n");
						log_INFO("swarm:cross_over: child=", result.str() + "\n");
					}
					return result;
				}
			
				template <int SIZE>
				std::array<Dynamic_Param, SIZE> cross_over(
					const std::array<Dynamic_Param, SIZE>& param1,
					const std::array<Dynamic_Param, SIZE>& param2,
					const float mutation_rate,
					std::default_random_engine& random_engine)
				{
					std::array<Dynamic_Param, SIZE> results;
					for (int i = 0; i < SIZE; ++i)
					{
						results[i] = priv::ga::cross_over(param1[i], param2[i], mutation_rate, random_engine);
					}
					return results;
				}
			}

			void random_param(
				Dynamic_Param& param,
				std::default_random_engine& engine)
			{
				std::uniform_real_distribution<float> SP_LOCAL_AREA_DENSITY_dist(0.0f, 0.1f);

				std::uniform_int_distribution<int> SP_PD_PERMANENCE_THRESHOLD_dist(0, 30);
				std::uniform_int_distribution<int> SP_PD_PERMANENCE_INIT_dist(0, 30);
				std::uniform_int_distribution<int> SP_PD_PERMANENCE_INC_dist(0, 30);
				std::uniform_int_distribution<int> SP_PD_PERMANENCE_DEC_dist(0, 30);
				std::uniform_int_distribution<int> SP_PD_PERMANENCE_INC_WEAK_dist(0, 30);

				std::uniform_int_distribution<int> TP_DD_ACTIVATION_THRESHOLD_dist(1, 30);
				std::uniform_int_distribution<int> MIN_DD_ACTIVATION_THRESHOLD_dist(1, 30);
				std::uniform_int_distribution<int> TP_DD_MAX_NEW_SYNAPSE_COUNT_dist(1, 40);

				std::uniform_int_distribution<int> TP_DD_PERMANENCE_THRESHOLD_dist(1, 30);
				std::uniform_int_distribution<int> TP_DD_PERMANENCE_INC_dist(1, 30);
				std::uniform_int_distribution<int> TP_DD_PERMANENCE_DEC_dist(1, 30);
				std::uniform_int_distribution<int> TP_DD_PREDICTED_SEGMENT_DEC_dist(1, 30);

				param.SP_LOCAL_AREA_DENSITY = SP_LOCAL_AREA_DENSITY_dist(engine);

				param.SP_PD_PERMANENCE_INIT = static_cast<signed char>(SP_PD_PERMANENCE_INIT_dist(engine));
				param.SP_PD_PERMANENCE_INC = static_cast<signed char>(SP_PD_PERMANENCE_INC_dist(engine));
				param.SP_PD_PERMANENCE_DEC = static_cast<signed char>(SP_PD_PERMANENCE_DEC_dist(engine));
				param.SP_PD_PERMANENCE_INC_WEAK = static_cast<signed char>(SP_PD_PERMANENCE_INC_WEAK_dist(engine));

				param.TP_DD_SEGMENT_ACTIVE_THRESHOLD = TP_DD_ACTIVATION_THRESHOLD_dist(engine);
				param.TP_MIN_DD_ACTIVATION_THRESHOLD = MIN_DD_ACTIVATION_THRESHOLD_dist(engine);
				param.TP_DD_MAX_NEW_SYNAPSE_COUNT = TP_DD_MAX_NEW_SYNAPSE_COUNT_dist(engine);

				param.TP_DD_PERMANENCE_INC = static_cast<signed char>(TP_DD_PERMANENCE_INC_dist(engine));
				param.TP_DD_PERMANENCE_DEC = static_cast<signed char>(TP_DD_PERMANENCE_DEC_dist(engine));
				param.TP_DD_PREDICTED_SEGMENT_DEC = static_cast<signed char>(TP_DD_PREDICTED_SEGMENT_DEC_dist(engine));
			}

			template <int SIZE>
			void random_param(
				std::array<Dynamic_Param, SIZE>& param,
				std::default_random_engine& engine)
			{
				for (int i = 0; i < SIZE; ++i) random_param(param[i], engine);
			}

			template <int SIZE>
			Configuration<SIZE> random_config(
				const std::array<Dynamic_Param, SIZE>& param,
				std::default_random_engine& random_engine)
			{
				std::array<Dynamic_Param, SIZE> rand_param;
				for (int i = 0; i < SIZE; ++i)
				{
					Dynamic_Param new_param(param[i]);
					priv::random_param(new_param, random_engine);
					rand_param[i] = new_param;
				}
				return Configuration<SIZE>(rand_param);
			}

			template <typename P>
			void do_work(
				const DataStream<P>& datastream,
				Layer<P>& layer,
				Swarm_Options& options,
				Configuration<1>& result)
			{
				if (!result.flag)
				{
					result.flag = true; // using std::atomic_flag would be better...

					htm::layer::init(layer, result.param[0]);
					const int total_mismatch = htm::layer::run(datastream, layer, result.param[0]);

					if (result.mismatch != Configuration<1>::MISMATCH_INVALID)
					{
						log_WARNING("swarm:do_work: concurrency problem: but not serious, you just did some duplicate work.\n");
					}
					result.mismatch = static_cast<float>(total_mismatch) / (result.param[0].n_time_steps * result.param[0].n_times);
					if (!options.quiet) log_INFO(result.str(), "\n");
					options.writeline_outputfile(result.str());
				}
			}

			template <typename P1, typename P2>
			void do_work(
				const DataStream<P1>& datastream,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				const Swarm_Options& options,
				Configuration<2>& result)
			{
				if (!result.flag)
				{
					result.flag = true; // using std::atomic_flag would be better...

					htm::layer::init(layer1, result.get_param(0));
					htm::layer::init(layer2, result.get_param(1));
					const int total_mismatch = htm::network::run(datastream, result.param, layer1, layer2);

					if (result.mismatch != Configuration<2>::MISMATCH_INVALID)
					{
						log_WARNING("swarm:do_work: concurrency problem: but not serious, you just did some duplicate work.\n");
					}
					result.mismatch = static_cast<float>(total_mismatch) / (result.param[0].n_time_steps * result.param[0].n_times);
					if (!options.quiet) log_INFO(result.str(), "\n");
					options.writeline_outputfile(result.str());
				}
			}

			template <typename P1, typename P2, typename P3>
			void do_work(
				const DataStream<P1>& datastream,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				Layer<P3>& layer3,
				const Swarm_Options& options,
				Configuration<3>& result)
			{
				if (!result.flag)
				{
					result.flag = true; // using std::atomic_flag would be better...

					htm::layer::init(layer1, result.get_param(0));
					htm::layer::init(layer2, result.get_param(1));
					htm::layer::init(layer3, result.get_param(2));
					const int total_mismatch = htm::network::run(datastream, result.param, layer1, layer2, layer3);

					if (result.mismatch != Configuration<1>::MISMATCH_INVALID)
					{
						log_WARNING("swarm:do_work: concurrency problem: but not serious, you just did some duplicate work.\n");
					}
					result.mismatch = static_cast<float>(total_mismatch) / (result.param[0].n_time_steps * result.param[0].n_times);
					if (!options.quiet) log_INFO(result.str(), "\n");
					options.writeline_outputfile(result.str());
				}
			}

			template <typename P>
			void do_work_all(
				const DataStream<P>& datastream,
				const int start_pos,
				Swarm_Options& options,
				std::vector<Configuration<1>>& results)
			{
				Layer<P> layer;
				const int nResults = static_cast<int>(results.size());

				for (int i = start_pos; i < nResults; ++i) do_work(datastream, layer, options, results[i]);
				for (int i = 0; i < nResults; ++i) do_work(datastream, layer, options, results[i]);
			}

			template <typename P1, typename P2>
			void do_work_all(
				const DataStream<P1>& datastream,
				const int start_pos,
				Swarm_Options& options,
				std::vector<Configuration<2>>& results)
			{
				Layer<P1> layer1;
				Layer<P2> layer2;
				const int nResults = static_cast<int>(results.size());

				for (int i = start_pos; i < nResults; ++i) do_work(datastream, layer1, layer2, options, results[i]);
				for (int i = 0; i < nResults; ++i) do_work(datastream, layer1, layer2, options, results[i]);
			}

			template <typename P1, typename P2, typename P3>
			void do_work_all(
				const DataStream<P1>& datastream,
				const int start_pos,
				Swarm_Options& options,
				std::vector<Configuration<3>>& results)
			{
				Layer<P1> layer1;
				Layer<P2> layer2;
				Layer<P3> layer3;
				const int nResults = static_cast<int>(results.size());

				for (int i = start_pos; i < nResults; ++i) do_work(datastream, layer1, layer2, layer3, options, results[i]);
				for (int i = 0; i < nResults; ++i) do_work(datastream, layer1, layer2, layer3, options, results[i]);
			}
		}

		template <typename P>
		std::vector<Configuration<1>> run_random(
			const DataStream<P>& datastream,
			const std::array<Dynamic_Param, 1>& param,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (!options.quiet) log_INFO("swarm:run_random: running random swarm with ", nThreads, " threads.");
			if (!options.quiet) std::cout << Configuration<1>::header_str() << std::endl;

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::vector<Configuration<1>> results;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				results.push_back(priv::random_config(param, random_engine));
			}

			std::vector<std::thread> workers;
			for (int thread = 0; thread < nThreads; ++thread)
			{
				workers.push_back(std::thread(priv::do_work_all<P>, datastream, thread * 4, std::ref(options), std::ref(results)));
			}
			for (auto& t : workers) if (t.joinable()) t.join();
			return results;
		}

		template <typename P>
		std::vector<Configuration<1>> run_ga(
			const DataStream<P>& datastream,
			const std::array<Dynamic_Param, 1>& param,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (true) log_INFO("swarm:run_ga: running ga swarm with ", nThreads, " threads.\n");
			if (true) log_INFO(Configuration<1>::header_str(), "\n");

			options.writeline_outputfile(Configuration<1>::header_str());

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::uniform_int_distribution<int> individual_dist(0, options.population_size-1);

			std::vector<Configuration<1>> pool = run_random<P>(datastream, param, options);
			std::vector<Configuration<1>> results;
			std::vector<std::thread> workers;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				results.push_back(priv::random_config(param, random_engine));
			}

			for (int epoch = 0; epoch < options.n_epochs; ++epoch)
			{
				if (!options.quiet) log_INFO("swarm:run_ga: epoch ", epoch, "\n");

				results.clear();
				const auto best_indiduals = priv::ga::get_best(options.population_size, pool);
				for (int i = 0; i < options.population_size; ++i)
				{
					const auto parent1 = best_indiduals[individual_dist(random_engine)].param;
					const auto parent2 = best_indiduals[individual_dist(random_engine)].param;
					const auto child = priv::ga::cross_over(parent1, parent2, options.mutation_rate, random_engine);
					results.push_back(Configuration<1>(child));
				}

				workers.clear();
				for (int thread = 0; thread < nThreads; ++thread)
				{
					workers.push_back(std::thread(priv::do_work_all<P>, datastream, thread * 4, std::ref(options), std::ref(results)));
				}
				for (auto& t : workers) if (t.joinable()) t.join();
				for (const auto& config : results) pool.push_back(config);
			}
			return pool;
		}
		
		template <typename P1, typename P2>
		std::vector<Configuration<2>> run_ga(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 2>& param,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (true) log_INFO("swarm:run_ga: running ga swarm with ", nThreads, " threads.\n");
			if (true) log_INFO(Configuration<2>::header_str(), "\n");

			options.writeline_outputfile(Configuration<2>::header_str());

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::uniform_int_distribution<int> individual_dist(0, options.population_size - 1);

			std::vector<Configuration<2>> pool;// = run_random<P1>(input_filename, param, options);
			std::vector<Configuration<2>> results;
			std::vector<std::thread> workers;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				pool.push_back(priv::random_config(param, random_engine));
			}

			for (int epoch = 0; epoch < options.n_epochs; ++epoch)
			{
				if (!options.quiet) log_INFO("swarm:run_ga: epoch ", epoch, "\n");

				results.clear();
				const auto best_indiduals = priv::ga::get_best(options.population_size, pool);
				for (int i = 0; i < options.population_size; ++i)
				{
					const auto parent1 = best_indiduals[individual_dist(random_engine)].param;
					const auto parent2 = best_indiduals[individual_dist(random_engine)].param;
					const auto child = priv::ga::cross_over(parent1, parent2, options.mutation_rate, random_engine);
					results.push_back(Configuration<2>(child));
				}

				workers.clear();
				for (int thread = 0; thread < nThreads; ++thread)
				{
					workers.push_back(std::thread(priv::do_work_all<P1, P2>, datastream, thread * 4, std::ref(options), std::ref(results)));
				}
				for (auto& t : workers) if (t.joinable()) t.join();
				for (const auto& config : results) pool.push_back(config);
			}
			return pool;
		}
		
		template <typename P1, typename P2, typename P3>
		std::vector<Configuration<3>> run_ga(
			const DataStream<P1>& datastream,
			const std::array<Dynamic_Param, 3>& param, 
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (true) log_INFO("swarm:run_ga: running ga swarm with ", nThreads, " threads.\n");
			if (true) log_INFO(Configuration<3>::header_str(), "\n");

			options.writeline_outputfile(Configuration<3>::header_str());

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::uniform_int_distribution<int> individual_dist(0, options.population_size - 1);

			std::vector<Configuration<3>> pool;// = run_random<P1>(input_filename, param, options);
			std::vector<Configuration<3>> results;
			std::vector<std::thread> workers;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				pool.push_back(priv::random_config(param, random_engine));
			}

			for (int epoch = 0; epoch < options.n_epochs; ++epoch)
			{
				if (!options.quiet) log_INFO("swarm:run_ga: epoch ", epoch, "\n");

				results.clear();
				const auto best_indiduals = priv::ga::get_best(options.population_size, pool);
				for (int i = 0; i < options.population_size; ++i)
				{
					const auto parent1 = best_indiduals[individual_dist(random_engine)].param;
					const auto parent2 = best_indiduals[individual_dist(random_engine)].param;
					const auto child = priv::ga::cross_over(parent1, parent2, options.mutation_rate, random_engine);
					results.push_back(Configuration<3>(child));
				}

				workers.clear();
				for (int thread = 0; thread < nThreads; ++thread)
				{
					workers.push_back(std::thread(priv::do_work_all<P1, P2, P3>, datastream, thread * 4, std::ref(options), std::ref(results)));
				}
				for (auto& t : workers) if (t.joinable()) t.join();
				for (const auto& config : results) pool.push_back(config);
			}
			return pool;
		}
	}
}