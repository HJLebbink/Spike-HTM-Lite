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

#include "constants.ipp"
#include "types.ipp"
#include "layer.ipp"

namespace htm
{
	namespace swarm
	{
		using namespace htm::types;

		struct Swarm_Options
		{
			std::mutex write_file;

			bool quiet;
			int population_size;
			int n_epochs;
			float mutation_rate;
			std::ofstream outputfile;
			
			void writeline_outputfile(const std::string& str)
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

		struct Configuration
		{
			bool flag;
			float mismatch;
			Dynamic_Param param;

			static constexpr float MISMATCH_INVALID = 4533234.0f; // just a high random number

			Configuration(const Dynamic_Param& param):
				param(param)
			{
				this->mismatch = MISMATCH_INVALID;
				this->flag = false;
			}
			static std::string header_str()
			{
				std::ostringstream result;
				result << "mismatch\t" << Dynamic_Param::header_str();
				return result.str();
			}
			std::string str() const
			{
				std::ostringstream result;
				//result << std::setw(3) << std::min(6, static_cast<int>(this->mismatch)) << "\t" << param.str();
				result << std::fixed << std::setw(7) << std::setprecision(3) << this->mismatch << "\t" << param.str();
				return result.str();
			}
		};

		struct Configuration2
		{
			bool flag;
			float mismatch;
			Dynamic_Param param1;
			Dynamic_Param param2;

			static constexpr float MISMATCH_INVALID = 4533234.0f; // just a high random number

			Configuration2(const Dynamic_Param& param1, const Dynamic_Param& param2) :
				param1(param1),
				param2(param2)
			{
				this->mismatch = MISMATCH_INVALID;
				this->flag = false;
			}
			static std::string header_str()
			{
				std::ostringstream result;
				result << "mismatch\t" << Dynamic_Param::header_str() << "\t" << Dynamic_Param::header_str();
				return result.str();
			}
			std::string str() const
			{
				std::ostringstream result;
				//result << std::setw(3) << std::min(6, static_cast<int>(this->mismatch)) << "\t" << param.str();
				result << std::fixed << std::setw(7) << std::setprecision(3) << this->mismatch << "\t" << param1.str() << "\t" << param2.str();
				return result.str();
			}
		};

		namespace priv
		{
			namespace ga
			{
				template <typename T>
				struct Configuration_better
				{
					constexpr bool operator()(const T& left, const T& right) const
					{
						return (left.mismatch < right.mismatch);
					}
				};

				template <typename T> 
				std::vector<T> get_best(const int selection_size, std::vector<T>& pool)
				{
					if (pool.size() < selection_size) log_ERROR("swarm::get_best: pool is too small.\n");

					std::vector<T> results;
					if (pool.size() == selection_size)
					{
						for (auto& p : pool) results.push_back(p);
					}
					else
					{
						std::nth_element(pool.begin(), pool.begin() + selection_size, pool.end(), Configuration_better<T>());
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

					result.SP_PD_CONNECTED_THRESHOLD = (bool_dist(engine) == 1) ? param1.SP_PD_CONNECTED_THRESHOLD : param2.SP_PD_CONNECTED_THRESHOLD;
					result.SP_PD_PERMANENCE_INIT = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INIT : param2.SP_PD_PERMANENCE_INIT;
					result.SP_PD_PERMANENCE_INC = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INC : param2.SP_PD_PERMANENCE_INC;
					result.SP_PD_PERMANENCE_DEC = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_DEC : param2.SP_PD_PERMANENCE_DEC;
					result.SP_PD_PERMANENCE_INC_WEAK = (bool_dist(engine) == 1) ? param1.SP_PD_PERMANENCE_INC_WEAK : param2.SP_PD_PERMANENCE_INC_WEAK;
					result.SP_LOCAL_AREA_DENSITY = (bool_dist(engine) == 1) ? param1.SP_LOCAL_AREA_DENSITY : param2.SP_LOCAL_AREA_DENSITY;

					result.TP_DD_SEGMENT_ACTIVE_THRESHOLD = (bool_dist(engine) == 1) ? param1.TP_DD_SEGMENT_ACTIVE_THRESHOLD : param2.TP_DD_SEGMENT_ACTIVE_THRESHOLD;
					result.MIN_DD_ACTIVATION_THRESHOLD = (bool_dist(engine) == 1) ? param1.MIN_DD_ACTIVATION_THRESHOLD : param2.MIN_DD_ACTIVATION_THRESHOLD;
					result.TP_DD_MAX_NEW_SYNAPSE_COUNT = (bool_dist(engine) == 1) ? param1.TP_DD_MAX_NEW_SYNAPSE_COUNT : param2.TP_DD_MAX_NEW_SYNAPSE_COUNT;

					result.TP_DD_ACTIVE_THRESHOLD = (bool_dist(engine) == 1) ? param1.TP_DD_ACTIVE_THRESHOLD : param2.TP_DD_ACTIVE_THRESHOLD;
					result.TP_DD_PERMANENCE_INC = (bool_dist(engine) == 1) ? param1.TP_DD_PERMANENCE_INC : param2.TP_DD_PERMANENCE_INC;
					result.TP_DD_PERMANENCE_DEC = (bool_dist(engine) == 1) ? param1.TP_DD_PERMANENCE_DEC : param2.TP_DD_PERMANENCE_DEC;
					result.TP_DD_PREDICTED_SEGMENT_DEC = (bool_dist(engine) == 1) ? param1.TP_DD_PREDICTED_SEGMENT_DEC : param2.TP_DD_PREDICTED_SEGMENT_DEC;

					mutate_byte(result.SP_PD_CONNECTED_THRESHOLD, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_INIT, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_INC, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_DEC, 0, 128, mutation_rate, engine);
					mutate_byte(result.SP_PD_PERMANENCE_INC_WEAK, 0, 128, mutation_rate, engine);
					mutate_float(result.SP_LOCAL_AREA_DENSITY, 0, 0.1, mutation_rate, engine);

					mutate_int(result.TP_DD_SEGMENT_ACTIVE_THRESHOLD, 1, 100, mutation_rate, engine);
					mutate_int(result.MIN_DD_ACTIVATION_THRESHOLD, 1, 100, mutation_rate, engine);
					mutate_int(result.TP_DD_MAX_NEW_SYNAPSE_COUNT, 1, 100, mutation_rate, engine);

					mutate_byte(result.TP_DD_ACTIVE_THRESHOLD, 0, 127, mutation_rate, engine);
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

				param.SP_PD_CONNECTED_THRESHOLD = static_cast<signed char>(SP_PD_PERMANENCE_THRESHOLD_dist(engine));
				param.SP_PD_PERMANENCE_INIT = static_cast<signed char>(SP_PD_PERMANENCE_INIT_dist(engine));
				param.SP_PD_PERMANENCE_INC = static_cast<signed char>(SP_PD_PERMANENCE_INC_dist(engine));
				param.SP_PD_PERMANENCE_DEC = static_cast<signed char>(SP_PD_PERMANENCE_DEC_dist(engine));
				param.SP_PD_PERMANENCE_INC_WEAK = static_cast<signed char>(SP_PD_PERMANENCE_INC_WEAK_dist(engine));

				param.TP_DD_SEGMENT_ACTIVE_THRESHOLD = TP_DD_ACTIVATION_THRESHOLD_dist(engine);
				param.MIN_DD_ACTIVATION_THRESHOLD = MIN_DD_ACTIVATION_THRESHOLD_dist(engine);
				param.TP_DD_MAX_NEW_SYNAPSE_COUNT = TP_DD_MAX_NEW_SYNAPSE_COUNT_dist(engine);

				param.TP_DD_ACTIVE_THRESHOLD = static_cast<signed char>(TP_DD_PERMANENCE_THRESHOLD_dist(engine));
				param.TP_DD_PERMANENCE_INC = static_cast<signed char>(TP_DD_PERMANENCE_INC_dist(engine));
				param.TP_DD_PERMANENCE_DEC = static_cast<signed char>(TP_DD_PERMANENCE_DEC_dist(engine));
				param.TP_DD_PREDICTED_SEGMENT_DEC = static_cast<signed char>(TP_DD_PREDICTED_SEGMENT_DEC_dist(engine));
			}

			template <typename P>
			void do_work(
				const std::vector<Layer<P>::Active_Visible_Sensors>& data,
				Layer<P>& layer,
				Swarm_Options& options,
				Configuration& result)
			{
				if (!result.flag)
				{
					result.flag = true; // using std::atomic_flag would be better...
					result.param.progress = false;
					result.param.quiet = true;

					htm::layer::init(layer, result.param);
					const int total_mismatch = htm::layer::run(data, layer, result.param);

					if (result.mismatch != Configuration::MISMATCH_INVALID)
					{
						log_WARNING("swarm:do_work: concurrency problem: but not serious, you just did some duplicate work.\n");
					}
					result.mismatch = static_cast<float>(total_mismatch) / (result.param.n_time_steps * result.param.n_times);
					if (!options.quiet) log_INFO(result.str(), "\n");
					options.writeline_outputfile(result.str());
				}
			}

			template <typename P1, typename P2>
			void do_work(
				const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
				Layer<P1>& layer1,
				Layer<P2>& layer2,
				Swarm_Options& options,
				Configuration2& result)
			{
				if (!result.flag)
				{
					result.flag = true; // using std::atomic_flag would be better...
					result.param1.progress = false;
					result.param1.quiet = true;

					htm::layer::init(layer1, result.param1);
					htm::layer::init(layer2, result.param2);
					const int total_mismatch = htm::network::run(data, layer1, layer2, result.param1, result.param2);

					if (result.mismatch != Configuration::MISMATCH_INVALID)
					{
						log_WARNING("swarm:do_work: concurrency problem: but not serious, you just did some duplicate work.\n");
					}
					result.mismatch = static_cast<float>(total_mismatch) / (result.param1.n_time_steps * result.param1.n_times);
					if (!options.quiet) log_INFO(result.str(), "\n");
					options.writeline_outputfile(result.str());
				}
			}

			template <typename P>
			void do_work_all(
				const std::vector<Layer<P>::Active_Visible_Sensors>& data,
				const int start_pos,
				Swarm_Options& options,
				std::vector<Configuration>& results)
			{
				Layer<P> layer;
				const int nResults = static_cast<int>(results.size());

				for (int i = start_pos; i < nResults; ++i) do_work(data, layer, options, results[i]);
				for (int i = 0; i < nResults; ++i) do_work(data, layer, options, results[i]);
			}

			template <typename P1, typename P2>
			void do_work_all(
				const std::vector<Layer<P1>::Active_Visible_Sensors>& data,
				const int start_pos,
				Swarm_Options& options,
				std::vector<Configuration2>& results)
			{
				Layer<P1> layer1;
				Layer<P2> layer2;
				const int nResults = static_cast<int>(results.size());

				for (int i = start_pos; i < nResults; ++i) do_work(data, layer1, layer2, options, results[i]);
				for (int i = 0; i < nResults; ++i) do_work(data, layer1, layer2, options, results[i]);
			}
		}

		template <typename P>
		std::vector<Configuration> run_random(
			const std::string& input_filename,
			const Dynamic_Param& param,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (!options.quiet) log_INFO("swarm:run_random: running random swarm with ", nThreads, " threads.");
			if (!options.quiet) std::cout << Configuration::header_str() << std::endl;

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::vector<Configuration> results;

			for (int i = 0; i < options.population_size; ++i)
			{
				Dynamic_Param param_local(param);
				param_local.progress = false;
				param_local.quiet = true;
				priv::random_param(param_local, random_engine);
				results.push_back(Configuration(param_local));
			}

			const auto data = encoder::encode_pass_through<P>(input_filename, param);

			std::vector<std::thread> workers;
			for (int thread = 0; thread < nThreads; ++thread)
			{
				workers.push_back(std::thread(priv::do_work_all<P>, data, thread * 4, std::ref(options), std::ref(results)));
			}
			for (auto& t : workers) if (t.joinable()) t.join();
			return results;
		}

		template <typename P>
		std::vector<Configuration> run_ga(
			const std::string& input_filename,
			const Dynamic_Param& param,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (true) log_INFO("swarm:run_ga: running ga swarm with ", nThreads, " threads.\n");
			if (true) log_INFO(Configuration::header_str(), "\n");

			options.writeline_outputfile(Configuration::header_str());

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::uniform_int_distribution<int> individual_dist(0, options.population_size-1);

			const auto data = encoder::encode_pass_through<P>(input_filename, param);

			std::vector<Configuration> pool = run_random<P>(input_filename, param, options);
			std::vector<Configuration> results;
			std::vector<std::thread> workers;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				Dynamic_Param new_param(param);
				priv::random_param(new_param, random_engine);
				pool.push_back(Configuration(new_param));
			}

			for (int epoch = 0; epoch < options.n_epochs; ++epoch)
			{
				if (!options.quiet) log_INFO("swarm:run_ga: epoch ", epoch, "\n");

				results.clear();
				const std::vector<Configuration> best_indiduals = priv::ga::get_best(options.population_size, pool);
				for (int i = 0; i < options.population_size; ++i)
				{
					const Dynamic_Param parent1 = best_indiduals[individual_dist(random_engine)].param;
					const Dynamic_Param parent2 = best_indiduals[individual_dist(random_engine)].param;
					const Dynamic_Param child = priv::ga::cross_over(parent1, parent2, options.mutation_rate, random_engine);
					results.push_back(Configuration(child));
				}

				workers.clear();
				for (int thread = 0; thread < nThreads; ++thread)
				{
					workers.push_back(std::thread(priv::do_work_all<P>, data, thread * 4, std::ref(options), std::ref(results)));
				}
				for (auto& t : workers) if (t.joinable()) t.join();

				for (const auto& config : results)
				{
					pool.push_back(config);
				}
			}
			return pool;
		}
		
		template <typename P1, typename P2>
		std::vector<Configuration2> run_ga(
			const std::string& input_filename,
			const Dynamic_Param& param1,
			const Dynamic_Param& param2,
			Swarm_Options& options)
		{
			const int nThreads = static_cast<int>(std::thread::hardware_concurrency());
			if (true) log_INFO("swarm:run_ga: running ga swarm with ", nThreads, " threads.\n");
			if (true) log_INFO(Configuration2::header_str(), "\n");

			options.writeline_outputfile(Configuration2::header_str());

			std::random_device r;
			std::default_random_engine random_engine(r());
			std::uniform_int_distribution<int> individual_dist(0, options.population_size - 1);

			const auto data = encoder::encode_pass_through<P1>(input_filename, param1);

			std::vector<Configuration2> pool;// = run_random<P1>(input_filename, param, options);
			std::vector<Configuration2> results;
			std::vector<std::thread> workers;

			// initialize the pool with random individuals
			for (int i = 0; i < options.population_size; ++i)
			{
				Dynamic_Param new_param1(param1);
				Dynamic_Param new_param2(param2);
				priv::random_param(new_param1, random_engine);
				priv::random_param(new_param2, random_engine);
				pool.push_back(Configuration2(new_param1, new_param2));
			}

			for (int epoch = 0; epoch < options.n_epochs; ++epoch)
			{
				if (!options.quiet) log_INFO("swarm:run_ga: epoch ", epoch, "\n");

				results.clear();
				const std::vector<Configuration2> best_indiduals = priv::ga::get_best(options.population_size, pool);
				for (int i = 0; i < options.population_size; ++i)
				{
					const int individual_id1 = individual_dist(random_engine);
					const int individual_id2 = individual_dist(random_engine);

					const Dynamic_Param parent1_A= best_indiduals[individual_id1].param1;
					const Dynamic_Param parent2_A = best_indiduals[individual_id2].param1;

					const Dynamic_Param parent1_B = best_indiduals[individual_id1].param2;
					const Dynamic_Param parent2_B = best_indiduals[individual_id2].param2;

					const Dynamic_Param child_A = priv::ga::cross_over(parent1_A, parent2_A, options.mutation_rate, random_engine);
					const Dynamic_Param child_B = priv::ga::cross_over(parent1_B, parent2_B, options.mutation_rate, random_engine);

					results.push_back(Configuration2(child_A, child_B));
				}

				workers.clear();
				for (int thread = 0; thread < nThreads; ++thread)
				{
					workers.push_back(std::thread(priv::do_work_all<P1, P2>, data, thread * 4, std::ref(options), std::ref(results)));
				}
				for (auto& t : workers) if (t.joinable()) t.join();

				for (const auto& config : results)
				{
					pool.push_back(config);
				}
			}
			return pool;
		}
	}
}