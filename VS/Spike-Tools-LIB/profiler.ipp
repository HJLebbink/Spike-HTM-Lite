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
#include <type_traits>
#include <iostream>		// for cerr and cout
#include <array>

namespace tools
{
	namespace profiler
	{
		namespace priv
		{
			inline std::string elapsed_cycles_str(const unsigned long long cycles)
			{
				std::ostringstream os;
				os << cycles << " cycles = "
					<< cycles / 1000 << " kcycles = "
					<< cycles / 1000000 << " mcycles";
				return os.str();
			}
		}

		//WARNING: is NOT thread safe
		template <int SIZE_IN, bool ON_IN>
		struct Profiler
		{
			static const int SIZE = SIZE_IN;
			static const bool ON = ON_IN;

			static std::array<unsigned long long, ON_IN> startTime_;
			static std::array<unsigned long long, ON_IN> totalTime_;

			static inline void reset()
			{}

			template <int I>
			static inline void start()
			{
				static_assert(I < SIZE, "provided I is invalid");
			}

			// stop stopwatch for profiler identified with parameter i, and add the time since starting it to the total time.
			template <int I>
			static inline void stop()
			{
				static_assert(I < SIZE, "provided I is invalid");
			}

			static inline void print_elapsed_time()
			{};
		};

		template <int SIZE_In>
		struct Profiler <SIZE_In, true>
		{
			static const int SIZE = SIZE_In;
			static const bool ON = true;

			using ThisClass = Profiler <SIZE, true>;

			static std::array<unsigned long long, SIZE_In> startTime_;
			static std::array<unsigned long long, SIZE_In> totalTime_;

			// constructor
			Profiler()
			{
				ThisClass::totalTime_.fill(0);
			}

			static void reset()
			{
				ThisClass::totalTime_.fill(0);
			}

			// start stopwatch for profiler identified with parameter i
			template <int I>
			static void start()
			{
				static_assert(I < SIZE, "provided I is invalid");
				//::tools::assert::static_assert_msg(I < SIZE, "provided i is invalid");
				ThisClass::startTime_[I] = rdtsc();
			}

			// stop stopwatch for profiler identified with parameter i, and add the time since starting it to the total time.
			template <int I>
			static void stop()
			{
				static_assert(I < SIZE, "provided I is invalid");
				ThisClass::totalTime_[I] += (rdtsc() - ThisClass::startTime_[I]);
			}
			
			static void print_elapsed_cycles()
			{
				unsigned long long sum = 0;
				for (int i = 0; i < SIZE; ++i)
				{
					sum += ThisClass::totalTime_[i];
				}

				for (int i = 0; i < SIZE; ++i)
				{
					const double timeInMs = static_cast<double>(ThisClass::totalTime_[i]) / (1024. * 1024.);
					const double percent = 100 * (static_cast<double>(ThisClass::totalTime_[i]) / sum);
					printf("\tprofiler::elapsed cycles: i=%2zu: percent %5.2f; %s\n", i, percent, priv::elapsed_cycles_str(ThisClass::totalTime_[i]).c_str());
				}
			}
		};

		template <int SIZE>
		std::array<unsigned long long, SIZE> Profiler <SIZE, true>::startTime_;

		template <int SIZE>
		std::array<unsigned long long, SIZE> Profiler <SIZE, true>::totalTime_;
	}
}