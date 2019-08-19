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
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream

namespace tools
{
	namespace assert
	{
		#if _DEBUG
		constexpr bool DEBUG_ON = true;
		#else 
		constexpr bool DEBUG_ON = false;
		#endif


		namespace detail
		{
			void addToStream(std::ostringstream&)
			{
				//intentionally empty
			}

			template<typename T, typename... Args>
			void addToStream(std::ostringstream& a_stream, T&& a_value, Args&&... a_args)
			{
				a_stream << std::forward<T>(a_value);
				addToStream(a_stream, std::forward<Args>(a_args)...);
			}
			inline void assert_msg([[maybe_unused]] const std::string& message)
			{
				if constexpr (DEBUG_ON) {
					std::cerr << message;
					//std::cerr << std::flush; // no need to flush cerr
					if (true)
					{
						__debugbreak(); // generate int 3
					}
					else
					{
						char dummy;
						std::cin.get(dummy);
						//std::abort();
					}
				}
			}
		}

		inline void assert_msg([[maybe_unused]] const bool cond)
		{
			if constexpr (DEBUG_ON) {
				if (!cond)
				{
					detail::assert_msg("ASSERT:");
				}
			}
		}

		template<typename... Args>
		inline void assert_msg([[maybe_unused]] const bool cond, [[maybe_unused]] Args&&... a_args)
		{
			if constexpr (DEBUG_ON) {
				if (!cond)
				{
					std::ostringstream s;
					s << "----------------------------------" << std::endl << "ASSERT:";
					detail::addToStream(s, std::forward<Args>(a_args)...);
					detail::assert_msg(s.str());
				}
			}
		}
	}
}