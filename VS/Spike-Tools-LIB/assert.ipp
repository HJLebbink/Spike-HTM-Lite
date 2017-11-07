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
		namespace priv
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
			inline void assert_msg(const std::string& message)
			{
				#if _DEBUG
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
				#endif
			}
		}

		inline void assert_msg(const bool cond)
		{
			#if _DEBUG
			if (!cond)
			{
				priv::assert_msg("ASSERT:");
			}
			#endif
		}

		template<typename... Args>
		inline void assert_msg(const bool cond, Args&&... a_args)
		{
			#if _DEBUG
			if (!cond)
			{
				std::ostringstream s;
				s << "----------------------------------" << std::endl << "ASSERT:";
				priv::addToStream(s, std::forward<Args>(a_args)...);
				priv::assert_msg(s.str());
			}
			#endif
		}
	}
}