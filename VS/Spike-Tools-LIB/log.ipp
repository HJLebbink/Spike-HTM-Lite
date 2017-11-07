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
#include <mutex>

namespace tools
{
	namespace log
	{
		namespace priv
		{
			std::mutex log_screen;

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
		}
		inline void log_ERROR(const std::string message, const bool pause = true)
		{
			priv::log_screen.lock();
			std::cout << "ERROR:" << message << std::endl;
			priv::log_screen.unlock();

			if (pause)
			{
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
		template<bool PAUSE = true, typename... Args >
		inline void log_ERROR(Args&&... a_args)
		{
			std::ostringstream s;
			priv::addToStream(s, std::forward<Args>(a_args)...);
			log_ERROR(s.str(), PAUSE);
		}
		inline void WARNING(const std::string message, const bool pause = false)
		{
			std::cout << "WARNING:" << message << std::flush;
			if (pause)
			{
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
		template<bool PAUSE = false, typename... Args>
		inline void log_WARNING(Args&&... a_args)
		{
			std::ostringstream s;
			priv::addToStream(s, std::forward<Args>(a_args)...);
			WARNING(s.str(), PAUSE);
		}
		inline void log_INFO(const std::string message, const bool pause = false)
		{
			priv::log_screen.lock();
			std::cout << "INFO:" << message;
			priv::log_screen.unlock();

			if (pause)
			{
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
		template<bool PAUSE = false, typename... Args>
		inline void log_INFO(Args&&... a_args)
		{
			std::ostringstream s;
			priv::addToStream(s, std::forward<Args>(a_args)...);
			log_INFO(s.str(), PAUSE);
		}
		inline void log_INFO_DEBUG(const std::string message, const bool pause = false)
		{
			#if _DEBUG
			priv::log_screen.lock();
			std::cout << "INFO:" << message << std::flush;
			priv::log_screen.unlock();
			if (pause)
			{
				char dummy;
				std::cin.get(dummy);
			}
			#endif
		}
		template<bool PAUSE = false, typename... Args >
		inline void log_INFO_DEBUG(Args&&... a_args)
		{
			#if _DEBUG
			std::ostringstream s;
			priv::addToStream(s, std::forward<Args>(a_args)...);
			log_INFO_DEBUG(s.str(), PAUSE);
			#endif
		}
	}
}