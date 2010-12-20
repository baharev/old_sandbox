/*
 * Macros.hpp
 *
 *  Created on: Dec 19, 2010
 *      Author: ali
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_

#include <stdexcept>
#include "log4cxx/logger.h"

#define ASSERT(condition, message) { \
	if (!(condition)) { \
		LOG4CXX_FATAL(log, message); \
		throw std::logic_error("assertion failed"); \
	} \
}

#define DBG(message) LOG4CXX_DEBUG(log, message);

#define INIT_LOGGER( class_name ) \
		log4cxx::LoggerPtr \
		class_name::log(log4cxx::Logger::getLogger(#class_name));

#endif /* MACROS_HPP_ */
