
#include "log4cxx/logger.h"
#include "log4cxx/propertyconfigurator.h"
#include "LoggerConf.hpp"

using namespace log4cxx;

void init(const char* config_file_name) {

	PropertyConfigurator::configure(config_file_name);

}

