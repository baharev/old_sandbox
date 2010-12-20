#ifndef CAR_HPP_
#define CAR_HPP_

#include "DeclareLogger.hpp"

namespace ttext {

class car {

public:

	car();

	void drive();

private:

	int counter;

	LOGGER
};

}

#endif /* CAR_HPP_ */
