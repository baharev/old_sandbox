#include "car.hpp"
#include "Macros.hpp"

INIT_LOGGER(ttext::car)

namespace ttext {

car::car() {

	counter = 0;
}

void car::drive() {

	DBG(counter);

	for (int i=0; i<100; ++i) {
		counter += i;
	}

	DBG(counter);
}

}
