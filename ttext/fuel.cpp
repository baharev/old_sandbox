#include "fuel.hpp"
#include "Macros.hpp"

INIT_LOGGER(ttext::fuel)

namespace ttext {

fuel::fuel() {

	counter = 1;
}

void fuel::check() {

	DBG("Hello");

	ASSERT(counter==1, "counter = "<<counter);

	for (int i=1; i<=5; ++i) {

		counter *= i;
	}

	DBG("counter="<<counter);
}

}
