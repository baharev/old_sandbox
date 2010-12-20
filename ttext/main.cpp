
#include "LoggerConf.hpp"
#include "car.hpp"
#include "fuel.hpp"

using namespace ttext;

void run_car() {

	car Porsche;

	Porsche.drive();

}

void check_fuel() {

	fuel Diesel;

	Diesel.check();
}

int main(int argc, char* argv[]) {

	init(argv[1]);

	run_car();

	check_fuel();

    return 0;
}
