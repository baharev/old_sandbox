
#include "D.hpp"
#include "C.hpp"
#include "A.hpp"

int main() {

	C<A>* c = new C2<A>;

	c->f();

	delete c;

	D<A> d;

	d.f();

	return 0;
}
