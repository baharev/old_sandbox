#ifndef D_HPP_
#define D_HPP_

#include <vector>

template <typename T> struct B;

template <typename T>
class D {

public:

	D();

	void f();

	~D();

private:

	std::vector<T> v1;
	std::vector<B<T> > v2;
};

#endif // D_HPP_
