#ifndef C_HPP_
#define C_HPP_

#include <vector>

template <typename T> struct B;

template <typename T>
class C {

public:

	void virtual f() = 0;

	virtual ~C();

	static void set_vec1(std::vector<T>* vec);

	static void set_vec2(std::vector<B<T> >* vec);

protected:

	explicit C(int x);

	T value() { return vec1->at(x); }

	const int x;

	T dummy;

	static std::vector<T>* vec1;

	static std::vector<B<T> >* vec2;
};

template <typename T>
class C2 : public C<T> {

public:

	C2();

private:

	void virtual f();

	const int y;
};

#endif // C_HPP_
