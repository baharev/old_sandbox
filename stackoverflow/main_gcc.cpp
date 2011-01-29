class car {
	friend void drive(const car c);
};

void drive(const car c) { }

int main() {
	car porsche;
	drive(porsche);
	void (*f)(const car);
	f = drive;
	return 0;
}
