//==============================================================================
//
// This code is part of ASOL (nonlinear system solver using affine arithmetic)
//
// Copyright (C) 2011 Ali Baharev
// All rights reserved. E-mail: <my_first_name.my_last_name@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//==============================================================================

#include <cstdio>
#include <fstream>
#include <string>
#include "diagnostics.hpp"
#include "vector_dump.hpp"

using namespace std;

namespace asol {

const char* const DUMP_FILE_NAME = "v_dump.txt";

const char* const INTERVAL_HEXA_FMT = "%la\t%la\n";


void parse_line(const string& line, std::vector<interval>& v) {

	if (line.empty()) {

		return;
	}

	double lb = 0, ub = 0;

	int k = sscanf(line.c_str(), INTERVAL_HEXA_FMT, &lb, &ub);

	if (k==0) {

		ASSERT2(false, "reading line failed: "<<line);
	}

	v.push_back(interval(lb, ub));
}

void load(std::vector<interval>& v) {

	v.clear();

	ifstream in;

	in.exceptions(ios_base::failbit | ios_base::badbit | ios_base::eofbit);

	in.open(DUMP_FILE_NAME);

	string line;

	do {

		getline(in, line);

		parse_line(line, v);

	}
	while (!line.empty());

	ASSERT(!v.empty());
}

void write(FILE* file, const char* format, const std::vector<interval>& v) {

	const int n = static_cast<int>(v.size());

	for (int i=0; i<n; ++i) {

		int k = fprintf(file, format, v.at(i).unchecked_inf(), v.at(i).unchecked_sup());

		if (k==0) {
			ASSERT2(false, "writing failed; format, i: "<<format<<", "<<i);
		}
	}

	fprintf(file, "\n");
}

void dump(const std::vector<interval>& v) {

	FILE * file = fopen(DUMP_FILE_NAME,"w");

	if (file==NULL) {
		ASSERT2(false, "failed to open file for dump");
	}

	write(file, INTERVAL_HEXA_FMT, v);

	write(file, "%g\t%g\n", v);

	fclose(file);
}

}
