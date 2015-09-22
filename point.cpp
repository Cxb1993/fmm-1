/*
 * point.cpp
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#include "master.h"

Point::Point() :
		x(0.), y(0.), phi(0.), dxp(0.), dxm(0.), dyp(0.), dym(0.), bc(0) {
}

Point::Point(double x, double y, double phi, unsigned bc) :
		x(x), y(y), phi(phi), bc(bc) {
	dxp = dxm = dyp = dym = 0.;
}

void Point::set_x(double x) {
	Point::x = x;
}

void Point::set_y(double y) {
	Point::y = y;
}

void Point::set_phi(double phi) {
	Point::phi = phi;
}

void Point::set_bc(unsigned bc) {
	Point::bc = bc;
}


void Point::set_dxp(double dxp) {
	Point::dxp = dxp;
}

void Point::set_dxm(double dxm) {
	Point::dxm = dxm;
}

void Point::set_dyp(double dyp) {
	Point::dyp = dyp;
}

void Point::set_dym(double dym) {
	Point::dym = dym;
}


void Point::set_key(Key& key) {
	Point::key = key;
}

void Point::set_key(const unsigned i, const unsigned j) {
	Point::key.i = i;
	Point::key.j = j;
}

const double Point::get_x() const {
	return x;
}

const double Point::get_y() const {
	return y;
}

const double Point::get_phi() const {
	return phi;
}

const double Point::get_dxp() const {
	return dxp;
}
const double Point::get_dxm() const {
	return dxm;
}

const double Point::get_dyp() const {
	return dyp;
}

const double Point::get_dym() const {
	return dym;
}

const unsigned Point::get_bc() const {
	return bc;
}

const Key& Point::get_key() const {
	return key;
}

const std::array<Point*, 4> Point::get_neighbors() {
	return neighbors;
}


void Point::set_neighbors(Grid& grid) {

	unsigned i = key.i;
	unsigned j = key.j;
	Point a=grid(i , j);
	switch (bc) {
		case 0: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = &grid(i, j - 1);
		}
			break;
		case 1: {
			neighbors[0] = NULL;
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = &grid(i, j - 1);

		}
			break;
		case 2: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = NULL;
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = &grid(i, j - 1);

		}
			break;
		case 3: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = NULL;
			neighbors[3] = &grid(i, j - 1);

		}
			break;
		case 4: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = NULL;

		}
			break;
		case 5: {
			neighbors[0] = NULL;
			neighbors[1] = NULL;
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = &grid(i, j - 1);

		}
			break;
		case 6: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = NULL;
			neighbors[2] = NULL;
			neighbors[3] = &grid(i, j - 1);

		}
			break;
		case 7: {
			neighbors[0] = &grid(i - 1, j);
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = NULL;
			neighbors[3] = NULL;

		}
			break;
		case 8: {
			neighbors[0] = NULL;
			neighbors[1] = &grid(i, j + 1);
			neighbors[2] = &grid(i + 1, j);
			neighbors[3] = NULL;

		}
			break;

	}

}

Point::~Point(){}

