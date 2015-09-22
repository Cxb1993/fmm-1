/*
 * point.h
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#ifndef POINT_H_
#define POINT_H_

#include <stdlib.h>
#include <array>

class key;
class Grid;

class Point {

public:
	Point();

	Point(double x, double y, double phi, unsigned bc);

	void set_x(double x);
	void set_y(double y);
	void set_phi(double phi);
	void set_bc(unsigned bc);
	void set_dxp(double dxp);
	void set_dxm(double dxm);
	void set_dyp(double dyp);
	void set_dym(double dym);
	void set_key(Key& key);
	void set_key(const unsigned i, const unsigned j);
	const double get_x() const ;
	const double get_y() const ;
	const double get_phi() const ;
	const double get_dxp() const ;
	const double get_dxm() const ;
	const double get_dyp() const;
	const double get_dym() const ;
	const unsigned get_bc() const;
	const Key& get_key() const ;
	const std::array<Point*, 4> get_neighbors();
	void set_neighbors(Grid& grid);

	virtual ~Point();

private:
	double x, y, phi, dxp, dxm, dyp, dym;
	unsigned bc;
	Key key;
	std::array<Point*, 4> neighbors;
};



#endif /* POINT_H_ */
