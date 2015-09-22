/*
 * grid.h
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#ifndef GRID_H_
#define GRID_H_

#include "matrix.h"

class Key;
class Point;
//class Matrix<Point>;


struct PhiCompare {
	bool operator()(const Point& p_lhs, const Point& p_rhs) {
		if (fabs(p_lhs.get_phi()) < fabs(p_rhs.get_phi()))
			return true;

		return false;

	}
};

typedef double (*Pt2Path)(void* path, double x, double y);
typedef void (*Pt2Grad)(void* path, double x, double y, double* grad);
typedef std::map<Key, Point*> Map;
typedef std::map<Point&, Key, PhiCompare> OrderMap;
typedef std::vector<std::pair<Key, Point*> > VTPoint;

typedef std::map<Key, Point*> Map;

class Grid {
public:
	Grid(const int xsize, const int ysize, const double dx, const double dy);

	Point& operator()(const unsigned& row, const unsigned& col);

//	Point& operator()(const Key& key);

	const int get_xsize() const;

	const int get_ysize() const ;

	void print_result();
	void adjacent_to_interface(Map& accepted);
	void find_tentative(Map& accepted, Map& tentative);
	void compute_deivative();

	virtual ~Grid();

private:

	const unsigned xsize;
	const unsigned ysize;
	const double dx, dy;
	Matrix<Point> grid;

};

#endif /* GRID_H_ */
