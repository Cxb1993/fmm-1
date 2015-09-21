//============================================================================
// Name        : testfmm.cpp
// Author      : Hossein
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <map>
#include "matrix.h"
#include <math.h>

typedef double (*Pt2Path)(void* path, double x, double y);
typedef void (*Pt2Grad)(void* path, double x, double y, double* grad);

class Key {

public:
	Key() :
			i(-1), j(-1) {
	}

	Key(unsigned a, unsigned b) :
			i(a), j(b) {
	}

	unsigned get_i() const {
		return i;
	}

	unsigned get_j() const {
		return j;
	}
//	The map object uses this expression to determine both the order the elements
//	follow in the container and whether two element keys are equivalent
//	(by comparing them reflexively: they are equivalent if !comp(a,b) && !comp(b,a)).
	bool operator<(const Key& rkey) const {
		if (i < rkey.i || (i == rkey.i && j < rkey.j))
			return true;

		return false;
	}

	Key& operator=(const Key& rhs) {
		if (&rhs == this)
			return *this;

		i = rhs.i;
		j = rhs.j;

		return *this;
	}

	Key& operator()(const unsigned ii, const unsigned jj) {
		i = ii;
		j = jj;
		return *this;
	}

	bool operator==(const Key& key) {
		if (i == key.i && j == key.j)
			return true;
		return false;
	}

private:
	unsigned i, j;

};

class Point {
public:
	Point() :
			x(0.), y(0.), phi(0.), dxp(0.), dxm(0.), dyp(0.), dym(0.), bc(0) {
	}

	Point(double x, double y, double phi, unsigned bc) :
			x(x), y(y), phi(phi), bc(bc) {
		dxp = dxm = dyp = dym = 0.;
	}

	void set_x(double x) {
		Point::x = x;
	}
	void set_y(double y) {
		Point::y = y;
	}
	void set_phi(double phi) {
		Point::phi = phi;
	}

	void set_bc(unsigned bc) {
		Point::bc = bc;
	}

	void set_dxp(double dxp) {
		Point::dxp = dxp;
	}

	void set_dxm(double dxm) {
		Point::dxm = dxm;
	}

	void set_dyp(double dyp) {
		Point::dyp = dyp;
	}

	void set_dym(double dym) {
		Point::dym = dym;
	}

	const double get_x() const {
		return x;
	}
	const double get_y() const {
		return y;
	}
	const double get_phi() const {
		return phi;
	}
	const double get_dxp() const {
		return dxp;
	}
	const double get_dxm() const {
		return dxm;
	}
	const double get_dyp() const {
		return dyp;
	}
	const double get_dym() const {
		return dym;
	}
	const unsigned get_bc() const {
		return bc;
	}

private:
	double x, y, phi, dxp, dxm, dyp, dym;
	unsigned bc;
};

struct PointCompare {
	bool operator()(const Point* p_lhs, const Point* p_rhs) {
		if (fabs(p_lhs->get_phi()) < fabs(p_rhs->get_phi()))
			return true;

		return false;

	}
};
//typedef std::map<Key, Point*, KeyCompare> Map;
typedef std::map<Key, Point*> Map;
typedef std::map<Point*, Key, PointCompare> OrderMap;

class Grid {
public:
	Grid(const int xsize, const int ysize, const double dx, const double dy) :
			xsize(xsize), ysize(ysize), dx(dx), dy(dy), grid(xsize, ysize) {

		for (int i = 0; i < xsize; ++i)
			for (int j = 0; j < ysize; ++j) {
				grid(i, j).set_x(i * dx);
				grid(i, j).set_y(j * dy);

				if (i == 0 && j == 0)
					grid(i, j).set_bc(8);
				else if (i == 0 && j == ysize - 1)
					grid(i, j).set_bc(5);
				else if (i == xsize - 1 && j == ysize - 1)
					grid(i, j).set_bc(6);
				else if (i == xsize - 1 && j == 0)
					grid(i, j).set_bc(7);
				else if (i == 0)
					grid(i, j).set_bc(1);
				else if (i == xsize - 1)
					grid(i, j).set_bc(3);
				else if (j == ysize - 1)
					grid(i, j).set_bc(2);
				else if (j == 0)
					grid(i, j).set_bc(4);

			}
	}

	Point& operator()(const unsigned& row, const unsigned& col) {
		return this->grid(row, col);
	}

	const int get_xsize() const {
		return xsize;
	}

	const int get_ysize() const {
		return ysize;
	}

	void print_result();
	void adjacent_to_interface(Map& accepted);
	void find_tentative(Map& accepted, Map& tentative);
	void compute_deivative();

private:

	const unsigned xsize;
	const unsigned ysize;
	const double dx, dy;
	Matrix<Point> grid;

};

class Ellipse {
public:
	Ellipse(double x_center, double y_center, double x_radius, double y_radius) :
			x_center(x_center), y_center(y_center), x_radius(x_radius), y_radius(y_radius) {
	}

	const double get_x_center() const {
		return x_center;
	}

	const double get_y_center() const {
		return y_center;
	}

	const double get_x_radius() const {
		return x_radius;
	}

	const double get_y_radius() const {
		return y_radius;
	}

private:
	const double x_center;
	const double y_center;
	const double x_radius;
	const double y_radius;
};

void initialize_distance(Point& point, Pt2Path& p_value_func, Pt2Grad& p_grad_func, void* ctx);
double path_ellipse(void* path, double x, double y);
void grad_path_ellipse(void* path, double x, double y, double* grad);
void solve_Eikonal(Grid& grid, Map& tentative);

int main(void) {

	double radius_sq;

	const int xsize = 101;
	const int ysize = 101;
	const double dx = .01, dy = .01, sq_rad = .1, x0 = .5, y0 = .5;

	Grid grid(xsize, ysize, dx, dy);

	for (int i = 0; i < xsize; ++i)
		for (int j = 0; j < ysize; ++j) {
			radius_sq = (i * dx - x0) * (i * dx - x0) + (j * dy - y0) * (j * dy - y0);
			if (radius_sq < sq_rad)
				grid(i, j).set_phi(-1000000.);
			else
				grid(i, j).set_phi(1000000.);
		}

	Map accepted, tentative, distant;
//	OrderMap tentative;
	grid.adjacent_to_interface(accepted);

	double ellipse_radi = sqrt(sq_rad);
	Ellipse ellipse(x0, y0, ellipse_radi, ellipse_radi);

	void *ctx = (void*) &ellipse;
	Pt2Path p2path = path_ellipse;
	Pt2Grad p2grad = grad_path_ellipse;

	for (Map::iterator it = accepted.begin(); it != accepted.end(); ++it)
		initialize_distance(*it->second, p2path, p2grad, ctx);

//	for (unsigned i = 0; i < xsize; ++i)
//		for (unsigned j = 0; j < ysize; ++j)
//			initialize_distance(grid(i, j), p2path, p2grad, ctx);

	grid.find_tentative(accepted, tentative);

//	grid.compute_deivative();

	solve_Eikonal(grid, tentative);

	std::cout << "done! \n";

	grid.print_result();

	return EXIT_SUCCESS;
}

void Grid::print_result() {

	char filename[256];

	sprintf(filename, "data_file%d.tec", 1);

	FILE* fp;

	fp = fopen(filename, "w");

	fprintf(fp, "TITLE= \" similar_to_FD\"\n");

	fprintf(fp, "VARIABLES = \"X\", \"Y\", \"PHI\"\n");

	fprintf(fp, "ZONE I=%d, J=%d, DATAPACKING=POINT \n", xsize, ysize);

	for (unsigned i = 0; i < xsize; ++i)
		for (unsigned j = 0; j < ysize; ++j)
			fprintf(fp, "%e %e %e ", grid(i, j).get_x(), grid(i, j).get_y(), grid(i, j).get_phi());
	fprintf(fp, "\n");

	fprintf(fp, "\n");

	fclose(fp);
}
void initialize_distance(Point& point, Pt2Path& p_value_func, Pt2Grad& p_grad_func, void* ctx) {

	double d1[2] = { 0., 0. }, d2[2] = { 0., 0. }, grad[2], pvalue, grad_dot, x_new, y_new, x_old,
			y_old, x_half, y_half, epslon = .01;

	// initializing the solution
	x_new = x_old = point.get_x();
	y_new = y_old = point.get_y();

	if (p_value_func == path_ellipse && fabs(fabs(y_new) - fabs(x_new)) < epslon) {
		x_new = x_old + 2 * epslon;
		y_new = y_old - epslon;
	}

	int iter = 0;
	double sgn = 1.;

	do {
		iter++;

		p_grad_func(ctx, x_new, y_new, grad);
		pvalue = p_value_func(ctx, x_new, y_new);

		grad_dot = grad[0] * grad[0] + grad[1] * grad[1];

		d1[0] = -pvalue * grad[0] / grad_dot;
		d1[1] = -pvalue * grad[1] / grad_dot;
		x_half = x_new + d1[0];
		y_half = y_new + d1[1];

		d2[0] = (x_old - x_new) - (x_old - x_new) * grad[0] / grad_dot * grad[0];
		d2[1] = (y_old - y_new) - (y_old - y_new) * grad[1] / grad_dot * grad[1];

		x_new = x_half + d2[0];
		y_new = y_half + d2[1];

	} while (sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d2[0] * d2[0] + d2[1] * d2[1]) > 1e-7 && iter < 20);

	if (point.get_phi() < 0.)
		sgn = -1.;

	point.set_phi(sgn * sqrt((x_new - x_old) * (x_new - x_old) + (y_new - y_old) * (y_new - y_old)));
//	std::cout << point.get_phi() << std::endl;
}

double path_ellipse(void* path, double x, double y) {

	Ellipse* ellipse = (Ellipse*) path;

	const double x_radius = ellipse->get_x_radius();
	const double x_rad_sq = x_radius * x_radius;
	const double y_radius = ellipse->get_y_radius();
	const double y_rad_sq = y_radius * y_radius;
	const double x_center = ellipse->get_x_center();
	const double y_center = ellipse->get_y_center();

	return (x - x_center) * (x - x_center) / x_rad_sq + (y - y_center) * (y - y_center) / y_rad_sq - 1;
}

void grad_path_ellipse(void* path, double x, double y, double* grad) {

	Ellipse* ellipse = (Ellipse*) path;

	const double x_radius = ellipse->get_x_radius();
	const double x_rad_sq = x_radius * x_radius;
	const double y_radius = ellipse->get_y_radius();
	const double y_rad_sq = y_radius * y_radius;
	const double x_center = ellipse->get_x_center();
	const double y_center = ellipse->get_y_center();

	grad[0] = 2 * (x - x_center) / x_rad_sq;
	grad[1] = 2 * (y - y_center) / y_rad_sq;

}

void Grid::adjacent_to_interface(Map& accepted) {

	unsigned i, j;

//	for (i = 0; i < xsize; ++i)
//		for (j = 0; j < ysize; ++j)
//			if (!(grid(i, j).get_bc())
//					&& (grid(i, j).get_phi() * grid(i + 1, j).get_phi() < 0
//							|| grid(i, j).get_phi() * grid(i - 1, j).get_phi() < 0
//							|| grid(i, j).get_phi() * grid(i, j - 1).get_phi() < 0
//							|| grid(i, j).get_phi() * grid(i, j + 1).get_phi() < 0))
//				accepted[Key(i, j)] = &grid(i, j);
//
//			else if(grid(i, j).get_bc()){
//
//
//			}

	for (i = 1; i < xsize - 1; ++i)
		for (j = 1; j < ysize - 1; ++j)
			if (grid(i, j).get_phi() * grid(i + 1, j).get_phi() < 0
					|| grid(i, j).get_phi() * grid(i - 1, j).get_phi() < 0
					|| grid(i, j).get_phi() * grid(i, j - 1).get_phi() < 0
					|| grid(i, j).get_phi() * grid(i, j + 1).get_phi() < 0)
				accepted[Key(i, j)] = &grid(i, j);

	i = 0;
	for (j = 1; j < ysize - 1; ++j)
		if (grid(i, j).get_phi() * grid(i + 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j - 1).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j + 1).get_phi() < 0)
			accepted[Key(i, j)] = &grid(i, j);

	i = xsize - 1;
	for (j = 1; j < ysize - 1; ++j)
		if (grid(i, j).get_phi() * grid(i - 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j - 1).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j + 1).get_phi() < 0)
			accepted[Key(i, j)] = &grid(i, j);

	j = 0;
	for (i = 1; i < xsize - 1; ++i)
		if (grid(i, j).get_phi() * grid(i + 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i - 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j + 1).get_phi() < 0)
			accepted[Key(i, j)] = &grid(i, j);

	j = ysize - 1;
	for (i = 1; i < xsize - 1; ++i)
		if (grid(i, j).get_phi() * grid(i + 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i - 1, j).get_phi() < 0
				|| grid(i, j).get_phi() * grid(i, j - 1).get_phi() < 0)
			accepted[Key(i, j)] = &grid(i, j);

}

void Grid::compute_deivative() {

	double dxp = 0., dxm = 0., dyp = 0., dym = 0.;

	for (unsigned i = 0; i < xsize; ++i)
		for (unsigned j = 0; j < ysize; ++j) {
			if (!grid(i, j).get_bc()) {

				dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
				dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
				dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
				dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;

			} else {
				switch (grid(i, j).get_bc()) {
				case 1: {
					dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
					dxm = 0.;
					dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
					dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;
				}
					break;
				case 2: {
					dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
					dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
					dyp = 0.;
					dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;
				}
					break;
				case 3: {
					dxp = 0.;
					dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
					dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
					dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;
				}
					break;
				case 4: {

					dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
					dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
					dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
					dym = 0.;

				}
					break;
				case 5: {

					dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
					dxm = 0.;
					dyp = 0.;
					dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;

				}
					break;
				case 6: {

					dxp = 0.;
					dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
					dyp = 0.;
					dym = (grid(i, j).get_phi() - grid(i, j - 1).get_phi()) / dy;

				}
					break;
				case 7: {

					dxp = 0.;
					dxm = (grid(i, j).get_phi() - grid(i - 1, j).get_phi()) / dx;
					dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
					dym = 0.;

				}
					break;
				case 8: {

					dxp = (grid(i + 1, j).get_phi() - grid(i, j).get_phi()) / dx;
					dxm = 0.;
					dyp = (grid(i, j + 1).get_phi() - grid(i, j).get_phi()) / dy;
					dym = 0.;

				}
					break;
				}

			}
			grid(i, j).set_dxp(dxp);
			grid(i, j).set_dxm(dxm);
			grid(i, j).set_dyp(dyp);
			grid(i, j).set_dym(dym);
		}

}

void Grid::find_tentative(Map& accepted, Map& tentative) {

	Key key;
	unsigned i, j;

	for (Map::iterator it = accepted.begin(); it != accepted.end(); ++it) {

		key = it->first;
		i = key.get_i();
		j = key.get_j();

		if (i < xsize - 1 && i > 0 && j < ysize - 1 && j > 0) {
			tentative[Key(i + 1, j)] = &grid(i + 1, j);
			tentative[Key(i - 1, j)] = &grid(i - 1, j);
			tentative[Key(i, j + 1)] = &grid(i, j + 1);
			tentative[Key(i, j - 1)] = &grid(i, j - 1);
		} else if (i == xsize - 1) {
			tentative[Key(i - 1, j)] = &grid(i - 1, j);
			tentative[Key(i, j + 1)] = &grid(i, j + 1);
			tentative[Key(i, j - 1)] = &grid(i, j - 1);
		} else if (i == 0) {
			tentative[Key(i + 1, j)] = &grid(i + 1, j);
			tentative[Key(i, j + 1)] = &grid(i, j + 1);
			tentative[Key(i, j - 1)] = &grid(i, j - 1);
		} else if (j == 0) {
			tentative[Key(i + 1, j)] = &grid(i + 1, j);
			tentative[Key(i - 1, j)] = &grid(i - 1, j);
			tentative[Key(i, j + 1)] = &grid(i, j + 1);
		} else {
			tentative[Key(i + 1, j)] = &grid(i + 1, j);
			tentative[Key(i - 1, j)] = &grid(i - 1, j);
			tentative[Key(i, j - 1)] = &grid(i, j - 1);

		}
	}

	std::cout << "size of tentative   " << tentative.size() << "  size of accepted  "
			<< accepted.size() << std::endl;

	for (Map::iterator it = accepted.begin(); it != accepted.end(); ++it) {
		key = it->first;
//		Point* point = it->second;
//		OrderMap::iterator it2 = tentative.find(point);
//		if (it2 != tentative.end())
		if (tentative[key])
			tentative.erase(key);
	}

	std::cout << "size of tentative after  " << tentative.size() << std::endl;

}

void min_tentative(Map tentative, Key& key) {

	double phi_min = 1000000.;
	for (Map::iterator it = tentative.begin(); it != tentative.end(); ++it) {
		if ((*it->second).get_phi() < phi_min)
			key = it->first;
	}

}

void solve_Eikonal(Grid& grid, Map& tentative) {

	Key key;
	unsigned i, j;
	double a, b, phi, grid_size = .01;

	for (Map::iterator it = tentative.begin(); it != tentative.end(); ++it) {

		key = it->first;
		i = key.get_i();
		j = key.get_j();

		a = std::min(grid(i + 1, j).get_phi(), grid(i - 1, j).get_phi());
		b = std::min(grid(i, j + 1).get_phi(), grid(i, j - 1).get_phi());

		if (fabs(a - b) < grid_size)
			phi = .5 * (a + b + sqrt(2 * grid_size * grid_size - (a - b) * (a - b)));
		else
			phi = std::min(a, b) + grid_size;

		if (fabs(phi) < fabs((*it->second).get_phi())) {
			if ((*it->second).get_phi() < 0)
				(*it->second).set_phi(-phi);
			(*it->second).set_phi(phi);
		}

	}

}
