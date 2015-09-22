//============================================================================
// Name        : testfmm.cpp
// Author      : Hossein
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "master.h"


struct PairCompare {
	bool operator()(const std::pair<Key, Point*>& pair_l, const std::pair<Key, Point*>& pair_r) {
		if (fabs(pair_l.second->get_phi()) < fabs(pair_r.second->get_phi()))
			return true;
		return false;
	}
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
//	VTPoint tentative;
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

//	solve_Eikonal(grid, tentative);

	std::cout << "done! \n";

	grid.print_result();

	return EXIT_SUCCESS;
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
