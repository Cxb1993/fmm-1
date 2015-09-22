/*
 * grid.cpp
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#include "master.h"

Grid::Grid(const int xsize, const int ysize, const double dx, const double dy) :
		xsize(xsize), ysize(ysize), dx(dx), dy(dy), grid(xsize, ysize) {

	for (int i = 0; i < xsize; ++i)
		for (int j = 0; j < ysize; ++j) {
			// setting location
			grid(i, j).set_x(i * dx);
			grid(i, j).set_y(j * dy);

			// setting boundar condition
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

			// setting key
			grid(i, j).set_key(i, j);

		}
}

Point& Grid::operator()(const unsigned& row, const unsigned& col) {
	return this->grid(row, col);
}

//	Point& Grid::operator()(const Key& key) {
//		return this->grid(key.get_i(), key.get_j());
//	}

const int Grid::get_xsize() const {
	return xsize;
}

const int Grid::get_ysize() const {
	return ysize;
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

	//	std::cout << "size of tentative   " << temp_tent.size() << "  size of accepted  "
	//	    << accepted.size() << std::endl;
	//
	//	for (Map::iterator it = accepted.begin(); it != accepted.end(); ++it) {
	//		key = it->first;
	////		Point* point = it->second;
	////		OrderMap::iterator it2 = tentative.find(point);
	////		if (it2 != tentative.end())
	//		if (temp_tent[key])
	//			temp_tent.erase(key);
	//	}
	//
	//	copy(temp_tent.begin(), temp_tent.end(), back_inserter(tentative));
	//	sort(tentative.begin(), tentative.end(), PairCompare());
	//
	//	std::cout << "size of tentative after  " << temp_tent.size() << std::endl;
	//
	//	for (unsigned i = 0; i < tentative.size(); ++i)
	//		std::cout << "key  " << tentative[i].first.get_i() << " , " << tentative[i].first.get_i()
	//		    << "  phi  " << tentative[i].second->get_phi() << std::endl;

}

Grid::~Grid(){}

