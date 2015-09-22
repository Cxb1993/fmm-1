/*
 * key.cpp
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#include "key.h"

Key::Key() :
		i(-1), j(-1) {
}

Key::Key(unsigned a, unsigned b) :
		i(a), j(b) {
}

unsigned Key::get_i() const {
	return i;
}

unsigned Key::get_j() const {
	return j;
}
//	The map object uses this expression to determine both the order the elements
//	follow in the container and whether two element keys are equivalent
//	(by comparing them reflexively: they are equivalent if !comp(a,b) && !comp(b,a)).
bool Key::operator<(const Key& rkey) const {
	if (i < rkey.i || (i == rkey.i && j < rkey.j))
		return true;

	return false;
}

Key& Key::operator=(const Key& rhs) {
	if (&rhs == this)
		return *this;

	i = rhs.i;
	j = rhs.j;

	return *this;
}

Key& Key::operator()(const unsigned ii, const unsigned jj) {
	i = ii;
	j = jj;
	return *this;
}

bool Key::operator==(const Key& key) {
	if (i == key.i && j == key.j)
		return true;
	return false;
}

Key::~Key(){}
