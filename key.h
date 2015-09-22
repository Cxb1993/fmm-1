/*
 * key.h
 *
 *  Created on: Sep 21, 2015
 *      Author: haghakha
 */

#ifndef KEY_H_
#define KEY_H_

#include <stdio.h>
#include <stdlib.h>
#include <array>

class Key {

	friend class Point;

public:
	Key();
	Key(unsigned a, unsigned b);

	unsigned get_i() const;

	unsigned get_j() const;
//	The map object uses this expression to determine both the order the elements
//	follow in the container and whether two element keys are equivalent
//	(by comparing them reflexively: they are equivalent if !comp(a,b) && !comp(b,a)).
	bool operator<(const Key& rkey) const ;

	Key& operator=(const Key& rhs);

	Key& operator()(const unsigned ii, const unsigned jj) ;

	bool operator==(const Key& key) ;

	virtual ~Key();

private:
	unsigned i, j;

};


#endif /* KEY_H_ */
