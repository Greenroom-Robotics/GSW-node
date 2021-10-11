#ifndef MAP_STRUCT_H_INCLUDED
#define MAP_STRUCT_H_INCLUDED

/*************************
   Version 1.06
   by Randall Kent Whited
   rkwhited@gmail.com
**************************/

using namespace std;

/*************************************
(the value repeated in original)
"replaced" in VLA (very large arrays)
**************************************/
const double repeatSectionValue = 8.9999999999999998e+90;

/** used to map the zero based VLA (0-max) */
struct ARRAYINFO
{
/** array section beginning pos */
	unsigned begins;
/** array section ending (max pos) */
	unsigned ends;
/** number of values in the array section */
	unsigned elements;
/** was the array section composed of 'repeatSectionValue' elements */
	bool wasRepeatValue;

	/** initialize */
	void init()
	{
		begins = ends = elements = 0;
		wasRepeatValue = true;
	};
};

#endif // MAP_STRUCT_H_INCLUDED
