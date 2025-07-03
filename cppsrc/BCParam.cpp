/*
BCParam.cpp
struct stores parameters for BC and evaluates BC based thereon
Drew Murray
4/4/16
*/


//libraries
#include <math.h>
#include <iostream>
using std::cout;


//class includes
#include "BCParam.h"
#include "Grid.h"
#include <vector>
using std::vector;

#define pi 3.141592653589793238462643383279502884L

//default values, VERY important to keep these updated to prevent garbage memory from leaking in
extern const struct BCParam BCParamDefault = { {false, false, false}, {10,10,10,10,10,10}, -1.0F, false,false,false, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

