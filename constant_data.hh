#ifndef CONSTANT_DATA_HH_
#define CONSTANT_DATA_HH_

#include <vector>
// Useful typedefs
typedef std::vector<double> vec_dbl;
typedef std::vector<int>    vec_int;
// Base cross section data:  Fuel, Water, Fuel+Rod
const double D[3][2]   = {{1.3000, 0.5000}, {1.5000, 0.5000}, {1.3000, 0.5000}};
const double nsf[3][2] = {{0.0060, 0.1950}, {0.0000, 0.0000}, {0.0060, 0.1950}};
const double sr[3][2]  = {{0.0318, 0.1140}, {0.0322, 0.0100}, {0.0318, 0.1180}};
const double s21[3]    = {0.022, 0.022, 0.022};
// Number of precursors
const int npc = 8;
// Delayed fractions
const double beta[] = { 0.000218, 0.001023, 0.000605, 0.001310,
                        0.002200, 0.000600, 0.000540, 0.000152 };
// Total delayed fraction
const double beta_tot = 0.006648;
// Decay constants
const double lambda[] = { 0.012467, 0.028292, 0.042524, 0.133042,
                          0.292467, 0.666488, 1.634781, 3.554601 };
// Velocities
const double velocity[] = {4.373835248681492e+07, 4.373835248681493e+05};
// Base total fine mesh
const int nfm = 42;
// Base fine mesh per material region
const int fm[] = {3,  17,  2,  17,  3};
// Number of regions
const int number_regions = 5;
// Region material map
const int mt[] = {1, 0, 2, 0, 1};

#endif // CONSTANT_DATA_HH_ 
