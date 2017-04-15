#ifndef SRC_TEMPLATES_H_
#define SRC_TEMPLATES_H_

#define _USE_MATH_DEFINES
#include <cmath>

// refer http://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi
// Bring the 'difference' between two angles into [-pi; pi].
template<typename T>
T normalizeRadiansPiToMinusPi(T rad) {
  static const T PI2 = 2.*M_PI;
  // Copy the sign of the value in radians to the value of pi.
  T signed_pi = std::copysign(M_PI,rad);
  // Set the value of difference to the appropriate signed value between pi and -pi.
  rad = std::fmod(rad + signed_pi,(PI2)) - signed_pi;
  return rad;
}


#endif /* SRC_TEMPLATES_H_ */