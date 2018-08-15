#ifndef _num_int_h
#define _num_int_h

/******************************************************************/
/* FROM http://rosettacode.org/wiki/Numerical_integration#C.2B.2B */
/******************************************************************/

// the integration routine
template<typename Method, typename F, typename Float>
double integrate(F f, Float a, Float b, int steps, Method m)
{
  double s = 0;
  double h = (b-a)/steps;
  for (int i = 0; i < steps; ++i)
    s += m(f, a + h*i, h);
  return h*s;
}

// methods
class rectangular
{
public:
  enum position_type { left, middle, right };
  rectangular(position_type pos): position(pos) {}
  template<typename F, typename Float>
  double operator()(F f, Float x, Float h) const
  {
    switch(position)
      {
      case left:
	return f(x);
      case middle:
	return f(x+h/2);
      case right:
	return f(x+h);
      }
  }
private:
  const position_type position;
};
 
class trapezium
{
public:
  template<typename F, typename Float>
  double operator()(F f, Float x, Float h) const
  {
    return (f(x) + f(x+h))/2;
  }
};
 
class simpson
{
public:
  template<typename F, typename Float>
  double operator()(F f, Float x, Float h) const
  {
    return (f(x) + 4*f(x+h/2) + f(x+h))/6;
  }
};

#endif
