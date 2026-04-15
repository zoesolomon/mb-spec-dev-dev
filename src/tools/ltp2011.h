#ifndef LTP2011_H
#define LTP2011_H

namespace h2o {

    struct ltp2011 {

	void operator()(const double, const double, const double, 
		      double&, double&);

	void cartesian(const size_t, const double*,
		       double*);

    };

} // namespace h2o

#endif // LTP2011_H
