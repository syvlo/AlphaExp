#ifndef RAYLEIGH_2_VARS_HXX_
# define RAYLEIGH_2_VARS_HXX_

# include "Rayleigh2Vars.hh"
# include <cmath>
# include <limits>

template <typename Input>
double
Rayleigh2Vars<Input>::Compute(const Input& V,
							  const Input& Ubv,
							  const std::vector<Input>& gamma,
							  const double BetaS,
							  const double BetaBV)
{
    Input Us = Rayleigh2Vars<Input>::ComputeUs(V, Ubv, gamma, BetaS, BetaBV);

    if (Us == 0)
		return Rayleigh2Vars<Input>::ComputeDataTermOnly(V, Ubv, Us, BetaBV);

    return BetaS/BetaBV;
}


template <typename Input>
double
Rayleigh2Vars<Input>::ComputeDataTermOnly(const Input& V,
										  const Input& Ubv,
										  const Input& Us,
										  const double BetaBV)
{
	return (((double)V * (double)V) / (((double)Us + (double)Ubv) * ((double)Us + (double)Ubv)) + 2 * log((double)Ubv + (double)Us)) / BetaBV;
}

template <typename Input>
Input
Rayleigh2Vars<Input>::ComputeUs(const Input& V,
								const Input& Ubv,
								const std::vector<Input>& gamma,
								const double BetaS,
								const double BetaBV)
{
    double minValue = std::numeric_limits<double>::max();

    double tmp = Rayleigh2Vars<Input>::ComputeDataTermOnly(V, Ubv, 0, BetaBV);
    if (tmp < BetaS / BetaBV)
		return 0;

    Input min;
    for (typename std::vector<Input>::const_iterator it = gamma.begin();
		 it != gamma.end();
		 ++it)
    {
		double val = Rayleigh2Vars<Input>::ComputeDataTermOnly(V, Ubv,
															   *it, BetaBV);
		if (val < minValue)
		{
			minValue = val;
			min = *it;
		}
    }

    if (minValue + BetaS < ComputeDataTermOnly(V, Ubv, 0, BetaBV))
		return min;
    else
		return 0;
}

//From numerical recipies.
template <typename Input>
double
Rayleigh2Vars<Input>::i0(const Input& x)
{
    double ax,ans;
    double y;

    if ((ax=fabs(x)) < 3.75)
	{
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
    else
	{
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
    }
    return ans;
}

#endif /* !RAYLEIGH_2_VARS_HXX_ */
