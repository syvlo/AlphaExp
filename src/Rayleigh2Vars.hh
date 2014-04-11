#ifndef RAYLEIGH_2_VARS_HH_
# define RAYLEIGH_2_VARS_HH_

/**
 * \file Rayleigh2Vars.hh
 * \author Sylvain Lobry
 * \date 04/11/2014
 * \brief File that defines the class Rayleig2Vars
 */

# include <vector>

/**
 * \class Rayleigh2Vars
 * \brief Class to be used to compute the data term for a Rayleigh distribution
 *				with 2 variables (e.g.: Ubv, Us).
 */
template <typename Input>
class Rayleigh2Vars
{
public:
    /**
     * \brief	Method to compute the data term
     * \param	ObservedValue		Observed value (i.e. y_s).
     * \param	RegularizedValue	Value to be tested (i.e. x_s)
     * \return  The value for the data term.
     */
    static double Compute(const Input& V, const Input& Ubv, const std::vector<Input>& gamma, const double BetaS, const double BetaBv);

    static double ComputeDataTermOnly(const Input& ObservedValue, const Input& RegularizedValue, const Input& Us, const double BetaBV);

    static Input ComputeUs(const Input& V, const Input& Ubv, const std::vector<Input>& gamma, const double BetaS, const double BetaBV);

    static double i0(const Input& x);
};

#include "Rayleigh2Vars.hxx"

#endif /* !RAYLEIGH_2_VARS_HH_ */
