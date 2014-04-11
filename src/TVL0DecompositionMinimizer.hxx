#ifndef TV_L0_DECOMPOSITION_MINIMIZER_HXX_
# define TV_L0_DECOMPOSITION_MINIMIZER_HXX_

#include "TVL0DecompositionMinimizer.hh"
#include "Args.hh"

#include "maxflow/graph.h"

#include <stdexcept>
#include <iostream>

#define INFTY 1E20

template <typename DataTerm>
TVL0DecompositionMinimizer<DataTerm>::TVL0DecompositionMinimizer(const std::vector<unsigned>& alpha, const std::vector<unsigned>& gamma)
	: BetaBV_(DEFAULT_BBV), BetaS_(DEFAULT_BS), Alpha_(alpha), Gamma_(gamma)
{
}

template <typename DataTerm>
TVL0DecompositionMinimizer<DataTerm>::TVL0DecompositionMinimizer(const std::vector<unsigned>& alpha, const std::vector<unsigned>& gamma, double BetaBV, double BetaS)
    : BetaBV_(BetaBV), BetaS_(BetaS), Alpha_(alpha), Gamma_(gamma)
{
}

template <typename DataTerm>
TVL0DecompositionMinimizer<DataTerm>::~TVL0DecompositionMinimizer()
{
}

template <typename DataTerm>
bool
TVL0DecompositionMinimizer<DataTerm>::compute(const cv::Mat& input)
{

    return true;
}

template <typename DataTerm>
void
TVL0DecompositionMinimizer<DataTerm>::setBetaBV(const double BetaBV)
{
    BetaBV_ = BetaBV;
}

template <typename DataTerm>
void
TVL0DecompositionMinimizer<DataTerm>::setBetaS(const double BetaS)
{
    BetaS_ = BetaS;
}

template <typename DataTerm>
const cv::Mat&
TVL0DecompositionMinimizer<DataTerm>::getOutputBV() const
{
    return OutputBV_;
}

template <typename DataTerm>
const cv::Mat&
TVL0DecompositionMinimizer<DataTerm>::getOutputS() const
{
    return OutputS_;
}

template <typename DataTerm>
const cv::Mat&
TVL0DecompositionMinimizer<DataTerm>::getOutputComplete() const
{
    return OutputComplete_;
}


#endif /* !TV_L0_DECOMPOSITION_MINIMIZER_HXX_ */
