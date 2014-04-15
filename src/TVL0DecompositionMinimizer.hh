#ifndef TV_L0_DECOMPOSITION_MINIMIZER_HH_
# define TV_L0_DECOMPOSITION_MINIMIZER_HH_

/**
 * \file TVL0DecompositionMinimizer.hh
 * \author Sylvain Lobry
 * \date 04/11/2014
 * \brief File that defines the minimizer using alpha expansion.
 */


# include <cv.h>
# include <vector>

/**
 * \class TVL0DecompositionMinimizer
 * \brief Class doing the computations using the alpha expansion algorithm.
 */
template <typename DataTerm>
class TVL0DecompositionMinimizer
{
public:
    /**
     * \brief		Default constructor.
     * \details		Will use default values for betas.
     * \param		alpha Look-up table for the values of alpha_i.
     * \param		gamma Look-up table for the values of gamma_i.
     */
    TVL0DecompositionMinimizer(const std::vector<unsigned>& alpha, const std::vector<unsigned>& gamma);

    /**
     * \brief		Constructor.
     * \param		alpha Look-up table for the values of alpha_i.
     * \param		gamma Look-up table for the values of gamma_i.
     * \param		BetaBV	Value for BetaBV.
     * \param		BetaS	Value for BetaS.
     */
    TVL0DecompositionMinimizer(const std::vector<unsigned>& alpha, const std::vector<unsigned>& gamma, double BetaBV, double BetaS);

    /**
     * \brief		Destructor.
     * \details		Will destroy the output images.
     */
    ~TVL0DecompositionMinimizer();

    /**
     * \brief		Method that construct the graph, launch the mincut and
     *				output the images.
     * \param		input	The image to be denoised.
     */
    bool compute(const cv::Mat& input);

	/**
	 * \brief			Compute the energy of one configuration.
	 * \param			input	The noisy image.
	 * \param			image	The configuration
	 * \retrurn			The energy
	 */
	double computeEnergy(const cv::Mat& input, const cv::Mat& image);

    /**
     * \brief		Setter for BetaBV.
     * \param		BetaBV	New value for BetaBV.
     */
    void setBetaBV(const double BetaBV);

    /**
     * \brief		Setter for BetaS.
     * \param		BetaS	New value for BetaS.
     */
    void setBetaS(const double BetaS);

    /**
     * \brief		Getter for the Bounded variations image.
     * \return		The image.
     */
    const cv::Mat& getOutputBV() const;

    /**
     * \brief		Getter for the Scatterers image.
     * \return		The image.
     */
    const cv::Mat& getOutputS() const;

    /**
     * \brief		Getter for the complete image.
     * \return		The image.
     */
    const cv::Mat& getOutputComplete() const;

private:
    /**
     * Beta parameter for the bounded variation component.
     */
    double BetaBV_;

    /**
     * Beta parameter for the sparse component.
     */
    double BetaS_;

    /**
     * look-up table for alpha_i values, i.e. possible values for BV image.
     */
    const std::vector<unsigned> Alpha_;

    /**
     * look-up table for gamma_i values, i.e. possible values for S image.
     */
    const std::vector<unsigned> Gamma_;

    /**
     * Stores the last bounded variations image computed
     */
    cv::Mat OutputBV_;

    /**
     * Stores the last scatterers image computed.
     */
    cv::Mat OutputS_;

    /**
     * Stores the last complete image computed.
     */
    cv::Mat OutputComplete_;

    //Look-up tables must be provided, so we remove default constructor.
    TVL0DecompositionMinimizer();

	/**
	 * \brief Compute the best solution within one alpha expansion of x.
	 * \param y			The original image.
	 * \param x			The current solution.
	 * \param Alpha		The alpha to be expanded.
	 * \return The best solution within one alpha expansion of x.
	 */
	cv::Mat findBestAExp(const cv::Mat& y, const cv::Mat& x, const unsigned Alpha);
};

#include "TVL0DecompositionMinimizer.hxx"

#endif /* !TV_L0_DECOMPOSITION_MINIMIZER_HH_ */

