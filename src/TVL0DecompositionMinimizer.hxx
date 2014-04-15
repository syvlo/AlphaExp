#ifndef TV_L0_DECOMPOSITION_MINIMIZER_HXX_
# define TV_L0_DECOMPOSITION_MINIMIZER_HXX_

#include "TVL0DecompositionMinimizer.hh"
#include "Args.hh"
#include "ImwHelper.hh"
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
	//No need to free the images, OpenCV::Mat does it automatically.
}

template <typename DataTerm>
bool
TVL0DecompositionMinimizer<DataTerm>::compute(const cv::Mat& input)
{
	cv::Mat currentSolution = input;
	double currentSolutionEnergy = computeEnergy(input, currentSolution);

	std::clog << "Original energy = " << currentSolutionEnergy << std::endl;

	bool MadeProgress;
	do
	{
		std::clog << "Start iteration." << std::endl;

		MadeProgress = false;

		for (std::vector<unsigned>::const_iterator it = Alpha_.begin();
			 it != Alpha_.end();
			 ++it)
		{
			std::clog << "Find AExp..." << std::flush;
			cv::Mat contestant = findBestAExp(input, currentSolution,
											  *it);
			std::clog << "Done." << std::endl;
			std::clog << "Compute energy..." << std::flush;
			double contestantEnergy = computeEnergy(input, contestant);
			std::clog << "Done." << std::endl;
			std::clog << "Contestant (" << *it << ") = " << contestantEnergy << std::endl;
			if (contestantEnergy < currentSolutionEnergy)
			{
				currentSolution = contestant;
				currentSolutionEnergy = contestantEnergy;
				MadeProgress = true;
			}
		}

		if (MadeProgress)
		{
			std::clog << "Made progress, energy = " << currentSolutionEnergy
					  << std::endl;
		}

	} while (MadeProgress);

	std::clog << "Final energy = " << currentSolutionEnergy << std::endl;

	OutputBV_ = currentSolution;
    OutputS_ = cv::Mat(currentSolution.size(), input.type());
    OutputComplete_ = cv::Mat(currentSolution.size(), input.type());

	for (int i = 0; i < currentSolution.size().height; ++i)
		for (int j = 0; j < currentSolution.size().width; ++j)
		{
			OutputS_.at<unsigned short>(i, j) =
				DataTerm::ComputeUs(input.at<unsigned short>(i, j),
									currentSolution.at<unsigned short>(i, j),
									Gamma_, BetaS_, BetaBV_);
			OutputComplete_.at<unsigned short>(i, j) =
				OutputBV_.at<unsigned short>(i, j)
				+ OutputS_.at<unsigned short>(i, j);
		}

    return true;
}

template <typename DataTerm>
cv::Mat
TVL0DecompositionMinimizer<DataTerm>::findBestAExp(const cv::Mat& y, const cv::Mat& x, const unsigned Alpha)
{
    unsigned Height = y.size().height;
    unsigned Width = y.size().width;
    unsigned nbPix = Width * Height;

	//This part could be optimized regarding memory.
	unsigned maxNodes = 3 * nbPix;//NbPix + 2*NbPix (auxilary nodes in 4c)
	unsigned maxArcs = 6 * nbPix;//Not sure its true...

	Graph* g = new Graph(maxNodes, maxArcs);
    Graph::node_id* nodes = new Graph::node_id[nbPix];

	for (unsigned i = 0; i < Height; ++i)
		for (unsigned j = 0; j < Width; ++j)
		{
			unsigned index = i * Width + j;
			double linkToSource = DataTerm::Compute(y.at<unsigned short> (i, j),
													Alpha,
													Gamma_,
													BetaS_,
													BetaBV_);
			double linkToSink = ((x.at<unsigned short>(i, j) == Alpha) ?
								 INFTY
								 : DataTerm::Compute(y.at<unsigned short> (i, j),
													 x.at<unsigned short> (i, j),
													 Gamma_,
													 BetaS_,
													 BetaBV_));

			nodes[index] = g->add_node();
			g->set_tweights(nodes[index],
							linkToSource,
							linkToSink);
		}
	for (unsigned i = 0; i < Height; ++i)
		for (unsigned j = 0; j < Width; ++j)
		{
			if (i < Height - 1)
			{
				unsigned index1 = i * Width + j;
				unsigned index2 = (i + 1) * Width + j;
				if (x.at<unsigned short> (i, j)
					== x.at<unsigned short> (i + 1, j))
				{
					unsigned weight = abs(x.at<unsigned short> (i, j)
										  - Alpha);
					g->add_edge(nodes[index1], nodes[index2],
								weight, weight);
				}
				else
				{
					Graph::node_id aux = g->add_node();
					unsigned weight = abs(x.at<unsigned short> (i, j)
										  - Alpha);
					g->add_edge(nodes[index1], aux,
								weight, weight);
					weight = abs(x.at<unsigned short> (i + 1, j)
										  - Alpha);
					g->add_edge(nodes[index1], aux,
								weight, weight);
					g->set_tweights(aux, 0,
									abs(x.at<unsigned short> (i, j)
										- x.at<unsigned short> (i + 1, j)));
				}
			}
			if (j < Width - 1)
			{
				unsigned index1 = i * Width + j;
				unsigned index2 = i * Width + j + 1;
				if (x.at<unsigned short> (i, j)
					== x.at<unsigned short> (i, j + 1))
				{
					unsigned weight = abs(x.at<unsigned short> (i, j)
										  - Alpha);
					g->add_edge(nodes[index1], nodes[index2],
								weight, weight);
				}
				else
				{
					Graph::node_id aux = g->add_node();
					unsigned weight = abs(x.at<unsigned short> (i, j)
										  - Alpha);
					g->add_edge(nodes[index1], aux,
								weight, weight);
					weight = abs(x.at<unsigned short> (i, j + 1)
										  - Alpha);
					g->add_edge(nodes[index1], aux,
								weight, weight);
					g->add_tweights(aux, 0,
									abs(x.at<unsigned short> (i, j)
										- x.at<unsigned short> (i, j + 1)));
				}
			}
		}

	g->maxflow();
	cv::Mat newX(Height, Width, x.type());

	for (unsigned i = 0; i < Height; ++i)
		for (unsigned j = 0; j < Width; ++j)
		{
			if (g->what_segment(nodes[i * Width +j])
				!= Graph::SOURCE)
				newX.at<unsigned short>(i, j) = Alpha;
			else
				newX.at<unsigned short>(i, j) = x.at<unsigned short>(i, j);
		}

	delete g;
	delete[] nodes;

	return newX;
}

template <typename DataTerm>
double
TVL0DecompositionMinimizer<DataTerm>::computeEnergy(const cv::Mat& input, const cv::Mat& image)
{
	double Energy = 0;

	for (int i = 0; i < image.size().height; ++i)
	{
		for (int j = 0; j < image.size().width; ++j)
		{
			unsigned Us = DataTerm::ComputeUs(input.at<unsigned short>(i, j),
											  image.at<unsigned short>(i, j),
											  Gamma_, BetaS_, BetaBV_);

			if (image.at <unsigned short>(i, j) + Us)
				Energy += DataTerm::ComputeDataTermOnly(input.at<unsigned short>(i, j),
														image.at<unsigned short>(i, j),
														Us, BetaBV_);
			if (Us)
			{
				Energy += BetaS_;
			}

			if (i < image.size().height - 1)
			{
				Energy += BetaBV_ * abs(image.at<unsigned short>(i, j)
										- image.at<unsigned short>(i + 1, j));
			}
			if (j < image.size().width - 1)
			{
				Energy += BetaBV_ * abs(image.at<unsigned short>(i, j)
										- image.at<unsigned short>(i, j + 1));
			}
		}
	}

	return Energy;
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
