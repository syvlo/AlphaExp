#include <iostream>
#include "Args.hh"
#include "TVL0DecompositionMinimizer.hh"
#include "Rayleigh2Vars.hh"
#include "ImwHelper.hh"

typedef TVL0DecompositionMinimizer<Rayleigh2Vars<unsigned> > Minimizer;

int main(int argc, char* argv[])
{
	Args args(argc, argv);
	if (args.getHelp()) //Help was asked and has been printed, leave.
	{
		return (0);
	}

	if (!args.checkConsistency())
    {
		std::cerr << "You have to define input and output images." << std::endl;
		args.printHelp();
		return (1);
    }

	std::cout << args << std::endl;

	//Label definition
	std::vector <unsigned> alpha;
	alpha.push_back(1);
	// for (unsigned i = 2; i < 255; i+= 1)
	// {
	// 	alpha.push_back(i);
	// }
	for (unsigned i = 30; i < 1000; i+= 30)
	{
		alpha.push_back(i);
	}

	std::vector<unsigned> gamma;
	// gamma.push_back(1);
    for (unsigned i = 0; i < 10000; i += 100)
    {
		gamma.push_back(i);
    }
    for (unsigned i = 10000; i < 100000; i += 1000)
    {
		gamma.push_back(i);
    }

	Minimizer minimizer(alpha, gamma, args.getBetaBV(), args.getBetaS());

	cv::Mat input = ReadImw(args.getInputImage());

	if (minimizer.compute(input))
	{
		WriteImw(minimizer.getOutputBV(), args.getOutputImageBV());
		WriteImw(minimizer.getOutputS(), args.getOutputImageS());
		WriteImw(minimizer.getOutputComplete(), args.getOutputImageComplete());
	}
	else
	{
	    std::cerr << "Error during computation." << std::endl;
	    return (1);
	}
	return (0);
}
