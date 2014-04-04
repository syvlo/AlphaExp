#include <iostream>
#include "Args.hh"

int main(int argc, char* argv[])
{
	Args args(argc, argv);
	if (args.getHelp()) //Help was asked and has been printed, leave.
		return (0);

	if (!args.checkConsistency())
    {
		std::cerr << "You have to define input and output images." << std::endl;
		args.printHelp();
		return (1);
    }
	return (0);
}
