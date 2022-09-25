
#include <iostream>
#include <stdexcept>
#include "derivatives.h"

// entry point
auto main()->int
{
	try
	{
		auto xi{ 4.0 };
		auto f = [](double x)->double { return -2.0 * x * x * x + 1.0; };
		std::cout << "\nf(x) = -2x^3 + 1 ===> f'(4) = " << dfdx(f, xi) << '\n';

		return EXIT_SUCCESS;
	}
	catch (const std::exception& xxx)
	{
		std::cerr << xxx.what() << std::endl;
		return EXIT_FAILURE;
	}
}
