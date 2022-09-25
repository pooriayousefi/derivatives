
#pragma once
#include <algorithm>
#include <numeric>
#include <execution>
#include <concepts>

namespace
{
	// ------------------------------------------------
	//
	// 
	//        yet another concepts and types
	//
	// 
	// ------------------------------------------------
	template<typename T>
	concept arithmetic_value = std::integral<T> || std::floating_point<T>;

	template<arithmetic_value... T>
	using real_value = std::common_type_t<float, T...>;

	template<typename It>
	concept real_valued_iterator = std::input_or_output_iterator<It> && std::floating_point<std::iter_value_t<It> >;

	template<typename It>
	concept integral_valued_iterator = std::input_or_output_iterator<It> && std::integral<std::iter_value_t<It> >;

	template<typename It>
	concept arithmetic_valued_iterator = std::input_or_output_iterator<It> && arithmetic_value<std::iter_value_t<It> >;

	template<typename F, typename... Args>
	concept real_valued_invocable = std::invocable<F, Args...> && std::floating_point<std::invoke_result_t<F, Args...> > && std::floating_point<Args...>;

	// ------------------------------------------------
	//                     n
	//                    d f(x)
	//                   --------
	//                       n
	//                     dx
	// ------------------------------------------------

	// first derivative based on four-point central difference
	template<std::floating_point T, real_valued_invocable<T> F>
	auto dfdx(F f, T xi)->T
	{
		auto eps{ std::numeric_limits<T>::epsilon() };
		auto h{ sqrt(eps) * xi };
		if (h == T{ 0 })
			h = eps;

		/*
		*           f(xi - 2h) - 8f(xi - h) + 8f(xi + h) - f(xi + 2h)
		*  df/dx = ---------------------------------------------------
		*                                 12h
		*/

		return (std::invoke(f, xi - T{ 2 } *h) - T{ 8 } *std::invoke(f, xi - h) + T{ 8 } *std::invoke(f, xi + h) - std::invoke(f, xi + T{ 2 } *h)) / (T{ 12 } *h);
	}

	// second derivative based on five-point central difference
	template<std::floating_point T, real_valued_invocable<T> F>
	auto d2fdx2(F f, T xi)->T
	{
		auto eps{ std::numeric_limits<T>::epsilon() };
		auto h{ sqrt(eps) * xi };
		if (h == T{ 0 })
			h = eps;

		/*
		*            -f(xi - 2h) + 16f(xi - h) - 30f(xi) + 16f(xi + h) - f(xi + 2h)
		* d2f/dx2 = ----------------------------------------------------------------
		*                                          2
		*                                       12h
		*/

		return (-std::invoke(f, xi - T{ 2 } *h) + T{ 16 } *std::invoke(f, xi - h) - T{ 30 } *std::invoke(f, xi) + T{ 16 } *std::invoke(f, xi + h) - std::invoke(f, xi + T{ 2 } *h)) / (T{ 12 } *h * h);
	}

	// third derivative based on six-point central difference
	template<std::floating_point T, real_valued_invocable<T> F>
	auto d3fdx3(F f, T xi)->T
	{
		auto eps{ std::numeric_limits<T>::epsilon() };
		auto h{ sqrt(eps) * xi };
		if (h == T(0))
			h = eps;

		/*
		*            f(xi - 3h) - 8f(xi - 2h) + 13f(xi - h) - 13f(xi + h) + 8f(xi + 2h) - f(xi + 3h)
		* d3f/dx3 = --------------------------------------------------------------------------------
		*                                                 3
		*                                               8h
		*/

		return (std::invoke(f, xi - T{ 3 } *h) - T{ 8 } *std::invoke(f, xi - T{ 2 } *h) + T{ 13 } *std::invoke(f, xi - h) - T{ 13 } *std::invoke(f, xi + h) + T{ 8 } *std::invoke(f, xi + T{ 2 } *h) - std::invoke(f, xi + T{ 3 } *h)) / (T{ 8 } *h * h * h);
	}

	// fourth derivative based on seven-point central difference
	template<std::floating_point T, real_valued_invocable<T> F>
	auto d4fdx4(F f, T xi)->T
	{
		auto eps{ std::numeric_limits<T>::epsilon() };
		auto h{ sqrt(eps) * xi };
		if (h == T(0))
			h = eps;

		/*
		*            f(xi - 3h) + 12f(xi - 2h) - 39f(xi - h) + 56f(xi) + 39f(xi + h) + 12f(xi + 2h) - f(xi + 3h)
		* d4f/dx4 = ---------------------------------------------------------------------------------------------
		*                                                     4
		*                                                   6h
		*/

		return (std::invoke(f, xi - T{ 3 } *h) + T{ 12 } *std::invoke(f, xi - T{ 2 } *h) - T{ 39 } *std::invoke(f, xi - h) + T{ 56 } *std::invoke(f, xi) + T{ 39 } *std::invoke(f, xi + h) + T{ 12 } *std::invoke(f, xi + T{ 2 } *h) - std::invoke(f, xi + T{ 3 } *h)) / (T{ 6 } *h * h * h * h);
	}

}