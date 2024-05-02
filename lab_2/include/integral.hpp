#pragma once

#include <memory.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <mutex>
#include <thread>

inline double integrate_function(float x)
{
    return cos(1 / x);
}

inline double second_derivative_function(float x)
{
    return -(2 * x * sin(1 / x) + cos(1 / x)) / pow(x, 4);
}

void add_to_answer(double part_answer);

double get_answer();

void integrate_part(double a, double b, double error_rate);

double my_integrate(int processes, double a, double b, double error_rate);