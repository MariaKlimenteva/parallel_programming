#include "integral.hpp"
#include <iostream>

static double answer = 0.0;

void add_to_answer(double part_answer)
{
    std::lock_guard<std::mutex> lock_guard(std::mutex mtx);
    answer += part_answer;
}

double get_answer()
{
    return answer;
}

double get_max_abs_second_derivative(double a, double b, double n_points)
{
    double max_value = second_derivative_function(a);
    double step = (b - a) / n_points;
    for (int i = 1; i <= n_points; ++i) {
        double x = a + i * step;
        max_value = std::max(max_value, second_derivative_function(x));
    }
    return max_value;
}

inline double get_h(double a, double b, double error_rate)
{
    return sqrt(error_rate / get_max_abs_second_derivative(a, b, 100));
}

void integrate_part(double a, double b, double error_rate)
{
    auto h = get_h(a, b, error_rate);
    double answer = 0;
    auto left = a;
    for(; left + h < b; left += h)
    {
        answer += h * integrate_function(left + h / 2);
    }
    answer += (b - left) * integrate_function((b - left) / 2 + left);
    add_to_answer(answer);
}


double my_integrate(int processes, double a, double b, double error_rate)
{
    auto k1 = 1 / a / M_PI;
    auto k2 = 1 / b / M_PI;

    auto delta_k = (k2 - k1) / processes;

    auto thread_storage = std::unique_ptr<std::vector<std::thread>>(new std::vector<std::thread>());

    for(auto i = 0; i < processes; ++i)
    {
        thread_storage->push_back(std::thread(integrate_part, 1 / ((k1 + delta_k * i) * M_PI),
        1 / ((k1 + delta_k * (i + 1)) * M_PI), error_rate));
    }

    for (auto& thread : *thread_storage) {
        thread.join();
    }

    return get_answer();
}