#include <iostream>
#include <omp.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Неправильное количество аргументов. Использование: " 
                << argv[0] << " <N>" << std::endl;
        return 1; 
    }

    long long N = std::atoll(argv[1]);

    double sum = 0.0;
    #pragma omp parallel reduction(+:sum)
    {
        long long thread_id = omp_get_thread_num();
        long long num_threads = omp_get_num_threads();

        long long start = thread_id * N / num_threads;
        long long end = (thread_id + 1) * N / num_threads;

        for (long long i = start; i < end; ++i) {
            if (i + 1 == 0) { 
                #pragma omp critical
                {
                    std::cerr << "Ошибка: Деление на ноль!" << std::endl;
                }
                break;
            } 
            sum += 1.0 / (i + 1); 
        }
    }

    std::cout << "Сумма 1/N по N: " << sum << std::endl;

    return 0;
}