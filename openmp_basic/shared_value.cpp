#include <iostream>
#include <omp.h>

int main() {
    int shared_value = 0;

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        #pragma omp critical
        {
            std::cout << "Поток " << thread_id << ": Доступ получен. Значение: " << shared_value << std::endl;
            shared_value++;
            std::cout << "Поток " << thread_id << ": Значение изменено на: " << shared_value << std::endl;
        }
    }

    std::cout << "Итоговое значение: " << shared_value << std::endl;

    return 0;
}