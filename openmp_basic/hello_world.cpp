#include <iostream>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        #pragma omp critical
        {
            std::cout << "Поток " << thread_id << " из " << 
            num_threads << " потоков: Hello World!" << std::endl;
        }
    }

    return 0;
}