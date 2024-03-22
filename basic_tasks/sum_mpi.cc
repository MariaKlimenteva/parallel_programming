#include <iostream>
#include <mpi.h>

double calculatePartialSum(int start, int end) {
    double partialSum = 0.0;
    for (int i = start; i <= end; i++) {
        partialSum += 1.0 / i;
    }
    return partialSum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = std::stoi(argv[1]); // Получаем значение N из аргументов программы

    // Распределяем работу между процессами
    int chunkSize = N / size;
    int start = rank * chunkSize + 1;
    int end = (rank != size - 1) ? (rank + 1) * chunkSize : N;

    // Вычисляем частичную сумму для каждого процесса
    double partialSum = calculatePartialSum(start, end);

    // Собираем все частичные суммы на процессе с рангом 0
    double totalSum;
    MPI_Reduce(&partialSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Выводим результат на процессе с рангом 0
    if (rank == 0) {
        std::cout << "Сумма от 1 до " << N << " 1/n равна: " << totalSum << std::endl;
    }

    MPI_Finalize();
    return 0;
}
