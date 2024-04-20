#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <cmath>
using namespace std;

// Функция f(t, x)
double f(double t, double x) {
    return x + t;
}

// Функция φ(x) - начальное условие
double phi(double x) {
    return cos(M_PI * x);
}

// Функция ψ(t) - граничное условие
double psi(double t) {
    return exp(-t);
}

int main(int argc, char** argv) {
    // Инициализация MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Параметры задачи
    double T = 1.0;  // Время
    double X = 1.0;  // Пространство
    int K = 1000;    // Количество шагов по времени
    int M = 1000;    // Количество шагов по пространству
    double tau = T / K;  // Шаг по времени
    double h = X / M;    // Шаг по пространству
    double a = 1.0;      // Коэффициент переноса

    // Распределение данных между процессами
    int M_local = M / size;  // Количество узлов на процесс
    int start_m = rank * M_local;
    int end_m = (rank == size - 1) ? M : start_m + M_local - 1;

    // Выделение памяти для локальной части решения
    double** u_local = new double*[K + 1];
    for (int k = 0; k <= K; ++k) {
        u_local[k] = new double[M_local + 2]; // +2 для граничных значений
    }

    // Начальное условие (только для локальных узлов)
    for (int m = 1; m <= M_local; ++m) {
        u_local[0][m] = phi(start_m * h + (m - 1) * h);
    }

    // Граничное условие (только для процесса 0)
    if (rank == 0) {
        for (int k = 0; k <= K; ++k) {
            u_local[k][0] = psi(k * tau);
        }
    }

    // Численное решение (схема левый уголок)
    for (int k = 0; k < K; ++k) {
        // Обмен граничными значениями
        if (rank > 0) {
            MPI_Send(&u_local[k][1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (rank < size - 1) {
            MPI_Recv(&u_local[k][M_local + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Вычисление решения для локальных узлов
        for (int m = 1; m <= M_local; ++m) {
            u_local[k + 1][m] = u_local[k][m] - a * tau / h * (u_local[k][m] - u_local[k][m - 1]) + tau * f(k * tau, start_m * h + (m - 1) * h);
        }
    }

    // Сбор результатов на процессе 0
    if (rank == 0) {
        double** u = new double*[K + 1];
        for (int k = 0; k <= K; ++k) {
            u[k] = new double[M + 1];
        }

        // Копирование локальных результатов
        for (int k = 0; k <= K; ++k) {
            for (int m = 1; m <= M_local; ++m) {
                u[k][start_m + m - 1] = u_local[k][m]; 
            }
        }

        // Получение результатов от других процессов
        for (int p = 1; p < size; ++p) {
            start_m = p * M_local;
            end_m = (p == size - 1) ? M : start_m + M_local - 1;
            for (int k = 0; k <= K; ++k) {
                MPI_Recv(&u[k][start_m], end_m - start_m + 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // Вывод результата в файл .csv
        ofstream outfile("solution_1.csv");
        outfile << "t,x,u" << endl;
        for (int k = 0; k <= K; ++k) {
            for (int m = 0; m <= M; ++m) {
                outfile << k * tau << "," << m * h << "," << u[k][m] << endl;
            }
        }
        outfile.close();

        // Освобождение памяти
        for (int k = 0; k <= K; ++k) {
            delete[] u[k];
        }
        delete[] u;
    } else {
        // Отправка результатов на процесс 0
        for (int k = 0; k <= K; ++k) {
            MPI_Send(&u_local[k][1], M_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Освобождение памяти
    for (int k = 0; k <= K; ++k) {
        delete[] u_local[k];
    }
    delete[] u_local;

    MPI_Finalize();
    return 0;
}