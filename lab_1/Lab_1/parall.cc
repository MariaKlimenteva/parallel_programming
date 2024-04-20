#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <cmath>

class TransportEquation {
public:
    TransportEquation(double T, double X, int K, int M, double a) :
        T(T), X(X), K(K), M(M), a(a), tau(T / K), h(X / M) {}

    // Функция f(t, x)
    virtual double f(double t, double x) const {
        return x + t;
    }

    // Функция φ(x) - начальное условие
    virtual double phi(double x) const {
        return cos(M_PI * x);
    }

    // Функция ψ(t) - граничное условие
    virtual double psi(double t) const {
        return exp(-t);
    }

    // Решение уравнения переноса (схема левый уголок)
    void solve(double** u) const {
        // Начальное условие
        for (int m = 0; m <= M; ++m) {
            u[0][m] = phi(m * h);
        }

        // Граничное условие
        for (int k = 0; k <= K; ++k) {
            u[k][0] = psi(k * tau);
        }

        // Численное решение
        for (int k = 0; k < K; ++k) {
            for (int m = 1; m <= M; ++m) {
                u[k + 1][m] = u[k][m] - a * tau / h * (u[k][m] - u[k][m - 1]) + tau * f(k * tau, m * h);
            }
        }
    }

protected:
    double T, X;
    int K, M;
    double a, tau, h;
}; //класс уравнения переноса

// Класс для параллельного решения с использованием MPI
class ParallelTransportSolver {
public:
    ParallelTransportSolver(const TransportEquation& equation, int argc, char** argv) :
        equation(equation) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        M_local = equation.M / size;
        start_m = rank * M_local;
        end_m = (rank == size - 1) ? equation.M : start_m + M_local - 1;
    }

    ~ParallelTransportSolver() {
        MPI_Finalize();
    }

    // Решение уравнения переноса
    void solve(double** u) const {
        // Выделение памяти для локальной части решения
        double** u_local = new double*[equation.K + 1];
        for (int k = 0; k <= equation.K; ++k) {
            u_local[k] = new double[M_local + 2];
        }

        // Начальное условие (только для локальных узлов)
        for (int m = 1; m <= M_local; ++m) {
            u_local[0][m] = equation.phi(start_m * equation.h + (m - 1) * equation.h);
        }

        // Граничное условие (только для процесса 0)
        if (rank == 0) {
            for (int k = 0; k <= equation.K; ++k) {
                u_local[k][0] = equation.psi(k * equation.tau);
            }
        }

        // Численное решение (схема левый уголок)
        for (int k = 0; k < equation.K; ++k) {
            // Обмен граничными значениями
            if (rank > 0) {
                MPI_Send(&u_local[k][1], 1, MPI_LONG_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            }
            if (rank < size - 1) {
                MPI_Recv(&u_local[k][M_local + 1], 1, MPI_LONG_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Вычисление решения для локальных узлов
            for (int m = 1; m <= M_local; ++m) {
                u_local[k + 1][m] = u_local[k][m] - equation.a * equation.tau / equation.h * (u_local[k][m] - u_local[k][m - 1]) + equation.tau * equation.f(k * equation.tau, start_m * equation.h + (m - 1) * equation.h);
            }
        }

        // Сбор результатов на процессе 0
        if (rank == 0) {
            // Копирование локальных результатов
            for (int k = 0; k <= equation.K; ++k) {
                for (int m = 1; m <= M_local; ++m) {
                    u[k][start_m + m - 1] = u_local[k][m]; 
                }
            }

            // Получение результатов от других процессов
            for (int p = 1; p < size; ++p) {
                int p_start_m = p * M_local;
                int p_end_m = (p == size - 1) ? equation.M : p_start_m + M_local - 1;
                for (int k = 0; k <= equation.K; ++k) {
                    MPI_Recv(&u[k][p_start_m], p_end_m - p_start_m + 1, MPI_LONG_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        } else {
            // Отправка результатов на процесс 0
            for (int k = 0; k <= equation.K; ++k) {
                MPI_Send(&u_local[k][1], M_local, MPI_LONG_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }

        // Освобождение памяти
        for (int k = 0; k <= equation.K; ++k) {
            delete[] u_local[k];
        }
        delete[] u_local;
    }

private:
    const TransportEquation& equation;
    int rank, size;
    int M_local, start_m, end_m;
};

#include 
#include 

class SolutionOutput {
public:
    // Вывод результата в файл .csv
    static void writeToFile(const double** u, int K, int M, double tau, double h, const std::string& filename) {
        std::ofstream outfile(filename);
        outfile << "t,x,u" << std::endl;
        for (int k = 0; k <= K; ++k) {
            for (int m = 0; m <= M; ++m) {
                outfile << k * tau << "," << m * h << "," << u[k][m] << std::endl;
            }
        }
        outfile.close();
    }

    // Освобождение памяти
    static void freeMemory(double** u, int K) {
        for (int k = 0; k <= K; ++k) {
            delete[] u[k];
        }
        delete[] u;
    }
};

int main(int argc, char** argv) {
    // Параметры задачи
    double T = 1.0;  // Время
    double X = 1.0;  // Пространство
    int K = 100;    // Количество шагов по времени
    int M = 100;    // Количество шагов по пространству
    double a = 1.0;      // Коэффициент переноса

    // Создание объекта уравнения переноса
    TransportEquation equation(T, X, K, M, a);

    // Выделение памяти для решения
    double** u = new double*[K + 1];
    for (int k = 0; k <= K; ++k) {
        u[k] = new double[M + 1];
    }

    // Создание объекта параллельного решателя
    ParallelTransportSolver solver(equation, argc, argv);

    // Решение уравнения переноса
    solver.solve(u);

    // Вывод результата в файл
    SolutionOutput::writeToFile(u, equation.K, equation.M, equation.tau, equation.h, "solution_1.csv");

    // Освобождение памяти
    SolutionOutput::freeMemory(u, equation.K);  

    return 0;
}

