#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
// Функция f(t, x)
double f(double t, double x) {
  return x + t;
}
//cospix, exp(-t) , f = x+t
// Функция φ(x)
double phi(double x) {
  return cos(M_PI*x);
}

// Функция ψ(t)
double psi(double t) {
  return exp(-t); 
}

int main() {
  // Параметры задачи
  double T = 1.0;  // Время
  double X = 1.0;  // Пространство
  int K = 100;    // Количество шагов по времени
  int M = 100;    // Количество шагов по пространству
  double tau = T / K;  // Шаг по времени
  double h = X / M;    // Шаг по пространству
  double a = 1.0;      // Коэффициент переноса

  // Выделение памяти для решения
  double** u = new double*[K + 1];
  for (int k = 0; k <= K; ++k) {
      u[k] = new double[M + 1];
  }

// Начальное условие
  for (int m = 0; m <= M; ++m) {
      u[0][m] = phi(m * h);  // Используем функцию phi(x)
  }

  // Граничное условие
  for (int k = 0; k <= K; ++k) {
      u[k][0] = psi(k * tau);
  }

  // Численное решение (схема левый уголок)
  for (int k = 0; k < K; ++k) {
      for (int m = 1; m <= M; ++m) {
          u[k + 1][m] = u[k][m] - a * tau / h * (u[k][m] - u[k][m - 1]) + tau * f(k * tau, m * h);
      }
  }

  // Вывод результата в файл .csv чтобы потом визуалилировать
  ofstream outfile("solution.csv");
  outfile << "t,x,u" << endl;
  for (int k = 0; k <= K; ++k) {
    for (int m = 0; m <= M; ++m) {
      outfile << k * tau << "," << m * h << "," << u[k][m] << endl;
    }
  }
  outfile.close();

  return 0;
}