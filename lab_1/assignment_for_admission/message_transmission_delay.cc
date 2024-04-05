#include <mpi.h>
#include <iostream>
#include <chrono>

//берем третью вводную зацикливаем, полное время меряем и делим на число реальных пересылок
// 3 вводная - слайды семинар 3
//f = x+t
//cos x, exp.., склеивание функций в нуле

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 2) {
    std::cerr << "Программа требует ровно 2 процесса!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int message = 0;
  auto start = std::chrono::high_resolution_clock::now();

  for (auto i = 0; i < 100000; i++) {
    if (rank == 0) {
      MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/100000;

  if (rank == 0) {
    std::cout << "Задержка передачи сообщения: " << duration << " микросекунд" << std::endl;
  }

  MPI_Finalize();
  return 0;
}