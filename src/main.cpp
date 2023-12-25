#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

// Глобальные переменные
int rank, commsize, lb, ub, nrows, n = 1680;

void get_chunk(int *l, int *u) {
  int rows_per_process = n / commsize;
  int remaining_rows = n % commsize;
  if (rank < remaining_rows) {
    *l = rank * (rows_per_process + 1);
    *u = *l + rows_per_process;
  } else {
    *l = remaining_rows * (rows_per_process + 1) + (rank - remaining_rows) * rows_per_process;
    *u = *l + rows_per_process - 1;
  }
}

int get_proc(int idx) {
  int rows_per_process = n / commsize;
  int remaining_rows = n % commsize;
  int threshold = remaining_rows * (rows_per_process + 1);
  if (idx < threshold) {
    // Вычисляем индекс, чтобы он попал в диапазон с доп. строкой
    return idx / (rows_per_process + 1);
  } else {
    // То же самое, но без строки
    return remaining_rows + (idx - threshold) / rows_per_process;
  }
}

double *get_input_matrix() {
  double *matrix = (double*)malloc(nrows * n * sizeof(double));
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < n; j++) {
      matrix[i * n + j] = std::min(n - j, n - i - rank * nrows);
    }
  }
  return matrix;
}

double *get_connected_matrix() {
  double *x = (double*)malloc(nrows * n * sizeof(double));
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < n; j++) {
      if (i + rank * nrows == j) {
        x[i * n + j] = 1.0;
      } else {
        x[i * n + j] = 0.0;
      }
    }
  }
  return x;
}

// Функция для обнуления cur_col в матрице A и сохранения диагонали с 1
// Всё дублируем в матрицу b
bool inverse_matrix(double *a, double *b, int cur_col) {
  double local_max = 0.0;
  int local_index = -1;
  for (int i = 0; i < nrows; i++) {
    int global_index = i + lb;
    if (global_index >= cur_col) {
      double value = std::fabs(a[i * n + cur_col]);
      if (value > local_max) {
        local_max = value;
        local_index = global_index;
      }
    }
  }

  // Храним максимумы в структуре
  struct {
    double value;
    int index;
  } local_data = {local_max, local_index}, global_data;

  // Глобальный максимум
  MPI_Allreduce(&local_data, &global_data, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  if (global_data.value == 0.0) {
    return false;
  }

  int diag_p = get_proc(cur_col);
  int main_p = get_proc(global_data.index);

  //Заполняем новую матрицу нулями
  double *s1 = (double*)malloc(n * 4 * sizeof(double));
  for (int i = 0; i < n * 4; i++) {
    s1[i] = 0.0;
  }
  
  if (rank == diag_p) {
    for (int i = n; i < n * 2; ++i) {
      s1[i] = a[(cur_col - rank * nrows) * n + i - n];
      s1[i + n * 2] = b[(cur_col - rank * nrows) * n + i - n];
    }
  }

  if (rank == main_p) {
    for (int i = 0; i < n; ++i) {
      s1[i] = a[(global_data.index - rank * nrows) * n + i];
      s1[i + n * 2] = b[(global_data.index - rank * nrows) * n + i];
    }
  }

  double *s2 = (double*)malloc(n * 4 * sizeof(double));
  MPI_Allreduce(s1, s2, n * 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double c = s2[cur_col];
  for (int i = 0; i < n; i++) {
    s2[i] = s2[i] / c;
    s2[i + n * 2] = s2[i + n * 2] / c;
  }

  if (rank == main_p) {
    for (int i = n; i < n * 2; i++) {
      a[(global_data.index - rank * nrows) * n + i - n] = s2[i];
      b[(global_data.index - rank * nrows) * n + i - n] = s2[i + n * 2];
    }
  }

  if (rank == diag_p) {
    for (int i = 0; i < n; i++) {
      a[(cur_col - rank * nrows) * n + i] = s2[i];
      b[(cur_col - rank * nrows) * n + i] = s2[i + n * 2];
    }
  }


  for (int i = 0; i < nrows; i++) {
    if (i + rank * nrows != cur_col) {
      c = a[i * n + cur_col];
      for (int j = 0; j < n; j++) {
        a[i * n + j] = a[i * n + j] - s2[j] * c;
        b[i * n + j] = b[i * n + j] - s2[j + n * 2] * c;
      }
    }
  }

  free(s1);
  free(s2);

  return true;
}

int main(int argc, char **argv) {
  double tstart = MPI_Wtime();

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc > 1) {
    n = std::atoi(argv[1]);
  }

  // По верхней и нижней границы, вычисляем сколько поток будет обрабатывать строк
  get_chunk(&lb, &ub);
  nrows = ub - lb + 1;

  double *a = get_input_matrix();
  double *b = get_connected_matrix();

  // Операции с i-м столбцом
  for (int i = 0; i < n; ++i) {
    if (!inverse_matrix(a, b, i)) {
      std::cerr << "No inverse matrix\n";
      return 1;
    }
  }

  double *recvbuf = nullptr;
  int *recvcounts = (int*)malloc(commsize * sizeof(int));
  int *displs = (int*)malloc(commsize * sizeof(int));
  if (rank == 0) {
    recvbuf = (double*)malloc(n * n * sizeof(double));
    for (int i = 0; i < commsize; ++i) {
      int l, u;
      get_chunk(&l, &u);
      recvcounts[i] = (u - l + 1) * n;
      displs[i] = l * n;
    }
  }

  // Соединяем всё в процессе 0
  MPI_Gatherv(a, nrows * n, MPI_DOUBLE, recvbuf, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Окончание время выполнения
  int t = MPI_Wtime() - tstart;

  // Имеем полученную матрицу в recvbuf

  if (rank == 0) {
    std::cout << commsize << " procs, n = " << n << ", t = " << t << " sec\n";
  }

  free(a);
  free(b);
  free(recvcounts);
  free(displs);
  if (rank == 0) {
    free(recvbuf);
  }
  MPI_Finalize();
  return 0;
}
