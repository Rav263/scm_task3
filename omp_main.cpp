#include <iostream>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define THREADS 4


class Mat2 {
    std::vector<double> data;
    int size[2];
public:
    Mat2(int size_x, int size_y) {
        this->data = std::vector<double>(size_x * size_y);
        size[0] = size_x;
        size[1] = size_y;
    }

    double& operator()(int i, int j) {
        return data[i * size[1] + j];
    }

    double *get_data() {
        return this->data.data();
    }

    int get_size() {
        return this->data.size();
    }
};

class Mat3 {
    std::vector<double> data;
    int size[3]={};
    public:
    Mat3(int size_x, int size_y, int size_z) {
        this->data = std::vector<double>(size_y * size_z * size_x);
        size[0] = size_x;
        size[1] = size_y;
        size[2] = size_z;
    }

    Mat3(int size[3]) {
        this->data = std::vector<double>(size[0] * size[1] * size[2]);
        this->size[0] = size[0];
        this->size[1] = size[1];
        this->size[2] = size[2];
    }
    double& operator()(int i, int j, int k) {
        return data[(i * size[1] + j) * size[2] + k];
    }
    void get_facet(int i, int val, Mat2 &facet) {
        if (i == 0) {
#pragma omp parallel for collapse(2) num_threads(THREADS)
            for(int j = 0; j < size[1]; j++) {
                for(int k = 0; k < size[2]; k++) {
                    facet(j, k) = (*this)(val, j, k);
                }
            }
        }
        if (i == 1) {
#pragma omp parallel for collapse(2) num_threads(THREADS)
            for(int j = 0; j < size[0]; j++) {
                for(int k = 0; k < size[2]; k++) {
                    facet(j, k) = (*this)(j, val, k);
                }
            }
        }

        if (i == 2) {
#pragma omp parallel for collapse(2) num_threads(THREADS)
            for(int j = 0; j < size[0]; j++) {
                for(int k = 0; k < size[1]; k++) {
                    facet(j, k) = (*this)(j, k, val);
                }
            }
        }
    }

    void fill_zero() {
#pragma omp parallel for collapse(1) num_threads(THREADS)
        for(int i = 0; i < size[0] * size[1] * size[2]; i++) {
            data[i] = 0;
        }
    }
};

class Block{
    Mat3 *data[3];
    Mat2 *facets[6];
    Mat2 *buffs[6];
    int dims[3];
    int coords[3];
    int bsize[3];
    int shift[3];
    int start[3];
    int end[3];
    int grid_size;
    int rank;
    int size;

    double L[3];
    double a_t;
    double h;
    double tau;
    MPI_Comm grid;
public:
    Block(int rank, int size, MPI_Comm grid, int grid_size, int dims[3], double L[3]) :
        rank{rank}, size{size}, grid{grid}, grid_size{grid_size} {
        MPI_Cart_coords(this->grid, this->rank, 3, this->coords);
        for (int i = 0; i < 3; i++) {
            this->dims[i] = dims[i];
            this->L[i] = L[i];
            this->bsize[i] = this->grid_size / this->dims[i];
            this->shift[i] = this->coords[i] * (this->bsize[i]);
            this->start[i] = coords[i] == 0;
            this->end[i] = coords[i] == dims[i] - 1;
            //std::cout << this->end[i] << " ";
        }
       //std::cout <<  std::endl;
        for (int i = 0; i < 3; i++) { 
            this->data[i] = new Mat3(this->bsize);
            this->data[i]->fill_zero();
        }
        this->h = L[0] / (this->grid_size - 1);
        this->a_t = M_PI * std::sqrt(3.0/(L[0]*L[0]));
        this->tau = 0.0001;//1 / (2*std::sqrt(3.0 / (this->h * this->h)));
        //std::cout << this->h << " " << this->a_t << " " << this->tau << std::endl;
        
        buffs[0] = new Mat2(bsize[1], bsize[2]);
        buffs[1] = new Mat2(bsize[1], bsize[2]);
        buffs[2] = new Mat2(bsize[0], bsize[2]);
        buffs[3] = new Mat2(bsize[0], bsize[2]);
        buffs[4] = new Mat2(bsize[0], bsize[1]);
        buffs[5] = new Mat2(bsize[0], bsize[1]);
        
        facets[0] = new Mat2(bsize[1], bsize[2]);
        facets[1] = new Mat2(bsize[1], bsize[2]);
        facets[2] = new Mat2(bsize[0], bsize[2]);
        facets[3] = new Mat2(bsize[0], bsize[2]);
        facets[4] = new Mat2(bsize[0], bsize[1]);
        facets[5] = new Mat2(bsize[0], bsize[1]);
    }

    double get_point(int i, int dim) {
        return (double)(shift[dim] + i) / (double)(grid_size - 1);
    }
    Mat3 *get_data(int i) {
        return this->data[i];
    }
    double an_sol(int i, int j, int k, double t) {
        return std::sin((M_PI * this->get_point(i, 0)) / L[0]) * 
            std::sin((M_PI * this->get_point(j, 1)) / L[1]) * 
            std::sin((M_PI * this->get_point(k, 2)) / L[2]) * std::cos(a_t * t * tau);
    }

    double val(int i, int j, int k, Mat3 &data) {
        if (i < 0)         return (*facets[0])(j, k);
        if (i >= bsize[0]) return (*facets[1])(j, k);
        if (j < 0)         return (*facets[2])(i, k);
        if (j >= bsize[1]) return (*facets[3])(i, k);
        if (k < 0)         return (*facets[4])(i, j);
        if (k >= bsize[2]) return (*facets[5])(i, j);

        return data(i, j, k);
    }

    void update_facets(int ind) {
        for (int i = 0; i < 3; i++) {
            data[ind]->get_facet(i, 0, *facets[i*2]);
            data[ind]->get_facet(i, bsize[i] - 1, *facets[i * 2 + 1]);
        }
    }

    void exchange_facets() {
        int rank_src, rank_dst;
     

        for (int i = 0; i < 3; i++) {
            MPI_Request request1[1] = {nullptr};
            MPI_Request request2[1] = {nullptr};
            MPI_Cart_shift(grid, i, 1, &rank_src, &rank_dst);
            if (rank_src != MPI_PROC_NULL) {
                MPI_Isend(facets[i * 2]->get_data(), facets[i * 2]->get_size(),
                          MPI_DOUBLE, rank_src, i * 2, grid, &request1[0]);
            }
            if (rank_dst != MPI_PROC_NULL) {
                MPI_Isend(facets[i * 2 + 1]->get_data(), facets[i * 2 + 1]->get_size(),
                          MPI_DOUBLE, rank_dst, i * 2 + 1, grid, &request2[0]);
            }
            if (rank_src != MPI_PROC_NULL) {
                MPI_Recv(buffs[i * 2]->get_data(), buffs[i * 2]->get_size(),
                         MPI_DOUBLE, rank_src, i * 2 + 1, grid, MPI_STATUS_IGNORE);
            }
            if (rank_dst != MPI_PROC_NULL) {
                MPI_Recv(buffs[i * 2 + 1]->get_data(), buffs[i * 2 + 1]->get_size(),
                         MPI_DOUBLE, rank_dst, i * 2, grid, MPI_STATUS_IGNORE);
            }
            if (request1[0] != nullptr) {
                MPI_Waitall(1, request1, MPI_STATUS_IGNORE);
            }
            if (request2[0] != nullptr) {
                MPI_Waitall(1, request2, MPI_STATUS_IGNORE);
            }
            if (rank_src != MPI_PROC_NULL) {
                std::swap(facets[i * 2], buffs[i * 2]);
            }
            if (rank_dst != MPI_PROC_NULL) {
                std::swap(facets[i * 2 + 1], buffs[i * 2 + 1]);
            }
        }
    }

    double func_lap(int i, int j, int k, Mat3 &data) {
        double tmp_1 = val(i - 1, j, k, data) - 2 * data(i, j, k) + val(i + 1, j, k, data);
        tmp_1 /= h * h;

        double tmp_2 = val(i, j - 1, k, data) - 2 * data(i, j, k) + val(i, j + 1, k, data);
        tmp_2 /= h * h;

        double tmp_3 = val(i, j, k - 1, data) - 2 * data(i, j, k) + val(i, j, k + 1, data);
        tmp_3 /= h * h;

        return tmp_1 + tmp_2 + tmp_3;
    }
    
    void print(Mat3 &data) {
        std::cout << std::fixed << std::setprecision(2) << std::setfill('0');
            for (int j = 0; j < bsize[1]; j++) {
            std::cout << "index j: " << j << " -------------------------" << std::endl;
        for (int i = 0; i < bsize[0]; i++) {
                    std::cout << "i: " << i << " ";
                for (int k = 0; k < bsize[2]; k++) {
                    std::cout  << "\t"<< data(i, j, k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

    }

    void calc_zero(int ind) {
#pragma omp parallel for collapse(3) num_threads(THREADS)
        for (int i = 0; i < bsize[0]; i++) {
            for (int j = 0; j < bsize[1]; j++) {
                for (int k = 0; k < bsize[2]; k++) {
                    (*data[ind])(i, j, k) = an_sol(i, j, k, 0);
                }
            }
        }
    }

    double calc_error(int ind, int ptime) {
        double error = 0;
#pragma omp parallel for num_threads(THREADS) reduction(max:error)
        for (int i = start[0]; i < bsize[0] - end[0]; i++) {
            for (int j = start[1]; j < bsize[1] - end[1]; j++) {
                for (int k = start[2]; k < bsize[2] - end[2]; k++) {
                    double tmp = (*data[ind])(i, j, k) - an_sol(i, j, k, ptime);
                    tmp = std::abs(tmp);
                    if (tmp > error) error = tmp;
                }
            }
        }
        return error;
    }
    void calc_first(int ind_0, int ind_1) {
#pragma omp parallel for collapse(3) num_threads(THREADS)
        for (int i = start[0]; i < bsize[0] - end[0]; i++) {
            for (int j = start[1]; j < bsize[1] - end[1]; j++) {
                for (int k = start[2]; k < bsize[2] - end[2]; k++) {
                    double lap = func_lap(i, j, k, *data[ind_0]);
                    lap *= tau*tau / 2.0;
                    (*data[ind_1])(i, j, k) = (*data[ind_0])(i, j, k) + lap;
                }
            }
        }
    }

    void calc_next() {
#pragma omp parallel for collapse(3) num_threads(THREADS)
        for (int i = start[0]; i < bsize[0] - end[0]; i++) {
            for (int j = start[1]; j < bsize[1] - end[1]; j++) {
                for (int k = start[2]; k < bsize[2] - end[2]; k++) {
                    double lap = func_lap(i, j, k, *data[1]);
                    (*data[2])(i, j, k) = tau*tau*lap + 2 * (*data[1])(i, j, k) - 
                                          (*data[0])(i, j, k);
                }
            }
        } 
    }

    void change_buffs() {
        //0->0; 1->1; 2->2;
        std::swap(data[1], data[2]);
        //0->0; 1->2; 2->1;
        std::swap(data[0], data[2]);
        //0->1; 1->2; 2->0;
    }

    ~Block() {
        for (int i = 0; i < 3; i++) {
            delete this->data[i];
        }
        for (int i = 0; i < 6; i++) {
            delete this->facets[i];
            delete this->buffs[i];
        }
    }
};

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cout << "Wrong arguments: ";
        std::cout << "<var 0 or 1> <N> <Stamps>" << std::endl;
        return 0;
    }
    
    int var = std::atoi(argv[1]);
    int N = std::atoi(argv[2]);
    int stamps = std::atoi(argv[3]);
    int size, rank;
    int dim[3] = {0,0, 0};
    double L[3];
    for (int i = 0; i < 3; i++) {
        if (var == 0) L[i] = 1.0;
        else L[i] = M_PI;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Dims_create(size, 3, dim);
    int period[3] = {0, 0, 0};
    int reorder = 0;
    MPI_Comm grid;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &grid);
    
    double start_time = MPI_Wtime();
    Block block(rank, size, grid, N, dim, L);
    block.calc_zero(0);
    block.update_facets(0);
    block.exchange_facets();
    double error = block.calc_error(0, 0);
    double global_error = 0;
    MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, grid);
    if (!rank) {
        std::cout << "\t err: " << global_error  << std::endl;
	}

    
    block.calc_first(0, 1);
    block.update_facets(1);
    block.exchange_facets();
    error = block.calc_error(1, 1);
    global_error = 0;
    MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, grid);
	if (!rank) {
        std::cout << "\t err: " << global_error << std::endl;
	}

    //block.print(*block.get_data(1));
    for (int i = 0; i < stamps; i++) {
        block.calc_next();
        block.update_facets(2);
        block.exchange_facets();
        block.change_buffs();
        double error = block.calc_error(1, i + 2);
        double global_error = 0;
        MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, grid);
	    if (!rank) {
            std::cout << "it: " << i << "\t err: " << global_error << std::endl;
	    }
    }
    double el_time = MPI_Wtime() - start_time;
    if (!rank) {
        std::cout << "Time       elapsed:  " << el_time << "sec." << std::endl;
    }
    MPI_Finalize();
} 
