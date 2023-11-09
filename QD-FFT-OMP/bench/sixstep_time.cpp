#include <omp.h>

#include <chrono>
#include <iostream>

#include "Timer.hpp"
#include "fft.hpp"
#include "qd.hpp"

int main(int argc, char *argv[]) {
    if (argc < 3)
        return 1;

    qd base_cos_table[1 << atoi(argv[1])];
    qd base_sin_table[1 << atoi(argv[1])];
    make_cos_table(1 << atoi(argv[1]), base_cos_table);
    make_sin_table(1 << atoi(argv[1]), base_sin_table, base_cos_table);
    
   
    for(int p = 1;p <= atoi(argv[1]);p++) {
        int n = 1 << p;
        
        qd cos_table[n];
        qd sin_table[n];
        for(int i = 0;i < n;i++) {
            copy(base_cos_table[i * (1 << atoi(argv[1])) / n], cos_table[i]);
            copy(base_sin_table[i * (1 << atoi(argv[1])) / n], sin_table[i]);
        }
        qd x[n];
        qd ix[n];
        qd y[n];
        qd iy[n];
        for (int i = 0; i < n; i++) {
            rand(x[i]);
            rand(ix[i]);
        }

        for (int argi = 2;argi < argc;argi++) {
            Timer timer;
            omp_set_num_threads(atoi(argv[argi]));
            for (int i = 0; i < 10; i++) {
                timer.start();
                sixstep_fft(n, p, x, ix, y, iy, cos_table, sin_table);
                timer.stop();
            }
            std::cout << argv[argi] << "," <<  n << "," << timer.calc_ave_microsec() << std::endl;
        }
    }
    return 0;
}