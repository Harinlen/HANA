#ifndef HMR_PARALLEL_H
#define HMR_PARALLEL_H

#include <thread>
#include <cstdio>

#include "hmr_global.hpp"

typedef struct HMR_PFOR_BLOCK
{
    int32_t start, end;
} HMR_PFOR_BLOCK;

#define HMR_PFOR_FUNC(x, ...)   void x(int32_t idx, __VA_ARGS__)

template <class Function, class... Args>
void hmr_parallel_for(const int32_t threads, const int32_t range, Function &&task, Args&&... args)
{
    //Calculate the range for each thread.
    int32_t work_per_thread = (range + threads - 1) / threads;
    std::thread* workers = new std::thread[threads];
    for (int32_t i = 0; i < threads; ++i)
    {
        int32_t item_start = i * work_per_thread;
        int32_t item_end = hMin(item_start + work_per_thread, range);
        workers[i] = std::thread([=] {
            for (int32_t i = item_start; i < item_end; ++i)
            {
                task(i, args...);
            }
            });
    }
    for (int32_t i = 0; i < threads; ++i)
    {
        workers[i].join();
    }
    delete[] workers;
}

#endif // HMR_PARALLEL_H