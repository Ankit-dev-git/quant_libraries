#include "parallel_for.h"
#include <thread>
#include <vector>

void parallel_for(int start, int end, const std::function<void(int, int)>& worker)
{
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4;

    int total = end - start;
    int chunkSize = total / numThreads;
    int remainder = total % numThreads;

    std::vector<std::jthread> threads;
    int chunkStart = start;

    for (unsigned int t = 0; t < numThreads; ++t) {
        int chunkEnd = chunkStart + chunkSize + (t < static_cast<unsigned int>(remainder) ? 1 : 0);
        threads.emplace_back(worker, chunkStart, chunkEnd);
        chunkStart = chunkEnd;
    }

    threads.clear();
}
