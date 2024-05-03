#include <iostream>
#include <vector>
#include <thread>

const int NUM_THREADS = 5;

void helloWorld(int thread_id) {
    std::cout << "Hello World! Thread ID: " << 
    thread_id << ", Total Threads: " << NUM_THREADS << std::endl;
}

int main() {
    std::vector<std::thread> threads(NUM_THREADS);

    for (int t = 0; t < NUM_THREADS; ++t) {
        threads[t] = std::thread(helloWorld, t);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return 0;
}