/*********************************************************
*
*  Copyright (C) 2014 by Vitaliy Vitsentiy
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*********************************************************/

#pragma once

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <exception>
#include <future>
#include <mutex>
#include <queue>



// thread pool to run user's functors with signature
//      ret func(int id, other_params)
// where id is the index of the thread that runs the functor
// ret is some return type


namespace ctpl {


class Worker;

class ThreadPool {
public:
    using size_type = std::uint_fast64_t;
    using Task = std::function<void(size_type)>;
    using TaskQueue = std::queue<Task>;

    explicit ThreadPool(size_type nThreads = static_cast<size_type>(std::thread::hardware_concurrency())) : nThreads(nThreads) {
        threads.resize(nThreads);
        flags.resize(nThreads);
        for(size_type i = 0; i < nThreads; ++i) {
            initializeThread(i);
        }
    }

    ~ThreadPool() {
        wait();
    }

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool(ThreadPool &&) = delete;
    ThreadPool & operator=(const ThreadPool &) = delete;
    ThreadPool & operator=(ThreadPool &&) = delete;

    void wait() const {
        while (true) {
            if (!paused) {
                if (tasks_total == 0) {
                    break;
                }
            } else {
                if (get_tasks_running() == 0) {
                    break;
                }
            }
            std::this_thread::yield();
        }
    }

    [[nodiscard]] TaskQueue::size_type get_tasks_running() const {
        return tasks_total - get_tasks_queued();
    }

    [[nodiscard]] TaskQueue::size_type get_tasks_queued() const {
        std::scoped_lock lock {queueMutex};
        return tasks.size();
    }

    [[nodiscard]] TaskQueue::size_type get_tasks_total() const {
        return tasks_total;
    }

    // get the number of running threads in the pool
    size_type size() {
        return static_cast<size_type>(threads.size());
    }

    // number of idle threads
    int n_idle() { return this->nWaiting; }
    std::jthread & get_thread(int i) { return *this->threads[i]; }

    // empty the queue
    void clearQueue() {
        std::scoped_lock lock {queueMutex};
        tasks = {};
    }

    // wait for all computing threads to finish and stop all threads
    // may be called asynchronously to not pause the calling thread while waiting
    // if isWait == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
    void stop(bool isWait = false) {
        if (!isWait) {
            if (this->isStop)
                return;
            this->isStop = true;
            for (size_type i = 0, n = this->size(); i < n; ++i) {
                *this->flags[i] = true;  // command the threads to stop
            }
            this->clearQueue();  // empty the queue
        }
        else {
            if (this->isDone || this->isStop)
                return;
            this->isDone = true;  // give the waiting threads a command to finish
        }
        {
            std::unique_lock<std::mutex> lock(this->mutex);
            this->cv.notify_all();  // stop all waiting threads
        }
        for (auto & thread : this->threads) {  // wait for the computing threads to finish
            if (thread->joinable())
                thread->join();
        }
        // if there were no threads in the pool but some functors in the queue, the functors are not deleted by the threads
        // therefore delete them here
        this->clearQueue();
        this->threads.clear();
        this->flags.clear();
    }

    template <typename F, typename... A, typename R = std::invoke_result_t<std::decay_t<F>, size_type, std::decay_t<A>...>>
    std::future<R> push(F &&f, A&&... args) {
        ++tasks_total;
        auto task = std::make_shared<std::packaged_task<R(size_type)>>([&f, &args...](size_type threadIndex) {
            return f(threadIndex, std::forward<A>(args)...);
        });
        {
            std::scoped_lock lock {queueMutex};
            tasks.push([task](size_type id) { (*task)(id);} );
        }
        {
            std::unique_lock<std::mutex> lock(this->mutex);
            this->cv.notify_one();
        }
        return task->get_future();
    }

private:

    void initializeThread(size_type i) {
        auto worker = [this, i, flag = flags.at(i)]() {
            while (true) {
                Task task;
                auto tasksAvailable = pop_task(task);
                while (tasksAvailable) {  // if there is anything in the queue
                    task(i);
                    if (flag) {
                        return;  // the thread is wanted to stop, return even if the queue is not empty yet
                    } else {
                        tasksAvailable = pop_task(task);
                    }
                }
                // the queue is empty here, wait for the next command
                std::unique_lock<std::mutex> lock(this->mutex);
                ++this->nWaiting;
                this->cv.wait(lock, [this, &task, &tasksAvailable, &flag](){
                    tasksAvailable = pop_task(task);
                    return tasksAvailable || this->isDone || flag;
                });
                --this->nWaiting;
                if (!tasksAvailable) {
                    return;  // if the queue is empty and this->isDone == true or *flag then return
                }
            }
        };
        threads[i] = std::make_unique<std::jthread>(worker);
    }


    bool pop_task(Task &task)
    {
        std::scoped_lock lock {queueMutex};
        if (tasks.empty()) {
            return false;
        } else {
            task = std::move(tasks.front());
            tasks.pop();
            return true;
        }
    }

    void init() {
        this->nWaiting = 0;
        this->isStop = false;
        this->isDone = false;
    }

    mutable std::mutex queueMutex {};
    TaskQueue tasks = {};

    std::atomic<bool> isRunning = true;
    std::atomic<bool> isDone {false};
    std::atomic<bool> isStop {false};
    std::atomic<size_type> nWaiting {};  // how many threads are waiting
    size_type nThreads;

    std::vector<std::shared_ptr<std::atomic<bool>>> flags;
    std::vector<std::unique_ptr<std::jthread>> threads;

    std::atomic<size_type> nTasks;

    std::mutex mutex;
    std::condition_variable cv;

    std::atomic<bool> paused = false;
    std::atomic<size_type> tasks_total = 0;
};

}
