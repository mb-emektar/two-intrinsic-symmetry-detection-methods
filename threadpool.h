#ifndef THREADPOOL_H
#define THREADPOOL_H

// Version macros.
#define TASK_THREAD_POOL_VERSION_MAJOR 1
#define TASK_THREAD_POOL_VERSION_MINOR 0
#define TASK_THREAD_POOL_VERSION_PATCH 10

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>

// MSVC does not correctly set the __cplusplus macro by default, so we must read it from _MSVC_LANG
// See https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
#if __cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)
#define TTP_CXX17 1
#else
#define TTP_CXX17 0
#endif

#if TTP_CXX17
#define TTP_NODISCARD [[nodiscard]]
#else
#define TTP_NODISCARD
#endif

namespace task_thread_pool {

#if !TTP_CXX17
    /**
     * A reimplementation of std::decay_t, which is only available since C++14.
     */
    template <class T>
    using decay_t = typename std::decay<T>::type;
#endif

    /**
     * A fast and lightweight thread pool that uses C++11 threads.
     */
    class task_thread_pool {
    public:
        /**
         * Create a task_thread_pool and start worker threads.
         *
         * @param num_threads Number of worker threads. If 0 then number of threads is equal to the
         *                    number of physical cores on the machine, as given by std::thread::hardware_concurrency().
         */
        explicit task_thread_pool(unsigned int num_threads = 0) {
            if (num_threads < 1) {
                num_threads = std::thread::hardware_concurrency();
                if (num_threads < 1) { num_threads = 1; }
            }
            start_threads(num_threads);
        }

        /**
         * Finish all tasks left in the queue then shut down worker threads.
         * If the pool is currently paused then it is resumed.
         */
        ~task_thread_pool() {
            unpause();
            wait_for_queued_tasks();
            stop_all_threads();
        }

        /**
         * Drop all tasks that have been submitted but not yet started by a worker.
         *
         * Tasks already in progress continue executing.
         */
        void clear_task_queue() {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            tasks = {};
        }

        /**
         * Get number of enqueued tasks.
         *
         * @return Number of tasks that have been enqueued but not yet started.
         */
        TTP_NODISCARD size_t get_num_queued_tasks() const {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            return tasks.size();
        }

        /**
         * Get number of in-progress tasks.
         *
         * @return Approximate number of tasks currently being processed by worker threads.
         */
        TTP_NODISCARD size_t get_num_running_tasks() const {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            return num_inflight_tasks;
        }

        /**
         * Get total number of tasks in the pool.
         *
         * @return Approximate number of tasks both enqueued and running.
         */
        TTP_NODISCARD size_t get_num_tasks() const {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            return tasks.size() + num_inflight_tasks;
        }

        /**
         * Get number of worker threads.
         *
         * @return Number of worker threads.
         */
        TTP_NODISCARD unsigned int get_num_threads() const {
            const std::lock_guard<std::recursive_mutex> threads_lock(thread_mutex);
            return static_cast<unsigned int>(threads.size());
        }

        /**
         * Set number of worker threads. Will start or stop worker threads as necessary.
         *
         * @param num_threads Number of worker threads. If 0 then number of threads is equal to the
         *                    number of physical cores on the machine, as given by std::thread::hardware_concurrency().
         * @return Previous number of worker threads.
         */
        unsigned int set_num_threads(unsigned int num_threads) {
            const std::lock_guard<std::recursive_mutex> threads_lock(thread_mutex);
            unsigned int previous_num_threads = get_num_threads();

            if (num_threads < 1) {
                num_threads = std::thread::hardware_concurrency();
                if (num_threads < 1) { num_threads = 1; }
            }

            if (previous_num_threads <= num_threads) {
                // expanding the thread pool
                start_threads(num_threads - previous_num_threads);
            } else {
                // contracting the thread pool
                stop_all_threads();
                {
                    const std::lock_guard<std::mutex> tasks_lock(task_mutex);
                    pool_running = true;
                }
                start_threads(num_threads);
            }

            return previous_num_threads;
        }

        /**
         * Stop executing queued tasks. Use `unpause()` to resume. Note: Destroying the pool will implicitly unpause.
         *
         * Any in-progress tasks continue executing.
         */
        void pause() {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            pool_paused = true;
        }

        /**
         * Resume executing queued tasks.
         */
        void unpause() {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            pool_paused = false;
            task_cv.notify_all();
        }

        /**
         * Check whether the pool is paused.
         *
         * @return true if pause() has been called without an intervening unpause().
         */
        TTP_NODISCARD bool is_paused() const {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            return pool_paused;
        }

        /**
         * Submit a Callable for the pool to execute and return a std::future.
         *
         * @param func The Callable to execute. Can be a function, a lambda, std::packaged_task, std::function, etc.
         * @param args Arguments for func. Optional.
         * @return std::future that can be used to get func's return value or thrown exception.
         */
        template <typename F, typename... A,
#if TTP_CXX17
            typename R = std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>
#else
            typename R = typename std::result_of<decay_t<F>(decay_t<A>...)>::type
#endif
            >
        TTP_NODISCARD std::future<R> submit(F&& func, A&&... args) {
#if defined(_MSC_VER)
            // MSVC's packaged_task is not movable even though it should be.
            // Discussion about this bug and its future fix:
            // https://developercommunity.visualstudio.com/t/unable-to-move-stdpackaged-task-into-any-stl-conta/108672
            std::shared_ptr<std::packaged_task<R()>> ptask =
                std::make_shared<std::packaged_task<R()>>(std::bind(std::forward<F>(func), std::forward<A>(args)...));
            submit_detach([ptask] { (*ptask)(); });
            return ptask->get_future();
#else
            std::packaged_task<R()> task(std::bind(std::forward<F>(func), std::forward<A>(args)...));
            auto ret = task.get_future();
            submit_detach(std::move(task));
            return ret;
#endif
        }

        /**
         * Submit a zero-argument Callable for the pool to execute.
         *
         * @param func The Callable to execute. Can be a function, a lambda, std::packaged_task, std::function, etc.
         */
        template <typename F>
        void submit_detach(F&& func) {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            tasks.emplace(std::forward<F>(func));
            task_cv.notify_one();
        }

        /**
         * Submit a Callable with arguments for the pool to execute.
         *
         * @param func The Callable to execute. Can be a function, a lambda, std::packaged_task, std::function, etc.
         */
        template <typename F, typename... A>
        void submit_detach(F&& func, A&&... args) {
            const std::lock_guard<std::mutex> tasks_lock(task_mutex);
            tasks.emplace(std::bind(std::forward<F>(func), std::forward<A>(args)...));
            task_cv.notify_one();
        }

        /**
         * Block until the task queue is empty. Some tasks may be in-progress when this method returns.
         */
        void wait_for_queued_tasks() {
            std::unique_lock<std::mutex> tasks_lock(task_mutex);
            notify_task_finish = true;
            task_finished_cv.wait(tasks_lock, [&] { return tasks.empty(); });
            notify_task_finish = false;
        }

        /**
         * Block until all tasks have finished.
         */
        void wait_for_tasks() {
            std::unique_lock<std::mutex> tasks_lock(task_mutex);
            notify_task_finish = true;
            task_finished_cv.wait(tasks_lock, [&] { return tasks.empty() && num_inflight_tasks == 0; });
            notify_task_finish = false;
        }

    protected:

        /**
         * Main function for worker threads.
         */
        void worker_main() {
            bool finished_task = false;

            while (true) {
                std::unique_lock<std::mutex> tasks_lock(task_mutex);

                if (finished_task) {
                    --num_inflight_tasks;
                    if (notify_task_finish) {
                        task_finished_cv.notify_all();
                    }
                }

                task_cv.wait(tasks_lock, [&]() { return !pool_running || (!pool_paused && !tasks.empty()); });

                if (!pool_running) {
                    break;
                }

                // Must mean that (!pool_paused && !tasks.empty()) is true

                std::packaged_task<void()> task{std::move(tasks.front())};
                tasks.pop();
                ++num_inflight_tasks;
                tasks_lock.unlock();

                try {
                    task();
                } catch (...) {
                    // std::packaged_task::operator() may throw in some error conditions, such as if the task
                    // had already been run. Nothing that the pool can do anything about.
                }

                finished_task = true;
            }
        }

        /**
         * Start worker threads.
         *
         * @param num_threads How many threads to start.
         */
        void start_threads(const unsigned int num_threads) {
            const std::lock_guard<std::recursive_mutex> threads_lock(thread_mutex);

            for (unsigned int i = 0; i < num_threads; ++i) {
                threads.emplace_back(&task_thread_pool::worker_main, this);
            }
        }

        /**
         * Stop, join, and destroy all worker threads.
         */
        void stop_all_threads() {
            const std::lock_guard<std::recursive_mutex> threads_lock(thread_mutex);

            {
                const std::lock_guard<std::mutex> tasks_lock(task_mutex);
                pool_running = false;
                task_cv.notify_all();
            }

            for (auto& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
            threads.clear();
        }

        /**
         * The worker threads.
         *
         * Access protected by thread_mutex
         */
        std::vector<std::thread> threads;

        /**
         * A mutex for methods that start/stop threads.
         */
        mutable std::recursive_mutex thread_mutex;

        /**
         * The task queue.
         *
         * Access protected by task_mutex.
         */
        std::queue<std::packaged_task<void()>> tasks = {};

        /**
         * A mutex for all variables related to tasks.
         */
        mutable std::mutex task_mutex;

        /**
         * Used to notify changes to the task queue, such as a new task added, pause/unpause, etc.
         */
        std::condition_variable task_cv;

        /**
         * Used to notify of finished tasks.
         */
        std::condition_variable task_finished_cv;

        /**
         * A signal for worker threads that the pool is either running or shutting down.
         *
         * Access protected by task_mutex.
         */
        bool pool_running = true;

        /**
         * A signal for worker threads to not pull new tasks from the queue.
         *
         * Access protected by task_mutex.
         */
        bool pool_paused = false;

        /**
         * A signal for worker threads that they should notify task_finished_cv when they finish a task.
         *
         * Access protected by task_mutex.
         */
        bool notify_task_finish = false;

        /**
         * A counter of the number of tasks in-progress by worker threads.
         * Incremented when a task is popped off the task queue and decremented when that task is complete.
         *
         * Access protected by task_mutex.
         */
        int num_inflight_tasks = 0;
    };
}

// clean up
#undef TTP_NODISCARD
#undef TTP_CXX17


#endif // COLOR_H


