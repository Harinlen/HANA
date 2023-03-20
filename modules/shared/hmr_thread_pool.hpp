#ifndef HMR_THREAD_POOL_H
#define HMR_THREAD_POOL_H

#include <cstdio>

#include <thread>
#include <atomic>
#include <functional>
#include <condition_variable>

#include "hmr_ui.hpp"

namespace hmr
{
    typedef std::int_fast32_t   hmr_i32;
    typedef std::uint_fast32_t  hmr_ui32;
    typedef std::uint_fast64_t  hmr_ui64;

    template <typename T>
    class thread_pool_queue
    {
    public:
        thread_pool_queue() :
            m_finish(false),
            m_items(NULL),
            m_size(0),
            m_head(0),
            m_tail(0)
        {
        }

        ~thread_pool_queue()
        {
            delete[] m_items;
        }

        void initialize(hmr_i32 size)
        {
            m_items = new T[size];
            if (!m_items)
            {
                time_error(-1, "Failed to create thread pool queue for size %d", size);
            }
            m_size = size;
        }

        bool empty() const
        {
            return m_head == m_tail;
        }

        bool full() const
        {
            return next_pos(m_tail) == m_head;
        }

        bool pop(T& item)
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_popCv.wait(lock, [this] { return m_finish || (!empty()); });
            //Copy the data to item.
            if (!empty())
            {
                item = std::move(m_items[m_head]);
                m_head = next_pos(m_head);
                //Notify the push conditional variable.
                m_pushCv.notify_one();
                return true;
            }
            else
            {
                return false;
            }
        }

        void push(const T& item)
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_pushCv.wait(lock, [this] {return !full(); });
            //Copy the item data to the queue.
            m_items[m_tail] = std::move(item);
            m_tail = next_pos(m_tail);
            //Notify the pop conditional variable.
            m_popCv.notify_one();
        }

        void close()
        {
            m_finish = true;
            //Notify all the pop condition variable.
            m_popCv.notify_all();
        }

    private:
        inline const int32_t next_pos(int current_pos) const
        {
            return (current_pos == m_size - 1) ? 0 : (current_pos + 1);
        }
        bool m_finish;
        T* m_items;
        hmr_i32 m_size, m_head, m_tail;
        mutable std::mutex m_mutex;
        std::condition_variable m_pushCv, m_popCv;
    };

    template<typename T>
    class thread_pool
    {
    public:
        thread_pool(void (*task)(const T &), const hmr_i32 &queue_depth, const hmr_ui32 &thread_count = std::thread::hardware_concurrency()):
            m_task(task),
            m_threads(new std::thread[thread_count]),
            m_tasksTotal(0),
            m_threadCount(thread_count)
        {
            m_tasks.initialize(queue_depth);
            //Create threads.
            if (!m_threads)
            {
                time_error(-1, "Not enough memory to create thread array.");
            }
            for (hmr_ui32 i = 0; i < thread_count; i++)
            {
                m_threads[i] = std::thread(&thread_pool::worker, this);
            }
        }

        ~thread_pool()
        {
            //Wait for all the task finished.
            wait_for_tasks();
            //Close the thread pool.
            m_running = false;
            //Quit the queue.
            m_tasks.close();
            //Destory all the tasks.
            for (hmr_ui32 i = 0; i < m_threadCount; i++)
            {
                m_threads[i].join();
            }
        }

        void push_task(const T& task)
        {
            //Increase the task counter.
            {
                std::unique_lock<std::mutex> lock(m_taskMutex);
                ++m_tasksTotal;
            }
            //Push the task to queue.
            m_tasks.push(task);
        }

        void wait_for_tasks()
        {
            //Wait until the task total is finished.
            while (m_tasksTotal)
            {
                //Sleep for a while.
                std::this_thread::sleep_for(std::chrono::microseconds(1000));
            }
        }

    private:
        void worker()
        {
            while (m_running)
            {
                //Pop the task from the queue.
                T param;
                if (m_tasks.pop(param))
                {
                    m_task(param);
                    {
                        std::unique_lock<std::mutex> lock(m_taskMutex);
                        --m_tasksTotal;
                    }
                }
            }
        }

        thread_pool_queue<T> m_tasks;
        void (*m_task)(const T&);
        std::unique_ptr<std::thread[]> m_threads;
        std::atomic<bool> m_running{ true };
        std::mutex m_taskMutex;
        hmr_ui32 m_tasksTotal;
        hmr_ui32 m_threadCount;
    };
}

#endif // HMR_THREAD_POOL_H
