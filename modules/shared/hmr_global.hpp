#ifndef HMR_GLOBAL_H
#define HMR_GLOBAL_H

#include <deque>
#include <list>
#include <vector>
#include <unordered_set>

typedef struct THREAD_BLOCK
{
    int idx;
    int total;
} THREAD_BLOCK;

#define HMR_UNUSED(x) (void)x;

template <typename T>
constexpr inline const T &hMin(const T &a, const T &b) { return (a < b) ? a : b; }

template <typename T>
constexpr inline const T &hMax(const T &a, const T &b) { return (a > b) ? a : b; }

template <typename T>
constexpr inline const T hAbs(const T &x) { return (x < 0) ? (-1 * x) : x; }

template <typename T>
constexpr inline bool hInSet(const T &x, const std::unordered_set<T> &s) { return s.find(x) != s.end(); }

template <typename T>
constexpr inline const T hSquare(const T &x) { return x * x; };

template <typename T>
inline void hMoveListToVector(std::list<T>& t_list, std::vector<T>& t_vector)
{
    t_vector.reserve(t_list.size());
    while (!t_list.empty())
    {
        t_vector.emplace_back(t_list.front());
        t_list.pop_front();
    }
}

template <typename T>
inline void hDequeListToVector(std::deque<T>& t_list, std::vector<T>& t_vector)
{
    t_vector.reserve(t_list.size());
    while (!t_list.empty())
    {
        t_vector.emplace_back(t_list.front());
        t_list.pop_front();
    }
}

#endif // HMR_GLOBAL_H
