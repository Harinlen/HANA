#ifndef HMR_GLOBAL_H
#define HMR_GLOBAL_H

#include <unordered_set>

typedef struct THREAD_BLOCK
{
    int idx;
    int total;
} THREAD_BLOCK;

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

#endif // HMR_GLOBAL_H
