#ifndef HMR_ALGORITHM_H
#define HMR_ALGORITHM_H

#include <cstdint>
#include <vector>

template <typename T>
inline bool hmr_in_ordered_vector(const T& x, const std::vector<T>& x_vector)
{
    //Binary search inside vector.
    int32_t left = 0, right = static_cast<int32_t>(x_vector.size() - 1);
    while (left <= right)
    {
        int32_t mid = left + (right - left) / 2;
        if (x_vector[mid] == x)
        {
            return true;
        }
        else if (x_vector[mid] < x)
        {
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
        }
    }
    return false;
}

#endif // HMR_ALGORITHM_H
