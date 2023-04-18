#ifndef HMR_ALGORITHM_H
#define HMR_ALGORITHM_H

#include <cstdint>
#include <cstdlib>
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

template <typename T>
inline void hmr_remove_common(std::vector<T> &vec_a, std::vector<T> &vec_b)
{
    //Loop and wait until all goes to the beginning.
    size_t pa = 0, pb = 0, a_len = vec_a.size(), b_len = vec_b.size(),
            a_off = 0, b_off = 0;
    //Goes for all vector a and b.
    while(pa < a_len && pb < b_len)
    {
        //Check whether two vectors are the same.
        if(vec_a[pa] == vec_b[pb])
        {
            //Increase the offset.
            ++a_off; ++b_off;
            ++pa; ++pb;
            continue;
        }
        if(vec_a[pa] < vec_b[pb])
        {
            //Increase a to get closer to b.
            if(a_off > 0)
            {
                vec_a[pa - a_off] = vec_a[pa];
            }
            ++pa;
        }
        else
        {
            //Increase b to get closer to a.
            if(b_off > 0)
            {
                vec_b[pb - b_off] = vec_b[pb];
            }
            ++pb;
        }
    }
    //Check the reset of vector.
    if(a_off > 0)
    {
        //Move the reset of vector when necessary.
        if(pa < a_len)
        {
            while(pa < a_len)
            {
                vec_a[pa - a_off] = vec_a[pa];
                ++pa;
            }
        }
        //Shorten vector a.
        vec_a.resize(a_len - a_off);
    }
    if(b_off > 0)
    {
        if(pb < b_len)
        {
            while(pb < b_len)
            {
                vec_b[pb - b_off] = vec_b[pb];
                ++pb;
            }
        }
        vec_b.resize(b_len - b_off);
    }
}

template <typename T>
inline void hmr_remove_one(std::vector<T> &vec, T value)
{
    for(size_t i=0; i<vec.size(); ++i)
    {
        if(vec[i] == value)
        {
            for(size_t j=i+1; j<vec.size(); ++j)
            {
                vec[j-1] = vec[j];
            }
            vec.resize(vec.size()-1);
            return;
        }
    }
}

#endif // HMR_ALGORITHM_H
