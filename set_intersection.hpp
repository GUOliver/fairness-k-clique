#include <cstdio>
#include <x86intrin.h>
#include <unistd.h>
#include <cstdint>
#include <sys/time.h>
#include "immintrin.h"

constexpr int cyclic_shift1 = _MM_SHUFFLE(0,3,2,1); // rotating right
constexpr int cyclic_shift2 = _MM_SHUFFLE(2,1,0,3); // rotating left
constexpr int cyclic_shift3 = _MM_SHUFFLE(1,0,3,2); // between

static const uint8_t shuffle_pi8_array[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 255, 255, 255, 255, 255, 255, 255, 255, 
    8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 255, 255, 255, 255, 
    12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 255, 255, 255, 255, 
    8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 255, 255, 255, 255, 
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 255, 255, 255, 255, 
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
};

static const __m128i *shuffle_mask = (__m128i*)(shuffle_pi8_array);

auto intersect(const int *set_a, int size_a, const int *set_b, int size_b, int *set_c) -> int;
auto intersect_count(const int *set_a, int size_a, const int *set_b, int size_b) -> int;
auto intersect_simd4x(const int *set_a, int size_a,const int *set_b, int size_b, int *set_c) -> int;
auto intersect_simd4x_count(const int* set_a, int size_a, const int* set_b, int size_b) -> int;
