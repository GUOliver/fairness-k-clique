/*
To execute:
clang++ -std=c++20 -march=native set_intersection.cpp -o set_intersection
*/

#include "set_intersection.hpp"

auto intersect(const int *set_a, int size_a, const int *set_b, int size_b, int *set_c) -> int {
    int i = 0, j = 0, size_c = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;   
}

auto intersect_count(const int *set_a, int size_a, const int *set_b, int size_b) -> int {
    auto i = 0, j = 0, res = 0;
    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return res;
}

auto intersect_simd4x(const int *set_a, int size_a,const int *set_b, int size_b, int *set_c) -> int {
    auto i = 0, j = 0, size_c = 0;
    auto qs_a = size_a - (size_a & 3);
    auto qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
            //_mm_prefetch((char*) (set_b + j), _MM_HINT_T0);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
            //_mm_prefetch((char*) (set_b + j), _MM_HINT_T0);
        }

        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3));

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        __m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
        _mm_storeu_si128((__m128i*)(set_c + size_c), p);

        size_c += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            set_c[size_c++] = set_a[i];
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }

    return size_c;    
}

auto intersect_simd4x_count(const int* set_a, int size_a, const int* set_b, int size_b) -> int {
    auto i = 0, j = 0, res = 0;
    auto qs_a = size_a - (size_a & 3);
    auto qs_b = size_b - (size_b & 3);

    while (i < qs_a && j < qs_b) {
        __m128i v_a = _mm_lddqu_si128((__m128i*)(set_a + i));
        __m128i v_b = _mm_lddqu_si128((__m128i*)(set_b + j));

        int a_max = set_a[i + 3];
        int b_max = set_b[j + 3];
        // i += (a_max <= b_max) * 4;
        // j += (b_max <= a_max) * 4;
        if (a_max == b_max) {
            i += 4;
            j += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
            //_mm_prefetch((char*) (set_b + j), _MM_HINT_T0);
        } else if (a_max < b_max) {
            i += 4;
            _mm_prefetch((char*) (set_a + i), _MM_HINT_NTA);
        } else {
            j += 4;
            _mm_prefetch((char*) (set_b + j), _MM_HINT_NTA);
            //_mm_prefetch((char*) (set_b + j), _MM_HINT_T0);
        }

        __m128i cmp_mask0 = _mm_cmpeq_epi32(v_a, v_b); // pairwise comparison
        __m128i rot1 = _mm_shuffle_epi32(v_b, cyclic_shift1);   // shuffling
        __m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, cyclic_shift2);
        __m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, cyclic_shift3);
        __m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, rot3);
        __m128i cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask0, cmp_mask1),
                _mm_or_si128(cmp_mask2, cmp_mask3));

        int mask = _mm_movemask_ps((__m128)cmp_mask);
        res += _mm_popcnt_u32(mask);
    }

    while (i < size_a && j < size_b) {
        if (set_a[i] == set_b[j]) {
            res++;
            i++; j++;
        } else if (set_a[i] < set_b[j]) {
            i++;
        } else {
            j++;
        }
    }
    return res;
}


auto main() -> int {
    int x[] = {1,2,3,4};
    int y[] = {3,4,5,6};
    int res = intersect_simd4x_count(x, 4, y, 4);

    int z[] = {};
    int res2 = intersect_simd4x(x, 4, y, 4, z);
    printf("res2: %d\n", res2);
    printf("z[0]: %d\n", z[0]);
    printf("z[1]: %d\n", z[1]);
    printf("z[2]: %d\n", z[2]);
    return 0;
}
