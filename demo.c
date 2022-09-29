#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
int main(int argc, char** argv) {
    int* arr = malloc(sizeof(int) * 5);

    printf("%d %d %d %d %d\n", arr[0], arr[1], arr[2], arr[3], arr[4]);
    printf("%d\n", arr[100]);
    // size_t n = sizeof(*arr) / sizeof(arr[0]);
    int size = (&arr)[1] - arr;
    printf("Size of the arr: %d\n", size);

    int a[10];
    int* p = a;

    assert(sizeof(a) / sizeof(a[0]) == 10);
    printf("Size of int* (ptr to int): %lu\n", sizeof(int*)); // size of 8
    assert(sizeof(p) == sizeof(int*));
    assert(sizeof(*p) == sizeof(int));

    return 0;
}