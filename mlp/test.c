#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "rand.c"
// merge sort
inline void merge_key_value (int *a, double *v,  int n, int m) {
    int i, j, k;
    int *x = malloc(n * sizeof(int));
    for (i = 0, j = m, k = 0; k < n; k++) {
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : v[a[j]] < v[a[i]] ? a[j++]
             :               a[i++];
    }
    for (i = n; i--;) {
        a[i] = x[i];
    }
    free(x);
}

inline void merge_sort_key_value (int *a, double *v, int n) {
    if (n < 2)
        return;
    int m = n>>1;   // divide by 2
    merge_sort_key_value(a, v, m);
    merge_sort_key_value(a + m, v, n - m);
    merge_key_value(a, v, n, m);
} 


int main()
{
  double *v = malloc(10 * sizeof(double));
  int *a = malloc(10*sizeof(int));
  int i;
  initrand(1234);
  printf("\n");
  for(i = 0; i < 10; i++){
    a[i] = i;
    v[i] = U01();
    printf("%.4lf  ", v[i]);
  }
  printf("\n sorted array: \n");
  merge_sort_key_value(a, v, 10);
  for(i = 0; i < 10; i++){    
    printf("%d  ", a[i]);    
  }
  printf("\n");
  for(i = 0; i < 10; i++){    
    printf("%.4lf  ", v[a[i]]);
  }
  printf("\n");
  free(a); free(v);
  return 0;
}