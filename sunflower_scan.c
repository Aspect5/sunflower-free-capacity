/* Plain triple scan: stdin "N m1 ... mN"; prints count of distinct triples i<j<k with
   equal pairwise intersections (C&(A|B)) == (A&B) and (A&B)==(A&C) as a double check,
   plus the first witness. Written independently of check_c49.py. */
#include <stdio.h>
#include <stdlib.h>
int main(void) {
    int N; if (scanf("%d", &N) != 1) return 2;
    unsigned *m = malloc(sizeof(unsigned) * N);
    for (int i = 0; i < N; i++) if (scanf("%u", &m[i]) != 1) return 2;
    long long cnt = 0; int wi = -1, wj = -1, wk = -1;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            unsigned K = m[i] & m[j], U = m[i] | m[j];
            for (int k = j + 1; k < N; k++)
                if ((m[k] & U) == K) {
                    /* equal pairwise intersections; distinctness from indices + distinct input */
                    if ((m[i] & m[k]) == K && (m[j] & m[k]) == K) {
                        cnt++; if (wi < 0) { wi = i; wj = j; wk = k; }
                    }
                }
        }
    printf("triples=%lld", cnt);
    if (wi >= 0) printf(" witness=%u,%u,%u", m[wi], m[wj], m[wk]);
    printf("\n");
    return 0;
}
