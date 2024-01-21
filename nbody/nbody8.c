#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

typedef struct particle_s {
    __m256 *x, *y, *z;     
    __m256 *vx, *vy, *vz;  
} particle_t;

static inline float randf(int use_sign) {
    float sign = (rand() > (RAND_MAX / 2)) ? 1.0f : -1.0f;
    return (use_sign ? sign : 1.0f) * ((float)rand() / (float)RAND_MAX);
}

void init(particle_t *p, u64 n) {
    
    size_t size = (n / 8) * sizeof(__m256); 
    p->x = aligned_alloc(32, size);
    p->y = aligned_alloc(32, size);
    p->z = aligned_alloc(32, size);
    p->vx = aligned_alloc(32, size);
    p->vy = aligned_alloc(32, size);
    p->vz = aligned_alloc(32, size);

    
    for (u64 i = 0; i < n / 8; i++) {
        float *xi = (float *)(p->x + i);
        float *yi = (float *)(p->y + i);
        float *zi = (float *)(p->z + i);
        float *vxi = (float *)(p->vx + i);
        float *vyi = (float *)(p->vy + i);
        float *vzi = (float *)(p->vz + i);

        for (int j = 0; j < 8; j++) {
            xi[j] = randf(1);
            yi[j] = randf(0);
            zi[j] = randf(1);
            vxi[j] = randf(0);
            vyi[j] = randf(1);
            vzi[j] = randf(0);
        }
    }
}


void move_particles(particle_t *p, const float dt, u64 n) {
    const float softening = 1e-20f;

    #pragma omp parallel for
    for (u64 i = 0; i < n / 8; i++) {  
        __m256 fx_sum = _mm256_setzero_ps(), fy_sum = _mm256_setzero_ps(), fz_sum = _mm256_setzero_ps();

        for (u64 j = 0; j < n / 8; j++) {
            __m256 x_i = p->x[i];
            __m256 y_i = p->y[i];
            __m256 z_i = p->z[i];

            __m256 dx = _mm256_sub_ps(p->x[j], x_i);
            __m256 dy = _mm256_sub_ps(p->y[j], y_i);
            __m256 dz = _mm256_sub_ps(p->z[j], z_i);

            __m256 d_2 = _mm256_fmadd_ps(dx, dx, _mm256_fmadd_ps(dy, dy, _mm256_fmadd_ps(dz, dz, _mm256_set1_ps(softening))));
            __m256 inv_d_3 = _mm256_rsqrt_ps(_mm256_mul_ps(d_2, _mm256_mul_ps(d_2, d_2)));

            fx_sum = _mm256_fmadd_ps(dx, inv_d_3, fx_sum);
            fy_sum = _mm256_fmadd_ps(dy, inv_d_3, fy_sum);
            fz_sum = _mm256_fmadd_ps(dz, inv_d_3, fz_sum);
        }

        
        float fx[8], fy[8], fz[8];
        _mm256_store_ps(fx, fx_sum);
        _mm256_store_ps(fy, fy_sum);
        _mm256_store_ps(fz, fz_sum);

        for (int k = 0; k < 8; k++) {
            p->vx[i][k] += dt * fx[k];
            p->vy[i][k] += dt * fy[k];
            p->vz[i][k] += dt * fz[k];
        }
    }
}

f64 compute_delta(f64 *p_ref, f32 *p, u64 n) {
  f64 delta = 0.0;
  for (u64 i = 0; i < n; i++)
    delta += (p_ref[i] - p[i]);
  delta /= (f64)n;
  return delta;
}


int main(int argc, char **argv)
{
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 steps = 13;
  const f32 dt = 0.01;

  f64 rate = 0.0, drate = 0.0;
  const u64 warmup = 3;

  FILE *file_ref = fopen("particle_references/particle_positions.dat", "r");
  f64 *p_ref_x = malloc(n * sizeof(f64));
  f64 *p_ref_y = malloc(n * sizeof(f64));
  f64 *p_ref_z = malloc(n * sizeof(f64));

  if (file_ref) {
    for (u64 i = 0; i < n; i++) {
        fscanf(file_ref, "%lf %lf %lf", &p_ref_x[i], &p_ref_y[i], &p_ref_z[i]);
    }
    fclose(file_ref);
  }

  
  particle_t p;
  init(&p, n);

  const u64 s = sizeof(f32) * n * 6; 
  printf("\nTotal memory size: %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);

  printf("%5s %10s %10s %8s\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  for (u64 i = 0; i < steps; i++)
    {
      const f64 start = omp_get_wtime();

      move_particles(&p, dt, n);

      const f64 end = omp_get_wtime();
      const f32 h1 = (f32)(n) * (f32)(n);
      const f32 h2 = (17.0 * h1 + 6.0 * (f32)n + 6.0 * (f32)n) * 1e-9;

      if (i >= warmup)
      {
        rate += h2 / (f32)(end - start);
        drate += (h2 * h2) / (f32)((end - start) * (end - start));
      }

      printf("%5llu %10.3e %10.3e %8.1f %s\n",
             i,
             (end - start),
             h1 / (end - start),
             h2 / (end - start),
             (i < warmup) ? "(warm up)" : "");
      
      fflush(stdout);
    }

  
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

  f64 delta_x = compute_delta(p_ref_x, p.x, n);
  f64 delta_y = compute_delta(p_ref_y, p.y, n);
  f64 delta_z = compute_delta(p_ref_z, p.z, n);

  printf("-----------------------------------------------------\n");
  printf("%s %4s %10.1lf +- %.1lf GFLOP/s\n", "Average performance:", "", rate, drate);
  printf("Delta X: %lf | Delta Y: %lf | Delta Z: %lf\n", delta_x, delta_y, delta_z);
  printf("-----------------------------------------------------\n");

  
  free(p.x);
  free(p.y);
  free(p.z);
  free(p.vx);
  free(p.vy);
  free(p.vz);

  return 0;
}
