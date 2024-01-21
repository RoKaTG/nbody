#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

  f32 *x, *y, *z;
  f32 *vx, *vy, *vz;
  
} particle_t;

static inline f32 randrealf(int use_sign) {
    f32 sign = (rand() > (RAND_MAX / 2)) ? 1.0f : -1.0f;
    return (use_sign ? sign : 1.0f) * ((f32)rand() / (f32)RAND_MAX);
}

void init(particle_t *p, u64 n)
{
  p->x = aligned_alloc(32, n * sizeof(f32));
  p->y = aligned_alloc(32, n * sizeof(f32));
  p->z = aligned_alloc(32, n * sizeof(f32));
  p->vx = aligned_alloc(32, n * sizeof(f32));
  p->vy = aligned_alloc(32, n * sizeof(f32));
  p->vz = aligned_alloc(32, n * sizeof(f32));

  for (u64 i = 0; i < n; i++)
  {
    p->x[i] = randrealf(1);
    p->y[i] = randrealf(0);
    p->z[i] = randrealf(1);
    p->vx[i] = randrealf(0);
    p->vy[i] = randrealf(1);
    p->vz[i] = randrealf(0);
  }
}

void move_particles(particle_t *p, const f32 dt, u64 n)
{
    const f32 softening = 1e-20;

    #pragma omp parallel for
    for (u64 i = 0; i < n; i++)
    {
        __m128 fx_v = _mm_setzero_ps();
        __m128 fy_v = _mm_setzero_ps();
        __m128 fz_v = _mm_setzero_ps();
        __m128 softening_v = _mm_set1_ps(softening);

        u64 j;
        for (j = 0; j <= n - 4; j += 4)
        {
            __m128 x_i = _mm_load1_ps(&p->x[i]);
            __m128 y_i = _mm_load1_ps(&p->y[i]);
            __m128 z_i = _mm_load1_ps(&p->z[i]);

            __m128 x_j = _mm_loadu_ps(&p->x[j]);
            __m128 y_j = _mm_loadu_ps(&p->y[j]);
            __m128 z_j = _mm_loadu_ps(&p->z[j]);

            __m128 dx = _mm_sub_ps(x_j, x_i);
            __m128 dy = _mm_sub_ps(y_j, y_i);
            __m128 dz = _mm_sub_ps(z_j, z_i);

            __m128 d2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
            d2 = _mm_add_ps(d2, softening_v);

            __m128 d2_3 = _mm_mul_ps(d2, _mm_mul_ps(d2, d2));
            __m128 inv_d3 = _mm_div_ps(_mm_set1_ps(1.0f), _mm_sqrt_ps(d2_3));

            fx_v = _mm_add_ps(fx_v, _mm_mul_ps(dx, inv_d3));
            fy_v = _mm_add_ps(fy_v, _mm_mul_ps(dy, inv_d3));
            fz_v = _mm_add_ps(fz_v, _mm_mul_ps(dz, inv_d3));
        }

        // Handle remaining particles
        for (; j < n; j++)
        {
            f32 dx = p->x[j] - p->x[i];
            f32 dy = p->y[j] - p->y[i];
            f32 dz = p->z[j] - p->z[i];

            f32 d2 = dx*dx + dy*dy + dz*dz + softening;
            f32 inv_d3 = 1.0f / sqrt(d2 * d2 * d2);

            fx_v = _mm_add_ss(fx_v, _mm_set_ss(dx * inv_d3));
            fy_v = _mm_add_ss(fy_v, _mm_set_ss(dy * inv_d3));
            fz_v = _mm_add_ss(fz_v, _mm_set_ss(dz * inv_d3));
        }

        f32 fx, fy, fz;
        _mm_store_ss(&fx, fx_v);
        _mm_store_ss(&fy, fy_v);
        _mm_store_ss(&fz, fz_v);

        p->vx[i] += dt * fx;
        p->vy[i] += dt * fy;
        p->vz[i] += dt * fz;
    }
}




f64 compute_delta(f64 *p_ref, f32 *p, u64 n) {
  f64 delta = 0.0;
  for (u64 i = 0; i < n; i++)
    delta += (p_ref[i] - p[i]);
  delta /= (f64)n;
  return delta;
}

//
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

  const u64 s = sizeof(f32) * n * 6; // 6: x, y, z, vx, vy, vz
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

