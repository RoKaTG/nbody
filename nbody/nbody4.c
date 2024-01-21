#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;


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
        f32 fx = 0.0, fy = 0.0, fz = 0.0;
        
        #pragma omp parallel for
        for (u64 j = 0; j < n; j++)
        {
            f32 dx, dy, dz, d_2, inv_d_3;

            
            __asm__ volatile (
                
                "movss (%[pxj]), %%xmm0\n\t"  // xmm0 = p->x[j]
                "movss (%[pxi]), %%xmm1\n\t"  // xmm1 = p->x[i]
                "subss %%xmm1, %%xmm0\n\t"  // xmm0 = xmm0 - xmm1 (dx)
                "movss %%xmm0, %[dx]\n\t"  // dx = xmm0

                "movss (%[pyj]), %%xmm2\n\t"  // xmm2 = p->y[j]
                "movss (%[pyi]), %%xmm3\n\t"  // xmm3 = p->y[i]
                "subss %%xmm3, %%xmm2\n\t"  // xmm2 = xmm2 - xmm3 (dy)
                "movss %%xmm2, %[dy]\n\t"  // dy = xmm2

                "movss (%[pzj]), %%xmm4\n\t"  // xmm4 = p->z[j]
                "movss (%[pzi]), %%xmm5\n\t"  // xmm5 = p->z[i]
                "subss %%xmm5, %%xmm4\n\t"  // xmm4 = xmm4 - xmm5 (dz)
                "movss %%xmm4, %[dz]\n\t"  // dz = xmm4

                
                "movss %[dx], %%xmm0\n\t"
                "movss %[dy], %%xmm1\n\t"
                "movss %[dz], %%xmm2\n\t"
                "mulss %%xmm0, %%xmm0\n\t"  // xmm0 = dx*dx
                "mulss %%xmm1, %%xmm1\n\t"  // xmm1 = dy*dy
                "mulss %%xmm2, %%xmm2\n\t"  // xmm2 = dz*dz
                "addss %%xmm1, %%xmm0\n\t"  // xmm0 = dx*dx + dy*dy
                "addss %%xmm2, %%xmm0\n\t"  // xmm0 = dx*dx + dy*dy + dz*dz
                "addss %[softening], %%xmm0\n\t"  // xmm0 += softening
                "movss %%xmm0, %[d2]\n\t"  // d_2 = xmm0

                
                "mulss %[d2], %%xmm0\n\t"  // xmm0 = d_2 * d_2
                "mulss %[d2], %%xmm0\n\t"  // xmm0 = d_2 * d_2 * d_2
                "sqrtss %%xmm0, %%xmm0\n\t"  // xmm0 = sqrt(d_2 * d_2 * d_2)
                "rcpss %%xmm0, %%xmm0\n\t"  // xmm0 = 1.0 / sqrt(d_2 * d_2 * d_2)
                "movss %%xmm0, %[invd3]\n\t"  // inv_d_3 = xmm0
                : [dx] "=m" (dx), [dy] "=m" (dy), [dz] "=m" (dz), [d2] "=m" (d_2), [invd3] "=m" (inv_d_3)
                : [pxj] "r" (p->x + j), [pxi] "r" (p->x + i), [pyj] "r" (p->y + j), [pyi] "r" (p->y + i), [pzj] "r" (p->z + j), [pzi] "r" (p->z + i), [softening] "m" (softening)
                : "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5"
            );

            fx += dx * inv_d_3;
            fy += dy * inv_d_3;
            fz += dz * inv_d_3;
        }

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

