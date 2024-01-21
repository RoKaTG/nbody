#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

void move_particles(particle_t *p, const f32 dt, u64 n) {
    const f32 softening = 1e-20;

    #pragma omp parallel for
    for (u64 i = 0; i < n; i++) {
        f32 fx = 0.0, fy = 0.0, fz = 0.0;

        for (u64 j = 0; j < n; j++) {
            float dx[8], dy[8], dz[8], d_2[8], inv_d_3[8];

            
            __asm__ volatile (
                
                "vbroadcastss (%[pxi]), %%ymm0\n\t"  // ymm0 = p->x[i]
                "vbroadcastss (%[pyi]), %%ymm1\n\t"  // ymm1 = p->y[i]
                "vbroadcastss (%[pzi]), %%ymm2\n\t"  // ymm2 = p->z[i]

                
                "vmovups (%[pxj]), %%ymm3\n\t"       // ymm3 = p->x[j]
                "vmovups (%[pyj]), %%ymm4\n\t"       // ymm4 = p->y[j]
                "vmovups (%[pzj]), %%ymm5\n\t"       // ymm5 = p->z[j]

                
                "vsubps %%ymm0, %%ymm3, %%ymm6\n\t"  // ymm6 = ymm3 - ymm0 (dx)
                "vsubps %%ymm1, %%ymm4, %%ymm7\n\t"  // ymm7 = ymm4 - ymm1 (dy)
                "vsubps %%ymm2, %%ymm5, %%ymm8\n\t"  // ymm8 = ymm5 - ymm2 (dz)

                // Stocker les résultats intermédiaires dans dx, dy, dz
                "vmovups %%ymm6, %[dx]\n\t"
                "vmovups %%ymm7, %[dy]\n\t"
                "vmovups %%ymm8, %[dz]\n\t"
                : [dx] "=m" (dx), [dy] "=m" (dy), [dz] "=m" (dz)
                : [pxi] "r" (p->x + i), [pyi] "r" (p->y + i), [pzi] "r" (p->z + i), 
                  [pxj] "r" (p->x + j), [pyj] "r" (p->y + j), [pzj] "r" (p->z + j)
                : "%ymm0", "%ymm1", "%ymm2", "%ymm3", "%ymm4", "%ymm5", "%ymm6", "%ymm7", "%ymm8"
            );

            
            __asm__ volatile (
                "vmovups %[dx], %%ymm6\n\t"  // ymm6 = dx
                "vmovups %[dy], %%ymm7\n\t"  // ymm7 = dy
                "vmovups %[dz], %%ymm8\n\t"  // ymm8 = dz

                "vmulps %%ymm6, %%ymm6, %%ymm6\n\t"  // ymm6 = dx*dx
                "vmulps %%ymm7, %%ymm7, %%ymm7\n\t"  // ymm7 = dy*dy
                "vmulps %%ymm8, %%ymm8, %%ymm8\n\t"  // ymm8 = dz*dz

                "vaddps %%ymm7, %%ymm6, %%ymm6\n\t"  // ymm6 = dx*dx + dy*dy
                "vaddps %%ymm8, %%ymm6, %%ymm6\n\t"  // ymm6 = dx*dx + dy*dy + dz*dz

                "vbroadcastss %[softening], %%ymm7\n\t"
                "vaddps %%ymm7, %%ymm6, %%ymm6\n\t"  // ymm6 = dx*dx + dy*dy + dz*dz + softening

                "vmovups %%ymm6, %[d2]\n\t"
                : [d2] "=m" (d_2)
                : [dx] "m" (dx), [dy] "m" (dy), [dz] "m" (dz), [softening] "m" (softening)
                : "%ymm6", "%ymm7", "%ymm8"
            );

            
            __asm__ volatile (
                "vmovups %[d2], %%ymm0\n\t"  // ymm0 = d_2
                "vmulps %%ymm0, %%ymm0, %%ymm0\n\t"  // ymm0 = d_2 * d_2
                "vmulps %%ymm0, %%ymm0, %%ymm0\n\t"  // ymm0 = d_2 * d_2 * d_2
                "vsqrtps %%ymm0, %%ymm0\n\t"  // ymm0 = sqrt(d_2 * d_2 * d_2)
                "vrcpps %%ymm0, %%ymm0\n\t"  // ymm0 = 1.0 / sqrt(d_2 * d_2 * d_2)
                "vmovups %%ymm0, %[invd3]\n\t"
                : [invd3] "=m" (inv_d_3)
                : [d2] "m" (d_2)
                : "%ymm0"
            );

            
            for (int k = 0; k < 8; ++k) {
                fx += dx[k] * inv_d_3[k];
                fy += dy[k] * inv_d_3[k];
                fz += dz[k] * inv_d_3[k];
            }
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

