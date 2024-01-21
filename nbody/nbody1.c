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

//
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
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;

      p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->y[i] = (f32)rand() / (f32)RAND_MAX;
      p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->vx[i] = (f32)rand() / (f32)RAND_MAX;
      p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t *p, const f32 dt, u64 n)
{
  const f32 softening = 1e-20;

  for (u64 i = 0; i < n; i++)
  {
    f32 fx = 0.0, fy = 0.0, fz = 0.0;

    u64 j = 0;
    
    for (; j < n - 3; j += 4)
    {
      
      f32 dx1 = p->x[j] - p->x[i];
      f32 dy1 = p->y[j] - p->y[i];
      f32 dz1 = p->z[j] - p->z[i];
      f32 d_21 = (dx1 * dx1) + (dy1 * dy1) + (dz1 * dz1) + softening;
      f32 d_3_over_21 = pow(d_21, 3.0 / 2.0);
      fx += dx1 / d_3_over_21;
      fy += dy1 / d_3_over_21;
      fz += dz1 / d_3_over_21;

      
      f32 dx2 = p->x[j + 1] - p->x[i];
      f32 dy2 = p->y[j + 1] - p->y[i];
      f32 dz2 = p->z[j + 1] - p->z[i];
      f32 d_22 = (dx2 * dx2) + (dy2 * dy2) + (dz2 * dz2) + softening;
      f32 d_3_over_22 = pow(d_22, 3.0 / 2.0);
      fx += dx2 / d_3_over_22;
      fy += dy2 / d_3_over_22;
      fz += dz2 / d_3_over_22;

      
      f32 dx3 = p->x[j + 2] - p->x[i];
      f32 dy3 = p->y[j + 2] - p->y[i];
      f32 dz3 = p->z[j + 2] - p->z[i];
      f32 d_23 = (dx3 * dx3) + (dy3 * dy3) + (dz3 * dz3) + softening;
      f32 d_3_over_23 = pow(d_23, 3.0 / 2.0);
      fx += dx3 / d_3_over_23;
      fy += dy3 / d_3_over_23;
      fz += dz3 / d_3_over_23;

      
      f32 dx4 = p->x[j + 3] - p->x[i];
      f32 dy4 = p->y[j + 3] - p->y[i];
      f32 dz4 = p->z[j + 3] - p->z[i];
      f32 d_24 = (dx4 * dx4) + (dy4 * dy4) + (dz4 * dz4) + softening;
      f32 d_3_over_24 = pow(d_24, 3.0 / 2.0);
      fx += dx4 / d_3_over_24;
      fy += dy4 / d_3_over_24;
      fz += dz4 / d_3_over_24;
    }

    
    for (; j < n; j++)
    {
      f32 dx = p->x[j] - p->x[i];
      f32 dy = p->y[j] - p->y[i];
      f32 dz = p->z[j] - p->z[i];
      f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening;
      f32 d_3_over_2 = pow(d_2, 3.0 / 2.0);
      fx += dx / d_3_over_2;
      fy += dy / d_3_over_2;
      fz += dz / d_3_over_2;
    }

    p->vx[i] += dt * fx;
    p->vy[i] += dt * fy;
    p->vz[i] += dt * fz;
  }

  for (u64 i = 0; i < n; i++)
  {
    p->x[i] += dt * p->vx[i];
    p->y[i] += dt * p->vy[i];
    p->z[i] += dt * p->vz[i];
  }
}

//
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

