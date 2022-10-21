#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmath.h"

# define n 2
double abserr = 0.0001;
double relerr = 0.0001;

int diffs(int ndf, double t, double *x, double *dx)
{
  dx[0] = -44 * x[0] - 160 * x[1] + cos(t + 1);
  dx[1] = x[0] + atan(1 + t * t);
  return (0);
}

double adamsdx0(double t, double x1, double x2)
{
  return -44 * x1 - 160 * x2 + cos(t + 1);
}

double adamsdx1(double t, double x1, double x2)
{
  return x1 + atan(1 + t * t);
}

void adams(double h_adams, double hprint)
{
  double x[n], xp[n], arr1[4], arr2[4], t, tout, h;
  double z1[2000], z2[2000];
  int fail, nfe, flag = 1, maxnfe = 5000;

  rkfinit(2, &fail);
  if (fail)
  {
    printf("RKF initialization failed.\n");
    exit(-1);
  }

  t = 0;
  tout = h_adams;
  x[0] = 2.0;
  x[1] = 0.5;
  z1[0] = x[0];
  z2[0] = x[1];

  for (int i = 1; i <= 3; i++)
  {
    rkf45(diffs, n, x, xp, &t, tout, &relerr, abserr, &h, &nfe, maxnfe, &flag);
    z1[i] = x[0];
    z2[i] = x[1];
    t = tout;
    tout += h_adams;
  }
  rkfend();

  printf("---------------------------------------------------------------------------------\n");
  printf("|                             Adams, h(int) = %.3f                             |\n", h_adams);
  printf("---------------------------------------------------------------------------------\n");
  printf("|            t            |            z1            |            z2            |\n");
  printf("---------------------------------------------------------------------------------\n");
  for (int step = 4; step <= (1.6 / h_adams); step++)
  {
    tout = h_adams * step;
    for (int k = 0; k < 4; k++)
    {
      arr1[k] = adamsdx0(tout - h_adams * k, z1[step - (k + 1)], z2[step - (k + 1)]);
      arr2[k] = adamsdx1(tout - h_adams * k, z1[step - (k + 1)], z2[step - (k + 1)]);
    }
    z1[step] = z1[step - 1] + h_adams / 24 * (55 * arr1[0] - 59 * arr1[1] + 37 * arr1[2] - 9 * arr1[3]);
    z2[step] = z2[step - 1] + h_adams / 24 * (55 * arr2[0] - 59 * arr2[1] + 37 * arr2[2] - 9 * arr2[3]);
  }

  int m = hprint / h_adams;
  for (int i = 1; i <= (1.6 / hprint); i++)
  {
    t = hprint * i;
    printf("|%14.2f           |%18.6f        |%18.6f        |\n", t, z1[i * m], z2[i * m]);
  }
  printf("---------------------------------------------------------------------------------\n\n");
}

void rkfer(double hprint, int maxnfe)
{
  int fail, nfe, flag = 1;
  double x[n], xp[n], h, t, tout;
  x[0] = 2.0;
  x[1] = 0.5;
  rkfinit(2, &fail);
  if (fail)
  {
    printf("RKF initialization failed.\n");
    exit(-1);
  }
  printf("--------------------------------------------------------------------------------------------------------------\n");
  printf("|                                          RKF45, MAXNFE = %d                                              |\n", maxnfe);
  printf("--------------------------------------------------------------------------------------------------------------\n");
  printf("|            t            |            x1            |            x2            |            flag            |\n");
  printf("--------------------------------------------------------------------------------------------------------------\n");
  for (int step = 1; step <= 1.6 / hprint; step++)
  {
    tout = hprint * step;
    t = tout - hprint;
    rkf45(diffs, n, x, xp, &t, tout, &relerr, abserr, &h, &nfe, maxnfe, &flag);
    printf("|%14.2f           |%17.6f         |%17.6f         |%14d              |\n", t, x[0], x[1], flag);
  }
  rkfend();
  printf("--------------------------------------------------------------------------------------------------------------\n\n");
}

int main(void)
{
  int maxnfe = 5000;
  double h_adams = 0.008, h2_adams = 0.001, hprint = 0.08;

  // RKF45 part
  rkfer(hprint, maxnfe);

  // Adams part
  adams(h_adams, hprint);
  adams(h2_adams, hprint);
  adams(hprint, hprint);

  return (0);
}