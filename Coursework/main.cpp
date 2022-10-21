#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Forsythe.h"

const Float g = 9.81;
const Float M = 1;
Float ptout[50], px0[50], px[50];
Float Kglob, L = 0;
# define n 4

Float funcL(Float x)
{
  Float temp;
  if (x == 0.0) temp = 1.0;
  if (x != 0.0) temp = exp(x * x);

  return temp;
}

void diffs(Float t, Float* y, Float* dy)
{
  dy[0] = y[1];
  dy[1] = -Kglob * y[0] / M - g * (1 - cos(y[2])) + (L + y[0]) * y[3] * y[3];
  dy[2] = y[3];
  dy[3] = -g * sin(y[2]) / (L + y[0]) - 2 * y[1] * y[3] / (L + y[0]);
}

Float quadkrit(Float K)
{
  Kglob = K;
  Float x0[] = { 0.0, 0.0, 0.0, 4.0 };
  Float x[] = { 0, 0.303, -0.465, 0.592, -0.409, 0.164, 0.180 };
  Float krit = 0;
  unsigned char work[6 * (4 * sizeof(Float)) + sizeof(rkf_inside)];

  rkf rkfinit;
  rkfinit.f = diffs, rkfinit.Y = x0, rkfinit.t = 0, rkfinit.tout = 0;
  rkfinit.ae = 1e-9, rkfinit.re = 1e-9, rkfinit.neqn = n, rkfinit.flag = 1, rkfinit.work = work;
  rkf_inside* iwork = (rkf_inside*)rkfinit.work;

  for (int i = 0; rkfinit.tout <= 2.4; i++)
  {
    rkf45(&rkfinit);
    ptout[i] = rkfinit.tout;
    px0[i] = x0[0];
    px[i] = x[i];
    rkfinit.tout += 0.4;
    krit += ((x0[0] - x[i]) * (x0[0] - x[i]));
  }
  return krit;
}

void quadkritPrint()
{
  printf("---------------------------------------------------------------------------------\n");
  printf("|           tout          |           x(0)           |             x            |\n");
  printf("---------------------------------------------------------------------------------\n");
  Float k = 0;
  for (int i = 0; k <= 2.4; i++)
  {
    printf("| %16.6f        | %16.6f         | %15.3f          |\n", ptout[i], px0[i], px[i]);
    k += 0.4;
  }
  printf("---------------------------------------------------------------------------------\n\n");
}

int main()
{
  Float abserr = 1e-14, relerr = 1e-14, errest, flag;
  int nfe;
  Float result = Quanc8(funcL, 0.0, 1.0, abserr, relerr, &errest, &nfe, &flag);
  L = result * 0.6836896;

  std::cout << "Coursework, #15\n\n";
  std::cout << "Parameter L (the length of spring) = " << L << "\n\n";
  std::cout << "Starting FMin...\n";
  Float K = FMin(quadkrit, 36, 46, 1e-3);
  std::cout << "Parameter K (the stiffness of spring) = " << K << "\n\n";
  quadkritPrint();

  return 0;
}