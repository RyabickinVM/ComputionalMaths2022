#include <math.h>
#include <stdio.h>
#include "cmath.h"

// SPLINE values
#define ndim 11
double x[ndim], y[ndim], bspl[ndim], cspl[ndim], dspl[ndim], u[ndim], fspline[ndim];

// Lagrange values
double temp[ndim];

// Total comparison values
double quancArr[ndim], splineArr[ndim], lagrangeArr[ndim], splineDiff, lagrangeDiff;

double fquanc(double t) // integrand (fun)
{
  double temp;
  if (t == 0.0) temp = 0.0;
  if (t != 0.0) temp = (1 - cos(t)) / t;
  return (temp);
}

void quancer(double a, double b, double bx, double step) // QUANC8 function for interval
{
  double result, errest, posn;
  int nfe, flag;
  double epsrel = 1.0e-10; // relative epsilon
  double epsabs = 0.0; // absolute epsilon
  int i = 0;
  while(b < bx + step)
  {
    quanc8(fquanc, a, b, epsabs, epsrel, &result, &errest, &nfe, &posn, &flag);
    x[i] = b;
    y[i] = result;
    b += step;
    i++;
  }
}

void splinerXk(double n, double final, double step) // U array creation for SPLINE
{
  int i = 0;
  while (n < final + step)
  {
    u[i] = n;
    n += step;
    i++;
  }
}

void spliner() // SPLINE function for interval, uses QUANC8 calculations
{
  int last, flagspline;
  spline(ndim, 1, 1, 1, 1, x, y, bspl, cspl, dspl, &flagspline);
  for (int i = 0; i < ndim; ++i)
  {
    fspline[i] = seval(ndim, u[i], x, y, bspl, cspl, dspl, &last);
  }
}

double lagranger(double xlag) // Lagrange function for interval, uses QUANC8 calculations
{
  double resultlag = 0;
  double numer[ndim], denom[ndim];
  for (int i = 0; i < ndim; i++)
  {
    for (int j = 0; j < ndim; j++)
    {
      if (i != j)
      {
        numer[i] *= xlag - x[j];
        denom[i] *= x[i] - x[j];
      }
    }
    if (denom[i] != 0)
    {
      resultlag += (numer[i] / denom[i]) * y[i];;
    }
  }
  return resultlag;
}

void totalComp(double a, double b, double bx, double step, double bp, double bxp, double stepp) // the overall comparison
{
  printf("-----------------------------------------------------------------\n");
  printf("|                    %.2f <= x <= %.2f                          |\n", b, bx);
  printf("-----------------------------------------------------------------\n");
  printf("|               X               |             QUANC8            |\n");
  printf("-----------------------------------------------------------------\n");
  quancer(a, b, bx, step);
  for (int i = 0; i < ndim; i++)
  {
    quancArr[i] = y[i];
    printf("|%17.2f              |%22.9f         |\n", x[i], quancArr[i]);
  }
  splinerXk(2.05, 2.95, 0.1);
  spliner();
  for (int z = 0; z < ndim - 1; z++)
  {
    temp[z] = lagranger(u[z]);
  }
  quancer(0, bp, bxp, stepp);
  printf("-----------------------------------------------------------------\n");
  printf("|                    %.2f <= x <= %.2f                          |\n", 2.05, 2.95);
  printf("-----------------------------------------------------------------\n");
  printf("|    X     |     QUANC8     |     SPLINE     |     Lagrange     |\n");
  printf("-----------------------------------------------------------------\n");
  for (int i = 0; i < ndim - 1; i++)
  {
    quancArr[i] = y[i];
    splineArr[i] = fspline[i];
    lagrangeArr[i] = temp[i];
    printf("|%7.2f   |%13.9f   |%13.9f   |%14.9f    |\n", u[i], quancArr[i], splineArr[i], lagrangeArr[i]);
  }
  printf("-----------------------------------------------------------------\n");
  printf("|    X     |     QUANC8     |     Q8 - S     |     Q8 - Lgr     |\n");
  printf("-----------------------------------------------------------------\n");
  for (int i = 0; i < ndim - 1; i++)
  {
    splineDiff = fabs(quancArr[i] - splineArr[i]);
    lagrangeDiff = fabs(quancArr[i] - lagrangeArr[i]);
    printf("|%7.2f   |%13.9f   |%14.6f  |%15.6f   |\n", x[i], quancArr[i], splineDiff, lagrangeDiff);
  }
  printf("-----------------------------------------------------------------\n");
  printf("|                                                               |\n");
  printf("-----------------------------------------------------------------\n\n\n");
}

int main(void)
{
  totalComp(0, 2.0, 3.0, 0.1, 2.05, 2.95, 0.1);
  return (0);
}
