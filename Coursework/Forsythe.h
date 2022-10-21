//-------------------------------------------------------------------
//		ForsytheLib v 0.1 beta 5
//	����� �� �������������� ������, ����������� � Fortran'� �� �++
//
//	Copyright (c) 2002 Eugene S.B.
//
//	If you find any bugs in this program, please, contact me.
//	E-Mail: <esb@hotbox.ru>
//-------------------------------------------------------------------

#ifndef __FORSYTHE_H__
#define __FORSYTHE_H__

#include <float.h>

#define FT_FLOAT	0	//'float'
#define FT_DOUBLE	1	//'double'
#define FT_LDOUBLE	2	//'long double'

//����� ����� �������������� ������������ ���, ������� ����� �������������� � �����������.
//���� ����� �������� - ����������� 'float' (�.�. FT_FLOAT), ���� ����� �������� - ����������� 'double' (FT_DOUBLE)
//��� 'long double' (FT_LDOUBLE).
#define FLOATTYPE	FT_DOUBLE

#if FLOATTYPE == FT_FLOAT
	typedef float Float;
	#define EPSILON FLT_EPSILON
	#define MAXIMUM FLT_MAX
#elif FLOATTYPE == FT_DOUBLE
	typedef double Float;
	#define EPSILON DBL_EPSILON
	#define MAXIMUM DBL_MAX
#elif FLOATTYPE == FT_LDOUBLE
	typedef long double Float;
	#define EPSILON LDBL_EPSILON
	#define MAXIMUM LDBL_MAX
#else
	#error Invalid Float Type
#endif

//��������� �������� �������
#define ABS(x) (((x) >= 0) ? (x) : -(x))
#define MAX(x, y) (((x) >= (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define SIGN0(x) (((x) > 0) ? 1 : -1)
#define SIGN(x) (((x) == 0) ? 0 : SIGN0(x))
#define SIGN2(a, b) (SIGN(b)*ABS(a))

//��������� ��� ������� � �������� ���������� rkf45
struct rkf
{
	void (*f)(Float t, Float *Y, Float *dY);
	Float *Y, t, tout, re, ae;
	unsigned int neqn;
	int flag;
	unsigned char *work;//�� ����� ��������� ���������� �������� '6*neqn*sizeof(Float) + sizeof(struct rkf_inside)' ����
};

//��������� ��� ���������� ������� rkf45

//	REMin - ��� ����������� ���������� �������� ��� ������������� ����������� (re). ������� ��������� �
//������� ����� ������� �������� ������ ������������, � ������� ������ ����������.
#define REMIN (1.0e-12) //Relative error minimum

//	������������ ���������� �������� ����������� ���������� ����������� (�.�. ����������� ���������� ������� f)
//�������� ����� �������� ������������� �������� 5000 �����
#define MAXNFE 30000 //Maximum number of Function E...?

/* ��������� ��� ����������� ������������� � rkf45 */
//		E� �������� ���������� ����:
// h	-	�������������� ������ ���� ��� ���������� �����
// nfe	-	������� ����� ���������� �������
struct rkf_inside
{
	unsigned long nfe, kop, init;
	int JFlag, KFlag;
	Float h, SaveRE, SaveAE;
};
//		���� ��� ���������� �������� ����� ������ ��������� ����� ������ rkf45,
// �� �� ����� �������� ����� ��������� work � 'struct rkf', ��������:
/*
	struct rkf_inside *p;
	struct rkf RKF;
	unsigned char work[6*(NEQN*sizeof(Float)) + sizeof(struct rkf_inside)];
	...
	RKF.f = ...
	RKF.t = ...
	RKF.work = work;
	...
	rkf45(&RKF);
	p = (struct rkf_inside *)RKF.work;
	... = p->nfe;
	... = p->h;
	...
 */

//****************** ����� �����-����� ��������� 4-5 ������� ******************
//	�� ��������� � ������������ ������� ���� ��������� (�������������) ������� ��������� ���������:
//		- ��� ������������ ������� rkf45 �� �������� ������ ������ ����������� ���������� �������,
// ������� ����� �������� ������ 1. ��� �� ����� ��� ������������ ��������� ��������� � ��������� struct rkf.
//		- ������������� ��������� WORK � IWORK ���������� � 1 (work; ��. �������� struct rkf)
//		- ��� ���������� ������ ������� ���������� false. � ������ �� ��������� ������ ����� ���������� true.
//		!!!!!!  ������� ���������, ��� ���� ���� rkf45 ������� false, �� �� ����� �� �������, ��� �� ��������� -> �����
// ����������� ��������� �������� ����� �� ������, ������ ����� ����� ����� ���������, ��� �������������� ����������� �������.
// ������ �������� ������� ��� ����� ������, ��� ��� ������� ��������� ��������������� �������� (� ������ ������
// ��� ������������� ��������� ������ ����������� ��� ���������, �� ����� ������������ ����������� ���������� ������ ��������)

//Copyright (c) 1976 H.A.Watts, L.F.Shampine (Sandia Laboratories, Albuquerque, New Mexico)
bool rkf45(struct rkf *p);

//Copyright (c) 1973 Ricard Brent
Float FMin(Float (*F)(Float), Float a, Float b, Float tol);
Float Zeroin(Float (*F)(Float), Float a, Float b, Float tol);

Float Quanc8(Float (*F)(Float), Float a, Float b, Float ae, Float re, Float *errest, int *nofun, Float *flag);
void Decomp(unsigned int n, Float *A, Float *cond, int *ipvt);
void Solve(unsigned int n, Float *A, Float *b, int *ipvt);
void Spline(unsigned int n, Float *X, Float *Y, Float *B, Float *C, Float *D);
Float SEval(unsigned int n, Float u, Float *X, Float *Y, Float *B, Float *C, Float *D);

#endif //__FORSYTHE_H__
