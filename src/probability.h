

#include "stdio.h"

//////////////////////////////////////////////////////////////////////////
//���ʼ��㺯���ļ�

// ٤�꺯��
extern double gam(int n);

//  ��̬�ֲ��ķֲ�����ֵ: p(-��,u)
extern double gaos(double u);

//  ��̬�ֲ��ķ�����, p(-��,u) ; ��֪p, ����u
extern double re_gaos(double q);

//  chi2�ֲ��ķֲ�������p(0,x)
extern double chi2(int n,double x,double *f);

//  chi2�ֲ��ķ�������p=chi2(0,x),��֪x������p
extern double re_chi2(int n,double p);

//B�ֲ�����ֵ��F�ֲ���t�ֲ�����ֵ���㽫���ñ�����
extern double B(int n1,int n2, double x,double *Ux);

//  F�ֲ�����ֵ:����(0,x)�ϵĸ���p, f-�ܶ�ֵ
extern double F(int n1,int n2, double x,double *f);

//  F�ֲ��ķ�������p=F(0,x), ��֪p,����x
extern double re_F(int n1,int n2,double p);

//  t�ֲ��ķֲ�����ֵ
extern double student(int n, double x,double *f);

//  t�ֲ��ķ�����: p(-��,u)=p; ��֪p, ����u
extern double re_student(int n,double p); 


