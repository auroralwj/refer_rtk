

#include "stdio.h"

//////////////////////////////////////////////////////////////////////////
//概率计算函数文件

// 伽玛函数
extern double gam(int n);

//  正态分布的分布函数值: p(-∞,u)
extern double gaos(double u);

//  正态分布的反函数, p(-∞,u) ; 已知p, 返回u
extern double re_gaos(double q);

//  chi2分布的分布函数：p(0,x)
extern double chi2(int n,double x,double *f);

//  chi2分布的反函数：p=chi2(0,x),已知x，反求p
extern double re_chi2(int n,double p);

//B分布函数值，F分布和t分布函数值计算将调用本函数
extern double B(int n1,int n2, double x,double *Ux);

//  F分布函数值:区间(0,x)上的概率p, f-密度值
extern double F(int n1,int n2, double x,double *f);

//  F分布的反函数：p=F(0,x), 已知p,反求x
extern double re_F(int n1,int n2,double p);

//  t分布的分布函数值
extern double student(int n, double x,double *f);

//  t分布的反函数: p(-∞,u)=p; 已知p, 返回u
extern double re_student(int n,double p); 


