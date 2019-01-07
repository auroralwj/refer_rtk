//
#ifndef AMBIGUITY_H
#define AMBIGUITY_H

#include "constant.h"
#include "matrix.h"
#include <cmath>
//#include "Lambda.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif
typedef matrix<double> Matrix;


//取整法固定整周模糊度
short FixAmbiguity(double *N,short count);


/**************************************************************************
********* 函数功能:快速模糊度解算    **********
********* 参数说明:in&out:Amb->输入模糊度浮点解，输出模糊度固定解             *********
*********          in:sigma0->浮点解中误差            **********
*********          in:AmbCount-> 模糊度的个数             ********
*********          in:Qxx-> 浮点解的协因素阵
*********          in:DfCount->双差观测方程的个数
*********          out:ratio->ratio值                    ********
****************************************************************************/
short FixAmbFARA(double *Amb,double sigma0, 
									 short AmbCount, Matrix Qxx, 
									 long DfCount,double* ratio);


/**************************************************************************
********* 函数功能:快速模糊度解算（扩充函数）    **********
********* 参数说明:in&out:Amb->输入模糊度浮点解，输出模糊度固定解             *********
*********          in:Omiga0->浮点解参差平方和            **********
*********          in:AmbCount-> 模糊度的个数             ********
*********          in:Qnn-> 浮点解的模糊度的协方差矩阵
*********          out:ratio->ratio值                    ********
****************************************************************************/
  short FixAmbFARAEx(double *Amb, double Omiga0, 
								   short AmbCount, Matrix Qnn,double* ratio);


//函数功能:用递推公式计算模糊度向量固定后求出的双差解的单位权中误差
//参数说明:AmbF->模糊度浮点解
//         sigmaF->浮点解的单位权中误差
//		   AmbN->被选的模糊度向量
//         AmbCount->模糊度的个数 
//         Qxx->浮点解的协因数阵
//         DfCount->双差观测方程的个数
//返回值：用被选的模糊度向量AmbN作为固定解得出的单位权中误差
  	double GetSigma(double* AmbF,double SigmaF,
									double* AmbN,short AmbCount,
									Matrix Qxx,long DfCount);
/**************************************************************************/
/********* 函数功能:用一下规则过滤模糊度向量，以便加速模糊度的搜索    **********/
//规则为：如果某个模糊度向量中的任意两个模糊度之差N_JK不满足以下条件
//则从备选模糊度向量中删除该模糊度向量
//N_jk-t(a/2)*mN_jk<=N_JK<=N_jk+t(a/2)*mN_jk
//N_jk:为这两个模糊度实数解之差
//t(a/2):置信度为（1-a）的t分布值
//mN_jk:整周模糊度的实数解的差值的验后方差
//参数说明:AmbVector:备选的模糊度向量组
//         VectorCount:备选的模糊度向量组的个数
//         Amb:浮点解模糊度
//         AmbCount:浮点解模糊度个数
//         Qxx:浮点解协因素阵
  	short FilterAmbVector(double** AmbVector, 
										  long *VectorCount,double* Amb,
										  short AmbCount, Matrix Qxx);

//函数功能:用区间法获取模糊度向量组
//参数说明
/*in*/
//Amb:浮点解模糊度
//AmbCount:浮点解模糊度个数
//Qxx:浮点解协因素阵
/*out*/
//AmbVector:备选的模糊度向量组
//VectorCount:备选的模糊度向量组的个数
//(注：调用该函数时，先使参数AmbVector为NULL，从而获得VectorCount,即模糊度向量的个数
// 然后根据VectorCount的大小为AmbVector分配内存，之后在调用该函数 ）

  	short GetAmbVector(double* Amb,short AmbCount,
									   double** AmbVector,long* VectorCount,
									   Matrix Qxx );

/*
//函数功能:区间法固定整周模糊度
//返回值：1.模糊度固定成功
//        0.模糊度部分固定成功，一般还需迭代固定其他的模糊度
            具体做法为：把部分以固定的模糊度当成已知值，重新求解其他的
			  模糊度，然后再用该方法固定模糊度。
		  -1.迭代失败，该方法固定模糊度失败
	*/
  	short FixAmbSection(double* Amb,short count,Matrix Qxx );


/*************************************************************************************
函数功能：
       用Lambda方法固定模糊度(不输出Ratio值)
参数说明：
       AmbL1（i/o）:L1双差模糊度,输入为浮点解模糊度，输出为固定后的整周模糊度
	   AmbL2（i/o）:L2双差模糊度,输入为浮点解模糊度，输出为固定后的整周模糊度（输入为NULL时，表示是单频）
       AmbCount（i）:单频模糊度个数
	   Qx:浮点解模糊度协方差矩阵（双频时为AmbCount*2 行 AmbCount*2 列）
	                             （单频时为AmbCount 行 AmbCount 列 或 AmbCount*2 行 AmbCount*2 列）
       
*************************************************************************************/

  short FixAmbLAMBDAEx(double *AmbL1,double *AmbL2,
									 short AmbCount,Matrix Qx);


/*************************************************************************************
函数功能：
       用Lambda方法固定模糊度(输出Ratio值)
参数说明：
       AmbL1（i/o）:L1双差模糊度,输入为浮点解模糊度，输出为固定后的整周模糊度
	   AmbL2（i/o）:L2双差模糊度,输入为浮点解模糊度，输出为固定后的整周模糊度（输入为NULL时，表示是单频）
       AmbCount（i）:单频模糊度个数
	   Qx:浮点解模糊度协方差矩阵（双频时为AmbCount*2 行 AmbCount*2 列）
	                             （单频时为AmbCount 行 AmbCount 列 或 AmbCount*2 行 AmbCount*2 列）
       
*************************************************************************************/

  short FixAmbLAMBDA(double *AmbL1,double *AmbL2,
								   short AmbCount,Matrix Qnn,double omiga0,
								   double *ratio);

//功能：用递推公式计算模糊度向量固定后求出的双差解的参差平方和
  double GetOmiga(double *AmbFloat,double Omiga0,
								double *AmbFix, short AmbCount, 
								Matrix Qnn);

//功能：用递推公式计算模糊度向量固定后求出的双差解的参差平方和
  double GetOmigaEx(double *AmbF,double Omiga0,
								  double *AmbN,short AmbCount, 
								  Matrix InvQnn);


#endif