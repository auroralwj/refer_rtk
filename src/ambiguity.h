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


//ȡ�����̶�����ģ����
short FixAmbiguity(double *N,short count);


/**************************************************************************
********* ��������:����ģ���Ƚ���    **********
********* ����˵��:in&out:Amb->����ģ���ȸ���⣬���ģ���ȹ̶���             *********
*********          in:sigma0->����������            **********
*********          in:AmbCount-> ģ���ȵĸ���             ********
*********          in:Qxx-> ������Э������
*********          in:DfCount->˫��۲ⷽ�̵ĸ���
*********          out:ratio->ratioֵ                    ********
****************************************************************************/
short FixAmbFARA(double *Amb,double sigma0, 
									 short AmbCount, Matrix Qxx, 
									 long DfCount,double* ratio);


/**************************************************************************
********* ��������:����ģ���Ƚ��㣨���亯����    **********
********* ����˵��:in&out:Amb->����ģ���ȸ���⣬���ģ���ȹ̶���             *********
*********          in:Omiga0->�����β�ƽ����            **********
*********          in:AmbCount-> ģ���ȵĸ���             ********
*********          in:Qnn-> ������ģ���ȵ�Э�������
*********          out:ratio->ratioֵ                    ********
****************************************************************************/
  short FixAmbFARAEx(double *Amb, double Omiga0, 
								   short AmbCount, Matrix Qnn,double* ratio);


//��������:�õ��ƹ�ʽ����ģ���������̶��������˫���ĵ�λȨ�����
//����˵��:AmbF->ģ���ȸ����
//         sigmaF->�����ĵ�λȨ�����
//		   AmbN->��ѡ��ģ��������
//         AmbCount->ģ���ȵĸ��� 
//         Qxx->������Э������
//         DfCount->˫��۲ⷽ�̵ĸ���
//����ֵ���ñ�ѡ��ģ��������AmbN��Ϊ�̶���ó��ĵ�λȨ�����
  	double GetSigma(double* AmbF,double SigmaF,
									double* AmbN,short AmbCount,
									Matrix Qxx,long DfCount);
/**************************************************************************/
/********* ��������:��һ�¹������ģ�����������Ա����ģ���ȵ�����    **********/
//����Ϊ�����ĳ��ģ���������е���������ģ����֮��N_JK��������������
//��ӱ�ѡģ����������ɾ����ģ��������
//N_jk-t(a/2)*mN_jk<=N_JK<=N_jk+t(a/2)*mN_jk
//N_jk:Ϊ������ģ����ʵ����֮��
//t(a/2):���Ŷ�Ϊ��1-a����t�ֲ�ֵ
//mN_jk:����ģ���ȵ�ʵ����Ĳ�ֵ����󷽲�
//����˵��:AmbVector:��ѡ��ģ����������
//         VectorCount:��ѡ��ģ����������ĸ���
//         Amb:�����ģ����
//         AmbCount:�����ģ���ȸ���
//         Qxx:�����Э������
  	short FilterAmbVector(double** AmbVector, 
										  long *VectorCount,double* Amb,
										  short AmbCount, Matrix Qxx);

//��������:�����䷨��ȡģ����������
//����˵��
/*in*/
//Amb:�����ģ����
//AmbCount:�����ģ���ȸ���
//Qxx:�����Э������
/*out*/
//AmbVector:��ѡ��ģ����������
//VectorCount:��ѡ��ģ����������ĸ���
//(ע�����øú���ʱ����ʹ����AmbVectorΪNULL���Ӷ����VectorCount,��ģ���������ĸ���
// Ȼ�����VectorCount�Ĵ�СΪAmbVector�����ڴ棬֮���ڵ��øú��� ��

  	short GetAmbVector(double* Amb,short AmbCount,
									   double** AmbVector,long* VectorCount,
									   Matrix Qxx );

/*
//��������:���䷨�̶�����ģ����
//����ֵ��1.ģ���ȹ̶��ɹ�
//        0.ģ���Ȳ��̶ֹ��ɹ���һ�㻹������̶�������ģ����
            ��������Ϊ���Ѳ����Թ̶���ģ���ȵ�����ֵ֪���������������
			  ģ���ȣ�Ȼ�����ø÷����̶�ģ���ȡ�
		  -1.����ʧ�ܣ��÷����̶�ģ����ʧ��
	*/
  	short FixAmbSection(double* Amb,short count,Matrix Qxx );


/*************************************************************************************
�������ܣ�
       ��Lambda�����̶�ģ����(�����Ratioֵ)
����˵����
       AmbL1��i/o��:L1˫��ģ����,����Ϊ�����ģ���ȣ����Ϊ�̶��������ģ����
	   AmbL2��i/o��:L2˫��ģ����,����Ϊ�����ģ���ȣ����Ϊ�̶��������ģ���ȣ�����ΪNULLʱ����ʾ�ǵ�Ƶ��
       AmbCount��i��:��Ƶģ���ȸ���
	   Qx:�����ģ����Э�������˫ƵʱΪAmbCount*2 �� AmbCount*2 �У�
	                             ����ƵʱΪAmbCount �� AmbCount �� �� AmbCount*2 �� AmbCount*2 �У�
       
*************************************************************************************/

  short FixAmbLAMBDAEx(double *AmbL1,double *AmbL2,
									 short AmbCount,Matrix Qx);


/*************************************************************************************
�������ܣ�
       ��Lambda�����̶�ģ����(���Ratioֵ)
����˵����
       AmbL1��i/o��:L1˫��ģ����,����Ϊ�����ģ���ȣ����Ϊ�̶��������ģ����
	   AmbL2��i/o��:L2˫��ģ����,����Ϊ�����ģ���ȣ����Ϊ�̶��������ģ���ȣ�����ΪNULLʱ����ʾ�ǵ�Ƶ��
       AmbCount��i��:��Ƶģ���ȸ���
	   Qx:�����ģ����Э�������˫ƵʱΪAmbCount*2 �� AmbCount*2 �У�
	                             ����ƵʱΪAmbCount �� AmbCount �� �� AmbCount*2 �� AmbCount*2 �У�
       
*************************************************************************************/

  short FixAmbLAMBDA(double *AmbL1,double *AmbL2,
								   short AmbCount,Matrix Qnn,double omiga0,
								   double *ratio);

//���ܣ��õ��ƹ�ʽ����ģ���������̶��������˫���Ĳβ�ƽ����
  double GetOmiga(double *AmbFloat,double Omiga0,
								double *AmbFix, short AmbCount, 
								Matrix Qnn);

//���ܣ��õ��ƹ�ʽ����ģ���������̶��������˫���Ĳβ�ƽ����
  double GetOmigaEx(double *AmbF,double Omiga0,
								  double *AmbN,short AmbCount, 
								  Matrix InvQnn);


#endif