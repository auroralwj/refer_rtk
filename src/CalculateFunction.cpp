//CalculateFunction.cpp
#include "stdafx.h"
#include "..\INCLUDE\CalculateFunction.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//==============================================
//函数说明：对数据取整
//如：Round(2.5) = 3; Round(-2.5) = -3;
//    Round(2.1) = 2; Round(2.9) = 3
//    Round(-2.1) = -2; Round(-2.9) = -3
//=============================================
/*GPSCLASSDLL_API    double Round(double x)
{
	double y = 0;
	if (x>0)
	{
		y = floor(x+0.5);
	}
	else
	{
		y = ceil(x-0.5);
	}
	return y;
}*/
//changed by wkyu from zxw---------------------------------
//==============================================
//函数说明：对数据取整
//如：Round(2.5, -1) = 2; Round(2.5,1) = 3; Round(2.5, 0) = 3;
//    Round(-2.5, -1) = -3; Round(-2.5,1) = -2; Round(-2.5,0) = -3;
//    Round(2.1, -1) = 2; Round(2.1,1) = 3; Round(2.1,0) = 2;
//    Round(-2.1, -1) = -3; Round(-2.1,1)=-2; Round(-2.1,0) = -2;
//    Round(2.9, -1) = 2; Round(2.9,1) = 3; Round(2.9,0) = 3;
//    Round(-2.9,-1) = -3;Round(-2.9,1) = -2;Round(-2.9,0) = -3;
//=============================================
GPSCLASSDLL_API    double Round(double x, short type)
{
	double y = 0;
	switch (type)
	{
	case -1:   //往小的值
		y = floor(x);
		break;
	case 1: //取大值
		y = ceil(x);
		break;
	default:  //按照四舍五入取值
		y = x>0 ? floor(x+0.5) : ceil(x-0.5);
		break;
	}
	return y;
}
//----------------------------------------------------------
//-
//模糊度解算函数(适合单频双频等各种模糊度解算方法)
//   Amb->输入模糊度浮点解，输出模糊度固定解
//   ambCount-> 模糊度的个数
//   Cxx-> 浮点解的协因素阵
//   sigma0->浮点解中误差
//   ratio->ratio值      
//   DfCount->双差观测方程的个数
 /*
	AMB_INTEGER = 0,
	AMB_FARA,
	AMB_FARAEX,
	AMB_SECTION,
	AMB_LAMBDA,
	AMB_LAMBDAEX,
//
	AMB_FAST,	
	AMB_ARCE,
	AMB_LSSSC,
	AMB_LSSAMC,*/
GPSCLASSDLL_API short ResolveAmbiguity(short ambType, short obsType,
									   double* Amb, short ambCount, Matrix Cxx,
									   double sigma0, double* ratio, long DfCount)
{
	short rn = 0;
	switch(ambType)
	{
	case AMB_INTEGER:
		FixAmbiguity(Amb, ambCount); 
		break;
	case AMB_FARA:
		rn = FixAmbFARA(Amb, sigma0, ambCount, Cxx, DfCount, ratio);
		break;
	case AMB_FARAEX:
		rn = FixAmbFARAEx(Amb, sigma0, ambCount, Cxx, ratio);//yu:未完成
		break;
	case AMB_SECTION:
		rn = FixAmbSection(Amb, ambCount, Cxx);//yu:有问题
		break;
	case AMB_LAMBDA:
		{
			if(obsType==OBS_L1L2)
			{
				int i=0, n=ambCount/2;
				double *AmbL1 = new double[n], *AmbL2 = new double[n];
				for(i=0; i<n; i++)
				{
					AmbL1[i] = Amb[i];
					AmbL2[i] = Amb[i+n];
				}
				
				rn = FixAmbLAMBDA(AmbL1, AmbL2, n, Cxx, sigma0, ratio);
				
				for(i=0; i<n; i++)
				{
					Amb[i] = AmbL1[i];
					Amb[i+n] = AmbL2[i];
				}
				delete[] AmbL1; 
				delete[] AmbL2;
			}
			else if(obsType==OBS_L1)  //单频
			{
				rn = FixAmbLAMBDA(Amb, NULL, ambCount, Cxx, sigma0, ratio);
			}
		}
		break;
	default:
		return 0;
	}
	return rn;
}
//-
//组双差观测方程
GPSCLASSDLL_API short	FormDblDifference(CGPSStationObsSet* pFStObsSet, CGPSStationObsSet* pUStObsSet,
		              CGPSStationNavSet * pFNavSet,
					  CStationInfo* pFStationInfo, CStationInfo* pUStationInfo,
					  CGPSDDCoefSet* pDDCoefSet, CCalcSolutionBase* pSolSetting)
{
	CGPSEpochObsSet* pFEpochObsSet=NULL;
	CGPSEpochObsSet* pUEpochObsSet=NULL;

	CDDCoef* pDDCoef=NULL;//双差观测方程系数矩阵
    CGPSEpochDDCoef* pEpochDDCoef=NULL;
	
	POSITION fPos = pFStObsSet->GetHeadPosition();
	POSITION uPos = pUStObsSet->GetHeadPosition();
	
	while (NULL!=fPos && uPos!=NULL)
	{
		pFEpochObsSet = pFStObsSet->GetNext(fPos);
		pUEpochObsSet = pUStObsSet->GetNext(uPos);
		
		//计算双差观测方程系数
		pEpochDDCoef = pDDCoefSet->NewObjItem();		
		
		ComDDCoef(pFEpochObsSet, pFStationInfo,
					pUEpochObsSet, pUStationInfo, pFNavSet,
					pSolSetting->m_BaseSat, pEpochDDCoef, pSolSetting);
	}	
	return 1;
}

//功能:用拟合的方法探测周跳并存入周跳数据集中
//pDDCoefSet:双差观测数据集
//pCSSet: 周跳数据集
GPSCLASSDLL_API short DetectCycleslipWithPolynomial(CGPSDDCoefSet* pDDCoefSet, CGPSCircleSlipSet* pCSSet)
{
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;

	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	CGPSEpochDDCoef* pEpochDDcoef = NULL;
	CDDCoef* pObs = NULL;
	
	short i = 0, j = 0;
	double** DDOMC=NULL;
	long EpochCount=0;
    long satEphCount = 0;
    long ephCount = 0;
	for(i=0; i<satNum; i++)
	{
		ephCount = 0;
		satEphCount = pDDCoefSet->GetSatEpochCount(satPrn[i]);
		DDOMC=new double*[satEphCount];//二维数组，用来存储双差观测值序列
		for(j=0;j<satEphCount;j++)
		{
			DDOMC[j]=new double[4];
		}
		
		POSITION pPos = pDDCoefSet->GetHeadPosition();
		while(pPos != NULL)
		{
            pEpochDDcoef = pDDCoefSet->GetNext(pPos);
			//获取指定卫星的观测数据
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);
			
			if(pObs == NULL)  
			{
				continue;
			}
			
			DDOMC[ephCount][0] = pObs->EpochGpsWeek;
			DDOMC[ephCount][1] = pObs->EpochGpsSecond;
			//DDOMC[ephCount][2] = pObs->dW;
			DDOMC[ephCount][2] = pObs->dObsL1;
			
			if(pObs->IsL2Exist)
			{
				//DDOMC[ephCount][3] = pObs->dWL2;
				DDOMC[ephCount][3] = pObs->dObsL2;
			}
			else
			{
				DDOMC[ephCount][3] = 999;
			}
			ephCount++;
		}
		//探测周跳
		DetectCircleSlip(DDOMC,ephCount, satPrn[i], pCSSet, OBS_L1);
		if(DDOMC[0][2]!=999)
		{
			DetectCircleSlip(DDOMC,ephCount, satPrn[i], pCSSet, OBS_L1L2);
		}

		for(j=0; j<satEphCount; j++)
		{
			delete DDOMC[j];
		}
		if(DDOMC) delete DDOMC;
	}
	return 1;
}

//函数功能:从双差观测闭合差序列中DDOMC中探测周跳
//参数说明：DDOMC双差观测值序列（DDCount行4列的二维数组）
//          DDCount 双差观测值个数
//          pCSSet周跳数据集(用于存放探出的周跳)
GPSCLASSDLL_API short DetectCircleSlip(double** DDOMC,long& DDCount,
									   short SatID,CGPSCircleSlipSet* pCSSet, short Opt)
{
	if(Opt!=OBS_L1L2 && Opt!=OBS_L1) return 0;

	long j=0, k=0;
	double sigma=0;
	short mark=0;
	short ranknum=6;
	Matrix A(DDCount-1,ranknum);
	Matrix AT(ranknum,DDCount-1);
	Matrix w(DDCount-1,1);
	Matrix P(DDCount-1,DDCount-1);
	Matrix N(ranknum,ranknum);
	Matrix InvN(ranknum,ranknum);
	Matrix U(ranknum,1);
	Matrix x(ranknum,1);
	Matrix v(DDCount-1,1);
	P.Unit();

	//把时间标准化
	double t0 = DDOMC[0][1];
	
	do{	
		for(j=1; j<DDCount; j++)
		{			
			for(k=0; k<ranknum; k++)
			{
				A(j-1,k)=pow((DDOMC[j][1]-t0), k+1) - pow((DDOMC[j-1][1]-t0),k+1);
			}
			if(Opt==OBS_L1)
			{
				w(j-1,0)=DDOMC[j][2]-DDOMC[j-1][2];
			}
			else if(Opt==OBS_L1L2)
			{
				w(j-1,0)=DDOMC[j][3]-DDOMC[j-1][3];
			}
		}
		AT=~A;
		AT=AT*P;
		N=AT*A;
	try
	{
		InvN=!N;
	}
	catch (...)
	{
//		TRACEMatrix(N, "N");
		return 0;
	}
		U=AT*w;
		x=InvN*U;		
		v=A*x;
		v=v-w;			
		sigma=0;
		for(j=0;j<DDCount-1;j++)
		{
			sigma=sigma+v(j,0)*P(j,j)*v(j,0);
		}
		long delDDOMC=0;//记录被剔除的数据个数
		for(j=0;j<DDCount-1;j++)
		{
			if(P(j,j)==0) delDDOMC++;
		}
		sigma = sqrt(sigma/(DDCount-1-ranknum-delDDOMC));
		mark=0;
		for(j=0;j<DDCount-1;j++)
		{
			if(fabs(v(j,0))>3*sigma && P(j,j)!=0) 
			{
				P(j,j)=0; 
				mark=1;
			}
		}
	}while(mark==1);
	/////////////////////////////////////
	if(Opt==OBS_L1)
	{
		long usedDDCount=DDCount;
		for(j=0;j<DDCount-1;j++)
		{
			//理论上如果残差大于3倍中误差，且残差大于0.5则认为有周跳,
			if(fabs(v(j,0))>3*sigma && fabs(Round(v(j,0)))>0)
			{
				CCircleSlip* pcs = pCSSet->NewObjItem();            
				pcs->cs = Round(v(j,0));
				pcs->ObsType = 1;
				pcs->SatID = SatID;
				pcs->GpsWeek = DDOMC[j+1][0];
				pcs->GpsSecond = DDOMC[j+1][1];
				pcs->PreEphSecond = DDOMC[j][1];
				if(j==DDCount-2)
				{
					pcs->NextEphSecond = 0;
				}
				else 
				{
					pcs->NextEphSecond=DDOMC[j+2][1];
				}
			}
		}		
	}
	else if(Opt==OBS_L1L2)
	{
		long usedDDCount=DDCount;
		for(j=0;j<DDCount-1;j++)
		{
			//理论上如果残差大于3倍中误差，且残差大于0.5则认为有周跳,
			if(fabs(v(j,0))>3*sigma && fabs(Round(v(j,0)))>0)
			{
				CCircleSlip* pcs = pCSSet->NewObjItem();             
				pcs->cs = Round(v(j,0)-0.5);
				pcs->ObsType = 2;
				pcs->SatID = SatID;
				pcs->GpsWeek = DDOMC[j+1][0];
				pcs->GpsSecond = DDOMC[j+1][1];
				pcs->PreEphSecond = DDOMC[j][1];
				if(j==DDCount-2)
				{
					pcs->NextEphSecond=0;
				}
				else 
				{
					pcs->NextEphSecond=DDOMC[j+2][1];
				}
			}
		}	
	}
	return 1;
} 

//此函数暂时不用
GPSCLASSDLL_API void RemoveSamllandGetAmbCount(CGPSDDCoefSet *pStSet, short minContinue, short maxGap, CSatAmbSet* pSatAmbSet)
{
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
/*	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //观测记录的个数
	short ambCount = 0; //卫星的模糊度数
	long gapCount = 0;  //记录数据断开的长度

	CGPSEpochDDCoef* pEpochDDcoef = NULL;
	CGPSEpochDDCoef* pDelDDcoef = NULL;
	CDDCoef* pObs = NULL;
    CSatAmb* pSatAmb = NULL;

	for(i=0; i<satNum-1; i++)
	{
		POSITION pPos = pStSet->GetHeadPosition();
		ambCount = 0;
		obsCount = 0;
		gapCount = 0;

		while(pPos != NULL)
		{
            pEpochDDcoef = pStSet->GetNext(pPos);
			//获取指定卫星的观测数据
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);
			if(pObs != NULL)  //如果有观测数据，记录观测数据的个数
			{
				obsCount++;
			}
			else   //如果无观测数据，有以下情况
			{
				gapCount++;

				//数据先有，后无，说明数据中断 
				if(obsCount>0)
				{
					if(obsCount<minContinue)  //连续的数据记录小于指定数据记录长度，则删除该段数据
					{
						POSITION dPos = pPos;
						for(j=obsCount; j>0; j--)
						{
							pDelDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
							pDelDDcoef->Remove(satPrn[i]);   //删除该卫星该段的数据
						}
						gapCount = gapCount + obsCount; //该段数据删除后，相当于间隔加上数据长度
						obsCount = 0;
					}
					else// if(obsCount>=minContinue)
					{
						if(ambCount == 0) //第一个模糊度
						{
							//卫星增加一个新的模糊度
							ambCount = 1; 
							
							//记录模糊度信息
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = satPrn[i];       //卫星号
							pSatAmb->m_ambNum = ambCount;       //模糊度个数
							pSatAmb->m_GpsWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsSecond = pEpochDDcoef->m_GpsSecond;  
							
							obsCount = 0;
							gapCount = 1;
						}
						else //第2，3，4，...个模糊度
						{
							if(gapCount>maxGap)   //添加新的模糊度
							{
								ambCount++;
								//记录模糊度信息
								pSatAmb = pSatAmbSet->NewObjItem();
								pSatAmb->m_satID = satPrn[i];       //卫星号
								pSatAmb->m_ambNum = ambCount;       //模糊度个数
								pSatAmb->m_GpsWeek = pEpochDDcoef->m_GpsWeek;
								pSatAmb->m_GpsSecond = pEpochDDcoef->m_GpsSecond;  
								
								obsCount = 0;
								gapCount = 1;
							}
							else  //将上一个模糊度的时间改变
							{
								CSatAmb* plSat = pSatAmbSet->GetSatAmb(satPrn[i],ambCount);
								plSat->m_GpsWeek = pEpochDDcoef->m_GpsWeek;
								plSat->m_GpsSecond = pEpochDDcoef->m_GpsSecond;
							}
						}
					}
				}
				else
				{
					continue;
				}
			}
		}
	}*/
}

//-
//此函数用于编辑并获得删除短数据短后的双差系数集的模糊度信息获取
//对于断开时间长的数据，赋予新的模糊度，同时将卫星编号加50，相当于增加一个新卫星
BOOL EditSatAmb(CGPSDDCoefSet *pStSet, short maxGap, CSatAmbSet* pSatAmbSet)
{
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //观测记录的个数
	short ambCount = 0; //卫星的模糊度数
	long gapCount = 0;  //记录数据断开的长度

	CGPSEpochDDCoef* pEpochDDcoef = NULL;
	CGPSEpochDDCoef* pPrevEpochDDcoef = NULL;
	CGPSEpochDDCoef* pDelDDcoef = NULL;
	CGPSEpochDDCoef* pEditDDcoef = NULL;
	CDDCoef* pObs = NULL;
    CSatAmb* pSatAmb = NULL;
	POSITION tmpPos = NULL;
	POSITION pPos = NULL;
	POSITION dPos = NULL;

    short newSatid = 0;
	for(i=0; i<satNum; i++)
	{
		pPos = pStSet->GetHeadPosition();
		ambCount = 0;
		obsCount = 0;
		gapCount = 0;

		while(pPos != NULL)
		{
            pEpochDDcoef = pStSet->GetNext(pPos);
			//获取指定卫星的观测数据
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);
			if(pObs != NULL)  //如果有观测数据，记录观测数据的个数
			{
				obsCount++;

				if(pPos == NULL)  //已经到了最后一个历元
				{
					if(ambCount==0)   //此卫星数据一直连续，只记录卫星模糊度信息，不用附加处理
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						//卫星号变为原来的  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数
						
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = pEpochDDcoef->m_GpsSecond; 

						obsCount = 0;
						gapCount = 0;
					}
					else
					{
						ambCount++;
						//记录模糊度信息
						pSatAmb = pSatAmbSet->NewObjItem();
						//卫星号变为原来的  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = pEpochDDcoef->m_GpsSecond; 
					
						dPos = pStSet->GetTailPosition();
						//将新的数据段的卫星编号改变
						for(j=obsCount; j>0; j--)
						{
							pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
							pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
							if (NULL == pObs)
								return 0;
							pObs->SatID = newSatid;
						}

						obsCount = 0;
						gapCount = 0;
					}
				}
			}
			else //(pObs==NULL)  //如果无观测数据，有以下情况
			{
				gapCount++;

				//数据先有，后无，说明数据中断 
				if(obsCount>0)
				{
					if(ambCount == 0) //第一个模糊度
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						//卫星号变为原来的  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数						
						
						if(pPos == NULL) //倒数第一个历元
						{
							tmpPos = pStSet->GetTailPosition();
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
							pSatAmb->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond; 
							obsCount = 0;
							gapCount = 0;
						}
						else 
						{
							tmpPos = pPos;
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
							pSatAmb->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond; 
							obsCount = 0;
							gapCount = 1;
						}
					}
					else //第2，3，4，...个模糊度
					{
						if(gapCount>maxGap)   //添加新的模糊度
						{
							ambCount++;
							pSatAmb = pSatAmbSet->NewObjItem();
							//卫星号变为原来的  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							pSatAmb->m_satID = newSatid;       //卫星号
							pSatAmb->m_ambNum = ambCount;       //模糊度个数

							if(pPos == NULL)  //已经到最后一个历元
							{
								tmpPos = pStSet->GetTailPosition();
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								pSatAmb->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								pSatAmb->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
								//将此段的所有数据的卫星编号为改成新的卫星编号
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
									pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
									if (NULL == pObs)
										return 0;
									pObs->SatID = newSatid;
								}
								obsCount = 0;
								gapCount = 0;
							}		
							else //(pPos != NULL)
							{
								tmpPos = pPos;
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								pSatAmb->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								pSatAmb->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
															
								//将此段的所有数据的卫星编号为改成新的卫星编号
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
									pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
									if (NULL == pObs)
										return 0;
									pObs->SatID = newSatid;
								}
								
								obsCount = 0;
								gapCount = 1;
							}
						}
						else //(gapCount<maxGap)  //不添加新的模糊度，只将上一个模糊度的最后时间改变
						{
							newSatid = satPrn[i] + (ambCount-1) * 50;
							CSatAmb* plSat = pSatAmbSet->GetSatAmb(newSatid,ambCount);

							if(pPos == NULL)  //已经到最后一个历元
							{
								tmpPos = pStSet->GetTailPosition();
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								plSat->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								plSat->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
								//将此段的所有数据的卫星编号为改成新的卫星编号
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
									pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
									if (NULL == pObs)
										return 0;
									pObs->SatID = newSatid;
								}
								obsCount = 0;
								gapCount = 0;
							}		
							else //(pPos != NULL)
							{
								tmpPos = pPos;
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								plSat->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								plSat->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
															
								//将此段的所有数据的卫星编号为改成新的卫星编号
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
									pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
									if (NULL == pObs)
										return 0;
									pObs->SatID = newSatid;
								}
								
								obsCount = 0;
								gapCount = 1;
							}
						}
					}
				}
			}
		}
	}
	return 1;
}
//-
//功能:双频三差相对定位
//pDDCoefSet:  双差方程系数集
//XF0,YF0,ZF0：已知点的初始位置
//XU0,YU0,ZU0: 未知点的初始位置（I&O）
//sigma0:单位权中误差
//Cxx:协方差阵
//DfCount:组成的观测方程个数
//Mion:两历元间电离层观测值最大变化值（一般为4*WAVE1）
//SigmaL1:L1观测值的先验中误差
//SigmaL2:L2观测值的先验中误差
//dn1:L1周跳搜索的的宽度(一般为5)
//dn5:L5(宽巷)周跳搜索的的宽度(一般为2)
void TDRelativePosDF(CGPSDDCoefSet* pDDCoefSet,
					 double XF0, double YF0, double ZF0,
					 double *XU0, double *YU0, double *ZU0, 
					 double* sigma0, Matrix* Cxx, long* DfCount,
					 double Mion,double SigmaL1,double SigmaL2,short dn1,short dn5)
{

	//step:1--------------------矩阵准备与限差设置
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	short nRow = satNum * 2; //一个历元的最大双差观测方程个数 L1 & L2
    short nCol = 3;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow,1);
    Matrix Ui(nCol,1);
	Matrix U(nCol,1);
    Matrix x(nCol,1); 
	Matrix vv(nRow,1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);
	Matrix Qv(nCol, nCol);//协方差距阵
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//三差观测方程总数
	long   DeleteCount = 0 ;  //删除的历元
	double delta=0;//中误差
	double vL1=0;
	double vL2=0;
	double vL3=0;
	double X0=0, Y0=0, Z0=0;

	double Bata1=F1*F1/(F1*F1-F2*F2);
	double Bata2=-(F2*F2)/(F1*F1-F2*F2);

	//无电离层组合观测值先验中误差
	//yu:3倍的三差组合观测值中误差，即三差组合观测限差
	double SigmaL3=3*sqrt(8.0)*sqrt(Bata1*Bata1*SigmaL1*SigmaL1+Bata2*Bata2*SigmaL2*SigmaL2);
	short iteration=0;//迭代次数
	
	//step:2----------------------迭代求解,//yu:按历元列三差方程,把所有历元的法方程叠加在一起算一个值
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		vL1=0;
		vL2=0;
		delta=0;//中误差
		DbFuncCount=0;
		N.Null();
		U.Null();

       
        //step:2.1-----------------------yu:组成方程
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:按历元,把所有历元的法方程叠加在一起
		{             
			//获得第i历元双差观测方程系数
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//组成三差观测方程系数矩阵A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)//yu:逐双差方程
			{

				pDDCoef=pEpochDDCoef->GetNext(epochCoef);//当两个历元间组成三差时，观测卫星必须保持一一对应
				//获得第i+1历元双差观测方程系数 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 				
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//如果出现间断，判断此处是周跳或者是新增加模糊度
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //有数据了，退出循环
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//说明该卫星的观测数据只到第i历元为止 
				{
					continue;
				}

				if(iteration>0)
				{
					ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
					pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
					ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
						pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);					
				}
				
				A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;//yu:所有卫星都建立三差方程
				A(k,1)=pNextDDCoef->dEy-pDDCoef->dEy ;
				A(k,2)=pNextDDCoef->dEz-pDDCoef->dEz ;
				
				A(k+satNum,0)=pNextDDCoef->dEx-pDDCoef->dEx ;
				A(k+satNum,1)=pNextDDCoef->dEy-pDDCoef->dEy ;
				A(k+satNum,2)=pNextDDCoef->dEz-pDDCoef->dEz ;
				
				w(k,0)=pNextDDCoef->dW-pDDCoef->dW;  
				w(k+satNum,0)=pNextDDCoef->dWL2-pDDCoef->dWL2;  
				
				if(iteration>0)
				{                 
					//   vL1=A(k,0)*x(0,0)+A(k,1)*x(1,0)+A(k,2)*x(2,0)-w(k,0);
					//   w(k,0)=-vL1;
					//   vL2=A(k+satnum-1,0)*x(0,0)+A(k+satnum-1,1)*x(1,0)+A(k+satnum-1,2)*x(2,0)-w(k+satnum-1,0);
					//   w(k+satnum-1,0)=-vL2;
					vL1=-w(k,0);
					vL2=-w(k+satNum,0);
					vL3=Bata1*vL1+Bata2*vL2;
				}
	
				if(fabs(w(k,0)/(WAVE1))>10 || fabs(w(k+satNum,0)/(WAVE2))>10||
					0.5*fabs(vL1+F2*F2/(F1*F1)*vL2)>Mion || fabs(vL3)>SigmaL3 ||
					fabs(vL1)> Mion || fabs((F2*F2)/(F1*F1)*vL2) > Mion ) 
				{
					P(k,k)=0;
					P(k+satNum,k+satNum)=0;
					DeleteCount += 2;
				}
				if(P(k,k) > 0)
				{
					DbFuncCount = DbFuncCount + 2;//记录三差观测方程总数
					delta = delta+vL1*vL1+vL2*vL2;
				}
				
				k ++;
			}
            //yu:把所有历元的法方程叠加在一起
			//   A 行是L1，L2所在的不同卫星 , 列为dx dy dz
			AT=~A;
			ATP=AT*P;
			Ni=ATP*A;
			Ui=ATP*w;
			N=N+Ni;
			U=U+Ui;	
		}
		//step:2.2-----------------------yu:计算
		InvN=!N;
		x=InvN*U;
	
		*XU0=*XU0+x(0,0);
		*YU0=*YU0+x(1,0);
		*ZU0=*ZU0+x(2,0);
		
		iteration++;
		//若迭代超过10次还不收敛则认为是发散的,终止迭代
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);//yu:----迭代


	//step:3--------------------yu:结果整理
	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//单位权中误差
	Qv = delta * InvN;
	TRACE("%f\t%f\t%f\n", Qv(0,0), Qv(0,1), Qv(0,2));
	TRACE("%f\t%f\n", Qv(1,1), Qv(1,2));
	TRACE("%f\n", Qv(2,2));

	for(m=0; m<3; m++)
	{
		for(j=0;j<3;j++)
		{
			(*Cxx)(m,j)=Qv(m,j);
		}
	}
	return;
}
void TDRelativePosDF_IonoFree(CGPSDDCoefSet* pDDCoefSet,
					 double XF0, double YF0, double ZF0,
					 double *XU0, double *YU0, double *ZU0, 
					 double* sigma0, Matrix* Cxx, long* DfCount,
					 double Mion,double SigmaL1,double SigmaL2,short dn1,short dn5)
{

	//step:1--------------------矩阵准备与限差设置
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	short nRow = satNum; 
    short nCol = 3;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow,1);
    Matrix Ui(nCol,1);
	Matrix U(nCol,1);
    Matrix x(nCol,1); 
	Matrix vv(nRow,1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);
	Matrix Qv(nCol, nCol);//协方差距阵
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//三差观测方程总数
	long   DeleteCount = 0 ;  //删除的历元
	double delta=0;//中误差
	double X0=0, Y0=0, Z0=0;

	short iteration=0;//迭代次数
	
	//step:2----------------------迭代求解,//yu:按历元列三差方程,把所有历元的法方程叠加在一起算一个值
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		delta=0;//中误差
		DbFuncCount=0;
		N.Null();
		U.Null();

       
        //step:2.1-----------------------yu:组成方程
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:按历元,把所有历元的法方程叠加在一起
		{             
			//获得第i历元双差观测方程系数
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//组成三差观测方程系数矩阵A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)//yu:逐双差方程
			{

				pDDCoef=pEpochDDCoef->GetNext(epochCoef);//当两个历元间组成三差时，观测卫星必须保持一一对应
				//获得第i+1历元双差观测方程系数 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 				
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//如果出现间断，判断此处是周跳或者是新增加模糊度
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //有数据了，退出循环
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//说明该卫星的观测数据只到第i历元为止 
				{
					continue;
				}

				if(iteration>0)
				{
					ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
					pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
					ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
						pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);					
				}
				
				A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;//yu:所有卫星都建立三差方程
				A(k,1)=pNextDDCoef->dEy-pDDCoef->dEy ;
				A(k,2)=pNextDDCoef->dEz-pDDCoef->dEz ;

				w(k,0)=pNextDDCoef->dWL3-pDDCoef->dWL3;

				
				if(fabs(w(k,0)) > Mion ) 
				{
					P(k,k)=0;
					P(k+satNum,k+satNum)=0;
					DeleteCount += 1;
				}
				if(P(k,k) > 0)
				{
					DbFuncCount = DbFuncCount + 1;//记录三差观测方程总数
					delta += w(k,0)*w(k,0);
				}
				
				k ++;
			}
            //yu:把所有历元的法方程叠加在一起
			//   A 行是L1，L2所在的不同卫星 , 列为dx dy dz
			AT=~A;
			ATP=AT*P;
			Ni=ATP*A;
			Ui=ATP*w;
			N=N+Ni;
			U=U+Ui;	
		}
		//step:2.2-----------------------yu:计算
		InvN=!N;
		x=InvN*U;
	
		*XU0=*XU0+x(0,0);
		*YU0=*YU0+x(1,0);
		*ZU0=*ZU0+x(2,0);
		
		iteration++;
		//若迭代超过10次还不收敛则认为是发散的,终止迭代
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);//yu:----迭代


	//step:3--------------------yu:结果整理
	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//单位权中误差
	Qv = delta * InvN;
	TRACE("%f\t%f\t%f\n", Qv(0,0), Qv(0,1), Qv(0,2));
	TRACE("%f\t%f\n", Qv(1,1), Qv(1,2));
	TRACE("%f\n", Qv(2,2));

	for(m=0; m<3; m++)
	{
		for(j=0;j<3;j++)
		{
			(*Cxx)(m,j)=Qv(m,j);
		}
	}
	return;
}
//-
//三差法探测周跳(双频)
void DetectCycleslipWithTD(CGPSDDCoefSet* pDDCoefSet, 
					double XF0, double YF0, double ZF0,
					double XU0, double YU0, double ZU0,
					double Mion,double SigmaL1,double SigmaL2,short dn1,short dn5,
					CGPSCircleSlipSet* pCSSet)
{
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;

	double delta=0;//中误差
	double vL1=0;
	double vL2=0;
	double vL3=0;
    double Bata1=F1*F1/(F1*F1-F2*F2);
	double Bata2=-(F2*F2)/(F1*F1-F2*F2);

	//无电离层组合观测值先验中误差
	//yu:极限误差3倍的组合观测值单位权中误差
	double SigmaL3=3*sqrt(8.0)*sqrt(Bata1*Bata1*SigmaL1*SigmaL1+Bata2*Bata2*SigmaL2*SigmaL2);

	//记录周跳
	long CSL1=0;
	long CSL2=0;

	if(pCSSet!=NULL)
	{
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:所有历元
		{             
			//获得第i历元双差观测方程系数
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();
			while (NULL!=epochCoef && ddPos!=NULL)//yu:所有方程
			{
				//获得第i+1历元双差观测方程系数 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				//当两个历元间组成三差时，观测卫星必须保持一一对应
				pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);

				if(NULL == pNextDDCoef)//出现间断，然后找到第i+s历的双差观测值来组成三差
				{		
					POSITION tempPos = ddPos;
					pNextEpochDDCoef=pDDCoefSet->GetNext(tempPos); 
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef)//yu:直到找到该卫星某历元双差
						{
							break;
						}
						pNextDDCoef = NULL;
					}					
				}
				if(!pNextDDCoef)//说明该卫星的观测数据只到第i历元为止 
				{
					continue;//yu:仅当pNextDDCoef一直为空时才到这，下一方程
					
				}
			   
			   ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,XU0,YU0,ZU0,
					   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
					   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
    		   ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,XU0,YU0,ZU0,
					   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
					   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);

			   
               vL1 = -(pNextDDCoef->dW - pDDCoef->dW);  
			   vL2 = -(pNextDDCoef->dWL2 - pDDCoef->dWL2);  
     		   vL3 = Bata1 * vL1 + Bata2 * vL2;//yu:无电离层影响组合观测值三差残差
 //
			   if(0.5*fabs(vL1+F2*F2/(F1*F1)*vL2)>Mion || fabs(vL3)>SigmaL3 ||
				   	fabs(vL1)> Mion || fabs((F2*F2)/(F1*F1)*vL2) > Mion) //满足这两个条件时认为有周跳发生
			   {
                   CCircleSlip * pcs = NULL;
				   //此处如果是粗差，将探测不出周跳，一般用下一历元探测，如果还不行，则将此处的数据删除。
				   CGPSCircleSlipSet::SearchCycleSlip(vL1,vL2,CSL1,CSL2,Mion,SigmaL3,dn1,dn5);
	                if(CSL1!=0)
					{
						pcs = pCSSet->NewObjItem();
						pcs->cs=CSL1; 
						pcs->ObsType=1;
						pcs->SatID=pNextDDCoef->SatID;
						pcs->GpsSecond=pNextEpochDDCoef->m_GpsSecond;
						pcs->GpsWeek=pNextEpochDDCoef->m_GpsWeek;
					}
					if(CSL2!=0)
					{
						pcs = pCSSet->NewObjItem();
						pcs->cs=CSL2; 
						pcs->ObsType=2;
						pcs->SatID=pNextDDCoef->SatID;
						pcs->GpsSecond=pNextEpochDDCoef->m_GpsSecond;
						pcs->GpsWeek=pNextEpochDDCoef->m_GpsWeek;
					}
			   }
			}
		}
	}
	return;
}

//-
//双频双差相对定位浮点解
void DDFloatRelativePosDF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double *XU0, double *YU0, double *ZU0, 
					double *Ambiguity,double* sigma0, Matrix* Cxx, long* DfCount, double* rms)
{
    //获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
    
	short nRow = 2 * satNum;
	short nCol = 2 * satNum + 3;
	Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol,1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);     P.Null();
	Matrix ATP(nCol, nRow);
   	Matrix Cv(nCol, nCol);//协方差距阵
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//中误差
	long   DbFuncCount=0;//双差观测方程总数

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		P.Null();
		//P.Unit();

		satCount = pEpochDDCoef->GetItemCount();		
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
			   P(ii+satNum, jj+satNum) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		   P(ii+satNum, ii+satNum) = satCount;
		}
		
		short pos = -1;
		POSITION ddCoefPos =  pEpochDDCoef->GetHeadPosition();
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos = GetPosInMat(satPrn, satNum, pDDCoef->SatID);
			if(pos==-1)
			{
				ASSERT(FALSE);
				continue;
			} 

		   A(pos, pos+3)=WAVE1;
		   A(satNum+pos, satNum+pos+3) = WAVE2;
		  
		   A(pos,0) = pDDCoef->dEx;
		   A(pos,1) = pDDCoef->dEy;
		   A(pos,2) = pDDCoef->dEz;
		   A(satNum+pos,0) = pDDCoef->dEx;
		   A(satNum+pos,1) = pDDCoef->dEy;
		   A(satNum+pos,2) = pDDCoef->dEz;

		   w(pos,0) = pDDCoef->dW;
		   w(satNum+pos,0) = pDDCoef->dWL2;

		   DbFuncCount += 2;
		}
 
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}
	short i=0, j=0;
    for(i=0; i<nCol; i++)
	{
		for(j=0;j<nCol; j++)
		{
			TRACE("%f\t", N(i,j));
		}
		TRACE("\n");
	}

	InvN = !N;
	x = InvN * U;
      
    *XU0 = *XU0 + x(0,0);
	*YU0 = *YU0 + x(1,0);
	*ZU0 = *ZU0 + x(2,0);

	//模糊度浮点解,这里是双差模糊度
	for(i=0; i<nCol-3; i++)
	{
        Ambiguity[i] = x(i+3,0);
	}
 	
	///计算残差平方和
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //求RMS用

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
		P.Null();
		//P.Unit();
		short satCount = pEpochDDCoef->GetItemCount();
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
			   P(ii+satNum, jj+satNum) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		   P(ii+satNum, ii+satNum) = satCount;
		}

		POSITION ddCoefPos = pEpochDDCoef->GetHeadPosition();
		short pos=-1;
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
           pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
		   if(pos==-1)
		   {
			   continue;
		   }
		   
		   ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
			   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
			   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, &ww);
		  
		   w(pos,0) = ww - WAVE1 * Ambiguity[pos];    
		   //w(pos,0) = pDDCoef->dW - WAVE1 * Ambiguity[pos];    
		   w(satNum+pos,0) = pDDCoef->dWL2 - WAVE2 * Ambiguity[satNum+pos];

		   v(pos,0) = -w(pos,0);
		   v(satNum+pos,0) = -w(satNum+pos,0);
		}
		
        vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//---------
		vTv = vT*v;
		vv = vv + vTv(0,0);
	}		
	
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv = delta * InvN;	
	for(int m=0; m<nCol; m++)
	{
		for(int j=0; j<nCol;j++)
		{
			(*Cxx)(m,j)=Cv(m,j);
		}
	}
	return ;
}
void DDFloatRelativePosDF_WL(CGPSDDCoefSet *pDDCoefSet, 
								   double XF0, double YF0, double ZF0, 
								   double *XU0, double *YU0, double *ZU0, 
								   double *WL_AMB,double* sigma0, Matrix* Cxx, long* DfCount, double* rms){
	short satPrn[40];
	short satNum = 0;
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
    
	short nRow = satNum;
	short nCol = satNum + 3;
	Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol,1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);     P.Null();
	Matrix ATP(nCol, nRow);
   	Matrix Cv(nCol, nCol);//协方差距阵
	
	double ww = 0;
	short  satCount=0;
	double aa[3];
	
	double delta=0;//中误差
	long   DbFuncCount=0;//双差观测方程总数
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	
    POSITION epochPos = pDDCoefSet->GetHeadPosition();
	
	//------------------- 1: wide ambiguity resolution -----------------
	
	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
		//组成观测方程系数矩阵A
		A.Null();
		w.Null();
		P.Null();
		
		satCount = pEpochDDCoef->GetItemCount();		
		
		for(short ii=0; ii<satNum; ii++)
		{
			for(short jj=0; jj<satNum; jj++)
			{
				P(ii,jj) = -1.0;
			}
			P(ii,ii) = satCount; 
		}
		
		short pos = -1;
		POSITION ddCoefPos =  pEpochDDCoef->GetHeadPosition();
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos = GetPosInMat(satPrn, satNum, pDDCoef->SatID);
			if(pos==-1)
			{
				ASSERT(FALSE);
				continue;
			} 
			
			A(pos, pos+3)=WAVEW;
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;
			w(pos,0) = pDDCoef->dWL5;
			DbFuncCount += 1;
		}
		
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}
	InvN = !N;
	x = InvN * U;
    *XU0 = *XU0 + x(0,0);
	*YU0 = *YU0 + x(1,0);
	*ZU0 = *ZU0 + x(2,0);
	
	//wide ambiguity
	short i=0, j=0;
	for(i=0; i<nCol-3; i++)
	{
        WL_AMB[i] = ceil(x(i+3,0)-0.5);
	}
	///计算残差平方和
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //求RMS用
	
	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
		//组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
		P.Null();
		//P.Unit();
		short satCount = pEpochDDCoef->GetItemCount();
		for(short ii=0; ii<satNum; ii++)
		{
			for(short jj=0; jj<satNum; jj++)
			{
				P(ii,jj) = -1.0;
			}
			P(ii,ii) = satCount; 
		}
		
		POSITION ddCoefPos = pEpochDDCoef->GetHeadPosition();
		short pos=-1;
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1)
			{
				continue;
			}
			double w5,ifw;
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, &ww,&ifw,&w5);
			
			w(pos,0) = w5 - WAVEW * WL_AMB[pos];
			v(pos,0) = -w(pos,0);
		}
		
        vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//---------
		vTv = vT*v;
		vv = vv + vTv(0,0);
	}		
	
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv = delta * InvN;	
	for(int m=0; m<nCol; m++)
	{
		for(int j=0; j<nCol;j++)
		{
			(*Cxx)(m,j)=Cv(m,j);
		}
	}
	return ;
}
void DDFloatRelativePosDF_IonoFree(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double *XU0, double *YU0, double *ZU0, 
					double *AmbL5,double *AmbL1,double* sigma0, Matrix* Cxx, long* DfCount, double* rms)
{
    //获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
    
	short nRow = satNum;
	short nCol = satNum + 3;
	Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol,1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);     P.Null();
	Matrix ATP(nCol, nRow);
   	Matrix Cv(nCol, nCol);//协方差距阵
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//中误差
	long   DbFuncCount=0;//双差观测方程总数

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	//------------------- 1: wide ambiguity resolution -----------------

	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		P.Null();

		satCount = pEpochDDCoef->GetItemCount();		
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		}
		
		short pos = -1;
		POSITION ddCoefPos =  pEpochDDCoef->GetHeadPosition();
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos = GetPosInMat(satPrn, satNum, pDDCoef->SatID);
			if(pos==-1)
			{
				ASSERT(FALSE);
				continue;
			} 

		   A(pos, pos+3)=WAVEN;
		   A(pos,0) = pDDCoef->dEx;
		   A(pos,1) = pDDCoef->dEy;
		   A(pos,2) = pDDCoef->dEz;
		   w(pos,0) = pDDCoef->dWL3-AmbL5[pos]*C*F2*1e-6/(F1*F1-F2*F2);//L3 & L5
		   DbFuncCount += 1;
		}
 
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}
	InvN = !N;
	x = InvN * U;
    *XU0 = *XU0 + x(0,0);
	*YU0 = *YU0 + x(1,0);
	*ZU0 = *ZU0 + x(2,0);

	//ambiguity
	short i=0, j=0;
	for(i=0; i<nCol-3; i++)
	{
        AmbL1[i] = x(i+3,0);//L1
	}

    for(i=0; i<nCol; i++)
	{
		for(j=0;j<nCol; j++)
		{
			TRACE("%f\t", N(i,j));
		}
		TRACE("\n");
	}
	
	//----------------------------------------
	
	///计算残差平方和
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //求RMS用

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
		P.Null();
		//P.Unit();
		short satCount = pEpochDDCoef->GetItemCount();
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		}

		POSITION ddCoefPos = pEpochDDCoef->GetHeadPosition();
		short pos=-1;
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
           pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
		   if(pos==-1)
		   {
			   continue;
		   }
		   double w5,ifw;
		   ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
			   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
			   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, &ww,&ifw,&w5);

		   ww=ifw-AmbL5[pos]*C*F2*1e-6/(F1*F1-F2*F2)-AmbL1[pos]*WAVEN;    
		   v(pos,0) = -ww;
		}
		
        vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//---------
		vTv = vT*v;
		vv = vv + vTv(0,0);
	}		
	
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv = delta * InvN;	
	for(int m=0; m<nCol; m++)
	{
		for(int j=0; j<nCol;j++)
		{
			(*Cxx)(m,j)=Cv(m,j);
		}
	}
	return ;
}
void DDFixedRelativePosDF_IonoFree(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double* XU0, double* YU0, double* ZU0, 
					double* AmbL5,double* AmbL1, double* sigma0, Matrix* Cxx, long* DfCount, double* rms, char* chFile)
{
	//打开残差数据文件并写入文件头数据
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[60];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 60);
	
	if(fResidul)
	{
       fprintf(fResidul,"RESIDUAL\n");
       fprintf(fResidul,"%d\n",satNum);
	   //fprintf(fResidul,"%d\n",pDDCoefSet->BaseSat);
	   for(short f=0; f<satNum; f++)
	   {
		   fprintf(fResidul,"%d\t",satPrn[f]);		   
	   }
	   fprintf(fResidul,"\n");
	   long totalepochcount=pDDCoefSet->GetItemCount();
	   fprintf(fResidul,"%d\n",totalepochcount);
    } 
///////////////////////////////////////////////////////

    short nRow = satNum;
	short nCol = 3;
	long i=0;
	long j=0;
	long k=0;
	short m=0;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol, 1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);

	Matrix Cv(nCol, nCol);//协方差距阵
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//中误差

	long   DbFuncCount=0;//双差观测方程总数
///////////////////////////////////////////
	long epochCount=pDDCoefSet->GetItemCount();
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

	POSITION ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(ddCoefPos);
        
		A.Null();
		w.Null();
		P.Null();
		
		satCount=pEpochDDCoef->GetItemCount();
		
	    for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		}
					
		 //组成观测方程系数矩阵A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			double IFw,ww5;
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww,&IFw,&ww5);
				
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;

			w(pos,0) = IFw - WAVEN * AmbL1[pos] - C*F2*1e-6/(F1*F1-F2*F2)* AmbL5[pos];    

			DbFuncCount = DbFuncCount+1;
		}
		
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}

    InvN=!N;
	x=InvN*U;
 
    *XU0=*XU0+x(0,0);
	*YU0=*YU0+x(1,0);
	*ZU0=*ZU0+x(2,0);

	Matrix v(nRow,1);
	Matrix vT(1, nRow);
	Matrix vTP(1, nRow);
	Matrix vTPv(1,1);

	ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef=pDDCoefSet->GetNext(ddCoefPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
        P.Null();
		//P.Unit();
		
		satCount=pEpochDDCoef->GetItemCount();
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		}
	
		short pos=-1;
		if(fResidul)
		{
			fprintf(fResidul,"%.1f\t%.2f\t%d\t",pEpochDDCoef->m_GpsWeek,pEpochDDCoef->m_GpsSecond,satCount);
			POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
			while (NULL != epochddCoefPos)
			{
               pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			   fprintf(fResidul,"%d\t",pDDCoef->SatID);
			}
			fprintf(fResidul,"\n");
		}

		POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			double ifw;
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww,&ifw);
				  
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;

			w(pos,0) = ifw - WAVEN * AmbL1[pos] - C*F2*1e-6/(F1*F1-F2*F2)* AmbL5[pos]; 

			v(pos,0)=-w(pos,0);
		}
		vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//--------------------
		vTv = vT * v;
		vv = vv + vTv(0,0);
		 //把残差观测值写入残差数据文件中
		if(fResidul)
		{
			short vvv=0;
			for(vvv=0; vvv<satNum; vvv++)
			{
				fprintf(fResidul,"%lf\t",v(vvv,0));				
			}
			fprintf(fResidul, "\n");
			for(vvv=0; vvv<satCount; vvv++)
			{
				fprintf(fResidul,"%lf\t",v(vvv+satNum,0));				
			}
			fprintf(fResidul, "\n");
		}
	}
	if(fResidul)
	{
	   fclose(fResidul);
	}
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//单位权中误差
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv=delta*InvN;

	for(m=0;m<3;m++)
	{
		for(j=0;j< 3  ;j++)
		{
			(*Cxx)(m,j) = Cv(m,j);
		}
	}
	return ;
}
void DDFixedRelativePosDF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double* XU0, double* YU0, double* ZU0, 
					double* Ambiguity, double* sigma0, Matrix* Cxx, long* DfCount, double* rms, char* chFile)
{
	//打开残差数据文件并写入文件头数据
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[60];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 60);
	
	if(fResidul)
	{
       fprintf(fResidul,"RESIDUAL\n");
       fprintf(fResidul,"%d\n",satNum);
	   //fprintf(fResidul,"%d\n",pDDCoefSet->BaseSat);
	   for(short f=0; f<satNum; f++)
	   {
		   fprintf(fResidul,"%d\t",satPrn[f]);		   
	   }
	   fprintf(fResidul,"\n");
	   long totalepochcount=pDDCoefSet->GetItemCount();
	   fprintf(fResidul,"%d\n",totalepochcount);
    } 
///////////////////////////////////////////////////////

    short nRow = 2 * satNum;
	short nCol = 3;
	long i=0;
	long j=0;
	long k=0;
	short m=0;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol, 1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);

	Matrix Cv(nCol, nCol);//协方差距阵
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//中误差

	long   DbFuncCount=0;//双差观测方程总数
///////////////////////////////////////////
	long epochCount=pDDCoefSet->GetItemCount();
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

	POSITION ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(ddCoefPos);
        
		A.Null();
		w.Null();
		P.Null();
	//	P.Unit();
		
		satCount=pEpochDDCoef->GetItemCount();
		
	    for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
			   P(ii+satNum, jj+satNum) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		   P(ii+satNum, ii+satNum) = satCount;
		}
					
		 //组成观测方程系数矩阵A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			double IFw,ww5;
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww,&IFw,&ww5);
				
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;
			A(satNum+pos,0) = pDDCoef->dEx;
			A(satNum+pos,1) = pDDCoef->dEy;
			A(satNum+pos,2) = pDDCoef->dEz;

			//w(pos,0) = pDDCoef->dW - WAVE1*Ambiguity[pos];
			w(pos,0) = ww - WAVE1 * Ambiguity[pos];    
			w(satNum+pos,0) = pDDCoef->dWL2 - WAVE2*Ambiguity[pos+satNum];
			
			DbFuncCount = DbFuncCount+2;
		}
		
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}

		InvN=!N;
	x=InvN*U;
 
    *XU0=*XU0+x(0,0);
	*YU0=*YU0+x(1,0);
	*ZU0=*ZU0+x(2,0);

	Matrix v(nRow,1);
	Matrix vT(1, nRow);
	Matrix vTP(1, nRow);
	Matrix vTPv(1,1);

	ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef=pDDCoefSet->GetNext(ddCoefPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
        P.Null();
		//P.Unit();
		
		satCount=pEpochDDCoef->GetItemCount();
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
			   P(ii+satNum, jj+satNum) = -1.0;
		   }
		   P(ii,ii) = satCount; 
		   P(ii+satNum, ii+satNum) = satCount;
		}
	
		short pos=-1;
		if(fResidul)
		{
			fprintf(fResidul,"%.1f\t%.2f\t%d\t",pEpochDDCoef->m_GpsWeek,pEpochDDCoef->m_GpsSecond,satCount);
			POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
			while (NULL != epochddCoefPos)
			{
               pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			   fprintf(fResidul,"%d\t",pDDCoef->SatID);
			}
			fprintf(fResidul,"\n");
		}

		POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww);
				  
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;
			A(satNum+pos,0) = pDDCoef->dEx;
			A(satNum+pos,1) = pDDCoef->dEy;
			A(satNum+pos,2) = pDDCoef->dEz;

			//w(pos,0) = pDDCoef->dW - WAVE1*Ambiguity[pos];
			w(pos,0) = ww - WAVE1*Ambiguity[pos];
			w(satNum+pos,0) = pDDCoef->dWL2 - WAVE2*Ambiguity[pos+satNum];

			v(pos,0)=-w(pos,0);
			v(satNum+pos,0) = -w(satNum+pos,0);
		}
		vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//--------------------
		vTv = vT * v;
		vv = vv + vTv(0,0);
		 //把残差观测值写入残差数据文件中
		if(fResidul)
		{
			short vvv=0;
			for(vvv=0; vvv<satNum; vvv++)
			{
				fprintf(fResidul,"%lf\t",v(vvv,0));				
			}
			fprintf(fResidul, "\n");
			for(vvv=0; vvv<satCount; vvv++)
			{
				fprintf(fResidul,"%lf\t",v(vvv+satNum,0));				
			}
			fprintf(fResidul, "\n");
		}
	}
	if(fResidul)
	{
	   fclose(fResidul);
	}
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//单位权中误差
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv=delta*InvN;

	for(m=0;m<3;m++)
	{
		for(j=0;j< 3  ;j++)
		{
			(*Cxx)(m,j) = Cv(m,j);
		}
	}
	return ;
}

//功能:单频三差相对定位
//pDDCoefSet:  双差方程系数集
//XF0,YF0,ZF0：已知点的初始位置
//XU0,YU0,ZU0: 未知点的初始位置（I&O）
//sigma0:单位权中误差
//Cxx:协方差阵
//DfCount:组成的观测方程个数
BOOL TDRelativePosSF(CGPSDDCoefSet* pDDCoefSet,
					 double XF0, double YF0, double ZF0,
					 double *XU0, double *YU0, double *ZU0, 
					 double* sigma0, Matrix* Cxx, long* DfCount)
{
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	short nRow = satNum; //一个历元的最大双差观测方程个数
    short nCol = 3;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow,1);
    Matrix Ui(nCol,1);
	Matrix U(nCol,1);
    Matrix x(nCol,1);
	Matrix vv(nRow,1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);
	Matrix Qv(nCol, nCol);//协方差距阵
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//三差观测方程总数
	long   DeleteCount = 0 ;  //删除的历元
	double delta=0;//中误差
	double vL1=0;
	double X0=0, Y0=0, Z0=0;

	short iteration=0;//迭代次数
	
	//迭代求解
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		vL1=0;
		delta=0;//中误差
		DbFuncCount=0;
		N.Null();
		U.Null();

		POSITION ddPos = pDDCoefSet->GetHeadPosition();
	
		while (NULL != ddPos)
		{             
			//获得第i历元双差观测方程系数
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//组成三差观测方程系数矩阵A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)
			{
				//获得第i+1历元双差观测方程系数 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				//当两个历元间组成三差时，观测卫星必须保持一一对应
				pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//如果出现间断，判断此处是周跳或者是新增加模糊度
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 

					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //有数据了，退出循环
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//说明该卫星的观测数据只到第i历元为止 
				{
					continue;
				}

				if(iteration>0)
				{
					ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
					pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
					ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
						pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);					
				}
				
				A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;
				A(k,1)=pNextDDCoef->dEy-pDDCoef->dEy ;
				A(k,2)=pNextDDCoef->dEz-pDDCoef->dEz ;
				
				w(k,0)=pNextDDCoef->dW-pDDCoef->dW;  
						
				if(iteration>0)
				{                 
					//   vL1=A(k,0)*x(0,0)+A(k,1)*x(1,0)+A(k,2)*x(2,0)-w(k,0);
					//   w(k,0)=-vL1;
					//   vL2=A(k+satnum-1,0)*x(0,0)+A(k+satnum-1,1)*x(1,0)+A(k+satnum-1,2)*x(2,0)-w(k+satnum-1,0);
					//   w(k+satnum-1,0)=-vL2;
					vL1=-w(k,0);
				}
	
				if(fabs(w(k,0)/(WAVE1))>10 || (vL1)>0.8) 
				{
					P(k,k)=0;
					DeleteCount ++;
				}
				if(P(k,k) > 0)
				{
					DbFuncCount = DbFuncCount + 1;//记录三差观测方程总数
					delta = delta+vL1*vL1;
				}
				
				k ++;
			}
			
			AT=~A;
			ATP=AT*P;
			Ni=ATP*A;
			Ui=ATP*w;
			N=N+Ni;
			U=U+Ui;	
		}
		
		if (1 != CalcXWithNandU(&N, &U, & x))
		{
			return 0;
		}
	
		*XU0=*XU0+x(0,0);
		*YU0=*YU0+x(1,0);
		*ZU0=*ZU0+x(2,0);
		
		iteration++;
		//若迭代超过10次还不收敛则认为是发散的,终止迭代
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);

	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//单位权中误差
	Qv = delta * InvN;
	TRACE("%f\t%f\t%f\n", Qv(0,0), Qv(0,1), Qv(0,2));
	TRACE("%f\t%f\n", Qv(1,1), Qv(1,2));
	TRACE("%f\n", Qv(2,2));

	for(m=0; m<3; m++)
	{
		for(j=0;j<3;j++)
		{
			(*Cxx)(m,j)=Qv(m,j);
		}
	}
	return 1;
}

/*
//三差法探测周跳(单频)
BOOL DetectCycleslipWithTD(CGPSDDCoefSet* pDDCoefSet, 
					double XF0, double YF0, double ZF0,
					double XU0, double YU0, double ZU0,
					CGPSCircleSlipSet* pCSSet)
{

//	FILE* fOut;
//	fOut=fopen("matrixsult.txt","w");
	////////////////////////////////////////////////////////////////////////
    long i=0;
	long j=0;
	long k=0;
	long l=0;
	short m=0;

	short satPrn[40];
	short satnum = 0;
	pDDCoefSet->GetAllObsSat(&satnum, satPrn, 40);
	//
    Matrix A(satnum-1,3);
	Matrix AT(3,satnum-1);
	Matrix Ni(3,3);
	Matrix N(3,3);
	Matrix InvN(3,3);
	Matrix w(satnum-1,1);
    Matrix Ui(3,1);
	Matrix U(3,1);
    Matrix x(3,1);
	Matrix vv(satnum-1,1);
	Matrix P(satnum-1,satnum-1);
	Matrix ATP(3,satnum-1);
	Matrix Qv(3,3);//协方差距阵
	
     
////////////////////////////////////////////////////

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    short satcount=0;
///////////////////////////////////////////

	long   DbFuncCount=0;//三差观测方程总数
	double delta=0;//中误差
	double v=0;
	double X0,Y0,Z0;
	double aa[3],aa0[3];
	double ww,ww0;
	short iteration=0;//迭代次数
	//迭代求解

    do
	{
		X0 = XU0;
  	    Y0 = YU0;
	    Z0 = ZU0;
		v=0;
		delta=0;//中误差
		DbFuncCount=0;
		N.Null();
		U.Null();

		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)
		{
			//获得第i历元双差观测方程系数
            pEpochDDCoef=pDDCoefSet->GetNext(pos); 
			if (NULL == pos)
			{//数据处理完毕,无后续数据
				break;
			}
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 

			//组成三差观测方程系数矩阵A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)
			{
				//获得第i+1历元双差观测方程系数 
				//当两个历元间组成三差时，观测卫星必须保持一一对应
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				if (NULL == pNextDDCoef)
					continue;
//				ASSERT(NULL != pNextDDCoef); //数据必须是已经对齐的
				
				ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
				ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
				
				
				A(k,0)=aa[0]-aa0[0];
				A(k,1)=aa[1]-aa0[1];
				A(k,2)=aa[2]-aa0[2];	  
				w(k,0)=ww-ww0; 
				
				if(iteration>0)
				{
					v=-w(k,0);
				}
				
				if(P(k,k)>0)//权为0表示观测值已被剔除 
				{
					DbFuncCount=DbFuncCount+1;//记录三差观测方程总数
				}
				if(abs(w(k,0)/(WAVE1))>10||abs(v/(WAVE1))>0.5) 
				{
					P(k,k)=0;
				}
				
			 }
			 
			 AT=~A;
			 ATP=AT*P;
			 Ni=ATP*A;
			 Ui=ATP*w;
			 N=N+Ni;
			 U=U+Ui;	
		}
		///////////////////////////////////////
		try
		{
			InvN=!N;
			x=InvN*U;
		}
		catch(...)
		{
#ifdef SAVE_ERR_LOG_INFO
			double year,month,day, hour,minute,second;
			if (epochCount > 0)
			{
				pos = pDDCoefSet->GetHeadPositon();
			    CGPSEpochDDCoef * pEpochDDCoef=pDDCoefSet->GetNext(pos); 
				GetYMDHMS(pEpochDDCoef->m_GpsWeek, pEpochDDCoef->m_GpsSecond,
					& year,& month,& day,& hour,& minute,& second);
			}
			FILE * hFile = fopen("calcLog.txt", "a+");
			fprintf(hFile, "%s:%d\t %.f-%.f-%.f-%.f-%.f-%.f \terror\n", 
				__FILE__, __LINE__,
				year,month,day, hour,minute,second);
			PrintMatrix(N, N.RowNo(), N.ColNo(), hFile);
			fclose(hFile);
#endif		
		return FALSE;
		}
		XU0=XU0+x(0,0);
		YU0=YU0+x(1,0);
		ZU0=ZU0+x(2,0);
		
		iteration++;
		//若迭代办处10次还不收敛则认为是发散的,终止迭代
		if(iteration>10)
			return FALSE;

	}while(abs(x(0,0))>0.02 && abs(x(1,0))>0.02);

	delta=0;
	//记录周跳	
	if(pCSSet!=NULL)
	{
		X0=XU0;
  	    Y0=YU0;
	    Z0=ZU0;
		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)
		{
			//获得第i历元双差观测方程系数
			pEpochDDCoef=pDDCoefSet->GetNext(pos);
			if (NULL == pos)
				break;
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 
			//组成三差观测方程系数矩阵A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)
			{
				//当两个历元间组成三差时，观测卫星必须保持一一对应
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				//
				if (NULL == pNextDDCoef)
					continue;
				
//				ASSERT(NULL != pNextDDCoef); //数据必须是已经对齐的
				///////////////
               ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
			   ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
			   v=-(ww-ww0);
			   if((abs(v)/(WAVE1))>0.5) 
			   {
                   CCircleSlip* pcs= pCSSet->NewObjItem();;
				 // 残差大于0.5周的被认为有周跳
					pcs->cs=ceil(v/(WAVE1)-0.5); 
					pcs->ObsType=1;
					pcs->SatID=pNextDDCoef->SatID;
					pcs->GpsSecond=pNextEpochDDCoef->m_GpsSecond;
					pcs->GpsWeek=pNextEpochDDCoef->m_GpsWeek;
					//
			   }
			   else
			   {
				    delta=delta+v*v;
			   }
			   
			}
            pDDCoefSet->MoveNext();
		}	    
	}		  

	delta=delta/(DbFuncCount-3);
//	if(DfCount!=NULL) * DfCount=DbFuncCount;
//	*sigma0=sqrt(delta);//单位权中误差
//	Qv=delta*InvN;
//   for(m=0;m<3;m++)
//   {
//	  for(j=0;j<3;j++)
//
  // 	  (*Qxx)(m,j)=Qv(m,j);
 //  }
   return TRUE;
}
*/

//-
//yu:三差求解实数解，利用三差观测残差大于0.5周判断为周跳
//三差法探测周跳(单频)
BOOL DetectCycleslipWithTD(CGPSDDCoefSet* pDDCoefSet, 
					double XF0, double YF0, double ZF0,
					double XU0, double YU0, double ZU0,
					CGPSCircleSlipSet* pCSSet)
{
	 long i=0, j=0, k=0, l=0;
	 short m = 0;
	 long epochCount = 0;
	 
	 epochCount = pDDCoefSet->GetItemCount();
	 
	 //获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	 short satPrn[40];
	 short satNum = 0;
	 
	 pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
	 
	 short nRow = satNum; //一个历元的最大双差观测方程个数
	 short nCol = 3;
	 
	 Matrix A(nRow, nCol);
	 Matrix AT(nCol, nRow);
	 Matrix Ni(nCol, nCol);
	 Matrix N(nCol, nCol);
	 Matrix InvN(nCol, nCol);
	 Matrix w(nRow,1);
	 Matrix Ui(nCol,1);
	 Matrix U(nCol,1);
	 Matrix x(nCol,1);
	 Matrix vv(nRow,1);
	 Matrix P(nRow, nRow);
	 Matrix ATP(nCol, nRow);
	 Matrix Qv(nCol, nCol);//协方差距阵
	 
	 CGPSEpochDDCoef* pEpochDDCoef=NULL;
	 CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	 CDDCoef* pDDCoef=NULL;
	 CDDCoef* pNextDDCoef=NULL;
	 
	 short satcount=0;
	 long   DbFuncCount=0;//三差观测方程总数
	 long   DeleteCount = 0 ;  //删除的历元
	 double delta=0;//中误差
	 double vL1=0;
	 double X0=0, Y0=0, Z0=0;
	 
	 short iteration=0;//迭代次数
	 double aa[3],aa0[3];
 	 double ww,ww0;

	 //迭代求解
	 do
	 {
		 X0=XU0;
		 Y0=YU0;
		 Z0=ZU0;
		 vL1=0;
		 delta=0;//中误差
		 DbFuncCount=0;
		 N.Null();
		 U.Null();
		 
		 POSITION ddPos = pDDCoefSet->GetHeadPosition();
		 
		 while (NULL != ddPos)
		 {             
			 //获得第i历元双差观测方程系数
			 pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			 //组成三差观测方程系数矩阵A
			 A.Null();
			 w.Null();
			 P.Unit();
			 
			 k = 0;
			 POSITION epochCoef = pEpochDDCoef->GetHeadPosition();
			 
			 while (NULL!=epochCoef && ddPos!=NULL)
			 {
				 //获得第i+1历元双差观测方程系数 
				 pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				 //当两个历元间组成三差时，观测卫星必须保持一一对应
				 pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				 pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				 
				 //如果出现间断，判断此处是周跳或者是新增加模糊度
				 if(NULL == pNextDDCoef)
				 {					
					 POSITION tempPos = ddPos;
					 pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					 
					 while (NULL != tempPos)
					 {
						 pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						 pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						 if(pNextDDCoef) //有数据了，退出循环
						 {
							 break;
						 }
						 pNextDDCoef = NULL;
					 }
					 
				 }
				 if(!pNextDDCoef)//说明该卫星的观测数据只到第i历元为止 
				 {
					 continue;
				 }

				 ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
				 ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
				 
				 
				 A(k,0)=aa[0]-aa0[0];
				 A(k,1)=aa[1]-aa0[1];
				 A(k,2)=aa[2]-aa0[2];	  
				 w(k,0)=ww-ww0; 
				 /*
				 if(iteration>0)
				 {
					 ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
						 pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						 pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
					 ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,
						 pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
						 pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
				 }
				 
				 A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;
				 A(k,1)=pNextDDCoef->dEy-pDDCoef->dEy ;
				 A(k,2)=pNextDDCoef->dEz-pDDCoef->dEz ;
				 
				 w(k,0)=pNextDDCoef->dW-pDDCoef->dW;  
				 */
				 if(iteration>0)
				 {                 
					 vL1=-w(k,0);
				 }
				 
				 if(fabs(w(k,0)/(WAVE1))>10 || (vL1)>0.8) 
				 {
					 P(k,k)=0;
					 DeleteCount ++;
				 }
				 if(P(k,k) > 0)
				 {
					 DbFuncCount = DbFuncCount + 1;//记录三差观测方程总数
					 delta = delta+vL1*vL1;
				 }
				 
				 k ++;
			 }
			 
			 AT=~A;
			 ATP=AT*P;
			 Ni=ATP*A;
			 Ui=ATP*w;
			 N=N+Ni;
			 U=U+Ui;	
		 }
		 
		 if (1 != CalcXWithNandU(&N, &U, & x))
		 {
			 return 0;
		 }
		 
		 XU0=XU0+x(0,0);
		 YU0=YU0+x(1,0);
		 ZU0=ZU0+x(2,0);
		 
		 iteration++;
		 //若迭代超过10次还不收敛则认为是发散的,终止迭代
		 if(iteration>10) break;

		 
//		 TRACEMatrix(N, "N");
//		 TRACEMatrix(U, "U");
//		 TRACEMatrix(x, "X");
		 
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02 );
	
	delta = delta / (DbFuncCount-3);

	//================6/3/2011 added WKYU=========================
	try
	{
		InvN=!N;
	}
	catch (...)
	{
//		TRACEMatrix(N, "N");
		return 0;
	}
	//============================end=============================
	Qv = delta * InvN;
	TRACE("%f\t%f\t%f\n", Qv(0,0), Qv(0,1), Qv(0,2));
	TRACE("%f\t%f\n", Qv(1,1), Qv(1,2));
	TRACE("%f\n", Qv(2,2));
	delta=0;
	//记录周跳	
	if(pCSSet!=NULL)
	{
		X0=XU0;
  	    Y0=YU0;
	    Z0=ZU0;
		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)//yu:本基线所有历元双差方程
		{
			//获得第i历元双差观测方程系数
			pEpochDDCoef=pDDCoefSet->GetNext(pos);
			if (NULL == pos)
				break;
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 
			//组成三差观测方程系数矩阵A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)//yu:前后历元所有卫星双差方程
			{
				//当两个历元间组成三差时，观测卫星必须保持一一对应
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				//
				if (NULL == pNextDDCoef)//yu:如果中断就下一个方程
					continue;
				///////////////
               ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
			   ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
			   double v=-(ww-ww0);//yu:即三差残差
			   if((abs(v)/(WAVE1))>0.5) ///////////////////////////////////////////////
			   {
                   CCircleSlip* pcs= pCSSet->NewObjItem();
				 // 残差大于0.5周的被认为有周跳
					pcs->cs=ceil(v/(WAVE1)-0.5); 
					pcs->ObsType=1;
					pcs->SatID=pNextDDCoef->SatID;					
					pcs->GpsSecond=pNextEpochDDCoef->m_GpsSecond;
					pcs->GpsWeek=pNextEpochDDCoef->m_GpsWeek;
					//
			   }
			   else
			   {
				    delta=delta+v*v;
			   }
			   
			}
            pDDCoefSet->MoveNext();
		}	    
	}		  
   return TRUE;
}

//单频双差相对定位浮点解
void DDFloatRelativePosSF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double *XU0, double *YU0, double *ZU0, 
					double *Ambiguity,double* sigma0, Matrix* Cxx, long* DfCount, double* rms)
{
    //获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
    
	short nRow = satNum;
	short nCol = satNum + 3;
	Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol,1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);     P.Null();
	Matrix ATP(nCol, nRow);
   	Matrix Cv(nCol, nCol);//协方差距阵
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//中误差
	long   DbFuncCount=0;//双差观测方程总数

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	while (NULL != epochPos)//yu:所有方程
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		P.Null();
		//P.Unit();

		satCount = pEpochDDCoef->GetItemCount();		
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;			   
		   }
		   P(ii,ii) = satCount; 
		}
		
		short pos = -1;
		POSITION ddCoefPos =  pEpochDDCoef->GetHeadPosition();
		while (NULL != ddCoefPos)//yu:所有方程
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos = GetPosInMat(satPrn, satNum, pDDCoef->SatID);
			if(pos==-1)
			{
				ASSERT(FALSE);
				continue;
			} 

		   A(pos, pos+3)=WAVE1;
				  
		   A(pos,0) = pDDCoef->dEx;
		   A(pos,1) = pDDCoef->dEy;
		   A(pos,2) = pDDCoef->dEz;
		
		   w(pos,0) = pDDCoef->dW;
		
		   DbFuncCount += 1;
		}
 
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}

//	TRACEMatrix(N, "N");

	InvN = !N;
	x = InvN * U;
      
    *XU0 = *XU0 + x(0,0);
	*YU0 = *YU0 + x(1,0);
	*ZU0 = *ZU0 + x(2,0);
//	TRACEMatrix(x, "X");

	//模糊度浮点解
	for(int i=0; i<nCol-3; i++)
	{
        Ambiguity[i] = x(i+3,0);
	}
 	
	///计算残差平方和
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //求RMS用

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
		P.Null();
		//P.Unit();
		short satCount = pEpochDDCoef->GetItemCount();
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
			   		   }
		   P(ii,ii) = satCount; 		   
		}

		POSITION ddCoefPos = pEpochDDCoef->GetHeadPosition();
		short pos=-1;
		while (NULL != ddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
           pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
		   if(pos==-1)
		   {
			   continue;
		   }
		   
		   ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
			   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
			   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, &ww);
		  
		   w(pos,0) = ww - WAVE1 * Ambiguity[pos];    
		   //w(pos,0) = pDDCoef->dW - WAVE1 * Ambiguity[pos];    
		   
		   v(pos,0) = -w(pos,0);		   
		}
		
        vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//---------
		vTv = vT*v;
		vv = vv + vTv(0,0);
	}		
	
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv = delta * InvN;	
	for(int m=0; m<nCol; m++)
	{
		for(int j=0; j<nCol;j++)
		{
			(*Cxx)(m,j)=Cv(m,j);
		}
	}
	return ;
}

void DDFixedRelativePosSF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double* XU0, double* YU0, double* ZU0, 
					double* Ambiguity, double* sigma0, Matrix* Cxx, long* DfCount, double* rms, char* chFile)
{
	//打开残差数据文件并写入文件头数据
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[60];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 60);
	
	if(fResidul)
	{
       fprintf(fResidul,"RESIDUAL\n");
       fprintf(fResidul,"%d\n",satNum);
	   //fprintf(fResidul,"%d\n",pDDCoefSet->BaseSat);
	   for(short f=0; f<satNum; f++)
	   {
		   fprintf(fResidul,"%d\t",satPrn[f]);		   
	   }
	   fprintf(fResidul,"\n");
	   long totalepochcount=pDDCoefSet->GetItemCount();
	   fprintf(fResidul,"%d\n",totalepochcount);
    } 
///////////////////////////////////////////////////////

    short nRow = satNum;
	short nCol = 3;
	long i=0;
	long j=0;
	long k=0;
	short m=0;

    Matrix A(nRow, nCol);
	Matrix AT(nCol, nRow);
	Matrix Ni(nCol, nCol);
	Matrix N(nCol, nCol);
	Matrix InvN(nCol, nCol);
	Matrix w(nRow, 1);
    Matrix Ui(nCol, 1);
	Matrix U(nCol, 1);
    Matrix x(nCol, 1);
	Matrix P(nRow, nRow);
	Matrix ATP(nCol, nRow);

	Matrix Cv(nCol, nCol);//协方差距阵
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//中误差

	long   DbFuncCount=0;//双差观测方程总数
///////////////////////////////////////////
	long epochCount=pDDCoefSet->GetItemCount();
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

	POSITION ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(ddCoefPos);
        
		A.Null();
		w.Null();
		P.Null();
	//	P.Unit();
		
		satCount=pEpochDDCoef->GetItemCount();
		
	    for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;			   
		   }
		   P(ii,ii) = satCount; 		   
		}
					
		 //组成观测方程系数矩阵A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			double IFw,ww5;
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww,&IFw,&ww5);
				
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;
			
			//w(pos,0) = pDDCoef->dW - WAVE1*Ambiguity[pos];
			w(pos,0) = ww - WAVE1 * Ambiguity[pos];    			
			
			DbFuncCount = DbFuncCount+1;
		}
		
        AT=~A;
		ATP=AT*P;
		Ni=ATP*A;
		Ui=ATP*w;
		N=N+Ni;
		U=U+Ui;
	}

		InvN=!N;
	x=InvN*U;
 
    *XU0=*XU0+x(0,0);
	*YU0=*YU0+x(1,0);
	*ZU0=*ZU0+x(2,0);

	Matrix v(nRow,1);
	Matrix vT(1, nRow);
	Matrix vTP(1, nRow);
	Matrix vTPv(1,1);

	ddCoefPos = pDDCoefSet->GetHeadPosition();
	while (NULL != ddCoefPos)
	{
		pEpochDDCoef=pDDCoefSet->GetNext(ddCoefPos);
         //组成观测方程系数矩阵A
		A.Null();
		w.Null();
		v.Null();
        P.Null();
		//P.Unit();
		
		satCount=pEpochDDCoef->GetItemCount();
		
		for(short ii=0; ii<satNum; ii++)
		{
		   for(short jj=0; jj<satNum; jj++)
		   {
			   P(ii,jj) = -1.0;
		   }
		   P(ii,ii) = satCount; 		   
		}
	
		short pos=-1;
		if(fResidul)
		{
			fprintf(fResidul,"%.1f\t%.2f\t%d\t",pEpochDDCoef->m_GpsWeek,pEpochDDCoef->m_GpsSecond,satCount);
			POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
			while (NULL != epochddCoefPos)
			{
               pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			   fprintf(fResidul,"%d\t",pDDCoef->SatID);
			}
			fprintf(fResidul,"\n");
		}

		POSITION epochddCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochddCoefPos)
		{
			pDDCoef=pEpochDDCoef->GetNext(epochddCoefPos);
			//satprn数组中的卫星编号需排序，该算法才有效
			pos=GetPosInMat(satPrn,satNum,pDDCoef->SatID);
			if(pos==-1) continue;
			
			ComputeDDA(pDDCoef,XF0,YF0,ZF0,*XU0,*YU0,*ZU0,aa,
				pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1 ,
				pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2, 
				&ww);
				  
			A(pos,0) = pDDCoef->dEx;
			A(pos,1) = pDDCoef->dEy;
			A(pos,2) = pDDCoef->dEz;
			
			//w(pos,0) = pDDCoef->dW - WAVE1*Ambiguity[pos];
			w(pos,0) = ww - WAVE1*Ambiguity[pos];
			
			v(pos,0)=-w(pos,0);			
		}
		vT=~v;
		vTP=vT*P;
		vTPv=vTP*v;
        delta=vTPv(0,0)+delta;
		//--------------------
		vTv = vT * v;
		vv = vv + vTv(0,0);
		 //把残差观测值写入残差数据文件中
		if(fResidul)
		{
			short vvv=0;
			for(vvv=0; vvv<satNum; vvv++)
			{
				fprintf(fResidul,"%lf\t",v(vvv,0));				
			}
			fprintf(fResidul, "\n");			
		}
	}
	if(fResidul)
	{
	   fclose(fResidul);
	}
	////////////////////////////////////////////////
	//中误差及协方差阵
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//单位权中误差
	*rms = sqrt(vv/(DbFuncCount-1));   //Rms值
	
	Cv=delta*InvN;

	for(m=0;m<3;m++)
	{
		for(j=0;j< 3  ;j++)
		{
			(*Cxx)(m,j) = Cv(m,j);
		}
	}
	return ;
}

//-
//大气改正
void CorrectTrop(CGPSStationObsSet * pStationSet, 
				 CGPSStationMetSet * pMetSet, double B, double H,
				 double dTemperature, double dPressure, double dRh, double H0,
				 short modelType, TropMode tropParam)
{
	POSITION pos = pStationSet->GetHeadPosition();
	while (NULL != pos)
	{
		CGPSEpochObsSet * pEphObs = pStationSet->GetNext(pos);
		CMet * pMet = pMetSet->GetMet(pEphObs->m_GpsWeek, pEphObs->m_GpsSecond);
		CorrectTrop(pEphObs, pMet, B, H, dTemperature, dPressure, dRh, H0, modelType, tropParam);
	}
}
//-
//modelType:对流层模型
void CorrectTrop(CGPSEpochObsSet * pEphObs, CMet * pMet, 
				 double B, double H,
				 double dTemperature, double dPressure, double dRh, double H0,
				 short modelType, TropMode tropParam)
{
	//P0:参考高度的大气压(milliBar)
	//T0:参考高度的干温度(K)
	//RH0:参考高度的相对湿度(%) 
	//H0:参考高度(m)
	//Ph:高度h处的气压(milliBar)
	//Th:高度h处的温度(K)
	//Pv:高度h处的湿气压(milliBar)
	//RHh:高度h处的相对湿度%
    //h:高度(m)
	double T0=291.15;        
	double P0=1013.25;      
	double RH0=50;          
	double Pa=0, Td=0, Pv=0, trop=0;
	
	if (NULL != pMet)   //用记录的气象元素
	{
		//double tpa=0, ttd=0, tpv=0, HR=0;
		//calmeteor(H,P0,T0,RH0,&tpa,&ttd,&tpv, & HR);
		Td = pMet->TD + 273.15;
		Pa = pMet->PR; //
		Pv = GetPvByTdAndHR(Td, pMet->RH); 
	}
	else if(dTemperature!=0 && dPressure!=0 && dRh!=0)  //用海平面测站平均气象元素
	{
		T0 = dTemperature + 273.15;
		P0 = dPressure;
		RH0 = dRh;
		//calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
		double RHh = 0;
		CalMeteorParameter(P0, T0, RH0, H0, H, &Pa, &Td, &Pv, & RHh);
	}
	else  //用标准气象元素 
	{
		//calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
		double RHh = 0;
		CalMeteorParameter(P0, T0, RH0, H0, H, &Pa, &Td, &Pv, & RHh);
	}
	
	//对流层改正
	POSITION pos = pEphObs->GetHeadPosition();
	double week, tow;
	week=pEphObs->m_GpsWeek;
	tow=pEphObs->m_GpsSecond;
    CGPSTime t(week,tow);
	int doy=t.GetDayIndex();
	while (NULL != pos)
	{
		CObservation * pObsData = pEphObs->GetNext(pos);
	    double WZD,mf_wet;//added wkyu 2014
		TropoCorrectByModel( modelType,Pa, Td, Pv, pObsData->elevation, H, B, doy, &trop,&WZD, &mf_wet);
		
		pObsData->vtrop = trop;
		pObsData->vWZD = WZD;
		pObsData->vMF = mf_wet;
		
		if (trop < 200)
		{//trop异常时,不考虑
			double dtrop_ = trop;
			pObsData->vL1 -= dtrop_/WAVE1;
			if(fabs(pObsData->vL2) >0.001)//表示有L2观测值时才作改正，以方便后面判断是否有L2观测值
			{
				pObsData->vL2 -= dtrop_/WAVE2; 
			}
		}
	}

}
//-
//----------------zxw1 2008-3-8------------------------------------------ 
//根据参考高度的气压，干温，相对湿度等计算指定高度的气压，水气压，温度
//P0:参考高度的大气压(milliBar)
//T0:参考高度的干温度(K)
//RH0:参考高度的相对湿度(%) 
//H0:参考高度(m)
//Ph:高度h处的气压(milliBar)
//Th:高度h处的温度(K)
//Pv:高度h处的湿气压(milliBar)
//RHh:高度h处的相对湿度%
//h:高度(m)
void CalMeteorParameter(double P0, double T0, double RH0, double H0, double Hh, double *Ph, double *Th, double *Pv , double *RHh)
{
	double aa=0, bb=0, cc=0, dH=0;
	dH = Hh - H0;
	
	//求h处的气压
	aa = 1.0 - 2.26e-5 * dH;
	bb = pow(aa, 5.225);
	*Ph = P0 * bb;
	
	//求h处的温度	
	*Th = T0 - 0.0065 * dH;

	//求h处的相对湿度
	double g=0, d=0;
	g = -6.396e-4 * dH;
	d = exp(g); 
	
	double relhum = RH0 * d; 
	if (NULL != RHh)
	{
		*RHh = relhum;
	}

	//求h处的水气压
	double tt = *Th - 273.15;
	if(tt> 0.0)
	{
		cc=(7.5*tt/(tt+237.3))+0.786;
	}
	if(tt<=0.0) 
	{
		cc=(9.5*tt/(tt+265.5))+0.786;
	}
	
	double es=pow(10.0,cc);
	
	*Pv = es*(relhum/100.0);
}
//大气改正
void CorrectTrop(CGPSStationObsSet * pStationSet, 
				 CGPSStationMetSet * pMetSet, double B, double H,
				 double dTemperature, double dPressure, double dRh,
				 short modelType, BOOL bTropParam)
{
	POSITION pos = pStationSet->GetHeadPosition();
	while (NULL != pos)
	{
		CGPSEpochObsSet * pEphObs = pStationSet->GetNext(pos);
		CMet * pMet = pMetSet->GetMet(pEphObs->m_GpsWeek, pEphObs->m_GpsSecond);

		CorrectTrop(pEphObs, pMet, B, H, dTemperature, dPressure, dRh, modelType, bTropParam);
	}
}

void CorrectTrop(CGPSEpochObsSet * pEphObs, CMet * pMet, 
				 double B, double H,
				 double dTemperature, double dPressure, double dRh,
				 short modelType, BOOL bTropParam)
{
	double T0=291.16;
	double P0=1013.25;
	double RH0=50;
	double Pa=0, Td=0, Pv=0, trop=0;
	
	if (NULL != pMet)   //用记录的气象元素
	{
		//double tpa=0, ttd=0, tpv=0, HR=0;
		//calmeteor(H,P0,T0,RH0,&tpa,&ttd,&tpv, & HR);
		Td = pMet->TD + 273.16;
		Pa = pMet->PR; //
		Pv = GetPvByTdAndHR(Td, pMet->RH); 
	}
	else if(dTemperature!=0 && dPressure!=0 && dRh!=0)  //用海平面测站平均气象元素
	{
		T0 = dTemperature + 273.16;
		P0 = dPressure;
		RH0 = dRh;
		calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
	}
	else  //用标准气象元素 
	{
		calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
	}
	
	//对流层改正
	POSITION pos = pEphObs->GetHeadPosition();
	double week, tow;
	week=pEphObs->m_GpsWeek;
	tow=pEphObs->m_GpsSecond;
    CGPSTime t(week,tow);
	int doy=t.GetDayIndex();

	while (NULL != pos)
	{
		CObservation * pObsData = pEphObs->GetNext(pos);
	    double WZD,mf_wet;//added wkyu 2014
		TropoCorrectByModel( modelType,Pa, Td, Pv, pObsData->elevation, H, B, doy, &trop,&WZD, &mf_wet);
		
		pObsData->vtrop = trop;
		pObsData->vWZD = WZD;
		pObsData->vMF = mf_wet;

		if (trop < 200)
		{//trop异常时,不考虑
			double dtrop_ = (bTropParam?(trop - mf_wet*WZD):trop);
			pObsData->vL1 -= dtrop_/WAVE1;
			if(fabs(pObsData->vL2) >0.001)//表示有L2观测值时才作改正，以方便后面判断是否有L2观测值
			{
				pObsData->vL2 -= dtrop_/WAVE2; 
			}
		}
	}
}
//-
//对流层模型改正
void TropoCorrectByModel(short tropModel, double Pa, double Td, double Pv, 
						 double elev, double Ht, double Bt, int doy,
						 double* dtrop, double *WZD, double *mf_wet)
{
	double trop_dry=0;
	switch (tropModel)
	{
	case TROPO_MAP_TYPE_SAAS_MIT:
		TropoBySaastamoinenWithMIT(Pa, Td, Pv, elev, Ht, Bt, dtrop,WZD, mf_wet);
		break;
	case TROPO_MAP_TYPE_SAAS_CHAO:
		TropoBySaastamoienWithChao(Pa, Td, Pv, elev, Ht, Bt, dtrop,WZD, mf_wet);
		break;	
	case TROPO_TYPE_HOPF:
		TropoByHopfield(Pa, Td, Pv, elev, dtrop,WZD,mf_wet);
		break;
	case TROPO_TYPE_BROWN:
		TropoByBrown(Pa, Td, Pv, elev, Ht, dtrop,WZD,mf_wet);    
		break;
	case TROPO_TYPE_LATM_HK:
		TropoBySaastamoinenWithMIT(Pa, Td, Pv, elev, Ht, Bt, dtrop,WZD,mf_wet);
        trop_dry=*dtrop-(*mf_wet)*(*WZD);
		*dtrop=TropoBy_LATM_HK(Ht,Bt,doy)*TropoDryMapFuncNiell(elev,Ht,b,doy);
		*WZD=(*dtrop-trop_dry)/(*mf_wet+EPS);
		break;
	default:
		TropoBySaastamoinenWithMIT(Pa, Td, Pv, elev, Ht, Bt, dtrop,WZD,mf_wet);
		break;
	}
}
//====================================================================
//Sasstamoinen模型来计算对流层改正    ***************
// 参数说明：Pa->大气压力；Td->大气温度；Pv->水汽压  *********
//           elev->卫星高度角                        
//           Ht->Height of the station
//           B->The latitude of the station 
//MIT映射系数
short TropoBySaastamoinenWithMIT(double Pa, double Td, double Pv,double elev, 
										  double Ht,double Bt, double *dtrop, double *WZD, double *mf_wet)
{
	double B = Bt * PI / 180.0;
	double fBh = 1.0 - 0.00266 * cos(2*B) - 0.00028 * Ht / 1000.0;

	//天顶方向的干延迟分量
	double rkDry = 0.002277 * Pa / fBh;
	//天顶方向的湿延迟分量
	double rkWet = 0.002277 * Pv / fBh * (1255/Td + 0.05);   //  很多文章都是用的这个公式

	double MF1 = 0, MF2 = 0;
	
	TropoDryMapFuncMIT(Td, Bt, Ht, elev, &MF1);
	TropoWetMapFuncMIT(Td, Bt, Ht, elev, &MF2);

	*dtrop = rkDry * MF1 + rkWet * MF2;
	*WZD = rkWet;
	*mf_wet=MF2;
	
	return 1;
}
//Chao映射系数
short TropoBySaastamoienWithChao(double Pa, double Td, double Pv,double elev, 
										  double Ht,double Bt, double *dtrop, double *WZD, double *mf_wet)
{
	double B = Bt * PI / 180.0;
	double fBh = 1.0 - 0.00266 * cos(2*B) - 0.00028 * Ht / 1000.0;

	//天顶方向的干延迟分量
	double rkDry = 0.002277 * Pa / fBh;
	//天顶方向的湿延迟分量
	double rkWet = 0.002277 * Pv / fBh * (1255/Td + 0.05);   //  很多文章都是用的这个公式

	double MF1 = 0, MF2 = 0;
	
	MF1 = TropoDryMapFuncChao(elev);
    MF2 = TropoWetMapFuncChao(elev);

	*dtrop = rkDry * MF1 + rkWet * MF2;
	*WZD = rkWet;
	*mf_wet=MF2;

	return 1;
}

//=============================================================
//计算干分量映射函数，依据温度和高度
//Td：温度, elew: 高度角, B: 纬度, Ht: 高度, drymap: 干分量映射系数
short TropoDryMapFuncMIT(double Td, double elev, double Bt, double Ht, double *drymap)
{
	double B = Bt * PI / 180.0;
	double cosB = cos(Bt);
	double Hkm = Ht / 1000.0;  //换成千米

	double Ac = 1.232E-3 + 0.01391E-3 * cosB - 0.02089E-3 * Hkm + 0.002154E-3 * (Td - 10.0);
	double Bc = 3.16116E-3 - 0.16004E-3 * cosB - 0.03306E-3 * Hkm + 0.002064E-3 * (Td - 10.0);
	double Cc = 71.24372E-3 - 4.29342E-3 * cosB - 0.14908E-3 * Hkm - 0.002098E-3 * (Td - 10.0);
	
	double sinE = sin(elev * PI / 180.0);
	double cosE = cos(elev * PI / 180.0);
	double Beta = Bc / (sinE + Cc);
	double Gamma = Ac / (sinE + Beta);
	double topcon = 1.0 + Ac / (1.0 + Bc/(1.0+Cc));
	
	*drymap = topcon / (sinE + Gamma);
	return 1;

}

//=============================================================
//计算干分量映射函数，依据温度和高度
//Td：温度, elew: 高度角, B: 纬度, Ht: 高度, wetmap: 湿分量映射系数
short TropoWetMapFuncMIT(double Td, double elev, double Bt, double Ht, double *wetmap)
{
	double B = Bt * PI / 180;
	double cosB = cos(B);
	double Hkm = Ht / 1000;  //换成千米

	double Ac = 0.58266E-3 - 0.01105E-3 * cosB - 0.05181E-3 * Hkm + 0.001442E-3 * (Td - 10.0);
	double Bc = 1.40218E-3 + 0.10249E-3 * cosB - 0.10128E-3 * Hkm + 0.002046E-3 * (Td - 10.0);
	double Cc = 45.85450E-3 - 1.91277E-3 * cosB - 1.28787E-3 * Hkm + 0.015136E-3 * (Td - 10.0);
	
	double sinE = sin(elev * PI / 180);
	double cosE = cos(elev * PI / 180);
	double Beta = Bc / (sinE + Cc);
	double Gamma = Ac / (sinE + Beta);
	double topcon = 1.0 + Ac / (1.0 + Bc/(1.0+Cc));
	
	*wetmap = topcon / (sinE + Gamma);
	
	return 1;
}

//yu:此函数方法欠佳，需要改写。直接根据前后历元的时间间隔来判断是否需要增加新模糊度cntAmb++
//  prn=(cntAmb-1)*50+prn 即可，无需回溯
//-
//yu:弧段间间隔：弧段末点处 依据弧段开始处的即本弧段前 (ephGapEnd-ephGapStart) > maxGap 来判断弧段是否属于新的模糊度 
//yu:历元间间隔：弧段中 依据前后历元 (ephEnd-ephend0) > maxGap 历元间隔是否太大来判断是否是新弧段
//此函数用于编辑并获得删除短数据后的双差系数集的模糊度信息获取
//对于断开时间长的数据，赋予新的模糊度，同时将卫星编号加50，相当于增加一个新卫星
BOOL EditSatAmb(CGPSDDCoefSet *pStSet, double maxGap, CSatAmbSet* pSatAmbSet)
{
	//获取双差观测方程系数阵里面的所有卫星数组（除去参考卫星）
	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //观测记录的个数
	short ambCount = 0; //卫星的模糊度数
	long gapCount = 0;  //记录数据断开的长度
	double ephStart=0, ephEnd0=0, ephEnd=0;
	double ephGapStart=0, ephGapEnd0=0, ephGapEnd=0;

	CGPSEpochDDCoef* pEpochDDcoef = NULL;
	CGPSEpochDDCoef* pPrevEpochDDcoef = NULL;
	CGPSEpochDDCoef* pDelDDcoef = NULL;
	CGPSEpochDDCoef* pEditDDcoef = NULL;
	CDDCoef* pObs = NULL;
    CSatAmb* pSatAmb = NULL;
	POSITION tmpPos = NULL;
	POSITION pPos = NULL;
	POSITION dPos = NULL;
 
    short newSatid = 0;
	for(i=0; i<satNum; i++)//yu:逐颗卫星建立模糊度信息
	{
		
		ambCount = 0;
		obsCount = 0;
		gapCount = 0;
		ephStart = 0;
		ephEnd = 0;
		ephEnd0 = 0;//yu:记录前次有数据的时间,即双差方程组前一元素时间
		ephGapStart = 0;
		ephGapEnd = 0;
		ephGapEnd0 = 0;

		pPos = pStSet->GetHeadPosition();

		//yu:所谓弧段末尾是同步数据的末尾7，即其他卫星有数据而此卫星无数据
		//   最大间隙是前后两个历元的数据间隔而言如bc，即同步数据历元间隔的限制或其他卫星无数据而此卫星有时的前后历元间隔限制
		//
		//-------------------------------------------------------------------
		//                  |  >maxGap |                      |>maxGap|
        //          |  <min |            |>maxGap |         
		//prn01:    * * * * *          * *        * * * * * * *       * * * *
		//prn##:* * * * * * * * *                     * * * * * * * * * * * * 
		//      + + + + + + + + +      + +        + + + + + + + + + + + + + + 有双差方程的历元
		//------------------------------------------------------------------- 等间隔时间
		//      1 2 3 4 5 6 7 8 9      a b        c d e f g h i j k l m n o p 
		//-------------------------------|ephEnd0-|ephEnd--------------------
		//-----------------------------|ephStart-----------------------------
		//-------ephGapStart| >maxGap  |ephGapend----------------------------
		//弧段间间隔新模糊度：i点判断此弧段与前一弧段间隙太大应加入新模糊度
        //历元间间隔新弧段随即新模糊度：c点数组连续但他是新弧段开始处（历元间间隔大），m点数组不连续是新弧段开始处

		while(pPos != NULL)//yu:遍历所有双差方程搜索此颗卫星
		{
            pEpochDDcoef = pStSet->GetNext(pPos);
			//获取指定卫星的观测数据
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);

			//yu:  3-7,a-i,m-p三个弧段
			if(pObs != NULL)  //如果该历元有次卫星观测数据，记录观测数据的个数
			{
				obsCount++;

				//-------记录观测数据的起止时间
				if(obsCount == 1)//yu: 3,a,m 即弧段起点
				{
					ephStart = pEpochDDcoef->m_GpsSecond;  //记录开始时间
					ephEnd = pEpochDDcoef->m_GpsSecond;  
					ephEnd0 = ephEnd;
				}
				else
				{
					ephEnd0 = ephEnd;//yu:前一历元时间
					ephEnd = pEpochDDcoef->m_GpsSecond;  //记录结束时间
				}

				if(pPos==NULL)  //已经到了最后一个历元 //yu: p
				{
					if(obsCount<=0){continue;  }
					
					//此前此卫星数据一直连续，只记录卫星模糊度信息
					if(ambCount==0)   //yu:从始至终一直连续，图上无
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = satPrn[i];       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数
						
						pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsStartSecond = ephStart;
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd;					
					}
					else  //已经存在模糊度信息 //yu: p
					{
						//yu:弧段末判断该段是否与前弧段间隔太大而属于新模糊度
						if((ephGapEnd-ephGapStart) > maxGap)  //断开时间过长，添加新的模糊度,将该段数据的卫星编号改变
						{
							ambCount++;
							//卫星号变为原来的  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = newSatid;       //卫星号
							pSatAmb->m_ambNum = ambCount;       //模糊度个数
							
							pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsStartSecond = ephStart;
							pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = ephEnd;	
						}
						else
						{
							newSatid = satPrn[i] + (ambCount-1) * 50;							
							pSatAmb = pSatAmbSet->GetSatAmb(newSatid,ambCount);
							
							pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = ephEnd;	
						}
						
						gapCount = 0;
						ephGapEnd = 0;
						ephGapStart = 0;
						ephGapEnd0 = 0;

						//将新的数据段的卫星编号改变
						if(ambCount != 1)
						{
							dPos = pStSet->GetTailPosition();
							
							for(j=obsCount; j>0; j--)
							{
								pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
								pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
								if (NULL == pObs)
										return 0;
								pObs->SatID = newSatid;
							}
						}
					}
					obsCount = 0;			
					ephStart = 0;
					ephEnd = 0;
					ephEnd0 = 0;
				}
				else if((ephEnd-ephEnd0) > maxGap)  //前后历元间隔太大  //yu:c点,是遇新弧段，再判断刚借宿的弧段的模糊度
				{ 
					//此卫星数据一直连续，只记录卫星模糊度信息
					if(ambCount==0)   
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = satPrn[i];       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数
						
						pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsStartSecond = ephStart;
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd0;	
						
						ephGapStart = ephEnd0;
						ephGapEnd = ephEnd;
						ephGapEnd0 = ephEnd;
					}
					else  //已经存在模糊度信息 //yu:c点
					{
                        //判断弧段是否属于新模糊度
						if((ephGapEnd-ephGapStart) > maxGap)  //断开时间过长，添加新的模糊度ambCount++,将该段数据的卫星编号改变
						{
							ambCount++;
							//卫星号变为原来的  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = newSatid;       //卫星号
							pSatAmb->m_ambNum = ambCount;       //模糊度个数
							
							pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsStartSecond = ephStart;
							pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = ephEnd;	
						}
						else
						{
							newSatid = satPrn[i] + (ambCount-1) * 50;							
							pSatAmb = pSatAmbSet->GetSatAmb(newSatid,ambCount);
							
							pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsEndSecond = ephEnd;	
						}
						
						gapCount = 0;
						ephGapEnd = 0;
						ephGapStart = 0;
						ephGapEnd0 = 0;
						
						//将新的数据段的卫星编号改变
						if(ambCount !=1 )
						{
							tmpPos = pPos;
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							
							dPos = tmpPos;
					
							for(j=obsCount; j>0; j--)
							{
								pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
								pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
								if (NULL == pObs)
										return 0;
								pObs->SatID = newSatid;
							}
						}
					}
					ephGapStart = ephEnd0;
					ephGapEnd = ephEnd;
					ephGapEnd0 = ephEnd;
					obsCount = 1;
					ephStart = ephEnd;				
				}
			}
			else   //(pObs==NULL)  //如果无该历元无此卫星观测数据，有以下情况 
			{
				//------------------
				//判断有观测数据记录
				if(obsCount == 0) //yu:1,2,8,9,j,k,l
				{ 					
					continue;
				}

				//此卫星数据一直连续，只记录卫星模糊度信息
				if(ambCount==0)   //yu:8
				{
					ambCount = 1;
					pSatAmb = pSatAmbSet->NewObjItem();
					pSatAmb->m_satID = satPrn[i];       //卫星号
					pSatAmb->m_ambNum = ambCount;       //模糊度个数
					
					pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
					pSatAmb->m_GpsStartSecond = ephStart;	//yu:弧段首
					pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
					pSatAmb->m_GpsEndSecond = ephEnd;	//yu:前一历元为弧段结束点
					
					ephGapStart = pEpochDDcoef->m_GpsSecond;//yu:中断开始处
					ephGapEnd = pEpochDDcoef->m_GpsSecond;
					ephGapEnd0 = pEpochDDcoef->m_GpsSecond;
				}
				else  //已经存在模糊度信息    //yu:j
				{
					if((ephGapEnd-ephGapStart) > maxGap)  //断开时间过长，添加新的模糊度,将该段数据的卫星编号改变
					{
						ambCount++;
						//卫星号变为原来的  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = newSatid;       //卫星号
						pSatAmb->m_ambNum = ambCount;       //模糊度个数
						
						pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsStartSecond = ephStart;
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd;	
					}
					else//yu:k,l
					{
						newSatid = satPrn[i] + (ambCount-1) * 50;							
						pSatAmb = pSatAmbSet->GetSatAmb(newSatid,ambCount);
						
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd;	
					}
					
					gapCount = 0;
					ephGapEnd = 0;
					ephGapStart = 0;
					ephGapEnd0 = 0;

					//------------------
					gapCount++;
					
					if(gapCount==1)
					{
						ephGapStart = pEpochDDcoef->m_GpsSecond;  //记录断开的开始时间
						ephGapEnd = pEpochDDcoef->m_GpsSecond;  
						ephGapEnd0 = ephGapEnd;
					}
					else
					{
						ephGapEnd0 = ephGapEnd;                   //记录断开的结束时间
						ephGapEnd = pEpochDDcoef->m_GpsSecond;
					}

					//将新的数据段的卫星编号改变
					if(ambCount !=1 )
					{
						if(pPos==NULL)
						{
							dPos = pStSet->GetTailPosition();
						}
						else
						{
							tmpPos = pPos;
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							dPos = tmpPos;
						}
						
						for(j=obsCount; j>0; j--)
						{
							pEditDDcoef = pStSet->GetPrev(dPos); //当前记录的前一个记录
							pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
							if (NULL == pObs)
								return 0;
							pObs->SatID = newSatid;
						}
					}
				}

				obsCount = 0;
				ephStart = 0;
				ephEnd = 0;
				ephEnd0 = 0;				
			}
		}
	}
	return 1;
}
					
//==========================================
// tr: 信号接收时刻
// pEph: 卫星星历数据
// Xr, Yr, Zr：接收机的近似位置
// Xk, Yk, Zk, Xdot, Ydot, Zdot: 卫星位置、速度，平近点角
//由接收机的近似位置，计算卫星的发射时刻
double GetSatSentTime(double tr, ephemeris pEph, 
					  double Xr, double Yr, double Zr,
					  double* Xk, double* Yk, double* Zk, double* Ek,
					  double* Xdot, double* Ydot, double* Zdot, double* range)
{
	double ts0 = 0;
	double ts = 0;  //信号发射时刻
	double deltT = 0;
	double satDt = 0; //卫星钟差

    ts = tr - 0.075;

	do{
		//根据卫星发射时间钟面时，计算卫星钟差
		deltT = ts - pEph.toe;
		if(deltT>302400.00) //考虑星期交叠的开始和结束时间
		{
			deltT -= 604800.0;
		}
		if(deltT<-302400.0) 
		{
			deltT += 604800.0;
		}
		
		//卫星钟差
		satDt = pEph.a0 + (pEph.a1 * deltT) + (pEph.a2 * deltT * deltT);
        
		//卫星发射时间
		ts = ts + satDt;
		
		//保存卫星发射时间
		ts0 = ts;

		//计算卫星位置
		coordsat(pEph, ts0, Xk, Yk, Zk, Ek, Xdot, Ydot, Zdot);
		*range = sqrt((*Xk-Xr)*(*Xk-Xr)+(*Yk-Yr)*(*Yk-Yr)+(*Zk-Zr)*(*Zk-Zr));
		
		//重新计算计算卫星信号发射时间
		ts = tr - *range / C;
	}while((ts-ts0)>1.0E-7);  //当两次计算的发射时间小于10-7方时，迭代结束
	
	return ts;
}

//-----------------------copied by wkyu 20110412 from zxw's gpssm----------------------------
					
GPSCLASSDLL_API  short FilterAmbVectorSF(double** AmbVector, long *VectorCount,double* Amb,short AmbCount, Matrix Dxx)
{
	long i=0;
	long j=0;
	long k=0;
	double deltaN=0, deltan=0;
	///////////////////////////////////////////////////
	long count=0;
	short factor=3;
	
	//对模糊度进行判断
	/*do
	{
		count = *VectorCount;
		for(i=0; i<*VectorCount; i++)
		{
			short mark = 0;
			for(j=0; j<AmbCount; j++)
			{
				if(mark == 1) break;
				for(k=j+1; k<AmbCount; k++)
				{
					deltan = Amb[j] - Amb[k];
					deltaN = AmbVector[i][j] - AmbVector[i][k];
					double m0 = factor * sqrt(fabs(Dxx(j,j)+Dxx(k,k)-2*Dxx(j,k)));   //modified by zxw 2007/1/20
					
					if(deltaN<(deltan-m0) ||deltaN>(deltan+m0))  
					{
						count--;
						mark=1;
						break;
					}
				}
			}
		}
		if(count<1) factor++;
	}while(count<1);
	*/
	///////////////////////////////////////////////////
    count = *VectorCount;
	for(i=0;i<*VectorCount;i++)
	{
		short mark = 0;
		for(j=0; j<AmbCount-1; j++)
		{
			if(mark == 1)
			{
				break;
			}
			
			for(k=j+1; k<AmbCount; k++)
			{
				deltan = Amb[j] - Amb[k];
				deltaN = AmbVector[i][j] - AmbVector[i][k];
				double m0 = factor * sqrt(fabs(Dxx(j,j)+Dxx(k,k)-2*Dxx(j,k)));
				
				double minDelta = Round((deltan-m0),-1);
				double maxDelta = Round((deltan+m0),1);				
				if(deltaN<minDelta || deltaN>maxDelta)
				//if(deltaN<(deltan-m0) || deltaN>(deltan+m0) )
				{
					delete[] AmbVector[i];
					AmbVector[i]=NULL;
					count--;
					mark=1;
					break;
				}
			}
		}
	}
    
	//对模糊度向量组AmbVector进行整理	
	long delVectorCount=0;
	for(i=0;i<*VectorCount-delVectorCount;i++)
	{
		if(AmbVector[i]==NULL)
		{
			delVectorCount++;
			for(j=i;j<((*VectorCount)-delVectorCount);j++)
			{
				AmbVector[j]=AmbVector[j+1];
			}
			i--;
			AmbVector[*VectorCount-delVectorCount]=NULL;
		}
	}
	*VectorCount=count;
	return 1;
}
//==================================================================
//函数功能:用以下规则过滤模糊度向量，以便加速模糊度的搜索
//规则为：如果同一对卫星的L1,L2模糊度向量中的两个模糊度之差N_JK不满足以下条件
//则从备选模糊度向量中删除该模糊度向量
//N_jk-t(a/2)*mN_jk<=N_JK<=N_jk+t(a/2)*mN_jk
//N_jk ＝ XN1 - (f1/f2)*XN2:
//t(a/2):置信度为（1-a）的t分布值
//mN_jk:整周模糊度的实数解的差值的验后方差
//参数说明:AmbVector:备选的模糊度向量组
//         VectorCount:备选的模糊度向量组的个数
//         Amb:浮点解模糊度
//         AmbCount:浮点解模糊度个数
//         Dxx:浮点解方差阵
//========================================================
GPSCLASSDLL_API 	short FilterAmbVectorDF(double** AmbVector, long *VectorCount,double* Amb, short AmbCount, Matrix Dx)
{
	long i=0;
	long j=0;
	double deltaN=0, deltan=0;
	///////////////////////////////////////////////////
	long count=0;
	short factor=3;
	short k = AmbCount/2;
	//对模糊度进行判断
	/*do
	{
		count = *VectorCount;
		for(i=0; i<*VectorCount; i++)
		{
			short mark = 0;
			for(j=0; j<AmbCount/2; j++)
			{
				deltan = Amb[j] - (77.0/60.0) * Amb[j+k];   //xn1 - f1/f2 *xn2  //浮点模糊度
				deltaN = AmbVector[i][j] - (77.0/60.0) * AmbVector[i][j+k];
				
				double m0 = factor * sqrt(fabs(Dx(j,j)+Dx(j+k,j+k)-2*Dx(j,j+k))); 
				double minDelta = Round((deltan-m0),-1);
				double maxDelta = Round((deltan+m0),1);
				if(deltaN<minDelta || deltaN>maxDelta)  
				{
					count--;
					mark=1;
					break;
				}
			}
		}
		if(count<1) factor++;
	}while(count<1); 
	*/
	///////////////////////////////////////////////////
    count = *VectorCount;
	for(i=0;i<*VectorCount;i++)
	{
		short mark = 0;
		for(j=0; j<AmbCount/2; j++)
		{
			if(mark == 1)
			{
				break;
			}
			deltan = Amb[j] - (77.0/60.0) * Amb[j+k];   //xn1 - f1/f2 *xn2  //浮点模糊度
			deltaN = AmbVector[i][j] - (77.0/60.0) * AmbVector[i][j+k];
			
			double m0 = factor * sqrt(fabs(Dx(j,j)+Dx(j+k,j+k)-2*Dx(j,j+k))); 
			
			double minDelta = Round((deltan-m0),-1);
			double maxDelta = Round((deltan+m0),1);
			if(deltaN<minDelta || deltaN>maxDelta)   
			{
				delete [] AmbVector[i];
				AmbVector[i]=NULL;
				count--;
				mark=1;
				break;
			}
		}
	}
    
	//对模糊度向量组AmbVector进行整理
	
	long delVectorCount=0;
	for(i=0;i<*VectorCount-delVectorCount;i++)
	{
		if(AmbVector[i]==NULL)
		{
			delVectorCount++;
			for(j=i;j<((*VectorCount)-delVectorCount);j++)
			{
				AmbVector[j]=AmbVector[j+1];
			}
			i--;
			AmbVector[*VectorCount-delVectorCount]=NULL;
		}
	}
	*VectorCount=count;
	return 1;
}	

//================================================================
// 函数功能:快速模糊度解算(Kalman双频模糊度解算) 
// 参数说明:in&out:Amb1, Amb2,->输入模糊度浮点解，输出模糊度固定解
//          in:sigma0->浮点解中误差            **********
//          in:AmbCount-> 模糊度的个数             ********
//          in:Dx1, Dx2-> 浮点解的协方差阵
//          out:ratio->ratio值                    ********
//================================================================
GPSCLASSDLL_API short FixAmbFARADF(double *Amb, double sigma0, short AmbCount, Matrix Dx, long DfCount, double* ratio)
{
	//找出所有可能的模糊度向量
	long i=0;
	double** AmbVector=NULL;
	long VectorCount=0;
    
    short rn = GetAmbVectorDF(Amb, AmbCount, NULL, &VectorCount, Dx); //首先计算所有模糊度量个数
	if(rn==0) return 0;
	if(VectorCount>2000000 || VectorCount<=0 ) return 0;
	
	AmbVector = new double* [VectorCount];
	
	for(i=0; i<VectorCount; i++)
	{
		AmbVector[i]=new double[AmbCount];
	}

	//获得所有的模糊度向量组
	GetAmbVectorDF(Amb, AmbCount, AmbVector, &VectorCount, Dx);
	
	/*if(VectorCount>1)
	{	//过滤向量组
		FilterAmbVectorSF(AmbVector, &VectorCount, Amb, AmbCount, Dx);
	}*/

	double* sigma=new double[VectorCount];
	for(i=0;i<VectorCount;i++)
	{
		sigma[i] = GetSigma(Amb, sigma0, AmbVector[i], AmbCount, Dx, DfCount); //计算每一组备选模糊度的单位权中误差
	}
   
	if(VectorCount>1)
	{
	
		long MinVector=0;
		double first=1000;
		double second=1000;
		for(i=0;i<VectorCount;i++)
		{
			if(sigma[i]<first) 
			{
				first=sigma[i];
				MinVector=i;
			}
		}
		for(i=0;i<VectorCount;i++)
		{
			if(sigma[i]<second && i!=MinVector) 
			{
				second=sigma[i];
			}
		}
		*ratio = second/first;
		for(i=0;i<AmbCount;i++)
		{
			Amb[i]=AmbVector[MinVector][i];
		}
	}
	else if(VectorCount==1)
	{
        for(i=0;i<AmbCount;i++)
		{
			Amb[i]=AmbVector[0][i];
		}
		*ratio=1000;
	}
	else
	{
		for(i=0;i<VectorCount;i++)
		{
			if(AmbVector[i]!=NULL)
			{
				delete[] AmbVector[i];
			}
		}
		if(AmbVector!=NULL)
		{
			delete[] AmbVector;
		}
		
		if(sigma!=NULL)
		{
			delete[] sigma;
		}

		return 0;
	}
	for(i=0;i<VectorCount;i++)
	{
		if(AmbVector[i]!=NULL)
		{
		   delete[] AmbVector[i];
		}
	}
	if(AmbVector!=NULL)
	{
	  delete[] AmbVector;
	}

	if(sigma!=NULL)
	{
	  delete[] sigma;
	}
	return 1;
}

GPSCLASSDLL_API short FixAmbFARASF(double *Amb, double sigma0, short AmbCount, Matrix Dx, long DfCount, double* ratio)
{
	//找出所有可能的模糊度向量
	long i=0;
	double** AmbVector=NULL;
	long VectorCount=0;
    
    short rn = GetAmbVectorSF(Amb, AmbCount, NULL, &VectorCount, Dx); //首先计算所有模糊度量个数
	if(rn==0) return 0;
	if(VectorCount>2000000 || VectorCount<=0 ) return 0;
	
	AmbVector = new double* [VectorCount];
	
	for(i=0; i<VectorCount; i++)
	{
		AmbVector[i]=new double[AmbCount];
	}

	//获得所有的模糊度向量组
	GetAmbVectorSF(Amb, AmbCount, AmbVector, &VectorCount, Dx);
	
	if(VectorCount>10000)
	{	//过滤向量组
		FilterAmbVectorSF(AmbVector, &VectorCount, Amb, AmbCount, Dx);
	}

	double* sigma=new double[VectorCount];
	for(i=0;i<VectorCount;i++)
	{
		sigma[i] = GetSigma(Amb, sigma0, AmbVector[i], AmbCount, Dx, DfCount); //计算每一组备选模糊度的单位权中误差
	}
   
	if(VectorCount>1)
	{
	
		long MinVector=0;
		double first=1000;
		double second=1000;
		for(i=0;i<VectorCount;i++)
		{
			if(sigma[i]<first) 
			{
				first=sigma[i];
				MinVector=i;
			}
		}
		for(i=0;i<VectorCount;i++)
		{
			if(sigma[i]<second && i!=MinVector) 
			{
				second=sigma[i];
			}
		}
		*ratio = second/first;
		for(i=0;i<AmbCount;i++)
		{
			Amb[i]=AmbVector[MinVector][i];
		}
	}
	else if(VectorCount==1)
	{
        for(i=0;i<AmbCount;i++)
		{
			Amb[i]=AmbVector[0][i];
		}
		*ratio=1000;
	}
	else
	{
		for(i=0;i<VectorCount;i++)
		{
			if(AmbVector[i]!=NULL)
			{
				delete[] AmbVector[i];
			}
		}
		if(AmbVector!=NULL)
		{
			delete[] AmbVector;
		}
		
		if(sigma!=NULL)
		{
			delete[] sigma;
		}

		return 0;
	}
	for(i=0;i<VectorCount;i++)
	{
		if(AmbVector[i]!=NULL)
		{
		   delete[] AmbVector[i];
		}
	}
	if(AmbVector!=NULL)
	{
	  delete[] AmbVector;
	}

	if(sigma!=NULL)
	{
	  delete[] sigma;
	}
	return 1;
}

//用区间法获取模糊度向量组
/*in*/
//Amb:浮点解模糊度
//AmbCount:浮点解模糊度个数
//Qxx  :浮点解模糊度方差
/*out*/
//AmbVector:备选的模糊度向量组
//VectorCount:备选的模糊度向量组的个数
//(注：调用该函数时，先使参数AmbVector为NULL，从而获得VectorCount,即模糊度向量的个数
// 然后根据VectorCount的大小为AmbVector分配内存，之后在调用该函数 ）
GPSCLASSDLL_API short GetAmbVectorSF(double *Amb,short AmbCount,double** AmbVector,long* VectorCount,Matrix Dxx)
{
	long i=0, j=0, k=0;
	double  amb1=0,amb2=0;
	short factor=4;
	double sigma=0;
	double deltaN=0;
	
	short *Ni=new short[AmbCount];//每个模糊度的备选模糊度个数
	double **Ambiguitys=NULL;//每个模糊度的备选模糊度
	Ambiguitys=new double* [AmbCount];
	for(i=0;i<AmbCount;i++) 
	{
		Ni[i]=0;
	    Ambiguitys[i]=NULL;
	}

	for(i=0; i<AmbCount; i++)
	{			 
		deltaN = factor*sqrt(Dxx(i,i));
		amb1 = Round(Amb[i]-deltaN, -1);
		amb2 = Round(Amb[i]+deltaN, 1);
		Ambiguitys[i]=new double[amb2-amb1+1];
		short count=0;		
		
		for(k=amb1; k<=amb2; k++)
		{
			Ambiguitys[i][count] = k;
			count++;
		}
		Ni[i] = count;
	 }
	
	//对各种可能的模糊度进行组合，得出所有的模糊度向量
    long cpn=1;
	long ii=0;
	*VectorCount=1;
	for(i=0;i<AmbCount;i++)
	{
		(*VectorCount)=(*VectorCount)*Ni[i];//得出所有的模糊度向量的个数
	}
    if(AmbVector==NULL) 
	{
		if(Ni!=NULL)
		{
			delete[] Ni;
			Ni = NULL;
		}
		for(i=0;i<AmbCount;i++)
		{
			if(Ambiguitys[i]!=NULL)
			{
				delete[] Ambiguitys[i];
				Ambiguitys[i] = NULL;
			}
		}
		if(Ambiguitys!=NULL) 
		{
			delete[] Ambiguitys;
			Ambiguitys = NULL;
		}
		return 1;
	}

	for(i=0;i<*VectorCount;i++)
	{
		cpn=*VectorCount;
		ii=i;
		for(j=0;j<AmbCount;j++)
		{
			cpn=cpn/Ni[j];
			k=ii/cpn;
			AmbVector[i][j]=Ambiguitys[j][k];
			ii=ii-cpn*k;
		}
	}

	//释放内存
	if(Ni!=NULL)
	{
	   delete[] Ni;
	   Ni = NULL;
	}
	for(i=0;i<AmbCount;i++)
	{
		if(Ambiguitys[i]!=NULL)
		{
			delete[] Ambiguitys[i];
			Ambiguitys[i]=NULL;
		}
	}
	if(Ambiguitys!=NULL) 
	{
		delete[] Ambiguitys;
		Ambiguitys=NULL;
	}
	return 1;
}
//用区间法获取模糊度向量组
/*in*/
//Amb:浮点解模糊度
//AmbCount:浮点解模糊度个数
//Qxx  :浮点解模糊度方差
/*out*/
//AmbVector:备选的模糊度向量组
//VectorCount:备选的模糊度向量组的个数
//(注：调用该函数时，先使参数AmbVector为NULL，从而获得VectorCount,即模糊度向量的个数
// 然后根据VectorCount的大小为AmbVector分配内存，之后在调用该函数 ）
GPSCLASSDLL_API short GetAmbVectorDF(double *Amb,short AmbCount,
				                	 double** AmbVector,long* VectorCount,Matrix Qxx)
{
	short  i=0;
	long  j=0;
	long  amb11=0, amb12=0;
	long  amb21=0, amb22=0;
	short factor=4;
	long  k=0;
	short nSat = AmbCount/2;
	
	short *Ni=new short[nSat];//每个卫星的备选模糊度个数	
	
	typedef struct tagDualAmb
	{
		double AmbV1;
		double AmbV2;
	}DualAmb;

	DualAmb **SatAmb = new DualAmb * [nSat];  //模糊度的备选模糊度
	
    
	for(i=0; i<nSat; i++)
	{
		Ni[i] = 0;
		SatAmb[i] = NULL;
	}

	double deltan=0, deltaN=0;  
    long Nl1=0, Nl2=0;
    double m0=0, minDelta=0, maxDelta=0;
	
	for(i=0; i<nSat; i++)
	{
		//L1模糊度
		amb11 = Round(Amb[i]-factor*sqrt(Qxx(i,i)), -1);
		amb12 = Round(Amb[i]+factor*sqrt(Qxx(i,i)), 1);
		//L2模糊度
		amb21 = Round(Amb[i+nSat]-factor*sqrt(Qxx(i+nSat,i+nSat)), -1);
		amb22 = Round(Amb[i+nSat]+factor*sqrt(Qxx(i+nSat,i+nSat)), 1);
        //过滤L1,和L2的模糊度
		deltan = Amb[i] - (77.0/60.0) * Amb[i+nSat];   //xn1 - f1/f2 *xn2  //浮点模糊度
	    
		m0 = factor * sqrt(fabs(Qxx(i,i)+Qxx(i+nSat,i+nSat)-2*Qxx(i,i+nSat))); 
		minDelta = Round((deltan-m0),-1);
		maxDelta = Round((deltan+m0),1);
		
		short count=0;
        
		for(Nl1=amb11; Nl1<=amb12; Nl1++)    //L1可能的模糊度值
		{
			for(Nl2=amb21; Nl2<=amb22; Nl2++)  //L2可能的模糊度值
			{
				deltaN = Nl1 - (77.0/60.0) * Nl2;
				
				if(deltaN >= minDelta && deltaN <= maxDelta)  //此组模糊度符合备选模糊度条件,保留
				{
					count++;
					Ni[i] = count;   //第i个卫星的组合个数				
				}
			}
		}
		
		//记录卫星模糊度
		SatAmb[i] = new DualAmb[count];
		short ijk=0;
		for(Nl1=amb11; Nl1<=amb12; Nl1++)    //L1可能的模糊度值
		{
			for(Nl2=amb21; Nl2<=amb22; Nl2++)  //L2可能的模糊度值
			{
				deltaN = Nl1 - (77.0/60.0) * Nl2;
				if(deltaN >= minDelta && deltaN <= maxDelta)  //此组模糊度符合备选模糊度条件,保留
				{
					SatAmb[i][ijk].AmbV1 = Nl1;
					SatAmb[i][ijk].AmbV2 = Nl2;
					ijk++;
				}
			}
		}
	}

	//对各种可能的模糊度进行组合，得出所有的模糊度向量
    *VectorCount=1;
	
	for(i=0; i<nSat; i++)
	{
		(*VectorCount) = (*VectorCount)*Ni[i];//得出所有的模糊度向量的个数
	}
	
	if(AmbVector==NULL) 
	{
		if(Ni != NULL)
		{
			delete[] Ni;
			Ni = NULL;
		}

		for(i=0; i<nSat; i++)
		{
			if(SatAmb[i] != NULL)
			{
				delete[] SatAmb[i];
				SatAmb[i] = NULL;
			}
		}
		if(SatAmb != NULL) 
		{
			delete[] SatAmb;
			SatAmb = NULL;
		}
		
		return 1;
	}
	
	long cpn=1;
	long ii=0;
	//将模糊度填入向量组
	if(AmbVector != NULL)
	{
		for(j=0; j<*VectorCount; j++)
		{
			cpn = *VectorCount;
			ii = j;
			for(i=0; i<nSat; i++)
			{
				cpn = cpn/Ni[i];
				k = ii/cpn;
				AmbVector[j][i] = SatAmb[i][k].AmbV1;
				AmbVector[j][i+nSat] = SatAmb[i][k].AmbV2;
				ii=ii-cpn*k;
			}
		}
	}
	///////////////
	//释放内存
	if(Ni!=NULL)
	{
		delete [] Ni;
		Ni = NULL;
	}
	for(i=0; i<nSat; i++)
	{
		if(SatAmb[i]!=NULL)
		{
			delete[] SatAmb[i];
			SatAmb[i] =  NULL;
		}
	}
	if(SatAmb != NULL) 
	{
		delete[] SatAmb;
		SatAmb = NULL;
	}
	return 1;
}
//-
//获得双差的协方差阵
int GetDDCoVariance(CGPSEpochObsSet *pFEpochObs, CGPSEpochObsSet *pUEpochObs, 
					CGPSEpochDDCoef *pEpochDDCoef, 
					Matrix &Dx, CCalcSolutionBase *pSolutinBase)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	Matrix Dl(nSatCount, nSatCount); //方差协方差阵

	switch (pSolutinBase->m_WeightModel)
	{
	case WEIGHT_WITH_SIGNAL_STRENGTH:  //信号强度求方差
		{
			GetDDCoVarianceSignalStrength(pFEpochObs, pUEpochObs, pEpochDDCoef, Dl);
			Dx = Dl;
		}
		break;
	case WEIGHT_WITH_ELEV:				//高度角定权
		{
			GetDDCoVarianceElev(pFEpochObs, pUEpochObs, pEpochDDCoef, Dl);
			Dx = Dl;
		}
		break;
	/*case WEIGHT_WITH_CNO:				//Cno定权
		{
			WeightWithCNo(pUEpochObs,pFEpochObs, EpochDDCoef,refsat, P, satcount); 
		}
		break;
	case WEIGHT_WITH_STRENGTH_DELTA:				//信噪比定权
		{
			CNoTempl FcnoTempl;
			FcnoTempl.minCNo = pSolutinBase->m_FminCNo;
			FcnoTempl.maxCNo = pSolutinBase->m_FmaxCNo;
			FcnoTempl.midElev = pSolutinBase->m_FmidElev;
			CNoTempl UcnoTempl;
			UcnoTempl.minCNo = pSolutinBase->m_UminCNo;
			UcnoTempl.maxCNo = pSolutinBase->m_UmaxCNo;
			UcnoTempl.midElev = pSolutinBase->m_UmidElev;
			WeightWithStrengthDelta(pUEpochObs,pFEpochObs, EpochDDCoef,refsat, P, satcount,
				FcnoTempl, UcnoTempl); 
		}
		break;
	case WEIGHT_WITH_SIGMA_DELTA:
		{
			CNoTempl FcnoTempl;
			FcnoTempl.minCNo = pSolutinBase->m_FminCNo;
			FcnoTempl.maxCNo = pSolutinBase->m_FmaxCNo;
			FcnoTempl.midElev = pSolutinBase->m_FmidElev;
			CNoTempl UcnoTempl;
			UcnoTempl.minCNo = pSolutinBase->m_UminCNo;
			UcnoTempl.maxCNo = pSolutinBase->m_UmaxCNo;
			UcnoTempl.midElev = pSolutinBase->m_UmidElev;
			WeightWithSigmaDelta(pUEpochObs,pFEpochObs, EpochDDCoef,refsat, P, satcount,
		                       FcnoTempl, UcnoTempl);
		}break;
	case WEIGHT_WITH_COMBINED_METHOD:
		WeightWithCombinedMethod(pUEpochObs,pFEpochObs, EpochDDCoef,refsat, P, satcount); 
		break;
		*/
	case WEIGHT_WITH_UNIT:
		{
			GetDDCoVarianceUnit(pFEpochObs, pUEpochObs, pEpochDDCoef, Dl);
			Dx = Dl;
		}
		break;	
	default:
		return 0;
	}	
	return 1;
}
//-
//利用单位权模型获得双差观测值的方差
void GetDDCoVarianceUnit(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	short i=0, j=0;
	double m02 = 0.01 * 0.01;  //观测值精度给1cm
	Matrix Ql(nSatCount,nSatCount); //协因数阵

	for(i=0; i<nSatCount; i++)
	{
		for(j=i; j<nSatCount; j++)
		{
			if(i==j)
			{
				Ql(i,j) = 4.0;
			}
			else
			{
				Ql(i,j) = 2.0;
				Ql(j,i) = 2.0;
			}
		}
	}	
	Dll = Ql * m02;  //转换成协方差阵
}
//-
//利用高度角度模型获得双差观测值的方差
void GetDDCoVarianceElev(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	double m02 = 0.01 * 0.01;  //
	short j=0, k=0;
	Matrix Ql(nSatCount, nSatCount); //协因数阵


	for(j=0; j<nSatCount; j++)
	{
		CDDCoef * pDDCoef = pEpochDDCoef->GetAt(j);
		for(k=j; k<nSatCount; k++)
		{
			if(j == k)
			{
				Ql(j,k)= 1.0/( sin(pDDCoef->SatElev*rad) * sin(pDDCoef->SatElev*rad))+
					     1.0/( sin(pDDCoef->BaseSatElev*rad) * sin(pDDCoef->BaseSatElev*rad)) +
						 1.0/( sin(pDDCoef->FBaseSatElev*rad) * sin(pDDCoef->FBaseSatElev*rad)) + 
						 1.0/( sin(pDDCoef->FSatElev*rad) * sin(pDDCoef->FSatElev*rad));
			}
			else
			{
				Ql(j,k)= 1.0/(sin(pDDCoef->BaseSatElev*rad)* sin(pDDCoef->BaseSatElev*rad)) + 
					     1.0/(sin(pDDCoef->FBaseSatElev*rad)* sin(pDDCoef->FBaseSatElev*rad));
				Ql(k,j) = Ql(j,k);
			}
		}
	}	
	Dll = Ql * m02;  //转换成协方差阵
}
//-
//利用信号强度模型获得双差观测值的方差
void GetDDCoVarianceSignalStrength(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short satCount = pEpochDDCoef->GetItemCount();
	Matrix Dl(satCount,satCount);
	double CL1=0.224; //0.224相当于2.7mm的精度
	short refsat = pEpochDDCoef->m_BaseSat;
	CObservation *pObs=NULL;
    double URefSL1=0,FRefSL1=0,USL1=0,FSL1=0;
	
	pObs=pUEpochObs->GetObs(refsat);
	URefSL1=pObs->sL1;
    pObs=pFEpochObs->GetObs(refsat);
	FRefSL1=pObs->sL1;
    
	double URefSigma=CL1*pow(10,-URefSL1/2.0);
	double FRefSigma=CL1*pow(10,-FRefSL1/2.0);
	double USigma=0,FSigma=0;
	
	POSITION pos = pEpochDDCoef->GetHeadPosition();
	CDDCoef * pDDCoef = NULL;
	int k = 0, m=0;
	while (NULL != pos && k < satCount)
	{
		pDDCoef = pEpochDDCoef->GetNext(pos);
		pObs=pUEpochObs->GetObs(pDDCoef->SatID);
		USL1=pObs->sL1;
		pObs=pFEpochObs->GetObs(pDDCoef->SatID);
		FSL1=pObs->sL1;
		USigma=CL1*pow(10,-USL1/2.0);
		FSigma=CL1*pow(10,-FSL1/2.0);
		
		for(m=0;m<satCount;m++)
		{
			if(k==m)
			{
				Dl(k,m)=URefSigma+FRefSigma+USigma+FSigma;
			}
			else
			{
				Dl(k,m)=URefSigma+FRefSigma;
			}
		}
		k++;
	}
	Dll = Dl;
}
//-
//站心(S)空间直角坐标向地心(E)空间直角坐标转换矩阵
void SxyztoExyzMatrix(double B, double L, Matrix &tranMatrix)
{
	tranMatrix(0,0) = -cos(L)*sin(B); tranMatrix(0,1) = -sin(L); tranMatrix(0,2) = cos(L)*cos(B);
	tranMatrix(1,0) = -sin(L)*sin(B); tranMatrix(1,1) = cos(L);  tranMatrix(1,2) = sin(L)*cos(B);
	tranMatrix(2,0) = cos(B);         tranMatrix(2,1) = 0;       tranMatrix(2,2) = sin(B);
}
//-
//地心(E)空间直角坐标向站心(S)空间坐标转换矩阵
void ExyztoSxyzMatrix(double B, double L, Matrix &tranMatrix)
{
	tranMatrix(0,0) = -cos(L)*sin(B); tranMatrix(0,1) = -sin(L)*sin(B); tranMatrix(0,2) = cos(B);
	tranMatrix(1,0) = -sin(L);        tranMatrix(1,1) = cos(L);         tranMatrix(1,2) = 0;
	tranMatrix(2,0) = cos(B)*cos(L);  tranMatrix(2,1) = cos(B)*sin(L);  tranMatrix(2,2) = sin(B);
}

//==========================================
// tr: 信号接收时刻
// pEph: 卫星星历数据
// Xr, Yr, Zr：接收机的近似位置
// Xk, Yk, Zk, Xdot, Ydot, Zdot: 卫星位置、速度，平近点角
//由接收机的近似位置，计算卫星的发射时刻
double GetSatSentTime(double tr, ephemeris *pEph, 
					  double Xr, double Yr, double Zr,
					  double* Xk, double* Yk, double* Zk, double* Ek,
					  double* Xdot, double* Ydot, double* Zdot, double* range)
{
	double ts0 = 0;
	double ts = 0;  //信号发射时刻
	double deltT = 0;
	double satDt = 0; //卫星钟差
    double dift = 0;

    ts = tr - 0.075;

	do{
		//根据卫星发射时间钟面时，计算卫星钟差
		deltT = ts - pEph->toe;
		if(deltT>302400.00) //考虑星期交叠的开始和结束时间
		{
			deltT -= 604800.0;
		}
		if(deltT<-302400.0) 
		{
			deltT += 604800.0;
		}
		
		//卫星钟差
		satDt = pEph->a0 + (pEph->a1 * deltT) + (pEph->a2 * deltT * deltT);
        
		//卫星发射时间
		ts0 = ts + satDt;
		
		//计算卫星位置
		coordsat(*pEph, ts0, Xk, Yk, Zk, Ek, Xdot, Ydot, Zdot);
		*range = sqrt((*Xk-Xr)*(*Xk-Xr)+(*Yk-Yr)*(*Yk-Yr)+(*Zk-Zr)*(*Zk-Zr));
		
		//重新计算计算卫星信号发射时间
		ts = tr - *range / C;
		dift = (ts+satDt) - ts0;
	}while(fabs(dift)>1.0E-7);  //当两次计算的发射时间小于10-7方时，迭代结束
	
	return ts;
}

void PosToPos(double antN, double antE, double antH, 
		          double B, double L, 
				  double X, double Y, double Z,
				  double &newX, double &newY, double &newZ, short type)
{
	Matrix refTran(3,3); //转换矩阵(站心－>空间)
	Matrix antfNEH(3,1); //天线偏移量(站心)
	Matrix antfXYZ(3,1);  //天线偏移量(空间)
	
	//将天线偏移量由站心转到地心
	SxyztoExyzMatrix(B, L, refTran);

	antfNEH(0,0) = antN;
	antfNEH(1,0) = antE;
	antfNEH(2,0) = antH;
	
	antfXYZ = refTran*antfNEH;  //空间直角坐标下的天线改正量

	if( type== 0 )
	{
		newX = X + antfXYZ(0,0);
		newY = Y + antfXYZ(1,0);
		newZ = Z + antfXYZ(2,0);
	}
	else
	{
		newX = X - antfXYZ(0,0);
		newY = Y - antfXYZ(1,0);
		newZ = Z - antfXYZ(2,0);
	}
}

//函数功能：卡方分布
//根据观测值个数及置信水平，返回卡方极限值

double ChiQuasiThresold(short n, double alfa)
{
	double Threshold = 0;
	/*double ChiQus[5][25]={{0, 0.01, 0.072, 0.207,0.412, 0.676, 0.989, 1.344, 1.735, 2.156, 2.603, 3.074,
							3.565,4.075,4.601,5.142,5.697,6.265,6.844,7.734,8.034,8.643,9.260,9.886,10.520},    // 0.995
							{0, 0.02, 0.115, 0.297, 0.554, 0.872, 1.239, 1.646, 2.088, 2.558, 3.053, 3.571,
							 4.107,4.660,5.229,5.812,6.408,7.015,7.633,8.260,8.897,9.542,10.196,10.856,11.524},   //0.99
							{0.001,0.051,0.216,0.484,0.831,1.237,1.690,2.180,2.700,3.247,3.816,4.404,
							 5.009,5.629,6.262,6.908,7.564,8.231,8.907,9.591,10.283,10.982,11.689,12.401,13.120},         //0.975
							{0.004,0.103,0.352,0.711,1.145,1.635,2.167,2.733,3.325,3.940,4.575,5.226,
							5.892,6.571,7.261,7.962,8.672,9.390,10.117,10.851,11.591,12.338,13.091,13.848,14.611},         //0.95
							{0.016,0.211,0.584,1.064,1.610,2.204,2.833,3.490,4.168,4.865,5.578,6.304,
							 7.042,7.790,8.547,9.312,10.085,10.865,11.651,12.443,13.240,14.848,15.659,16.473}};        //0.90
	*/
	double Confidence[5] = {0.995, 0.99, 0.975, 0.95, 0.90};
	double ChiQus[5][25]={{7.789, 10.597, 12.838, 14.860, 16.750, 18.548, 20.278, 21.955, 23.589, 25.188, 26.757, 28.299,
						   29.819, 31.319, 32.801, 34.267, 35.718, 37.156, 38.582, 39.997, 41.401, 42.796, 44.181, 45.559, 46.928},    // 0.995
							{6.635, 9.210, 11.345, 13.277, 15.086, 16.812, 18.475, 20.090, 21.666, 23.209, 24.725, 26.217,
							 27.688,29.141,30.578, 32.000, 33.409, 34.805, 36.191, 37.566, 38.932,40.289, 41.638, 42.980, 44.314},   //0.99
							{5.024, 7.378, 9.348, 11.143, 12.833, 14.449, 16.013, 17.535, 19.023, 20.483, 21.920, 23.337,
							 24.736, 26.119, 27.488, 28.845, 30.191, 31.526, 32.852, 34.170, 35.479, 36.781, 38.076, 39.464, 40.646},         //0.975
							{3.841, 5.991, 7.815, 9.488, 11.071, 12.592, 14.067, 15.507, 16.919, 18.307, 19.675, 21.026,
							22.362, 23.685, 24.996, 26.296, 27.587, 28.869, 30.144, 31.410, 32.671, 33.924, 35.172, 36.415, 37.652},         //0.95
							{2.706, 4.605, 6.251, 7.779, 9.236, 10.645, 12.017, 13.362, 14.684, 15.987, 17.275, 18.549,
							 19.812, 21.064, 22.307, 23.542, 24.769, 25.989, 27.204, 28.412, 29.615, 30.813, 32.007, 33.196,34.382}};        //0.90  
	
	short i=0;
	for(i=0; i<5; i++)
	{
		if(Confidence[i] == (1-alfa))
		{
			Threshold = ChiQus[i][n-1];
			break;
		}
	}
	return  Threshold;
}
//----------------------------------------------------------------------------------
//------------------------------------------------------------------
//Remark :To use local area(HK)tropospheric delay model without met elements,
//        and it uses Niell dry Mapping function for the time being
//Author :WK.YU
//Added  :2011-04-07
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%% Theory Based:  ZH.Chen,Regional precision tropospheric delay modeling[D],
//%%                       HuNan:Central South University Master Thesis,2009
//%% Author      :  WK.YU
//%% Remark      :  To know more please see the reference.        
//%% --------------------------------------------------
//%% Function    :  To calculate the zenith tropospheric delay of station using 
//%%                local area model.
//%% Parameters  :  h->height ; B->atitude ; doy-> day of year
//%% Edited Date :  2011-04-06
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TropoBy_LATM_HK(double h,double B,short doy)
{
    double ztd_priori,ztd_hb,ztd_ht;
	ztd_priori=ztd_hb=ztd_ht=0;
    double tw=doy*0.01833;
    ztd_priori=2.579-0.07844*cos(tw)  -0.03079*sin(tw)
		            -0.01039*cos(2*tw)-0.01128*sin(2*tw)
                    -0.009456*cos(3*tw)+0.006091*sin(3*tw)
                    +0.00392*cos(4*tw)-0.006597*sin(4*tw)
                    -0.01046*cos(5*tw)+0.00533*sin(5*tw);

	B=cos(2*B);
    ztd_hb=  2.641E-03          + (-5.049E-05         *B)+     3.894E-07      *B*B +      (-1.042E-09    *pow(B,3))+       8.331E-13    *pow(B,4)
		+  (-6.857E-04*h)       +   3.817E-06*h       *B +     (-2.245E-08*h  *B*B)+        5.911E-11*h  *pow(B,3) +     (-6.787E-14*h  *pow(B,4))
		+    4.117E-06*h*h      + (-3.086E-08*h*h     *B)+       3.453E-11*h*h*B*B +        1.798E-13*h*h*pow(B,3) +     (-1.179E-16*h*h*pow(B,4))
		+  (-1.666E-08*pow(h,3))+   9.499E-11*pow(h,3)*B +  2.261E-13*pow(h,3)*B*B + (-1.871E-15*pow(h,3)*pow(B,3))+  1.375E-18*pow(h,3)*pow(B,4)
		+    2.228E-11*pow(h,4) + (-9.642E-14*pow(h,4)*B)+(-6.714E-16*pow(h,4)*B*B)+   3.648E-18*pow(h,4)*pow(B,3) +(-2.738E-21*pow(h,4)*pow(B,4));
    
    ztd_ht=(-8.0688E-03)          +   1.2628E-04    *doy + (-3.2562E-07    *doy*doy) + (-9.0762E-10    *pow(doy,3))+  2.9315E-12    *pow(doy,4)
		+    3.2512E-04*h         + (-1.1315E-05*h  *doy)+   1.0442E-07*h  *doy*doy  + (-3.4783E-10*h  *pow(doy,3))+  3.7815E-13*h  *pow(doy,4)
		+  (-4.1300E-06*h*h)      +   1.5927E-07*h*h*doy + (-1.5670E-09*h*h*doy*doy) +   5.4466E-12*h*h*pow(doy,3) +(-6.1154E-15*h*h*pow(doy,4))
		+    1.8158E-08*pow(h,3)  + (-7.0619E-10*pow(h,3)*doy)+   6.9873E-12*pow(h,3)*doy*doy  + (-2.4421E-14*pow(h,3)*pow(doy,3))+  2.7591E-17*pow(h,3)*pow(doy,4)
		+  (-2.5053E-11*pow(h,4)) +   9.8135E-13*pow(h,4)*doy + (-9.7469E-15*pow(h,4)*doy*doy) +   3.4169E-17*pow(h,4)*pow(doy,3) +(-3.8719E-20*pow(h,4)*pow(doy,4));
	
	
	return ztd_priori+ztd_hb+ztd_ht;
} 


void CtrpCrc::WriteObsHeader(CGPSRinexFileHeader rinexHeader, FILE* fobsfile)
{
    if (!fobsfile)
    {
       return;
    }
	fprintf(fobsfile,"%9.2f%11s%-20s%-20s%-20s\n",
		atof(rinexHeader.m_RinexVer),"","OBSERVATION DATA",rinexHeader.m_RinexDataType,"RINEX VERSION / TYPE");
	fprintf(fobsfile,"%-20s%-20s%-20s%-20s\n",rinexHeader.m_PRGName,"MGPSLIB",rinexHeader.m_RunDate,"PGM / RUN BY / DATE");
	fprintf(fobsfile,"%-60s%-20s\n",rinexHeader.m_MarkerName,"MARKER NAME");
	fprintf(fobsfile,"%-60s%-20s\n",rinexHeader.m_MarkerNumber,"MARKER NUMBER");
	fprintf(fobsfile,"%-20s%-40s%-20s\n",rinexHeader.m_Observer, rinexHeader.m_Agency,"OBSERVER / AGENCY");
	fprintf(fobsfile,"%-20s%-20s%-20s%-20s\n",rinexHeader.m_ReceiverNo, rinexHeader.m_ReceiverType, rinexHeader.m_ReceiverVer,"REC # / TYPE / VERS");
	fprintf(fobsfile,"%-20s%-20s%-20s%-20s\n",rinexHeader.m_AntennaNo, rinexHeader.m_AntennaType,"","ANT # / TYPE");
	fprintf(fobsfile,"%14.4f%14.4f%14.4f%18s%-20s\n",rinexHeader.m_ApproxX, rinexHeader.m_ApproxY, rinexHeader.m_ApproxZ," ","APPROX POSITION XYZ");
	fprintf(fobsfile,"%14.4f%14.4f%14.4f%18s%-20s\n",rinexHeader.m_AntDH, rinexHeader.m_AntDX, rinexHeader.m_AntDY," ","ANTENNA: DELTA H/E/N");
    fprintf(fobsfile,"%6d",rinexHeader.m_ObsTypeCount);
    int i;
	if (rinexHeader.m_ObsTypeCount<=9)
	{
	 
		for (i=0;i<rinexHeader.m_ObsTypeCount;i++)
		{
          fprintf(fobsfile,"%6s",rinexHeader.m_ObsType[i]);
		} 
		for (i=0;i<9-rinexHeader.m_ObsTypeCount;i++)
		{
           fprintf(fobsfile,"%6s","");
		}
		fprintf(fobsfile,"%-20s\n","# / TYPES OF OBSERV");	
	}
    else
	{
		for (i=0;i<rinexHeader.m_ObsTypeCount;i++)
		{
			if (i==9)
			{
			   fprintf(fobsfile,"%-20s\n%6s","# / TYPES OF OBSERV","");
			}
			fprintf(fobsfile,"%6s",rinexHeader.m_ObsType[i]);	
		} 
        fprintf(fobsfile,"\n");	
	}
      CGPSTime FirstEpotime=CGPSTime(rinexHeader.m_firstobsweek,rinexHeader.m_firstobssec);
	  CGPSTime LastEpotime=CGPSTime(rinexHeader.m_lastobsweek,rinexHeader.m_lastobssec);

	fprintf(fobsfile,"%6d%6d%6d%6d%6d%13.7f%5s%-12s%-20s\n",
	 FirstEpotime.m_year,
	 FirstEpotime.m_month,
	 FirstEpotime.m_day,
	 FirstEpotime.m_hour,
	 FirstEpotime.m_min,
     FirstEpotime.m_sec,
		"",
		"GPS",
		"TIME OF FIRST OBS");
	fprintf(fobsfile,"%6d%6d%6d%6d%6d%13.7f%5s%-12s%-20s\n",
		LastEpotime.m_year,
		LastEpotime.m_month,
		LastEpotime.m_day,
		LastEpotime.m_hour,
		LastEpotime.m_min,
		LastEpotime.m_sec,
		"",
		"GPS",
		"TIME OF LAST OBS");
	
	fprintf(fobsfile,"%10.3f%-50s%-20s\n",rinexHeader.m_Interval,"","INTERVAL");
	fprintf(fobsfile,"%-60s%-20s\n", "Observations of no nav info were removed","COMMENT");
	fprintf(fobsfile,"%-60s%-20s\n", "Tropospheric delay was corrected using LATM_HK","COMMENT");
	fprintf(fobsfile,"%-60s%-20s\n","","END OF HEADER");

}
void CtrpCrc::writeonesat(CObservation& obs, FILE* fobsfile,CGPSRinexFileHeader& header)
{
	int i,j=0;
	int typecnt=header.m_ObsTypeCount;

	for (i=0;i<typecnt;i++)
	{

		if (j==5)
		{
            fprintf(fobsfile,"\n");
			j=0;
		}
		if (0==strcmpi(header.m_ObsType[i],"C1"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vCA);
			j++;
		}
		else if (0==strcmpi(header.m_ObsType[i],"L1"))
		{
			fprintf(fobsfile,"%14.3f %1d", obs.vL1,obs.sL1);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"L2"))
		{
			fprintf(fobsfile,"%14.3f %1d", obs.vL2,obs.sL2);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"P1"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vY1);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"P2"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vY2);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"S1"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vS1);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"S2"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vS2);
			j++;
		}
		else if  (0==strcmpi(header.m_ObsType[i],"D1"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vD1);
			j++;
		}
		else if (0==strcmpi(header.m_ObsType[i],"D2"))
		{
			fprintf(fobsfile,"%14.3f  ", obs.vD2);
			j++;
		}
	}	
	fprintf(fobsfile,"\n");
}
void CtrpCrc::writeEphObsToRinexFile(CGPSEpochObsSet& obs, FILE* fobsfile,CGPSRinexFileHeader& rinexHeader)
{
	//yu:如果文件为空就退出
	if (!fobsfile)
	{
		return;
	}
	
	CGPSTime gpsT(obs.m_GpsWeek,obs.m_GpsSecond);
	int yeari = gpsT.m_year % 100;
	fprintf(fobsfile," %02d%3d%3d%3d%3d%11.7f%3d%3d",
		yeari, gpsT.m_month,gpsT.m_day,
		gpsT.m_hour, gpsT.m_min, gpsT.m_sec, 
		0, obs.GetSatCount()); //年月日时分秒 ，事件标志,卫星数
	int i = 0;
	
	if(obs.GetSatCount()<=12)  //如果小于12颗卫星
	{
		CString sLine, s;
		for(i=0; i<obs.GetSatCount(); i++)
		{
			s.Format("%s%2d","G",obs.GetAt(i)->satID);
			sLine += s;
		}
		sLine += CString(_T(' '),36 - sLine.GetLength());
		s.Format(_T("%12.9lf"), 0);
		sLine += s;
		fprintf(fobsfile,"%s\n", sLine.GetBuffer(0));
		//fprintf(m_hRinexObsFile,"\n");
	}
	else        //大于12颗卫星
	{
		CString sLine, s;
		for(i=0; i<12; i++)
		{
			s.Format("%s%2d","G",obs.GetAt(i)->satID);
			sLine += s;
		}
		sLine += CString(_T(' '),36 - sLine.GetLength());
		s.Format(_T("%12.9lf"), 0);
		sLine += s;
		fprintf(fobsfile,"%s\n", sLine.GetBuffer(0));
		////					
		fprintf(fobsfile,"%32s", "");
		for(i=12; i<obs.GetSatCount(); i++)
		{
			fprintf(fobsfile,"%s%2d","G",obs.GetAt(i)->satID);
		}
		fprintf(fobsfile,"\n");
	}


    int satcnt=obs.GetSatCount();
	for(i=0; i<satcnt; i++)
	{
		writeonesat(*obs.GetAt(i), fobsfile,rinexHeader);
	}

  
}


void CtrpCrc::OutputRinexObsfile(char* output,CGPSStationObsSet& obsset)
{
	POSITION Pos = obsset.GetHeadPosition();
	CGPSEpochObsSet *pEpochObsSet=NULL;
	
    FILE* p=fopen(output,"w");
	if (!p)
    {
		return;
    }
    WriteObsHeader(obsset.m_RinexHeader,p);
	while(NULL != Pos)
	{
		pEpochObsSet = obsset.GetNext(Pos);
        writeEphObsToRinexFile(*pEpochObsSet,p,obsset.m_RinexHeader);
	}
	
    fclose(p);
}


int CtrpCrc::Trpcrc(CString obsfile)
{

    CString navfile="";
	CGPSStationObsSet _obsset;
	CGPSStationNavSet _navset;
	CStationInfo stationInfo;

	_obsset.OpenObsFile(obsfile.GetBuffer(0),FALSE);
	//cmp nav file
	navfile = obsfile.Left(obsfile.GetLength()-1)+"N";
	_navset.OpenNavFile(navfile.GetBuffer(0),FALSE);

	if (_navset.GetItemCount()<=0||_obsset.GetItemCount()<=0)
	{
		return 0;
	}
	_obsset.RemoveNoNavSat(&_navset);
	double X = 0;
	double Y = 0;
	double Z = 0;
	if (abs(_obsset.m_RinexHeader.m_ApproxX)<0.01 )
	{
		if (0 == CalcStationOrgPos2(&_obsset, &_navset, &X, &Y, & Z))
		{
			return 0;
		}
		stationInfo.SetStationNewPos(X,Y,Z);
		_obsset.m_RinexHeader.m_ApproxX=X;
		_obsset.m_RinexHeader.m_ApproxY=Y;
		_obsset.m_RinexHeader.m_ApproxZ=Z;

	}
	
	stationInfo.SetStationNewPos(_obsset.m_RinexHeader.m_ApproxX,
								 _obsset.m_RinexHeader.m_ApproxY,
								 _obsset.m_RinexHeader.m_ApproxZ);
	
	//计算卫星高度角
	CalcStationElevation(&_obsset, &_navset, &stationInfo);

	//由标准气象元素计算
	double B=0, L=0, H=0;
	CCoordSysType CStype;
	CStype.XYZtoBLH(stationInfo.GetX(), stationInfo.GetY(), stationInfo.GetZ(), &B, &L, &H);
	CorrectTrop(&_obsset, NULL, B, H, 0, 0, 0,TROPO_TYPE_LATM_HK, false);	


	_obsset.m_RinexHeader.m_firstobsweek=_obsset.GetHead()->m_GpsWeek;
	_obsset.m_RinexHeader.m_firstobssec=_obsset.GetHead()->m_GpsSecond;
	_obsset.m_RinexHeader.m_lastobsweek=_obsset.GetTail()->m_GpsWeek;
	_obsset.m_RinexHeader.m_lastobssec=_obsset.GetTail()->m_GpsSecond;

	OutputRinexObsfile((obsfile+"_TRPCRC").GetBuffer(0),_obsset);

    _obsset.RemoveAll();
	_navset.RemoveAll();
	return 1;
}

//----------------------------------------------------------------------------
	