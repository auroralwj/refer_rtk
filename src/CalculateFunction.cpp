//CalculateFunction.cpp
#include "stdafx.h"
#include "..\INCLUDE\CalculateFunction.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//==============================================
//����˵����������ȡ��
//�磺Round(2.5) = 3; Round(-2.5) = -3;
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
//����˵����������ȡ��
//�磺Round(2.5, -1) = 2; Round(2.5,1) = 3; Round(2.5, 0) = 3;
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
	case -1:   //��С��ֵ
		y = floor(x);
		break;
	case 1: //ȡ��ֵ
		y = ceil(x);
		break;
	default:  //������������ȡֵ
		y = x>0 ? floor(x+0.5) : ceil(x-0.5);
		break;
	}
	return y;
}
//----------------------------------------------------------
//-
//ģ���Ƚ��㺯��(�ʺϵ�Ƶ˫Ƶ�ȸ���ģ���Ƚ��㷽��)
//   Amb->����ģ���ȸ���⣬���ģ���ȹ̶���
//   ambCount-> ģ���ȵĸ���
//   Cxx-> ������Э������
//   sigma0->����������
//   ratio->ratioֵ      
//   DfCount->˫��۲ⷽ�̵ĸ���
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
		rn = FixAmbFARAEx(Amb, sigma0, ambCount, Cxx, ratio);//yu:δ���
		break;
	case AMB_SECTION:
		rn = FixAmbSection(Amb, ambCount, Cxx);//yu:������
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
			else if(obsType==OBS_L1)  //��Ƶ
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
//��˫��۲ⷽ��
GPSCLASSDLL_API short	FormDblDifference(CGPSStationObsSet* pFStObsSet, CGPSStationObsSet* pUStObsSet,
		              CGPSStationNavSet * pFNavSet,
					  CStationInfo* pFStationInfo, CStationInfo* pUStationInfo,
					  CGPSDDCoefSet* pDDCoefSet, CCalcSolutionBase* pSolSetting)
{
	CGPSEpochObsSet* pFEpochObsSet=NULL;
	CGPSEpochObsSet* pUEpochObsSet=NULL;

	CDDCoef* pDDCoef=NULL;//˫��۲ⷽ��ϵ������
    CGPSEpochDDCoef* pEpochDDCoef=NULL;
	
	POSITION fPos = pFStObsSet->GetHeadPosition();
	POSITION uPos = pUStObsSet->GetHeadPosition();
	
	while (NULL!=fPos && uPos!=NULL)
	{
		pFEpochObsSet = pFStObsSet->GetNext(fPos);
		pUEpochObsSet = pUStObsSet->GetNext(uPos);
		
		//����˫��۲ⷽ��ϵ��
		pEpochDDCoef = pDDCoefSet->NewObjItem();		
		
		ComDDCoef(pFEpochObsSet, pFStationInfo,
					pUEpochObsSet, pUStationInfo, pFNavSet,
					pSolSetting->m_BaseSat, pEpochDDCoef, pSolSetting);
	}	
	return 1;
}

//����:����ϵķ���̽�������������������ݼ���
//pDDCoefSet:˫��۲����ݼ�
//pCSSet: �������ݼ�
GPSCLASSDLL_API short DetectCycleslipWithPolynomial(CGPSDDCoefSet* pDDCoefSet, CGPSCircleSlipSet* pCSSet)
{
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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
		DDOMC=new double*[satEphCount];//��ά���飬�����洢˫��۲�ֵ����
		for(j=0;j<satEphCount;j++)
		{
			DDOMC[j]=new double[4];
		}
		
		POSITION pPos = pDDCoefSet->GetHeadPosition();
		while(pPos != NULL)
		{
            pEpochDDcoef = pDDCoefSet->GetNext(pPos);
			//��ȡָ�����ǵĹ۲�����
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
		//̽������
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

//��������:��˫��۲�պϲ�������DDOMC��̽������
//����˵����DDOMC˫��۲�ֵ���У�DDCount��4�еĶ�ά���飩
//          DDCount ˫��۲�ֵ����
//          pCSSet�������ݼ�(���ڴ��̽��������)
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

	//��ʱ���׼��
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
		long delDDOMC=0;//��¼���޳������ݸ���
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
			//����������в����3�������Ҳв����0.5����Ϊ������,
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
			//����������в����3�������Ҳв����0.5����Ϊ������,
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

//�˺�����ʱ����
GPSCLASSDLL_API void RemoveSamllandGetAmbCount(CGPSDDCoefSet *pStSet, short minContinue, short maxGap, CSatAmbSet* pSatAmbSet)
{
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
/*	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //�۲��¼�ĸ���
	short ambCount = 0; //���ǵ�ģ������
	long gapCount = 0;  //��¼���ݶϿ��ĳ���

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
			//��ȡָ�����ǵĹ۲�����
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);
			if(pObs != NULL)  //����й۲����ݣ���¼�۲����ݵĸ���
			{
				obsCount++;
			}
			else   //����޹۲����ݣ����������
			{
				gapCount++;

				//�������У����ޣ�˵�������ж� 
				if(obsCount>0)
				{
					if(obsCount<minContinue)  //���������ݼ�¼С��ָ�����ݼ�¼���ȣ���ɾ���ö�����
					{
						POSITION dPos = pPos;
						for(j=obsCount; j>0; j--)
						{
							pDelDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
							pDelDDcoef->Remove(satPrn[i]);   //ɾ�������Ǹöε�����
						}
						gapCount = gapCount + obsCount; //�ö�����ɾ�����൱�ڼ���������ݳ���
						obsCount = 0;
					}
					else// if(obsCount>=minContinue)
					{
						if(ambCount == 0) //��һ��ģ����
						{
							//��������һ���µ�ģ����
							ambCount = 1; 
							
							//��¼ģ������Ϣ
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = satPrn[i];       //���Ǻ�
							pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
							pSatAmb->m_GpsWeek = pEpochDDcoef->m_GpsWeek;
							pSatAmb->m_GpsSecond = pEpochDDcoef->m_GpsSecond;  
							
							obsCount = 0;
							gapCount = 1;
						}
						else //��2��3��4��...��ģ����
						{
							if(gapCount>maxGap)   //����µ�ģ����
							{
								ambCount++;
								//��¼ģ������Ϣ
								pSatAmb = pSatAmbSet->NewObjItem();
								pSatAmb->m_satID = satPrn[i];       //���Ǻ�
								pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
								pSatAmb->m_GpsWeek = pEpochDDcoef->m_GpsWeek;
								pSatAmb->m_GpsSecond = pEpochDDcoef->m_GpsSecond;  
								
								obsCount = 0;
								gapCount = 1;
							}
							else  //����һ��ģ���ȵ�ʱ��ı�
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
//�˺������ڱ༭�����ɾ�������ݶ̺��˫��ϵ������ģ������Ϣ��ȡ
//���ڶϿ�ʱ�䳤�����ݣ������µ�ģ���ȣ�ͬʱ�����Ǳ�ż�50���൱������һ��������
BOOL EditSatAmb(CGPSDDCoefSet *pStSet, short maxGap, CSatAmbSet* pSatAmbSet)
{
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //�۲��¼�ĸ���
	short ambCount = 0; //���ǵ�ģ������
	long gapCount = 0;  //��¼���ݶϿ��ĳ���

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
			//��ȡָ�����ǵĹ۲�����
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);
			if(pObs != NULL)  //����й۲����ݣ���¼�۲����ݵĸ���
			{
				obsCount++;

				if(pPos == NULL)  //�Ѿ��������һ����Ԫ
				{
					if(ambCount==0)   //����������һֱ������ֻ��¼����ģ������Ϣ�����ø��Ӵ���
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
						
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = pEpochDDcoef->m_GpsSecond; 

						obsCount = 0;
						gapCount = 0;
					}
					else
					{
						ambCount++;
						//��¼ģ������Ϣ
						pSatAmb = pSatAmbSet->NewObjItem();
						//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = pEpochDDcoef->m_GpsSecond; 
					
						dPos = pStSet->GetTailPosition();
						//���µ����ݶε����Ǳ�Ÿı�
						for(j=obsCount; j>0; j--)
						{
							pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
			else //(pObs==NULL)  //����޹۲����ݣ����������
			{
				gapCount++;

				//�������У����ޣ�˵�������ж� 
				if(obsCount>0)
				{
					if(ambCount == 0) //��һ��ģ����
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						pSatAmb->m_satID = newSatid;       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���						
						
						if(pPos == NULL) //������һ����Ԫ
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
					else //��2��3��4��...��ģ����
					{
						if(gapCount>maxGap)   //����µ�ģ����
						{
							ambCount++;
							pSatAmb = pSatAmbSet->NewObjItem();
							//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							pSatAmb->m_satID = newSatid;       //���Ǻ�
							pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���

							if(pPos == NULL)  //�Ѿ������һ����Ԫ
							{
								tmpPos = pStSet->GetTailPosition();
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								pSatAmb->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								pSatAmb->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
								//���˶ε��������ݵ����Ǳ��Ϊ�ĳ��µ����Ǳ��
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
															
								//���˶ε��������ݵ����Ǳ��Ϊ�ĳ��µ����Ǳ��
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
									pObs = pEditDDcoef->GetDDCoef(satPrn[i]);
									if (NULL == pObs)
										return 0;
									pObs->SatID = newSatid;
								}
								
								obsCount = 0;
								gapCount = 1;
							}
						}
						else //(gapCount<maxGap)  //������µ�ģ���ȣ�ֻ����һ��ģ���ȵ����ʱ��ı�
						{
							newSatid = satPrn[i] + (ambCount-1) * 50;
							CSatAmb* plSat = pSatAmbSet->GetSatAmb(newSatid,ambCount);

							if(pPos == NULL)  //�Ѿ������һ����Ԫ
							{
								tmpPos = pStSet->GetTailPosition();
								pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
								pPrevEpochDDcoef = pStSet->GetAt(tmpPos);
								plSat->m_GpsEndWeek = pPrevEpochDDcoef->m_GpsWeek;
								plSat->m_GpsEndSecond = pPrevEpochDDcoef->m_GpsSecond;  

								dPos = tmpPos;
								//���˶ε��������ݵ����Ǳ��Ϊ�ĳ��µ����Ǳ��
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
															
								//���˶ε��������ݵ����Ǳ��Ϊ�ĳ��µ����Ǳ��
								for(j=obsCount; j>0; j--)
								{
									pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
//����:˫Ƶ������Զ�λ
//pDDCoefSet:  ˫���ϵ����
//XF0,YF0,ZF0����֪��ĳ�ʼλ��
//XU0,YU0,ZU0: δ֪��ĳ�ʼλ�ã�I&O��
//sigma0:��λȨ�����
//Cxx:Э������
//DfCount:��ɵĹ۲ⷽ�̸���
//Mion:����Ԫ������۲�ֵ���仯ֵ��һ��Ϊ4*WAVE1��
//SigmaL1:L1�۲�ֵ�����������
//SigmaL2:L2�۲�ֵ�����������
//dn1:L1���������ĵĿ��(һ��Ϊ5)
//dn5:L5(����)���������ĵĿ��(һ��Ϊ2)
void TDRelativePosDF(CGPSDDCoefSet* pDDCoefSet,
					 double XF0, double YF0, double ZF0,
					 double *XU0, double *YU0, double *ZU0, 
					 double* sigma0, Matrix* Cxx, long* DfCount,
					 double Mion,double SigmaL1,double SigmaL2,short dn1,short dn5)
{

	//step:1--------------------����׼�����޲�����
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
	short satPrn[40];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	short nRow = satNum * 2; //һ����Ԫ�����˫��۲ⷽ�̸��� L1 & L2
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
	Matrix Qv(nCol, nCol);//Э�������
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//����۲ⷽ������
	long   DeleteCount = 0 ;  //ɾ������Ԫ
	double delta=0;//�����
	double vL1=0;
	double vL2=0;
	double vL3=0;
	double X0=0, Y0=0, Z0=0;

	double Bata1=F1*F1/(F1*F1-F2*F2);
	double Bata2=-(F2*F2)/(F1*F1-F2*F2);

	//�޵������Ϲ۲�ֵ���������
	//yu:3����������Ϲ۲�ֵ������������Ϲ۲��޲�
	double SigmaL3=3*sqrt(8.0)*sqrt(Bata1*Bata1*SigmaL1*SigmaL1+Bata2*Bata2*SigmaL2*SigmaL2);
	short iteration=0;//��������
	
	//step:2----------------------�������,//yu:����Ԫ�������,��������Ԫ�ķ����̵�����һ����һ��ֵ
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		vL1=0;
		vL2=0;
		delta=0;//�����
		DbFuncCount=0;
		N.Null();
		U.Null();

       
        //step:2.1-----------------------yu:��ɷ���
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:����Ԫ,��������Ԫ�ķ����̵�����һ��
		{             
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//�������۲ⷽ��ϵ������A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)//yu:��˫���
			{

				pDDCoef=pEpochDDCoef->GetNext(epochCoef);//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
				//��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 				
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//������ּ�ϣ��жϴ˴�������������������ģ����
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //�������ˣ��˳�ѭ��
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//˵�������ǵĹ۲�����ֻ����i��ԪΪֹ 
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
				
				A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;//yu:�������Ƕ����������
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
					DbFuncCount = DbFuncCount + 2;//��¼����۲ⷽ������
					delta = delta+vL1*vL1+vL2*vL2;
				}
				
				k ++;
			}
            //yu:��������Ԫ�ķ����̵�����һ��
			//   A ����L1��L2���ڵĲ�ͬ���� , ��Ϊdx dy dz
			AT=~A;
			ATP=AT*P;
			Ni=ATP*A;
			Ui=ATP*w;
			N=N+Ni;
			U=U+Ui;	
		}
		//step:2.2-----------------------yu:����
		InvN=!N;
		x=InvN*U;
	
		*XU0=*XU0+x(0,0);
		*YU0=*YU0+x(1,0);
		*ZU0=*ZU0+x(2,0);
		
		iteration++;
		//����������10�λ�����������Ϊ�Ƿ�ɢ��,��ֹ����
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);//yu:----����


	//step:3--------------------yu:�������
	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//��λȨ�����
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

	//step:1--------------------����׼�����޲�����
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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
	Matrix Qv(nCol, nCol);//Э�������
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//����۲ⷽ������
	long   DeleteCount = 0 ;  //ɾ������Ԫ
	double delta=0;//�����
	double X0=0, Y0=0, Z0=0;

	short iteration=0;//��������
	
	//step:2----------------------�������,//yu:����Ԫ�������,��������Ԫ�ķ����̵�����һ����һ��ֵ
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		delta=0;//�����
		DbFuncCount=0;
		N.Null();
		U.Null();

       
        //step:2.1-----------------------yu:��ɷ���
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:����Ԫ,��������Ԫ�ķ����̵�����һ��
		{             
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//�������۲ⷽ��ϵ������A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)//yu:��˫���
			{

				pDDCoef=pEpochDDCoef->GetNext(epochCoef);//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
				//��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 				
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//������ּ�ϣ��жϴ˴�������������������ģ����
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //�������ˣ��˳�ѭ��
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//˵�������ǵĹ۲�����ֻ����i��ԪΪֹ 
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
				
				A(k,0)=pNextDDCoef->dEx-pDDCoef->dEx ;//yu:�������Ƕ����������
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
					DbFuncCount = DbFuncCount + 1;//��¼����۲ⷽ������
					delta += w(k,0)*w(k,0);
				}
				
				k ++;
			}
            //yu:��������Ԫ�ķ����̵�����һ��
			//   A ����L1��L2���ڵĲ�ͬ���� , ��Ϊdx dy dz
			AT=~A;
			ATP=AT*P;
			Ni=ATP*A;
			Ui=ATP*w;
			N=N+Ni;
			U=U+Ui;	
		}
		//step:2.2-----------------------yu:����
		InvN=!N;
		x=InvN*U;
	
		*XU0=*XU0+x(0,0);
		*YU0=*YU0+x(1,0);
		*ZU0=*ZU0+x(2,0);
		
		iteration++;
		//����������10�λ�����������Ϊ�Ƿ�ɢ��,��ֹ����
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);//yu:----����


	//step:3--------------------yu:�������
	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//��λȨ�����
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
//���̽������(˫Ƶ)
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

	double delta=0;//�����
	double vL1=0;
	double vL2=0;
	double vL3=0;
    double Bata1=F1*F1/(F1*F1-F2*F2);
	double Bata2=-(F2*F2)/(F1*F1-F2*F2);

	//�޵������Ϲ۲�ֵ���������
	//yu:�������3������Ϲ۲�ֵ��λȨ�����
	double SigmaL3=3*sqrt(8.0)*sqrt(Bata1*Bata1*SigmaL1*SigmaL1+Bata2*Bata2*SigmaL2*SigmaL2);

	//��¼����
	long CSL1=0;
	long CSL2=0;

	if(pCSSet!=NULL)
	{
		POSITION ddPos = pDDCoefSet->GetHeadPosition();
		while (NULL != ddPos)//yu:������Ԫ
		{             
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();
			while (NULL!=epochCoef && ddPos!=NULL)//yu:���з���
			{
				//��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
				pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);

				if(NULL == pNextDDCoef)//���ּ�ϣ�Ȼ���ҵ���i+s����˫��۲�ֵ���������
				{		
					POSITION tempPos = ddPos;
					pNextEpochDDCoef=pDDCoefSet->GetNext(tempPos); 
					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef)//yu:ֱ���ҵ�������ĳ��Ԫ˫��
						{
							break;
						}
						pNextDDCoef = NULL;
					}					
				}
				if(!pNextDDCoef)//˵�������ǵĹ۲�����ֻ����i��ԪΪֹ 
				{
					continue;//yu:����pNextDDCoefһֱΪ��ʱ�ŵ��⣬��һ����
					
				}
			   
			   ComputeDDAEx(pDDCoef,XF0,YF0,ZF0,XU0,YU0,ZU0,
					   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
					   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);
    		   ComputeDDAEx(pNextDDCoef,XF0,YF0,ZF0,XU0,YU0,ZU0,
					   pDDCoefSet->m_dx1,pDDCoefSet->m_dy1,pDDCoefSet->m_dz1,
					   pDDCoefSet->m_dx2,pDDCoefSet->m_dy2,pDDCoefSet->m_dz2);

			   
               vL1 = -(pNextDDCoef->dW - pDDCoef->dW);  
			   vL2 = -(pNextDDCoef->dWL2 - pDDCoef->dWL2);  
     		   vL3 = Bata1 * vL1 + Bata2 * vL2;//yu:�޵����Ӱ����Ϲ۲�ֵ����в�
 //
			   if(0.5*fabs(vL1+F2*F2/(F1*F1)*vL2)>Mion || fabs(vL3)>SigmaL3 ||
				   	fabs(vL1)> Mion || fabs((F2*F2)/(F1*F1)*vL2) > Mion) //��������������ʱ��Ϊ����������
			   {
                   CCircleSlip * pcs = NULL;
				   //�˴�����Ǵֲ��̽�ⲻ��������һ������һ��Ԫ̽�⣬��������У��򽫴˴�������ɾ����
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
//˫Ƶ˫����Զ�λ�����
void DDFloatRelativePosDF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double *XU0, double *YU0, double *ZU0, 
					double *Ambiguity,double* sigma0, Matrix* Cxx, long* DfCount, double* rms)
{
    //��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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
   	Matrix Cv(nCol, nCol);//Э�������
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//�����
	long   DbFuncCount=0;//˫��۲ⷽ������

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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

	//ģ���ȸ����,������˫��ģ����
	for(i=0; i<nCol-3; i++)
	{
        Ambiguity[i] = x(i+3,0);
	}
 	
	///����в�ƽ����
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //��RMS��

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
   	Matrix Cv(nCol, nCol);//Э�������
	
	double ww = 0;
	short  satCount=0;
	double aa[3];
	
	double delta=0;//�����
	long   DbFuncCount=0;//˫��۲ⷽ������
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	
    POSITION epochPos = pDDCoefSet->GetHeadPosition();
	
	//------------------- 1: wide ambiguity resolution -----------------
	
	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
		//��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	///����в�ƽ����
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //��RMS��
	
	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
		//��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
    //��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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
   	Matrix Cv(nCol, nCol);//Э�������
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//�����
	long   DbFuncCount=0;//˫��۲ⷽ������

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	//------------------- 1: wide ambiguity resolution -----------------

	while (NULL != epochPos)
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	
	///����в�ƽ����
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //��RMS��

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
	//�򿪲в������ļ���д���ļ�ͷ����
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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

	Matrix Cv(nCol, nCol);//Э�������
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//�����

	long   DbFuncCount=0;//˫��۲ⷽ������
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
					
		 //��ɹ۲ⷽ��ϵ������A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
		 //�Ѳв�۲�ֵд��в������ļ���
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//��λȨ�����
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
	//�򿪲в������ļ���д���ļ�ͷ����
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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

	Matrix Cv(nCol, nCol);//Э�������
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//�����

	long   DbFuncCount=0;//˫��۲ⷽ������
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
					
		 //��ɹ۲ⷽ��ϵ������A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
		 //�Ѳв�۲�ֵд��в������ļ���
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//��λȨ�����
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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

//����:��Ƶ������Զ�λ
//pDDCoefSet:  ˫���ϵ����
//XF0,YF0,ZF0����֪��ĳ�ʼλ��
//XU0,YU0,ZU0: δ֪��ĳ�ʼλ�ã�I&O��
//sigma0:��λȨ�����
//Cxx:Э������
//DfCount:��ɵĹ۲ⷽ�̸���
BOOL TDRelativePosSF(CGPSDDCoefSet* pDDCoefSet,
					 double XF0, double YF0, double ZF0,
					 double *XU0, double *YU0, double *ZU0, 
					 double* sigma0, Matrix* Cxx, long* DfCount)
{
    long i=0, j=0, k=0, l=0;
	short m = 0;
	long epochCount = 0;
	
	epochCount = pDDCoefSet->GetItemCount();
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
	short satPrn[40];
	short satNum = 0;
	
	pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);

	short nRow = satNum; //һ����Ԫ�����˫��۲ⷽ�̸���
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
	Matrix Qv(nCol, nCol);//Э�������
	
	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    
	short satcount=0;
	long   DbFuncCount=0;//����۲ⷽ������
	long   DeleteCount = 0 ;  //ɾ������Ԫ
	double delta=0;//�����
	double vL1=0;
	double X0=0, Y0=0, Z0=0;

	short iteration=0;//��������
	
	//�������
    do
	{
		X0=*XU0;
  	    Y0=*YU0;
	    Z0=*ZU0;
		vL1=0;
		delta=0;//�����
		DbFuncCount=0;
		N.Null();
		U.Null();

		POSITION ddPos = pDDCoefSet->GetHeadPosition();
	
		while (NULL != ddPos)
		{             
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
            pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			//�������۲ⷽ��ϵ������A
			A.Null();
			w.Null();
			P.Unit();

			k = 0;
			POSITION epochCoef = pEpochDDCoef->GetHeadPosition();

			while (NULL!=epochCoef && ddPos!=NULL)
			{
				//��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
				pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				
				//������ּ�ϣ��жϴ˴�������������������ģ����
				if(NULL == pNextDDCoef)
				{					
					POSITION tempPos = ddPos;
					pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 

					while (NULL != tempPos)
					{
						pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						if(pNextDDCoef) //�������ˣ��˳�ѭ��
						{
							break;
						}
						pNextDDCoef = NULL;
					}
					
				}
				if(!pNextDDCoef)//˵�������ǵĹ۲�����ֻ����i��ԪΪֹ 
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
					DbFuncCount = DbFuncCount + 1;//��¼����۲ⷽ������
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
		//����������10�λ�����������Ϊ�Ƿ�ɢ��,��ֹ����
		if(iteration>10) break;
		
	}while(abs(x(0,0)) > 0.02 && abs(x(1,0)) > 0.02);

	delta = delta / (DbFuncCount-3);
	if(DfCount!=NULL)
	{
		*DfCount = DbFuncCount;
	}

	*sigma0 = sqrt(delta);//��λȨ�����
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
//���̽������(��Ƶ)
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
	Matrix Qv(3,3);//Э�������
	
     
////////////////////////////////////////////////////

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;
	CDDCoef* pNextDDCoef=NULL;
    short satcount=0;
///////////////////////////////////////////

	long   DbFuncCount=0;//����۲ⷽ������
	double delta=0;//�����
	double v=0;
	double X0,Y0,Z0;
	double aa[3],aa0[3];
	double ww,ww0;
	short iteration=0;//��������
	//�������

    do
	{
		X0 = XU0;
  	    Y0 = YU0;
	    Z0 = ZU0;
		v=0;
		delta=0;//�����
		DbFuncCount=0;
		N.Null();
		U.Null();

		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)
		{
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
            pEpochDDCoef=pDDCoefSet->GetNext(pos); 
			if (NULL == pos)
			{//���ݴ������,�޺�������
				break;
			}
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 

			//�������۲ⷽ��ϵ������A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)
			{
				//��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				if (NULL == pNextDDCoef)
					continue;
//				ASSERT(NULL != pNextDDCoef); //���ݱ������Ѿ������
				
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
				
				if(P(k,k)>0)//ȨΪ0��ʾ�۲�ֵ�ѱ��޳� 
				{
					DbFuncCount=DbFuncCount+1;//��¼����۲ⷽ������
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
		//�������촦10�λ�����������Ϊ�Ƿ�ɢ��,��ֹ����
		if(iteration>10)
			return FALSE;

	}while(abs(x(0,0))>0.02 && abs(x(1,0))>0.02);

	delta=0;
	//��¼����	
	if(pCSSet!=NULL)
	{
		X0=XU0;
  	    Y0=YU0;
	    Z0=ZU0;
		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)
		{
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
			pEpochDDCoef=pDDCoefSet->GetNext(pos);
			if (NULL == pos)
				break;
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 
			//�������۲ⷽ��ϵ������A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)
			{
				//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				//
				if (NULL == pNextDDCoef)
					continue;
				
//				ASSERT(NULL != pNextDDCoef); //���ݱ������Ѿ������
				///////////////
               ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
			   ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
			   v=-(ww-ww0);
			   if((abs(v)/(WAVE1))>0.5) 
			   {
                   CCircleSlip* pcs= pCSSet->NewObjItem();;
				 // �в����0.5�ܵı���Ϊ������
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
//	*sigma0=sqrt(delta);//��λȨ�����
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
//yu:�������ʵ���⣬��������۲�в����0.5���ж�Ϊ����
//���̽������(��Ƶ)
BOOL DetectCycleslipWithTD(CGPSDDCoefSet* pDDCoefSet, 
					double XF0, double YF0, double ZF0,
					double XU0, double YU0, double ZU0,
					CGPSCircleSlipSet* pCSSet)
{
	 long i=0, j=0, k=0, l=0;
	 short m = 0;
	 long epochCount = 0;
	 
	 epochCount = pDDCoefSet->GetItemCount();
	 
	 //��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
	 short satPrn[40];
	 short satNum = 0;
	 
	 pDDCoefSet->GetAllObsSat(&satNum, satPrn, 40);
	 
	 short nRow = satNum; //һ����Ԫ�����˫��۲ⷽ�̸���
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
	 Matrix Qv(nCol, nCol);//Э�������
	 
	 CGPSEpochDDCoef* pEpochDDCoef=NULL;
	 CGPSEpochDDCoef* pNextEpochDDCoef=NULL;
	 CDDCoef* pDDCoef=NULL;
	 CDDCoef* pNextDDCoef=NULL;
	 
	 short satcount=0;
	 long   DbFuncCount=0;//����۲ⷽ������
	 long   DeleteCount = 0 ;  //ɾ������Ԫ
	 double delta=0;//�����
	 double vL1=0;
	 double X0=0, Y0=0, Z0=0;
	 
	 short iteration=0;//��������
	 double aa[3],aa0[3];
 	 double ww,ww0;

	 //�������
	 do
	 {
		 X0=XU0;
		 Y0=YU0;
		 Z0=ZU0;
		 vL1=0;
		 delta=0;//�����
		 DbFuncCount=0;
		 N.Null();
		 U.Null();
		 
		 POSITION ddPos = pDDCoefSet->GetHeadPosition();
		 
		 while (NULL != ddPos)
		 {             
			 //��õ�i��Ԫ˫��۲ⷽ��ϵ��
			 pEpochDDCoef=pDDCoefSet->GetNext(ddPos);
			 //�������۲ⷽ��ϵ������A
			 A.Null();
			 w.Null();
			 P.Unit();
			 
			 k = 0;
			 POSITION epochCoef = pEpochDDCoef->GetHeadPosition();
			 
			 while (NULL!=epochCoef && ddPos!=NULL)
			 {
				 //��õ�i+1��Ԫ˫��۲ⷽ��ϵ�� 
				 pNextEpochDDCoef=pDDCoefSet->GetAt(ddPos); 
				 //��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
				 pDDCoef=pEpochDDCoef->GetNext(epochCoef);
				 pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				 
				 //������ּ�ϣ��жϴ˴�������������������ģ����
				 if(NULL == pNextDDCoef)
				 {					
					 POSITION tempPos = ddPos;
					 pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
					 
					 while (NULL != tempPos)
					 {
						 pNextEpochDDCoef = pDDCoefSet->GetNext(tempPos); 
						 pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
						 if(pNextDDCoef) //�������ˣ��˳�ѭ��
						 {
							 break;
						 }
						 pNextDDCoef = NULL;
					 }
					 
				 }
				 if(!pNextDDCoef)//˵�������ǵĹ۲�����ֻ����i��ԪΪֹ 
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
					 DbFuncCount = DbFuncCount + 1;//��¼����۲ⷽ������
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
		 //����������10�λ�����������Ϊ�Ƿ�ɢ��,��ֹ����
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
	//��¼����	
	if(pCSSet!=NULL)
	{
		X0=XU0;
  	    Y0=YU0;
	    Z0=ZU0;
		POSITION pos = pDDCoefSet->GetHeadPosition();
		while (NULL != pos)//yu:������������Ԫ˫���
		{
			//��õ�i��Ԫ˫��۲ⷽ��ϵ��
			pEpochDDCoef=pDDCoefSet->GetNext(pos);
			if (NULL == pos)
				break;
			pNextEpochDDCoef=pDDCoefSet->GetAt(pos); 
			//�������۲ⷽ��ϵ������A
			 A.Null();
			 w.Null();
			 P.Unit();
			 POSITION ephPos = pEpochDDCoef->GetHeadPosition();
			 while (NULL != ephPos)//yu:ǰ����Ԫ��������˫���
			{
				//��������Ԫ���������ʱ���۲����Ǳ��뱣��һһ��Ӧ
			    pDDCoef=pEpochDDCoef->GetNext(ephPos);
				pNextDDCoef=pNextEpochDDCoef->GetDDCoef(pDDCoef->SatID);
				//
				if (NULL == pNextDDCoef)//yu:����жϾ���һ������
					continue;
				///////////////
               ComputeDDA(pDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa0,&ww0);
			   ComputeDDA(pNextDDCoef,XF0,YF0,ZF0,X0,Y0,Z0,aa,&ww);
			   double v=-(ww-ww0);//yu:������в�
			   if((abs(v)/(WAVE1))>0.5) ///////////////////////////////////////////////
			   {
                   CCircleSlip* pcs= pCSSet->NewObjItem();
				 // �в����0.5�ܵı���Ϊ������
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

//��Ƶ˫����Զ�λ�����
void DDFloatRelativePosSF(CGPSDDCoefSet *pDDCoefSet, 
					double XF0, double YF0, double ZF0, 
					double *XU0, double *YU0, double *ZU0, 
					double *Ambiguity,double* sigma0, Matrix* Cxx, long* DfCount, double* rms)
{
    //��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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
   	Matrix Cv(nCol, nCol);//Э�������
 
	double ww = 0;
	short  satCount=0;
	double aa[3];
	 	
	double delta=0;//�����
	long   DbFuncCount=0;//˫��۲ⷽ������

	CGPSEpochDDCoef* pEpochDDCoef=NULL;
	CDDCoef* pDDCoef=NULL;

    POSITION epochPos = pDDCoefSet->GetHeadPosition();

	while (NULL != epochPos)//yu:���з���
	{
		pEpochDDCoef = pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
		while (NULL != ddCoefPos)//yu:���з���
		{
			pDDCoef=pEpochDDCoef->GetNext(ddCoefPos);
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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

	//ģ���ȸ����
	for(int i=0; i<nCol-3; i++)
	{
        Ambiguity[i] = x(i+3,0);
	}
 	
	///����в�ƽ����
	Matrix v(nRow,1);
	Matrix vT(1,nRow);
	Matrix vTP(1,nRow);
	Matrix vTPv(1,1);
	//-----------
	Matrix vTv(1,1);
	double vv = 0;  //��RMS��

	epochPos = pDDCoefSet->GetHeadPosition();
	while (NULL != epochPos)
	{
        pEpochDDCoef=pDDCoefSet->GetNext(epochPos);
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
    
	*sigma0 = sqrt(delta);
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
	//�򿪲в������ļ���д���ļ�ͷ����
	FILE* fResidul=NULL;
	if(chFile)
	{
	   fResidul=fopen(chFile,"w");
	}
	
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
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

	Matrix Cv(nCol, nCol);//Э�������
    

	short satCount=0;
	double ww = 0;
	double aa[3];
	Matrix vTv(1,1);
	double vv = 0;
	double delta=0;//�����

	long   DbFuncCount=0;//˫��۲ⷽ������
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
					
		 //��ɹ۲ⷽ��ϵ������A
		short pos=-1;
		POSITION epochDDCoefPos = pEpochDDCoef->GetHeadPosition();
		while (NULL != epochDDCoefPos)
		{
			pDDCoef = pEpochDDCoef->GetNext(epochDDCoefPos);
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
         //��ɹ۲ⷽ��ϵ������A
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
			//satprn�����е����Ǳ�������򣬸��㷨����Ч
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
		 //�Ѳв�۲�ֵд��в������ļ���
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
	//����Э������
	delta = delta/(DbFuncCount-nCol);
	if(DfCount!=NULL) *DfCount=DbFuncCount;
	*sigma0=sqrt(delta);//��λȨ�����
	*rms = sqrt(vv/(DbFuncCount-1));   //Rmsֵ
	
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
//��������
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
//modelType:������ģ��
void CorrectTrop(CGPSEpochObsSet * pEphObs, CMet * pMet, 
				 double B, double H,
				 double dTemperature, double dPressure, double dRh, double H0,
				 short modelType, TropMode tropParam)
{
	//P0:�ο��߶ȵĴ���ѹ(milliBar)
	//T0:�ο��߶ȵĸ��¶�(K)
	//RH0:�ο��߶ȵ����ʪ��(%) 
	//H0:�ο��߶�(m)
	//Ph:�߶�h������ѹ(milliBar)
	//Th:�߶�h�����¶�(K)
	//Pv:�߶�h����ʪ��ѹ(milliBar)
	//RHh:�߶�h�������ʪ��%
    //h:�߶�(m)
	double T0=291.15;        
	double P0=1013.25;      
	double RH0=50;          
	double Pa=0, Td=0, Pv=0, trop=0;
	
	if (NULL != pMet)   //�ü�¼������Ԫ��
	{
		//double tpa=0, ttd=0, tpv=0, HR=0;
		//calmeteor(H,P0,T0,RH0,&tpa,&ttd,&tpv, & HR);
		Td = pMet->TD + 273.15;
		Pa = pMet->PR; //
		Pv = GetPvByTdAndHR(Td, pMet->RH); 
	}
	else if(dTemperature!=0 && dPressure!=0 && dRh!=0)  //�ú�ƽ���վƽ������Ԫ��
	{
		T0 = dTemperature + 273.15;
		P0 = dPressure;
		RH0 = dRh;
		//calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
		double RHh = 0;
		CalMeteorParameter(P0, T0, RH0, H0, H, &Pa, &Td, &Pv, & RHh);
	}
	else  //�ñ�׼����Ԫ�� 
	{
		//calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
		double RHh = 0;
		CalMeteorParameter(P0, T0, RH0, H0, H, &Pa, &Td, &Pv, & RHh);
	}
	
	//���������
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
		{//trop�쳣ʱ,������
			double dtrop_ = trop;
			pObsData->vL1 -= dtrop_/WAVE1;
			if(fabs(pObsData->vL2) >0.001)//��ʾ��L2�۲�ֵʱ�����������Է�������ж��Ƿ���L2�۲�ֵ
			{
				pObsData->vL2 -= dtrop_/WAVE2; 
			}
		}
	}

}
//-
//----------------zxw1 2008-3-8------------------------------------------ 
//���ݲο��߶ȵ���ѹ�����£����ʪ�ȵȼ���ָ���߶ȵ���ѹ��ˮ��ѹ���¶�
//P0:�ο��߶ȵĴ���ѹ(milliBar)
//T0:�ο��߶ȵĸ��¶�(K)
//RH0:�ο��߶ȵ����ʪ��(%) 
//H0:�ο��߶�(m)
//Ph:�߶�h������ѹ(milliBar)
//Th:�߶�h�����¶�(K)
//Pv:�߶�h����ʪ��ѹ(milliBar)
//RHh:�߶�h�������ʪ��%
//h:�߶�(m)
void CalMeteorParameter(double P0, double T0, double RH0, double H0, double Hh, double *Ph, double *Th, double *Pv , double *RHh)
{
	double aa=0, bb=0, cc=0, dH=0;
	dH = Hh - H0;
	
	//��h������ѹ
	aa = 1.0 - 2.26e-5 * dH;
	bb = pow(aa, 5.225);
	*Ph = P0 * bb;
	
	//��h�����¶�	
	*Th = T0 - 0.0065 * dH;

	//��h�������ʪ��
	double g=0, d=0;
	g = -6.396e-4 * dH;
	d = exp(g); 
	
	double relhum = RH0 * d; 
	if (NULL != RHh)
	{
		*RHh = relhum;
	}

	//��h����ˮ��ѹ
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
//��������
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
	
	if (NULL != pMet)   //�ü�¼������Ԫ��
	{
		//double tpa=0, ttd=0, tpv=0, HR=0;
		//calmeteor(H,P0,T0,RH0,&tpa,&ttd,&tpv, & HR);
		Td = pMet->TD + 273.16;
		Pa = pMet->PR; //
		Pv = GetPvByTdAndHR(Td, pMet->RH); 
	}
	else if(dTemperature!=0 && dPressure!=0 && dRh!=0)  //�ú�ƽ���վƽ������Ԫ��
	{
		T0 = dTemperature + 273.16;
		P0 = dPressure;
		RH0 = dRh;
		calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
	}
	else  //�ñ�׼����Ԫ�� 
	{
		calmeteor(H,P0,T0,RH0,&Pa,&Td,&Pv);
	}
	
	//���������
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
		{//trop�쳣ʱ,������
			double dtrop_ = (bTropParam?(trop - mf_wet*WZD):trop);
			pObsData->vL1 -= dtrop_/WAVE1;
			if(fabs(pObsData->vL2) >0.001)//��ʾ��L2�۲�ֵʱ�����������Է�������ж��Ƿ���L2�۲�ֵ
			{
				pObsData->vL2 -= dtrop_/WAVE2; 
			}
		}
	}
}
//-
//������ģ�͸���
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
//Sasstamoinenģ����������������    ***************
// ����˵����Pa->����ѹ����Td->�����¶ȣ�Pv->ˮ��ѹ  *********
//           elev->���Ǹ߶Ƚ�                        
//           Ht->Height of the station
//           B->The latitude of the station 
//MITӳ��ϵ��
short TropoBySaastamoinenWithMIT(double Pa, double Td, double Pv,double elev, 
										  double Ht,double Bt, double *dtrop, double *WZD, double *mf_wet)
{
	double B = Bt * PI / 180.0;
	double fBh = 1.0 - 0.00266 * cos(2*B) - 0.00028 * Ht / 1000.0;

	//�춥����ĸ��ӳٷ���
	double rkDry = 0.002277 * Pa / fBh;
	//�춥�����ʪ�ӳٷ���
	double rkWet = 0.002277 * Pv / fBh * (1255/Td + 0.05);   //  �ܶ����¶����õ������ʽ

	double MF1 = 0, MF2 = 0;
	
	TropoDryMapFuncMIT(Td, Bt, Ht, elev, &MF1);
	TropoWetMapFuncMIT(Td, Bt, Ht, elev, &MF2);

	*dtrop = rkDry * MF1 + rkWet * MF2;
	*WZD = rkWet;
	*mf_wet=MF2;
	
	return 1;
}
//Chaoӳ��ϵ��
short TropoBySaastamoienWithChao(double Pa, double Td, double Pv,double elev, 
										  double Ht,double Bt, double *dtrop, double *WZD, double *mf_wet)
{
	double B = Bt * PI / 180.0;
	double fBh = 1.0 - 0.00266 * cos(2*B) - 0.00028 * Ht / 1000.0;

	//�춥����ĸ��ӳٷ���
	double rkDry = 0.002277 * Pa / fBh;
	//�춥�����ʪ�ӳٷ���
	double rkWet = 0.002277 * Pv / fBh * (1255/Td + 0.05);   //  �ܶ����¶����õ������ʽ

	double MF1 = 0, MF2 = 0;
	
	MF1 = TropoDryMapFuncChao(elev);
    MF2 = TropoWetMapFuncChao(elev);

	*dtrop = rkDry * MF1 + rkWet * MF2;
	*WZD = rkWet;
	*mf_wet=MF2;

	return 1;
}

//=============================================================
//����ɷ���ӳ�亯���������¶Ⱥ͸߶�
//Td���¶�, elew: �߶Ƚ�, B: γ��, Ht: �߶�, drymap: �ɷ���ӳ��ϵ��
short TropoDryMapFuncMIT(double Td, double elev, double Bt, double Ht, double *drymap)
{
	double B = Bt * PI / 180.0;
	double cosB = cos(Bt);
	double Hkm = Ht / 1000.0;  //����ǧ��

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
//����ɷ���ӳ�亯���������¶Ⱥ͸߶�
//Td���¶�, elew: �߶Ƚ�, B: γ��, Ht: �߶�, wetmap: ʪ����ӳ��ϵ��
short TropoWetMapFuncMIT(double Td, double elev, double Bt, double Ht, double *wetmap)
{
	double B = Bt * PI / 180;
	double cosB = cos(B);
	double Hkm = Ht / 1000;  //����ǧ��

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

//yu:�˺�������Ƿ�ѣ���Ҫ��д��ֱ�Ӹ���ǰ����Ԫ��ʱ�������ж��Ƿ���Ҫ������ģ����cntAmb++
//  prn=(cntAmb-1)*50+prn ���ɣ��������
//-
//yu:���μ���������ĩ�㴦 ���ݻ��ο�ʼ���ļ�������ǰ (ephGapEnd-ephGapStart) > maxGap ���жϻ����Ƿ������µ�ģ���� 
//yu:��Ԫ������������ ����ǰ����Ԫ (ephEnd-ephend0) > maxGap ��Ԫ����Ƿ�̫�����ж��Ƿ����»���
//�˺������ڱ༭�����ɾ�������ݺ��˫��ϵ������ģ������Ϣ��ȡ
//���ڶϿ�ʱ�䳤�����ݣ������µ�ģ���ȣ�ͬʱ�����Ǳ�ż�50���൱������һ��������
BOOL EditSatAmb(CGPSDDCoefSet *pStSet, double maxGap, CSatAmbSet* pSatAmbSet)
{
	//��ȡ˫��۲ⷽ��ϵ��������������������飨��ȥ�ο����ǣ�
	short satPrn[40];
	short satNum = 0;

	pStSet->GetAllObsSat(&satNum, satPrn, 40);

	short i = 0, j = 0;
	long obsCount = 0; //�۲��¼�ĸ���
	short ambCount = 0; //���ǵ�ģ������
	long gapCount = 0;  //��¼���ݶϿ��ĳ���
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
	for(i=0; i<satNum; i++)//yu:������ǽ���ģ������Ϣ
	{
		
		ambCount = 0;
		obsCount = 0;
		gapCount = 0;
		ephStart = 0;
		ephEnd = 0;
		ephEnd0 = 0;//yu:��¼ǰ�������ݵ�ʱ��,��˫�����ǰһԪ��ʱ��
		ephGapStart = 0;
		ephGapEnd = 0;
		ephGapEnd0 = 0;

		pPos = pStSet->GetHeadPosition();

		//yu:��ν����ĩβ��ͬ�����ݵ�ĩβ7�����������������ݶ�������������
		//   ����϶��ǰ��������Ԫ�����ݼ��������bc����ͬ��������Ԫ��������ƻ��������������ݶ���������ʱ��ǰ����Ԫ�������
		//
		//-------------------------------------------------------------------
		//                  |  >maxGap |                      |>maxGap|
        //          |  <min |            |>maxGap |         
		//prn01:    * * * * *          * *        * * * * * * *       * * * *
		//prn##:* * * * * * * * *                     * * * * * * * * * * * * 
		//      + + + + + + + + +      + +        + + + + + + + + + + + + + + ��˫��̵���Ԫ
		//------------------------------------------------------------------- �ȼ��ʱ��
		//      1 2 3 4 5 6 7 8 9      a b        c d e f g h i j k l m n o p 
		//-------------------------------|ephEnd0-|ephEnd--------------------
		//-----------------------------|ephStart-----------------------------
		//-------ephGapStart| >maxGap  |ephGapend----------------------------
		//���μ�����ģ���ȣ�i���жϴ˻�����ǰһ���μ�϶̫��Ӧ������ģ����
        //��Ԫ�����»����漴��ģ���ȣ�c�����������������»��ο�ʼ������Ԫ�����󣩣�m�����鲻�������»��ο�ʼ��

		while(pPos != NULL)//yu:��������˫��������˿�����
		{
            pEpochDDcoef = pStSet->GetNext(pPos);
			//��ȡָ�����ǵĹ۲�����
			pObs = pEpochDDcoef->GetDDCoef(satPrn[i]);

			//yu:  3-7,a-i,m-p��������
			if(pObs != NULL)  //�������Ԫ�д����ǹ۲����ݣ���¼�۲����ݵĸ���
			{
				obsCount++;

				//-------��¼�۲����ݵ���ֹʱ��
				if(obsCount == 1)//yu: 3,a,m ���������
				{
					ephStart = pEpochDDcoef->m_GpsSecond;  //��¼��ʼʱ��
					ephEnd = pEpochDDcoef->m_GpsSecond;  
					ephEnd0 = ephEnd;
				}
				else
				{
					ephEnd0 = ephEnd;//yu:ǰһ��Ԫʱ��
					ephEnd = pEpochDDcoef->m_GpsSecond;  //��¼����ʱ��
				}

				if(pPos==NULL)  //�Ѿ��������һ����Ԫ //yu: p
				{
					if(obsCount<=0){continue;  }
					
					//��ǰ����������һֱ������ֻ��¼����ģ������Ϣ
					if(ambCount==0)   //yu:��ʼ����һֱ������ͼ����
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = satPrn[i];       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
						
						pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsStartSecond = ephStart;
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd;					
					}
					else  //�Ѿ�����ģ������Ϣ //yu: p
					{
						//yu:����ĩ�жϸö��Ƿ���ǰ���μ��̫���������ģ����
						if((ephGapEnd-ephGapStart) > maxGap)  //�Ͽ�ʱ�����������µ�ģ����,���ö����ݵ����Ǳ�Ÿı�
						{
							ambCount++;
							//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = newSatid;       //���Ǻ�
							pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
							
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

						//���µ����ݶε����Ǳ�Ÿı�
						if(ambCount != 1)
						{
							dPos = pStSet->GetTailPosition();
							
							for(j=obsCount; j>0; j--)
							{
								pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
				else if((ephEnd-ephEnd0) > maxGap)  //ǰ����Ԫ���̫��  //yu:c��,�����»��Σ����жϸս��޵Ļ��ε�ģ����
				{ 
					//����������һֱ������ֻ��¼����ģ������Ϣ
					if(ambCount==0)   
					{
						ambCount = 1;
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = satPrn[i];       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
						
						pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsStartSecond = ephStart;
						pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
						pSatAmb->m_GpsEndSecond = ephEnd0;	
						
						ephGapStart = ephEnd0;
						ephGapEnd = ephEnd;
						ephGapEnd0 = ephEnd;
					}
					else  //�Ѿ�����ģ������Ϣ //yu:c��
					{
                        //�жϻ����Ƿ�������ģ����
						if((ephGapEnd-ephGapStart) > maxGap)  //�Ͽ�ʱ�����������µ�ģ����ambCount++,���ö����ݵ����Ǳ�Ÿı�
						{
							ambCount++;
							//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
							newSatid = satPrn[i] + (ambCount-1) * 50;
							
							pSatAmb = pSatAmbSet->NewObjItem();
							pSatAmb->m_satID = newSatid;       //���Ǻ�
							pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
							
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
						
						//���µ����ݶε����Ǳ�Ÿı�
						if(ambCount !=1 )
						{
							tmpPos = pPos;
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							pPrevEpochDDcoef = pStSet->GetPrev(tmpPos);
							
							dPos = tmpPos;
					
							for(j=obsCount; j>0; j--)
							{
								pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
			else   //(pObs==NULL)  //����޸���Ԫ�޴����ǹ۲����ݣ���������� 
			{
				//------------------
				//�ж��й۲����ݼ�¼
				if(obsCount == 0) //yu:1,2,8,9,j,k,l
				{ 					
					continue;
				}

				//����������һֱ������ֻ��¼����ģ������Ϣ
				if(ambCount==0)   //yu:8
				{
					ambCount = 1;
					pSatAmb = pSatAmbSet->NewObjItem();
					pSatAmb->m_satID = satPrn[i];       //���Ǻ�
					pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
					
					pSatAmb->m_GpsStartWeek = pEpochDDcoef->m_GpsWeek;
					pSatAmb->m_GpsStartSecond = ephStart;	//yu:������
					pSatAmb->m_GpsEndWeek = pEpochDDcoef->m_GpsWeek;
					pSatAmb->m_GpsEndSecond = ephEnd;	//yu:ǰһ��ԪΪ���ν�����
					
					ephGapStart = pEpochDDcoef->m_GpsSecond;//yu:�жϿ�ʼ��
					ephGapEnd = pEpochDDcoef->m_GpsSecond;
					ephGapEnd0 = pEpochDDcoef->m_GpsSecond;
				}
				else  //�Ѿ�����ģ������Ϣ    //yu:j
				{
					if((ephGapEnd-ephGapStart) > maxGap)  //�Ͽ�ʱ�����������µ�ģ����,���ö����ݵ����Ǳ�Ÿı�
					{
						ambCount++;
						//���Ǻű�Ϊԭ����  ambCount * 50 + satID;
						newSatid = satPrn[i] + (ambCount-1) * 50;
						
						pSatAmb = pSatAmbSet->NewObjItem();
						pSatAmb->m_satID = newSatid;       //���Ǻ�
						pSatAmb->m_ambNum = ambCount;       //ģ���ȸ���
						
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
						ephGapStart = pEpochDDcoef->m_GpsSecond;  //��¼�Ͽ��Ŀ�ʼʱ��
						ephGapEnd = pEpochDDcoef->m_GpsSecond;  
						ephGapEnd0 = ephGapEnd;
					}
					else
					{
						ephGapEnd0 = ephGapEnd;                   //��¼�Ͽ��Ľ���ʱ��
						ephGapEnd = pEpochDDcoef->m_GpsSecond;
					}

					//���µ����ݶε����Ǳ�Ÿı�
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
							pEditDDcoef = pStSet->GetPrev(dPos); //��ǰ��¼��ǰһ����¼
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
// tr: �źŽ���ʱ��
// pEph: ������������
// Xr, Yr, Zr�����ջ��Ľ���λ��
// Xk, Yk, Zk, Xdot, Ydot, Zdot: ����λ�á��ٶȣ�ƽ�����
//�ɽ��ջ��Ľ���λ�ã��������ǵķ���ʱ��
double GetSatSentTime(double tr, ephemeris pEph, 
					  double Xr, double Yr, double Zr,
					  double* Xk, double* Yk, double* Zk, double* Ek,
					  double* Xdot, double* Ydot, double* Zdot, double* range)
{
	double ts0 = 0;
	double ts = 0;  //�źŷ���ʱ��
	double deltT = 0;
	double satDt = 0; //�����Ӳ�

    ts = tr - 0.075;

	do{
		//�������Ƿ���ʱ������ʱ�����������Ӳ�
		deltT = ts - pEph.toe;
		if(deltT>302400.00) //�������ڽ����Ŀ�ʼ�ͽ���ʱ��
		{
			deltT -= 604800.0;
		}
		if(deltT<-302400.0) 
		{
			deltT += 604800.0;
		}
		
		//�����Ӳ�
		satDt = pEph.a0 + (pEph.a1 * deltT) + (pEph.a2 * deltT * deltT);
        
		//���Ƿ���ʱ��
		ts = ts + satDt;
		
		//�������Ƿ���ʱ��
		ts0 = ts;

		//��������λ��
		coordsat(pEph, ts0, Xk, Yk, Zk, Ek, Xdot, Ydot, Zdot);
		*range = sqrt((*Xk-Xr)*(*Xk-Xr)+(*Yk-Yr)*(*Yk-Yr)+(*Zk-Zr)*(*Zk-Zr));
		
		//���¼�����������źŷ���ʱ��
		ts = tr - *range / C;
	}while((ts-ts0)>1.0E-7);  //�����μ���ķ���ʱ��С��10-7��ʱ����������
	
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
	
	//��ģ���Ƚ����ж�
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
    
	//��ģ����������AmbVector��������	
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
//��������:�����¹������ģ�����������Ա����ģ���ȵ�����
//����Ϊ�����ͬһ�����ǵ�L1,L2ģ���������е�����ģ����֮��N_JK��������������
//��ӱ�ѡģ����������ɾ����ģ��������
//N_jk-t(a/2)*mN_jk<=N_JK<=N_jk+t(a/2)*mN_jk
//N_jk �� XN1 - (f1/f2)*XN2:
//t(a/2):���Ŷ�Ϊ��1-a����t�ֲ�ֵ
//mN_jk:����ģ���ȵ�ʵ����Ĳ�ֵ����󷽲�
//����˵��:AmbVector:��ѡ��ģ����������
//         VectorCount:��ѡ��ģ����������ĸ���
//         Amb:�����ģ����
//         AmbCount:�����ģ���ȸ���
//         Dxx:����ⷽ����
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
	//��ģ���Ƚ����ж�
	/*do
	{
		count = *VectorCount;
		for(i=0; i<*VectorCount; i++)
		{
			short mark = 0;
			for(j=0; j<AmbCount/2; j++)
			{
				deltan = Amb[j] - (77.0/60.0) * Amb[j+k];   //xn1 - f1/f2 *xn2  //����ģ����
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
			deltan = Amb[j] - (77.0/60.0) * Amb[j+k];   //xn1 - f1/f2 *xn2  //����ģ����
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
    
	//��ģ����������AmbVector��������
	
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
// ��������:����ģ���Ƚ���(Kalman˫Ƶģ���Ƚ���) 
// ����˵��:in&out:Amb1, Amb2,->����ģ���ȸ���⣬���ģ���ȹ̶���
//          in:sigma0->����������            **********
//          in:AmbCount-> ģ���ȵĸ���             ********
//          in:Dx1, Dx2-> ������Э������
//          out:ratio->ratioֵ                    ********
//================================================================
GPSCLASSDLL_API short FixAmbFARADF(double *Amb, double sigma0, short AmbCount, Matrix Dx, long DfCount, double* ratio)
{
	//�ҳ����п��ܵ�ģ��������
	long i=0;
	double** AmbVector=NULL;
	long VectorCount=0;
    
    short rn = GetAmbVectorDF(Amb, AmbCount, NULL, &VectorCount, Dx); //���ȼ�������ģ����������
	if(rn==0) return 0;
	if(VectorCount>2000000 || VectorCount<=0 ) return 0;
	
	AmbVector = new double* [VectorCount];
	
	for(i=0; i<VectorCount; i++)
	{
		AmbVector[i]=new double[AmbCount];
	}

	//������е�ģ����������
	GetAmbVectorDF(Amb, AmbCount, AmbVector, &VectorCount, Dx);
	
	/*if(VectorCount>1)
	{	//����������
		FilterAmbVectorSF(AmbVector, &VectorCount, Amb, AmbCount, Dx);
	}*/

	double* sigma=new double[VectorCount];
	for(i=0;i<VectorCount;i++)
	{
		sigma[i] = GetSigma(Amb, sigma0, AmbVector[i], AmbCount, Dx, DfCount); //����ÿһ�鱸ѡģ���ȵĵ�λȨ�����
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
	//�ҳ����п��ܵ�ģ��������
	long i=0;
	double** AmbVector=NULL;
	long VectorCount=0;
    
    short rn = GetAmbVectorSF(Amb, AmbCount, NULL, &VectorCount, Dx); //���ȼ�������ģ����������
	if(rn==0) return 0;
	if(VectorCount>2000000 || VectorCount<=0 ) return 0;
	
	AmbVector = new double* [VectorCount];
	
	for(i=0; i<VectorCount; i++)
	{
		AmbVector[i]=new double[AmbCount];
	}

	//������е�ģ����������
	GetAmbVectorSF(Amb, AmbCount, AmbVector, &VectorCount, Dx);
	
	if(VectorCount>10000)
	{	//����������
		FilterAmbVectorSF(AmbVector, &VectorCount, Amb, AmbCount, Dx);
	}

	double* sigma=new double[VectorCount];
	for(i=0;i<VectorCount;i++)
	{
		sigma[i] = GetSigma(Amb, sigma0, AmbVector[i], AmbCount, Dx, DfCount); //����ÿһ�鱸ѡģ���ȵĵ�λȨ�����
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

//�����䷨��ȡģ����������
/*in*/
//Amb:�����ģ����
//AmbCount:�����ģ���ȸ���
//Qxx  :�����ģ���ȷ���
/*out*/
//AmbVector:��ѡ��ģ����������
//VectorCount:��ѡ��ģ����������ĸ���
//(ע�����øú���ʱ����ʹ����AmbVectorΪNULL���Ӷ����VectorCount,��ģ���������ĸ���
// Ȼ�����VectorCount�Ĵ�СΪAmbVector�����ڴ棬֮���ڵ��øú��� ��
GPSCLASSDLL_API short GetAmbVectorSF(double *Amb,short AmbCount,double** AmbVector,long* VectorCount,Matrix Dxx)
{
	long i=0, j=0, k=0;
	double  amb1=0,amb2=0;
	short factor=4;
	double sigma=0;
	double deltaN=0;
	
	short *Ni=new short[AmbCount];//ÿ��ģ���ȵı�ѡģ���ȸ���
	double **Ambiguitys=NULL;//ÿ��ģ���ȵı�ѡģ����
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
	
	//�Ը��ֿ��ܵ�ģ���Ƚ�����ϣ��ó����е�ģ��������
    long cpn=1;
	long ii=0;
	*VectorCount=1;
	for(i=0;i<AmbCount;i++)
	{
		(*VectorCount)=(*VectorCount)*Ni[i];//�ó����е�ģ���������ĸ���
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

	//�ͷ��ڴ�
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
//�����䷨��ȡģ����������
/*in*/
//Amb:�����ģ����
//AmbCount:�����ģ���ȸ���
//Qxx  :�����ģ���ȷ���
/*out*/
//AmbVector:��ѡ��ģ����������
//VectorCount:��ѡ��ģ����������ĸ���
//(ע�����øú���ʱ����ʹ����AmbVectorΪNULL���Ӷ����VectorCount,��ģ���������ĸ���
// Ȼ�����VectorCount�Ĵ�СΪAmbVector�����ڴ棬֮���ڵ��øú��� ��
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
	
	short *Ni=new short[nSat];//ÿ�����ǵı�ѡģ���ȸ���	
	
	typedef struct tagDualAmb
	{
		double AmbV1;
		double AmbV2;
	}DualAmb;

	DualAmb **SatAmb = new DualAmb * [nSat];  //ģ���ȵı�ѡģ����
	
    
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
		//L1ģ����
		amb11 = Round(Amb[i]-factor*sqrt(Qxx(i,i)), -1);
		amb12 = Round(Amb[i]+factor*sqrt(Qxx(i,i)), 1);
		//L2ģ����
		amb21 = Round(Amb[i+nSat]-factor*sqrt(Qxx(i+nSat,i+nSat)), -1);
		amb22 = Round(Amb[i+nSat]+factor*sqrt(Qxx(i+nSat,i+nSat)), 1);
        //����L1,��L2��ģ����
		deltan = Amb[i] - (77.0/60.0) * Amb[i+nSat];   //xn1 - f1/f2 *xn2  //����ģ����
	    
		m0 = factor * sqrt(fabs(Qxx(i,i)+Qxx(i+nSat,i+nSat)-2*Qxx(i,i+nSat))); 
		minDelta = Round((deltan-m0),-1);
		maxDelta = Round((deltan+m0),1);
		
		short count=0;
        
		for(Nl1=amb11; Nl1<=amb12; Nl1++)    //L1���ܵ�ģ����ֵ
		{
			for(Nl2=amb21; Nl2<=amb22; Nl2++)  //L2���ܵ�ģ����ֵ
			{
				deltaN = Nl1 - (77.0/60.0) * Nl2;
				
				if(deltaN >= minDelta && deltaN <= maxDelta)  //����ģ���ȷ��ϱ�ѡģ��������,����
				{
					count++;
					Ni[i] = count;   //��i�����ǵ���ϸ���				
				}
			}
		}
		
		//��¼����ģ����
		SatAmb[i] = new DualAmb[count];
		short ijk=0;
		for(Nl1=amb11; Nl1<=amb12; Nl1++)    //L1���ܵ�ģ����ֵ
		{
			for(Nl2=amb21; Nl2<=amb22; Nl2++)  //L2���ܵ�ģ����ֵ
			{
				deltaN = Nl1 - (77.0/60.0) * Nl2;
				if(deltaN >= minDelta && deltaN <= maxDelta)  //����ģ���ȷ��ϱ�ѡģ��������,����
				{
					SatAmb[i][ijk].AmbV1 = Nl1;
					SatAmb[i][ijk].AmbV2 = Nl2;
					ijk++;
				}
			}
		}
	}

	//�Ը��ֿ��ܵ�ģ���Ƚ�����ϣ��ó����е�ģ��������
    *VectorCount=1;
	
	for(i=0; i<nSat; i++)
	{
		(*VectorCount) = (*VectorCount)*Ni[i];//�ó����е�ģ���������ĸ���
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
	//��ģ��������������
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
	//�ͷ��ڴ�
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
//���˫���Э������
int GetDDCoVariance(CGPSEpochObsSet *pFEpochObs, CGPSEpochObsSet *pUEpochObs, 
					CGPSEpochDDCoef *pEpochDDCoef, 
					Matrix &Dx, CCalcSolutionBase *pSolutinBase)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	Matrix Dl(nSatCount, nSatCount); //����Э������

	switch (pSolutinBase->m_WeightModel)
	{
	case WEIGHT_WITH_SIGNAL_STRENGTH:  //�ź�ǿ���󷽲�
		{
			GetDDCoVarianceSignalStrength(pFEpochObs, pUEpochObs, pEpochDDCoef, Dl);
			Dx = Dl;
		}
		break;
	case WEIGHT_WITH_ELEV:				//�߶ȽǶ�Ȩ
		{
			GetDDCoVarianceElev(pFEpochObs, pUEpochObs, pEpochDDCoef, Dl);
			Dx = Dl;
		}
		break;
	/*case WEIGHT_WITH_CNO:				//Cno��Ȩ
		{
			WeightWithCNo(pUEpochObs,pFEpochObs, EpochDDCoef,refsat, P, satcount); 
		}
		break;
	case WEIGHT_WITH_STRENGTH_DELTA:				//����ȶ�Ȩ
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
//���õ�λȨģ�ͻ��˫��۲�ֵ�ķ���
void GetDDCoVarianceUnit(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	short i=0, j=0;
	double m02 = 0.01 * 0.01;  //�۲�ֵ���ȸ�1cm
	Matrix Ql(nSatCount,nSatCount); //Э������

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
	Dll = Ql * m02;  //ת����Э������
}
//-
//���ø߶ȽǶ�ģ�ͻ��˫��۲�ֵ�ķ���
void GetDDCoVarianceElev(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short nSatCount = pEpochDDCoef->GetItemCount();
	double m02 = 0.01 * 0.01;  //
	short j=0, k=0;
	Matrix Ql(nSatCount, nSatCount); //Э������


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
	Dll = Ql * m02;  //ת����Э������
}
//-
//�����ź�ǿ��ģ�ͻ��˫��۲�ֵ�ķ���
void GetDDCoVarianceSignalStrength(CGPSEpochObsSet *pFEpochObs, 
						 CGPSEpochObsSet *pUEpochObs, 
						 CGPSEpochDDCoef *pEpochDDCoef, Matrix& Dll)
{
	short satCount = pEpochDDCoef->GetItemCount();
	Matrix Dl(satCount,satCount);
	double CL1=0.224; //0.224�൱��2.7mm�ľ���
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
//վ��(S)�ռ�ֱ�����������(E)�ռ�ֱ������ת������
void SxyztoExyzMatrix(double B, double L, Matrix &tranMatrix)
{
	tranMatrix(0,0) = -cos(L)*sin(B); tranMatrix(0,1) = -sin(L); tranMatrix(0,2) = cos(L)*cos(B);
	tranMatrix(1,0) = -sin(L)*sin(B); tranMatrix(1,1) = cos(L);  tranMatrix(1,2) = sin(L)*cos(B);
	tranMatrix(2,0) = cos(B);         tranMatrix(2,1) = 0;       tranMatrix(2,2) = sin(B);
}
//-
//����(E)�ռ�ֱ��������վ��(S)�ռ�����ת������
void ExyztoSxyzMatrix(double B, double L, Matrix &tranMatrix)
{
	tranMatrix(0,0) = -cos(L)*sin(B); tranMatrix(0,1) = -sin(L)*sin(B); tranMatrix(0,2) = cos(B);
	tranMatrix(1,0) = -sin(L);        tranMatrix(1,1) = cos(L);         tranMatrix(1,2) = 0;
	tranMatrix(2,0) = cos(B)*cos(L);  tranMatrix(2,1) = cos(B)*sin(L);  tranMatrix(2,2) = sin(B);
}

//==========================================
// tr: �źŽ���ʱ��
// pEph: ������������
// Xr, Yr, Zr�����ջ��Ľ���λ��
// Xk, Yk, Zk, Xdot, Ydot, Zdot: ����λ�á��ٶȣ�ƽ�����
//�ɽ��ջ��Ľ���λ�ã��������ǵķ���ʱ��
double GetSatSentTime(double tr, ephemeris *pEph, 
					  double Xr, double Yr, double Zr,
					  double* Xk, double* Yk, double* Zk, double* Ek,
					  double* Xdot, double* Ydot, double* Zdot, double* range)
{
	double ts0 = 0;
	double ts = 0;  //�źŷ���ʱ��
	double deltT = 0;
	double satDt = 0; //�����Ӳ�
    double dift = 0;

    ts = tr - 0.075;

	do{
		//�������Ƿ���ʱ������ʱ�����������Ӳ�
		deltT = ts - pEph->toe;
		if(deltT>302400.00) //�������ڽ����Ŀ�ʼ�ͽ���ʱ��
		{
			deltT -= 604800.0;
		}
		if(deltT<-302400.0) 
		{
			deltT += 604800.0;
		}
		
		//�����Ӳ�
		satDt = pEph->a0 + (pEph->a1 * deltT) + (pEph->a2 * deltT * deltT);
        
		//���Ƿ���ʱ��
		ts0 = ts + satDt;
		
		//��������λ��
		coordsat(*pEph, ts0, Xk, Yk, Zk, Ek, Xdot, Ydot, Zdot);
		*range = sqrt((*Xk-Xr)*(*Xk-Xr)+(*Yk-Yr)*(*Yk-Yr)+(*Zk-Zr)*(*Zk-Zr));
		
		//���¼�����������źŷ���ʱ��
		ts = tr - *range / C;
		dift = (ts+satDt) - ts0;
	}while(fabs(dift)>1.0E-7);  //�����μ���ķ���ʱ��С��10-7��ʱ����������
	
	return ts;
}

void PosToPos(double antN, double antE, double antH, 
		          double B, double L, 
				  double X, double Y, double Z,
				  double &newX, double &newY, double &newZ, short type)
{
	Matrix refTran(3,3); //ת������(վ�ģ�>�ռ�)
	Matrix antfNEH(3,1); //����ƫ����(վ��)
	Matrix antfXYZ(3,1);  //����ƫ����(�ռ�)
	
	//������ƫ������վ��ת������
	SxyztoExyzMatrix(B, L, refTran);

	antfNEH(0,0) = antN;
	antfNEH(1,0) = antE;
	antfNEH(2,0) = antH;
	
	antfXYZ = refTran*antfNEH;  //�ռ�ֱ�������µ����߸�����

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

//�������ܣ������ֲ�
//���ݹ۲�ֵ����������ˮƽ�����ؿ�������ֵ

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
	//yu:����ļ�Ϊ�վ��˳�
	if (!fobsfile)
	{
		return;
	}
	
	CGPSTime gpsT(obs.m_GpsWeek,obs.m_GpsSecond);
	int yeari = gpsT.m_year % 100;
	fprintf(fobsfile," %02d%3d%3d%3d%3d%11.7f%3d%3d",
		yeari, gpsT.m_month,gpsT.m_day,
		gpsT.m_hour, gpsT.m_min, gpsT.m_sec, 
		0, obs.GetSatCount()); //������ʱ���� ���¼���־,������
	int i = 0;
	
	if(obs.GetSatCount()<=12)  //���С��12������
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
	else        //����12������
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
	
	//�������Ǹ߶Ƚ�
	CalcStationElevation(&_obsset, &_navset, &stationInfo);

	//�ɱ�׼����Ԫ�ؼ���
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
	