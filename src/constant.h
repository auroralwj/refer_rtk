//: constants.h
//
/*
#ifdef CREATE_GPS_DLL
#define GPSCLASSDLL_API __declspec(dllexport)
#else
#define GPSCLASSDLL_API __declspec(dllimport)
#endif
*/
#define GPSCLASSDLL_API

#ifdef CREATE_GPS_DLL
#define GPSCLASSDLL_API2 __declspec(dllexport)
#else
#define GPSCLASSDLL_API2 __declspec(dllimport)
#endif

#ifndef CONSTANTS_H
#define CONSTANTS_H //Debug Begin

/***************************************************************************/
/*** ------------------------------------------------------------------- ***/
/***  constants.h: the header file to define some constant parameters    ***/
/***************************************************************************/

#include <cmath>
#include "..\GNSSPRODLL\src\rtklib.h"

typedef int (* RecvDataFun) (BYTE * pData, int lgh);
typedef int (* ToRinexBackFun) (char * pBuffer);

#define MAX_COMMAND_LEGHT	1025
#define MAX_BUFFER_SIZE 1025
#define MAX_RECV_BUFFER_SIZE 1025
//#define MAX_RECV_BUFFER_SIZE 400
#define MAX_RAWTORINEX_BUFFER_SIZE 102400
#define ARRAY_COUNT		2
#define MAX_SNV_COUNT		40
const double vC				=	2.99792458e+8;         //��������е��ٶ�(metres/sec)


typedef int (* RecvDataFun) (BYTE * pData, int lgh);

typedef BOOL (* DoSendCmdFun) (int cmdType);

enum{
	ADD_ERROR=0,
	ADD_OK,			
	NEED_SNV,
};

const double EPS		= 1e-13;			//small number
const double NaN		= -999;             //Not a number

/***** File Format *****/
const unsigned int BINARY  = 1;        //Binary
const unsigned int ASCI    = 2;        //Ascii

/***** Maximum Line Lenght ******/
const long MAXILINELEN = 150;
const long MAXIFILENAMELEN = 100;

/***** Project Mode *****/
const unsigned int POST		=	1;
const unsigned int AUTO		=	2;

/***** Process Mode *****/
const int STATIC  = 1;
const int SNGEPC  = 2;
const int RTK     = 3;

/****** Coordinate System  *******/
const unsigned int WGS84   =  0;
const unsigned int GAUSS   =  1;

/***** Ambiguity Method ***********/
/*
const unsigned int FARA		=	0;
const unsigned int LAMBDA	=	1;
const unsigned int ARCE		=	2;   //??? Different with Dai's
const unsigned int FAST		=	3;
const unsigned int INTG		=	4;
const unsigned int LSSSC	=	5;
const unsigned int LSSAMC	=	6;
*/

/***** Observation Type ************/
const unsigned int L1AL2	=	0;                    // L1 && L2
const unsigned int L1		=	1;                    // L1
const unsigned int L2		=	2;                    // L2
const unsigned int L1PL2	=	3;                    // L1+L2
const unsigned int L1ML2	=	4;                    // L1-L2
const unsigned int L5		=	5;                    // L5


/**** PI *****/
const double PI2			=	PI/2.0;
const double rad			=	PI/180.0;            //����ת��ϵ��

/*** earth const parameters ***/
const double GM				=	3.986004418e+14;          //������������
const double W				=	7.2921151467e-5;       //������ת���ٶ�(rad/sec)
const double C				=	2.99792458e+8;         //��������е��ٶ�(metres/sec)
const double a				=	6378137.0;             //WGS84��ת����ĳ��뾶(metres)
const double b				=	6356752.3142;          //WGS84��ת����Ķ̰뾶(metres) 
const double F				=	-4.442807633e-10;      //�����ЧӦ���� F=-SQU(u)/(C*C)(sec/metre_0.5)

/*** wave length ***/
const double F1				=	1575.42;               //L1��Ƶ��(MHz)
const double F2				=	1227.6;                //L2��Ƶ��(MHz)
const double F11DF2			=	((F1*F1)/(F2*F2));               //L1��Ƶ��(MHz)
const double WAVE1			=	((C*1.0e-6)/F1);        //f1�Ĳ���
const double WAVE2			=	((C*1.0e-6)/F2);        //f2�Ĳ���
const double WAVELC			=	WAVE1;					//LC�Ĳ���
const double WAVEN			=	((C*1.0e-6)/(F1+F2));   //(f1+f2)�Ĳ���
const double WAVEW			=	((C*1.0e-6)/(F1-F2));   //(f1-f2)�Ĳ���

//some const
const unsigned short	MAX_CHANCEL_COUNT	=	36;
const unsigned short	MAX_MULREF_COUNT	=	10;
const unsigned short	MAXEPOCHSATCOUNT	=	24;
const unsigned short	MAXPRNID			=	36;
const unsigned short	MAXGEOSTATIONARYID	=	99;
const unsigned short	MAXSATPEREPOCH		=	24;
const unsigned short	RINEXRECSIZE		=	83;   // 80 cols plus \r \n etc.
const unsigned short	MAXOBSTYPES			=	11;
const  long				MAXEPOCHINSET		=	9999999; //��¼���������ܵ���Ԫ��
const unsigned long		SEC_COUNT_IN_WEEK   =   604800;

enum{CALCFINISH = 0, CALCOK, CALCERROR, CALC_NEEDOBS, CALC_CANCEL};

enum{READOK = 1, 
	READEND,	//����
	READFINISH,	//�Ѿ����
	READERROR
};

enum{
	CALC_ERROR_1 = 1,   //����1
	CALC_NORMAL,
	CALC_CANCEL_D,			//ȡ������
	CALC_CALC_DATA_INFO,	//��������
	CALC_CIRCLESLIP_FOUND,	//��������
	CALC_EXIT_END,			//�������
	CALC_RESULT_SHOW,		//��������ʾ
};

enum{
	CHOOSESAT_HIGHEST_ELEV_EACH_EPH = -1,
	CHOOSESAT_LOGEST_EPH = 0,
	CHOOSESAT_1,
	CHOOSESAT_2,
	CHOOSESAT_3,
	CHOOSESAT_4,
	CHOOSESAT_5,
	CHOOSESAT_6,
	CHOOSESAT_7,
	CHOOSESAT_8,
	CHOOSESAT_9,
	CHOOSESAT_10,
	CHOOSESAT_11,
	CHOOSESAT_12,
	CHOOSESAT_13,
	CHOOSESAT_14,
	CHOOSESAT_15,
	CHOOSESAT_16,
	CHOOSESAT_17,
	CHOOSESAT_18,
	CHOOSESAT_19,
	CHOOSESAT_20,
	CHOOSESAT_21,
	CHOOSESAT_22,
	CHOOSESAT_23,
	CHOOSESAT_24,
	CHOOSESAT_25,
	CHOOSESAT_26,
	CHOOSESAT_27,
	CHOOSESAT_28,
	CHOOSESAT_29,
	CHOOSESAT_30,
	CHOOSESAT_31,
	CHOOSESAT_32	
}; //ѡ����Ƿ���


enum{
	AMB_INTEGER = 0,
	AMB_ARCE,
	AMB_FARA,
	AMB_LAMBDA,
////
	AMB_SECTION,
	AMB_QIF,
	AMB_FARAEX,
	AMB_LAMBDAEX,
//
	AMB_FAST,	

	AMB_LSSSC,
	AMB_LSSAMC,
};

enum{
	CLAC_USE_L1 =  1,	// for sgn, 
	CLAC_USE_L1L2,		// for sgn, 
}; //ѡ�񲨷�ʽ

//##ModelId=462C195302DE
typedef int (* CallBackMsgShow)(char * msgInfo, int infoType);
//##ModelId=462C195302E0
typedef int (* CallBackExitCalc)();		//�˳�����

#endif// Debug End
///:~
