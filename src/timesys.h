//
#ifndef TIMESYS_H_GGGGG
#define TIMESYS_H_GGGGG

#include "constant.h"

class CGPSTime
{
public:
	CGPSTime();
	CGPSTime(CString strTime);
	CGPSTime(short year, short month, short day, short hour, short min, double sec);
	CGPSTime(long gpsWeek, double gpsSec);
	CGPSTime(double gpsTime);
	virtual ~CGPSTime();
	double GetYearDays(long year); //���һ���ж�����
	void ToDayAndHour(long * dayIndex, char * hourIndex);
	int GetDayIndex();
	void ToGpsTime(long *Week, double *TOW);
	double ToGpsTime();
	BOOL SetTime(CString strTime);
	BOOL SetTime2(CString strTime);
	CString ToTimeStr();
	CString ToTimeFormatStr();
	short m_year; 
	short m_month;
	short m_day;
	short m_hour;
	short m_min; 
	double m_sec; 

	double ToReportDayTime(char * dateTimeStr);
	double ToReportDayHour(char * dateTimeStr);
private:
	double GpsTimeToDay(double gpsTime,short *year, short  *month,
				 short *day,short *hour,short *minute,double *second);
	double GpsTimeToDay(double gpsweek,double gpssecond,
						short *year, short  *month,short *day,
						short *hour,short *minute,double *second);
	
};
/*

 double GetYMDHMS2(double gpsweek,double gpssecond, char * timeStr);
//��GPSʱ��ת��Ϊ  ��,��,��,ʱ,��,��--4-2-2-2-2-other..
 double GetYMDHMS(double gpsweek,double gpssecond, char * timeStr);

 void GetYMDHMSByStr(char * timeStr, short *year,
								   short  *month,short *day,short *hour,short *minute,double *second);

  void DateTimeToWeekTOW(short year, short month, short day, short hour, short min, double sec, 
                             long * Week, double * TOW);

  void DateTimeToWeekTOW(char * timeStr, long *Week, double *TOW);

double GetYMDHMSTimeReportDayTime(double gpsTime, char * DateTimeStr);

// void DateTimeToWeekTOW(DATETIMEEX DateTime, long *Week, double *TOW);

//ʱ��ת��,�ӱ���ʱ�任�㵽gpsʱ��
// year month, day, hour ,min ,sec ����ʱ��
//hourOffset  ���ʱ��,����۱���ʱ����gpsʱ���8Сʱ hourOffset = 8
//short LocalTimeToGPSTime(short hourOffset, short year, short month, short day, 
//											  short hour, short min, double sec, 
//											  long &GpsWeek, double &GpsSecond);

// double GetYMDHMS(double gpsTime,double *year,
//								 double  *month,double *day,double *hour,double *minute,double *second);


 void GetGPSTimeByStr(char * timeStr, double * gpsTime);

 void DateTimeToDayAndHour(short year, short month, short day, short hour, short min, double sec, 
										long * dayIndex, char * hourIndex);
										*/
#endif