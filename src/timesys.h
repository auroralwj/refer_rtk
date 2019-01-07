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
	double GetYearDays(long year); //获得一年有多少天
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
//把GPS时间转换为  年,月,日,时,分,秒--4-2-2-2-2-other..
 double GetYMDHMS(double gpsweek,double gpssecond, char * timeStr);

 void GetYMDHMSByStr(char * timeStr, short *year,
								   short  *month,short *day,short *hour,short *minute,double *second);

  void DateTimeToWeekTOW(short year, short month, short day, short hour, short min, double sec, 
                             long * Week, double * TOW);

  void DateTimeToWeekTOW(char * timeStr, long *Week, double *TOW);

double GetYMDHMSTimeReportDayTime(double gpsTime, char * DateTimeStr);

// void DateTimeToWeekTOW(DATETIMEEX DateTime, long *Week, double *TOW);

//时区转换,从本地时间换算到gps时间
// year month, day, hour ,min ,sec 本地时间
//hourOffset  相减时间,入香港本地时间与gps时间差8小时 hourOffset = 8
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