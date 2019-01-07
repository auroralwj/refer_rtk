#include"../stdafx.h"
#include "timesys.h"
#include <math.h>

#include "constant.h"
//////////////////////////////////////////////////////////////////////
//计算一年的总天数

CGPSTime::CGPSTime()
{
	CTime curTime = CTime::GetCurrentTime();
	m_year =  curTime.GetYear();
	m_month =  curTime.GetMonth();
	m_day =  curTime.GetDay();
	m_hour =  curTime.GetHour();
	m_min =  curTime.GetMinute();
	m_sec =  curTime.GetSecond();
}

CGPSTime::CGPSTime(CString strTime)
{
	SetTime2(strTime);
}

CGPSTime::~CGPSTime()
{

}

CGPSTime::CGPSTime(short year, short month, short day, 
				   short hour, short min, double sec)
{
	m_year =  year;
	m_month =  month;
	m_day =  day;
	m_hour =  hour;
	m_min =  min;
	m_sec =  sec;
}

CGPSTime::CGPSTime(long gpsWeek, double gpsSec)
{
	GpsTimeToDay(gpsWeek, gpsSec,
		& m_year, & m_month, & m_day, & m_hour, & m_min , & m_sec);
}

CGPSTime::CGPSTime(double gpsTime)
{
	GpsTimeToDay(gpsTime,
		& m_year, & m_month, & m_day, & m_hour, & m_min , & m_sec);
}

double CGPSTime::GetYearDays(long year)
{
	double   total,check;	
	check =fmod (year,4.);
	if (check == 0.)
		  {
		check = fmod(year,100.);
		if (check == 0.)
		{
			check = fmod(year,400.);
			if (check == 0.)
				total = 366.;
			else
				total = 365.;
		}
		else
			total = 366.;
		  }
	else
		total = 365.;
	return(total);
}

//把GPS时间转换为  年,月,日,时,分,秒
//返回值为一天中的秒
double CGPSTime::GpsTimeToDay(double gpsTime,short *year, short  *month,
				 short *day,short *hour,short *minute,double *second)
{
	unsigned long   lnggpsTime = (unsigned long) gpsTime;
	long gpsWeek = (int)(lnggpsTime /  SEC_COUNT_IN_WEEK);
	double gpssecond = gpsTime - gpsWeek * SEC_COUNT_IN_WEEK;
	//
	return GpsTimeToDay(gpsWeek, gpssecond,  year, month, day,
					 hour, minute, second);
}


double CGPSTime::GpsTimeToDay(double gpsweek,double gpssecond,short *year,
								 short  *month,short *day,short *hour,short *minute,double *second)
{
	double totsec,totday,daysecond,hoursecond;
	
	totsec=gpsweek*SEC_COUNT_IN_WEEK+gpssecond;
	totday=(int ) (totsec/ 86400);
	daysecond= totsec - totday * 86400;
	
	*hour= (int ) (daysecond / 3600);
	hoursecond=daysecond - *hour * 3600;
	
	*minute= (int) (hoursecond / 60);
	*second=hoursecond - *minute * 60;
	
	if (*second > 59.9999)
	{
		*second = 0;
		*minute = *minute + 1;
		if (*minute> 59.9999)
		{
			*minute = 0;
			*hour = *hour + 1;
		}
	}
	
	double yeardays;
    double daycount=0;
	long yy;
	for(yy=1980;;yy++)
	{
		daycount+=GetYearDays(yy);//由年份确定的天数
		if(daycount>totday+6)//GPS时间从1980年1月6日协调世界时开始计算
		{
			yeardays=totday+6-daycount+GetYearDays(yy);
			break;
		}
	}
	*year=yy;
	double regu_month_day[13]={ 0, 31, 59, 90, 120, 151, 181, 212,
		243, 273, 304, 334, 365};
    double leap_month_day[13]={ 0, 31, 60, 91, 121, 152, 182, 213,
		244, 274, 305, 335, 366 };
	
    short guess=short(yeardays*0.032)+1;  
	if(fabs(GetYearDays((long)*year)-365)<0.0000000001)//考虑润年的2月份
	{
		if (guess < 13 && guess >=0)
		{
			if((yeardays-regu_month_day[guess])>0)
			{
				*month=guess+1;
			}
			else
			{
				*month=guess;
			}
		}
		else
		{
			
//	FILE * hFile = fopen("c:\\debug.txt", "a+");
//	fprintf(hFile, "%lf\t%lf\n", gpsweek, gpssecond);
//	fclose(hFile);
//	hFile = NULL;
			*month=guess;
		}
		*day = yeardays-regu_month_day[(short)*month-1];
		
	}
	else
    {
		if (guess < 13 && guess >=0)
		{
			if((yeardays-leap_month_day[guess])>0)
			{
				*month=guess+1;
			}
			else
			{
				*month=guess;
			}
		}
		else
		{
			
//	FILE * hFile = fopen("c:\\debug.txt", "a+");
//	fprintf(hFile, "%lf\t%lf\n", gpsweek, gpssecond);
//	fclose(hFile);
//	hFile = NULL;
			*month=guess;
		}
		*day=yeardays-leap_month_day[(short)*month-1];
	}
	
	/////////////////////////////	
	if (*second - floor(*second) < 0.01)
	{
		*second = floor(*second);
	}

	if (*second > 59.9999)
	{
		*second = 0;
		*minute ++;
	}
	
	if (ceil(*second) - *second < 0.01)
	{
		*second = ceil(*second);
	}
	
	if (*minute > 59.9999)
	{
		*minute = 0;
		*hour ++;
	}
	
	if (*hour> 23.999999)
	   {
		   *hour = 0;
		   *day = *day + 1;
	   }
	   
	   BOOL bcorrentFlag = FALSE;
	   if (*day == 0)
	   {
		   *month = *month - 1;
		   bcorrentFlag = TRUE;
	   }
	   
	   if (*month == 0 && bcorrentFlag)
	   {
		   *month = 12;
		   *year = *year - 1;
		   bcorrentFlag = TRUE;
	   }
	   
	   if (bcorrentFlag)
	   {
		   int iMonth = (int) *month;
		   if (iMonth < 8 )
		   {
			   if ((iMonth % 2) == 1)
				   *day = 31;
			   else
				   *day = 30;
		   }
		   else
		   {
			   if ((iMonth % 2) == 0)
				   *day = 31;
			   else
				   *day = 30;
		   }
	   }
	   return daysecond;
}

/*---------------------------------------------------------------------------*/
/*
* Input:   DateTime containing the Date and Time in GPS time system
* Output:
*   *Week: The corresponding GPS week number (from Jan 06,1980)
*
*   *TOW:  The corresponding GPS seconds elapsed since begin of the
*          week (TOW)
*/
double CGPSTime::ToGpsTime()
{
	long week;
	double tow;
	ToGpsTime(& week, & tow);
	return week * SEC_COUNT_IN_WEEK + tow;
}

void CGPSTime::ToGpsTime(long *Week, double *TOW)
{
	short DayCount[13] = {0,31,59,90,120,151,181,212,
		243,273,304,334,365};
	double j;
	
	if (m_year<80) m_year += 2000;
	else if (m_year<100) m_year += 1900;
	
	j = (m_year-1980)*365.25;
	
	if ((m_year%4) != 0) j += 1.0 - (m_year%4)*0.25;
	
	j += DayCount[m_month-1]+m_day-1;
	
	if (((m_year%4)==0) && (m_month>=3)) j++;
	
	j -= 5;
	
	*Week = (long)j/7;
	*TOW = fmod(j,7.0)*86400.0 
		+ 3600.0*(double)m_hour 
		+ 60.0*(double)m_min + m_sec;
}
/*---------------------------------------------------------------------------*/
/*
* Input:   DateTime containing the Date and Time in GPS time system
* Output:
*   *Week: The corresponding GPS week number (from Jan 06,1980)
*
*   *TOW:  The corresponding GPS seconds elapsed since begin of the
*          week (TOW)
*/
void CGPSTime::ToDayAndHour(long * dayIndex, char * hourIndex)
{
	
	short DayCount[13] = {0,31,59,90,120,151,181,212,
		243,273,304,334,365};
	double j = 0;
	
	if ((m_year%4) != 0) j += 1.0 - (m_year%4)*0.25;
	
	j += DayCount[m_month-1]+m_day-1;
	
	if (((m_year%4)==0) && (m_month>=3)) j++;
		
	* dayIndex = (int) (j + 1);
	*hourIndex = 'a' + m_hour;
}


int CGPSTime::GetDayIndex()
{
	
	short DayCount[13] = {0,31,59,90,120,151,181,212,
		243,273,304,334,365};
	double j = 0;
	
	if ((m_year%4) != 0) j += 1.0 - (m_year%4)*0.25;
	
	j += DayCount[m_month-1]+m_day-1;
	
	if (((m_year%4)==0) && (m_month>=3)) j++;
		
	return (int) (j + 1);
}

CString CGPSTime::ToTimeFormatStr()
{
	char timeStr[256];
	memset(timeStr, 0, 256);
	sprintf(timeStr, "%4d-%02d-%02d %02d:%02d:", m_year, m_month, m_day, m_hour, m_min);  //zxw 25/2/2008
	double fmData = fmod(m_sec, 1.0);
	if (m_sec >= 10)
		sprintf(timeStr + 17, "%2.f", m_sec);
	else
	{
		sprintf(timeStr + 17, "0%1.f", m_sec);
		if (strlen(timeStr) > 22)
		{
			sprintf(timeStr + 17, "%2.f", m_sec);
			timeStr[22] = 0;
		}
	}
	ASSERT(FALSE);
	return _T("");
}

double CGPSTime::ToReportDayHour(char * dateTimeStr)
{
	sprintf(dateTimeStr, "%4d/%02d/%02d %02d %02d:%02d", m_year, m_month, m_day, m_hour, 59, 59);
	///////////////////////
	return 0;
}

double CGPSTime::ToReportDayTime(char * dateTimeStr)
{
	sprintf(dateTimeStr, "%4d/%02d/%02d %02d %02d:%02.f", m_year, m_month, m_day, m_hour, m_min, m_sec);
	///////////////////////
	return 0;
}

BOOL CGPSTime::SetTime2(CString strTime)
{
	return FALSE;
}

BOOL CGPSTime::SetTime(CString strTime)
{
	return FALSE;
}

