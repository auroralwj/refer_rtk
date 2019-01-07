// FileSearcher.cpp: implementation of the CFileSearcher class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileSearcher.h"
#include "TimeSys.h"
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CFileSearcher::CFileSearcher()
{

}

CFileSearcher::~CFileSearcher()
{

}

//bSubDir 是否查询子目录
int CFileSearcher::SearchFile(CString strRstPath, CString fileType, BOOL bSubDir)
{
	if (strRstPath.IsEmpty())
		return 0;
	CString strSubPath;
	/*** set the path format in C++/C custom ***/
    if(strRstPath.Right(3) != "*.*")
	{
		if(strRstPath.Right(1) != "\\" )
		{
			strRstPath+="\\*.*";
		}
		else
		{
			 strRstPath+="*.*";
		}
	}
	
   CFileFind filesearch;
   CString strRstFile, strRstFileName,strStaID;
   BOOL bFind=filesearch.FindFile(strRstPath);
  
   int typeSize = fileType.GetLength();
   fileType.MakeUpper();
   while(bFind)
   {
	   bFind = filesearch.FindNextFile();	
	   
	   if(!filesearch.IsDots())  //
	   {
		   if(filesearch.IsDirectory () && bSubDir)
		   {
			   strSubPath=filesearch.GetFilePath();
			   SearchFile(strSubPath, fileType, TRUE);
		   }
		   else
		   {				
			   //the result file
			   strRstFile=filesearch.GetFilePath();
			   strRstFileName=filesearch.GetFileName();
			   strRstFileName.MakeUpper();
			   if (fileType == strRstFileName.Right(typeSize)
				   || fileType == _T(""))
			   {
				   Add(strRstFile);
			   }
		   }
	   }
   }
   filesearch.Close();
   SortByFileName();
   return GetSize();
}

BOOL CFileSearcher::CompareAndSwapByFileName( int pos )
{
   CString temp;   
   int posFirst = pos;   
   int posNext = pos + 1;
   CString pHead = GetAt(posFirst);
   CString pNext = GetAt(posNext);
   int idx = pHead.Find('\\');
   CString posHeadFile = pHead.Right(pHead.GetLength() - idx- 1);
   idx = pNext.Find('\\');
   CString posNextFile = pNext.Right(pNext.GetLength() - idx- 1);
   if (posHeadFile.CompareNoCase(posNextFile) > 0)   
   {
      temp = pHead;      
	  SetAt(posFirst, pNext);
      SetAt(posNext, temp);      
	  return TRUE;   
   }   
   return FALSE;
}

void CFileSearcher::SortByFileName()
{
	BOOL bNotDone = TRUE;   
	while (bNotDone)   
	{      bNotDone = FALSE;
      for(int pos = 0;pos < GetUpperBound();pos++)
        bNotDone |= CompareAndSwapByFileName(pos);   
	}
}

int CFileSearcher::GetFilesByTimeRange(double gpsStart, double gpsEnd, CStringArray & files)
{
	/*
	CGPSTime gpsStaT(gpsStart);
	CGPSTime gpsStaE(gpsEnd);
	long startDayIndex, endDayIndex;
	char startHour, endHour;
	int startYear = gpsStaT.m_year;
	int endYear = gpsStaE.m_year;
	gpsStaT.ToDayAndHour(& startDayIndex, & startHour);
	gpsStaE.ToDayAndHour(& endDayIndex, & endHour);
	int size = GetSize();
	for (int i = 0 ; i < size ; i ++)
	{
		CString obsFile = GetAt(i);
		int inx = obsFile.Find('.');
		CString sTmp = obsFile.Mid(inx + 1, 2);
		int iYear = atoi(sTmp);
		if (0 == iYear)
		{
			TRACE("File Name Error %s.", obsFile);
			return -1;
		}
		sTmp = obsFile.Mid(4, 3);
		int iDayInex = atoi(sTmp.GetBuffer(0));
		sTmp.ReleaseBuffer();
		if (0 == iDayInex)
		{
			TRACE("File Name Error %s.", obsFile);
			return -1;
		}
		char cHour = obsFile.GetAt(inx + 3);
		int iHour = cHour - 'A';
	}
	*/
	return 1;
}