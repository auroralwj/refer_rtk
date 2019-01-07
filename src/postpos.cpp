/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
*          Copyright (C) 2007-2014 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioing api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*-----------------------------------------------------------------------------*/

#include"../stdafx.h"
#include "rtklib.h"

static const char rcsid[]="$Id: postpos.c,v 1.1 2008/07/17 21:48:06 ttaka Exp $";

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */

/* constants/global variables ------------------------------------------------*/

static pcvs_t pcvss={0};        /* receiver antenna parameters */
static pcvs_t pcvsr={0};        /* satellite antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sbs_t sbss={0};          /* sbas messages */
static lex_t lexs={0};          /* lex messages */
static sta_t stas[MAXRCV];      /* station infomation */
static int nepoch=0;            /* number of observation epochs */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int isbs  =0;            /* current sbas message index */
static int ilex  =0;            /* current lex message index */
static int revs  =0;            /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */
static char rtcm_file[1024]=""; /* rtcm data file */
static char rtcm_path[1024]=""; /* rtcm data path */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE *fp_rtcm=NULL;      /* rtcm data file pointer */

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    //if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    //return showmsg(buff);
	return 0;
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];
    const char *sep=opt->sep;
    
    trace(3,"outrpos :\n");
    
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1);
            deg2dms(pos[1]*R2D,dms2);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];
    
    trace(3,"outheader: n=%d\n",n);
    
    if (sopt->posf==SOLF_NMEA) return;
    
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
		SYSTEMTIME st; 
		GetLocalTime(&st); 
		int year=st.wYear,month=st.wMonth,day=st.wDay,hour=st.wHour,minute=st.wMinute,sec=st.wSecond;
		fprintf(fp,"%s cal time  : %4d-%02d-%02d %02d:%02d:%02d\n",COMMENTH,year,month,day,hour ,minute,sec);
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss.data[i].time;
        te=obss.data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);
    
    outsolhead(fp,sopt);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}


//Modified by XGAO,2014/7/21
//构建共视卫星列表
static int sat_select(rtk_t *rtk,int *satprn,int *satprn1, int *satprn2, int ns1, int ns2)
{
	int i,j;
	int ns=0;
	//两历元共视卫星
	for(i=0;i<ns1;i++)
	{
		for(j=0;j<ns2;j++)
		{
			if(satprn1[i]==satprn2[j])
				satprn[ns++]=satprn1[i];
		}
	}
	//非共视历元变0
	int count=0;
	for(i=0;i<ns1;i++)
	{
		count=0;
		for(j=0;j<ns;j++)
		{
			if(satprn1[i]!=satprn[j])
				count++;
			else
				break;
		}
		if(count==ns)
			rtk->ssat[satprn1[i]-1].D1_resc[0]=rtk->ssat[satprn1[i]-1].D1_resc[1]=rtk->ssat[satprn1[i]-1].D1_resc[2]=0.0;
	}

	for(i=0;i<ns2;i++)
	{
		count=0;
		for(j=0;j<ns;j++)
		{
			if(satprn2[i]!=satprn[j])
				count++;
			else
				break;
		}
		if(count==ns)
			rtk->ssat[satprn2[i]-1].D2_resc[0]=rtk->ssat[satprn2[i]-1].D2_resc[1]=rtk->ssat[satprn2[i]-1].D2_resc[2]=0.0;
	}
	return ns;
}


//Modified by XGAO,2014/7/18
//三差方程
static void cycleslip_TripeD(rtk_t *rtk,int *sat,int nf,int ns)
{
	int i,j,k,f;
	double lami;
	double resi;
	trace (3,"slip detected using Tripe-Difference method \n");
	for(f=0;f<nf;f++)
	{
		for(i=0;i<ns;i++)
		{
			if(rtk->ssat[sat[i]-1].D1_resc[f]==0 || rtk->ssat[sat[i]-1].D2_resc[f]==0.0)
				continue;
			lami=rtk->ssat[sat[i]-1].lam[f];
			//三差
			resi=rtk->ssat[sat[i]-1].D2_resc[f]-rtk->ssat[sat[i]-1].D1_resc[f];
			rtk->ssat[sat[i]-1].TD_resc[f] = resi;
			//方差
			rtk->ssat[sat[i]-1].TD_Stoc[f] = rtk->ssat[sat[i]-1].D1_Stoc[f]+rtk->ssat[sat[i]-1].D2_Stoc[f];
			//设计矩阵
			for (j=0;j<3;j++) 
				rtk->ssat[sat[i]-1].HH[f][j] = rtk->ssat[sat[i]-1].H2[f][j]-rtk->ssat[sat[i]-1].H1[f][j];  
			//周跳探测
			if ((fabs(resi)>0.5))
			{
				rtk->ssat[sat[i]-1].slip[f]=1;
				//rtk->ssat[sat[i]-1].slip[f]
				rtk->ssat[sat[i]-1].slipc[f]++;
				trace (3,"slip detected with TD method (sat=%2d TD_resi=%10.4f   Slip-Count= %4d)\n",sat[i],fabs(resi),int(resi/lami));

				if(rtk->ssat[sat[i]-1].slipc[f]>15)  //某卫星超过100个历元被探测出周跳
				{
					rtk->opt.exsats[sat[i]-1]=1;      //忽略该卫星exclude
					trace (3,"Excluded satellite with TD method (sat=%2d \n",sat[i]);
				}
			}
			//将当前历元替换上一历元
			rtk->ssat[sat[i]-1].D1_resc[f] = rtk->ssat[sat[i]-1].D2_resc[f];
		}
	}
}

/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    
    return !*outfile?stdout:fopen(outfile,"a");
}


//Modified by XGAO,2014/7/18
//三差法定位计算
static int posi_TripeD(rtk_t *rtk,int *sat,int nf,int ns)
{
	int i,j,k,m,t=0;
	double *H,*P,*D,*D1,*D2,*w,*HP;
	double N[3*3]={0},W[3]={0},dx[3]={0};
	double temp[3]={0};
	int obs_n=0;
	int info;
	int stat=0;
	
    m=(ns-1)*nf;
	if(m<4) return stat;

	trace (3,"Relative positioning with TD method \n");

	//动态申请内存
	H   = zeros(m,3);
	P   = zeros(m,m);   //三差矩阵  P=C*D*Ct
	D   = zeros(m,m);   //三差方差矩阵（对角）
	D1  = zeros(m,m);   //D(t1)
	D2  = zeros(m,m);   //D(t2)
    w   = zeros(m,1);
	HP  = zeros(3,m);   //HP=H'*P
   
	//设计矩阵
	for (i=0;i<nf;i++)for(j=0,t=0;j<ns;j++,t++)    //Row
	{
		if(rtk->ssat[sat[j]-1].refsatflag[i]==1)
		{
			t--;
			continue;
		}
		for(k=0;k<3;k++)temp[k]=rtk->ssat[sat[j]-1].HH[i][k];
		if(norm(temp,3)>0)
		{
			obs_n++;
			for(k=0;k<3;k++) //Column
			{
				H[(t+i*(ns-1))+k*m]=rtk->ssat[sat[j]-1].HH[i][k];
			}
		}
		else
			continue;
	}

	if(obs_n<5)
	{	
		free(H); free(D); free(D1); free(D2); 
		free(HP);free(P); free(w); 
		return stat;
	}
    //trace(3,"H= \n");tracemat(4,H,m,3,12,8);

    //常数项
	for(i=0;i<nf;i++)for(t=0,j=0;j<ns;j++,t++)
	{
		if(rtk->ssat[sat[j]-1].refsatflag[i]==1)
		{
			t--;
			continue;
		}
		w[i*(ns-1)+t] = rtk->ssat[sat[j]-1].TD_resc[i];
	}
	//trace(3,"w= \n");tracemat(4,w,m,1,10,6);
	

	//构建D(t1)
	//双差方差矩阵
	double a1[NFREQ]={0.0};
	double a2[NFREQ]={0.0};
	for(i=0;i<nf;i++)for(t=0,j=0;j<ns;j++,t++)
	{
		if(rtk->ssat[sat[j]-1].refsatflag[i]==1)
		{
			a1[i]=rtk->ssat[sat[j]-1].D1_Stoc[i];
			a2[i]=rtk->ssat[sat[j]-1].D2_Stoc[i];
		}
	}

	//D1
	for(i=0;i<nf;i++)for(t=0,j=0;j<ns;j++,t++)
	{
		if(rtk->ssat[sat[j]-1].refsatflag[i]==1)
		{
			t--;
			continue;
		}
		else
		{
			for(k=0;k<m;k++)
			{
				if((t+i*(ns-1))==k) D1[t*m+i*(ns-1)*m+k] = a1[i]+rtk->ssat[sat[j]-1].D1_Stoc[i];
				else                D1[t*m+i*(ns-1)*m+k] = a1[i];
			}
		}
	}

	//D2
	for(i=0;i<nf;i++)for(t=0,j=0;j<ns;j++,t++)
	{
		if(rtk->ssat[sat[j]-1].refsatflag[i]==1)
		{
			t--;
			continue;
		}
		else
		{
			for(k=0;k<m;k++)
			{
				if((t+i*(ns-1))==k) D2[t*m+i*(ns-1)*m+k] = a2[i]+rtk->ssat[sat[j]-1].D2_Stoc[i];
				else                D2[t*m+i*(ns-1)*m+k] = a2[i];
			}
		}
	}

	//D =D1+D2
	for(i=0;i<m;i++)for(j=0;j<m;j++)
	{
		D[j+i*m]=D1[j+i*m]+D2[j+i*m];
	}
	//trace(3,"D= \n");tracemat(4,D,m,m,10,6);

	//计算权阵
	if (norm(D,m*m)>0 &&(!(info=matinv(D,m))))     /* P =D      */
	{
		for(i=0;i<m;i++)
			for(j=0;j<m;j++)
				P[j+i*m]=D[j+i*m];	
		//LSE
		matmul("TN",3,m,m,1.0,H,P,0.0,HP);         /* HP=H'*P   */
		matmul("NN",3,3,m,1.0,HP,H,0.0,N);         /* N =H'*P*H */
		//trace(3,"N= \n");tracemat(4,N,3,3,10,6);
		matmul("NN",3,1,m,1.0,HP,w,0.0,W);         /* W =H'*P*w */
		//trace(3,"W= \n");tracemat(4,W,1,3,8,4);
		if (norm(N,9)>0 && (!(info=matinv(N,3))))  /* Q =N^-1   */
		{
			//trace(3,"Q= \n");tracemat(4,N,3,3,10,6);
			matmul("NN",3,1,3,1.0,N,W,0.0,dx);     /* dx=(N^-1)*W*/
			//trace(3,"dx= \n");tracemat(4,dx,1,3,10,6);
			if(norm(dx,3)<3 && rtk->epoch_count<5){
				for(i=0;i<3;i++) rtk->x[i]+=dx[i];}
			//trace(4,"x(TD)="); tracemat(4,rtk->x,1,3,13,4);
			stat=1;
		}
	}

	//Free
	free(H); free(D); free(D1); free(D2); 
	free(HP);free(P); free(w); 
	return stat;
}


/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    char path[1024];
    int i,nu,nr,n=0;
    
    trace(3,"infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n",revs,iobsu,iobsr,isbs);
    
    if (0<=iobsu&&iobsu<obss.n) {
        //settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq))
		{
            aborts=1; 
			//showmsg("aborted"); 
			return -1;
        }
    }
    if (!revs) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        }
        nr=nextobsf(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=nu;
        
        /* update sbas corrections */
        while (isbs<sbss.n) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update lex corrections */
        while (ilex<lexs.n) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            }
            ilex++;
        }
        /* update rtcm corrections */
        if (*rtcm_file) {
            
            /* open or swap rtcm file */
            reppath(rtcm_file,path,obs[0].time,"","");
            
            if (strcmp(path,rtcm_path)) {
                strcpy(rtcm_path,path);
                
                if (fp_rtcm) fclose(fp_rtcm);
                fp_rtcm=fopen(path,"rb");
                if (fp_rtcm) {
                    rtcm.time=obs[0].time;
                    input_rtcm3f(&rtcm,fp_rtcm);
                    trace(2,"rtcm file open: %s\n",path);
                }
            }
            if (fp_rtcm) {
                while (timediff(rtcm.time,obs[0].time)<0.0) {
                    if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;
                }
                for (i=0;i<MAXSAT;i++) navs.ssr[i]=rtcm.ssr[i];
            }
        }
    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        }
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;
        
        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
        /* update lex corrections */
        while (ilex>=0) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            }
            ilex--;
        }
    }
    return n;
}
/* process positioning -------------------------------------------------------*/
/*//////////////////////////////////////////////////////////////////////////////
//函数：   static void procpos（）;
//功能：   定位计算函数
//参数：   FILE *fp         %   输出文件指针
//         prcopt_t *popt   %   解算参数
//         int mode         %   0：顺序或倒序；1：顺序+倒序
//XGAO, 2014/9
//////////////////////////////////////////////////////////////////////////////*/ 
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                    int mode)
{
    gtime_t time={0};
    sol_t sol={{0}};
	sol.sigma=1000;	//Modified by XGAO,2014/10/17
    rtk_t rtk;
    obsd_t obs[MAXOBS+20]; /* for rover and base */
    double rb[3]={0};
    int i,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6};

	//Modified by XGAO , 2014/10/13
	int count0=0,count1=0;        //计数器
	double sum_sol[3]={{0.0}}; //记录每固定解
	double sum_sol0[3]={0.0};  //固定解均值初值
	double sum_ratio=0.0;
    double pi=0.0;
    trace(3,"procpos : mode=%d\n",mode);
    
	//sopt->solstatic (0:all,1:single)
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
	//rtk初始化
    rtkinit(&rtk,popt);
    rtcm_path[0]='\0';
    
	//读取每个历元观测数据，保存在obs数组中
    while ((nobs=inputobs(obs,rtk.sol.stat,popt))>=0)
	{
        /* exclude satellites */
		//卫星剔除
        for (i=n=0;i<nobs;i++) 
		{
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }
		//观测数据不大于0，则进入下一历元
        if (n<=0) continue;
        
		//精密定位
		if (!rtkpos(&rtk,obs,n,&navs)) continue;


		//////加权定位
		//if(rtk.sol.stat==1 && rtk.sol.ratio>10)
		//{
		//	count0++;
		//	pi = SQRT(rtk.sol.ratio);
		//	sum_ratio += pi;
		//	for(i=0;i<3;i++) sum_sol[i]+=(pi*rtk.sol.rr[i]);
		//}

		//VTPV
		//if(rtk.sol.stat==1 && rtk.sol.ratio>10)
		//{
		//	count0++;
		//	pi=1.0/rtk.sol.sigma;
		//	sum_ratio += pi;
		//	for(i=0;i<3;i++) sum_sol[i]+=(pi*rtk.sol.rr[i]);
		//}

		//Modified by XGAO+CYJ
		//if(rtk.sol.stat==1)count0++;
		//else  count0=0;

		//if((rtk.tt*count0)>(10*60.0))//15min
		//{
		//	count1++;
		//	//记录固定解均值初值(由计数第一个历元赋值)
		//	if(count1==1) for(i=0;i<3;i++) sum_sol0[i]=rtk.sol.rr[i];  
		//	for(i=0;i<3;i++)
		//	{
		//		sum_sol[i]+=(rtk.sol.rr[i]-sum_sol0[i]);//(考虑数据长度超出浮点解有效数的影响)			
		//	}
		//}

		if (mode==0) { /* forward/backward */
			if (!solstatic) 
			{
				//结果输出
				if (rtk.sol.ratio>=3.0)
				{
					outsol(fp,&rtk.sol,rtk.rb,sopt);
				}
				
			}
			//Modified by XGAO
			else if ((time.time==0|| pri[rtk.sol.stat]<=pri[sol.stat]) && (rtk.sol.ratio >=sol.ratio))
			{
					sol=rtk.sol;
					for (i=0;i<3;i++) rb[i]=rtk.rb[i];
					if (time.time==0||timediff(rtk.sol.time,time)<0.0) time=rtk.sol.time;
			}
		}
		else if (!revs) { /* combined-forward */
			if (isolf>=nepoch) return;
			solf[isolf]=rtk.sol;
			for (i=0;i<3;i++) rbf[i+isolf*3]=rtk.rb[i];
			isolf++;
		}
		else { /* combined-backward */
			if (isolb>=nepoch) return;
            solb[isolb]=rtk.sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk.rb[i];
            isolb++;
        }
	}


	////取均值
	//if(count1>0)
	//{
	//	for(i=0;i<3;i++) sum_sol[i]/=count1;
	//}		
	//if(norm(sum_sol,3)>0.0)
	//	for(i=0;i<3;i++) sol.rr[i]=sum_sol[i]+sum_sol0[i];

	//if(count0>1)
	//	for(i=0;i<3;i++) sol.rr[i] = sum_sol[i]/sum_ratio;

    if (mode==0&&solstatic&&time.time!=0.0) 
	{
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }


	//释放rtk
    rtkfree(&rtk);
}

//Modified by XGAO,2014/7/18
/* process positioning -------------------------------------------------------*/
static void procpos_DDR(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
					int mode)
{
	gtime_t time={0};
	sol_t sol={{0}};
	rtk_t rtk;
	obsd_t obs[MAXOBS];
	double rb[3]={0};
	int i,j,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6};
	int satprn[MAXOBS]={0};          //存储卫星PRN
    int satprn1[MAXOBS];         //存储卫星PRN
	int satprn2[MAXOBS];         //存储卫星PRN
	int RefSat1[4][NFREQ]={0};   //参考卫星
	int RefSat2[4][NFREQ]={0};   //参考卫星
	int flag=0;                  //参考卫星标识符
	int ns,ns1=0,ns2=0;          //共视卫星数
	//trace(3,"procpos : mode=%d\n",mode);

	solstatic=sopt->solstatic&&
		(popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);

	rtkinit(&rtk,popt);
	rtcm_path[0]='\0';


	int EpochCount=0;
	//Slip初始化
	for(i=0;i<MAXSAT;i++)
		for(int f=0;f<NFREQ;f++)
			rtk.ssat[i].slip[f]&=0xFC;
	//RefSat初始化
	for(i=0;i<4;i++)
		for(int f=0;f<NFREQ;f++)
			rtk.Ref_Sat[i][f]=0;  //参考卫星

	while ((nobs=inputobs(obs,rtk.sol.stat,popt))>=0) 
	{		
		for (i=n=0;i<nobs;i++)
		{	
			/* exclude satellites */
			if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
				popt->exsats[obs[i].sat-1]!=1) 
			{obs[n++]      = obs[i];}
		}
		
		if (n<=0) continue;
		rtk.epoch_count++;

		//初始历元
		if(EpochCount==0)
		{
			if(rtkpos_DDR(&rtk,obs,n,&navs,1))
			{
				EpochCount++;
				ns1=rtk.ns;
                for(i=0;i<ns1;i++)
					satprn1[i]=rtk.satprn[i];
				for(i=0;i<4;i++)
					for(int f=0;f<NFREQ;f++)
						RefSat1[i][f]=rtk.Ref_Sat[i][f];  //参考卫星
				continue;
			}	
		}
		else
		{
			if(rtkpos_DDR(&rtk,obs,n,&navs,2))
			{
				//计数器累加
				if(EpochCount>0)
				{
					for(i=0;i<ns2;i++)
						satprn1[i]=satprn2[i];
					for(i=0;i<4;i++)
						for(int f=0;f<NFREQ;f++)
							RefSat1[i][f]=RefSat2[i][f];  //参考卫星

				}
				ns2=rtk.ns;
				for(i=0;i<ns2;i++)
					satprn2[i]=rtk.satprn[i];
				for(i=0;i<4;i++)
					for(int f=0;f<NFREQ;f++)
						RefSat2[i][f]=rtk.Ref_Sat[i][f];  //参考卫星
			}
		}

		if(EpochCount==0) continue;

		//判断参考卫星是否更替
		if(RefSat1[0][0]==RefSat2[0][0] && RefSat1[1][0]==RefSat2[1][0] && RefSat1[2][0]==RefSat2[2][0])
		{
			//构建共视卫星列表
            ns=sat_select(&rtk,satprn,satprn1,satprn2,ns1,ns2);

			//三差法周跳探测
			cycleslip_TripeD(&rtk,satprn,popt->nf,ns);

			////三差法定位计算
			if(!(posi_TripeD(&rtk,satprn,popt->nf,ns)))
				trace(3,"TD Positioning failed \n");

		}
		else   //参考卫星变更，重新进行三差法
		{
			EpochCount--;
		}

		//转回双差相对定位
		if (!rtkpos(&rtk,obs,n,&navs)) continue;


		if (mode==0) { /* forward/backward */
			if (!solstatic) 
			{
				outsol(fp,&rtk.sol,rtk.rb,sopt);
			}
			else if (time.time==0||pri[rtk.sol.stat]<=pri[sol.stat]) {
				sol=rtk.sol;
				for (i=0;i<3;i++) rb[i]=rtk.rb[i];
				if (time.time==0||timediff(rtk.sol.time,time)<0.0) {
					time=rtk.sol.time;
				}
			}
		}
		else if (!revs) { /* combined-forward */
			if (isolf>=nepoch) return;
			solf[isolf]=rtk.sol;
			for (i=0;i<3;i++) rbf[i+isolf*3]=rtk.rb[i];
			isolf++;
		}
		else { /* combined-backward */
			if (isolb>=nepoch) return;
			solb[isolb]=rtk.sol;
			for (i=0;i<3;i++) rbb[i+isolb*3]=rtk.rb[i];
			isolb++;
		}
	}
	if (mode==0&&solstatic&&time.time!=0.0) {
		sol.time=time;
		outsol(fp,&sol,rb,sopt);
	}
	rtkfree(&rtk);
}



/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb)
{
    double dr[3],var[3];
    int i;
    char tstr[32];
    
    trace(3,"valcomb :\n");
    
    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];
        var[i]=solf->qr[i]+solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */
        
        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0};
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};
    
    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);
    
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        
        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        else if (solf[i].stat<solb[j].stat) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (solf[i].stat>solb[j].stat) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);
            
            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {
                
                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j)) sols.stat=SOLQ_FLOAT;
            }
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];
            
            if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,sopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char **infile, int n, const prcopt_t *prcopt,
                        nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    seph_t seph0={0};
    int i;
    char *ext;
    
    trace(3,"readpreceph: n=%d\n",n);
    
    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    sbs->n =sbs->nmax =0;
    lex->n =lex->nmax =0;
    
    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readrnxc(infile[i],nav);
    }
    /* read sbas message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        sbsreadmsg(infile[i],prcopt->sbassatsel,sbs);
    }
    /* read lex message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        lexreadmsg(infile[i],0,lex);
    }
    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         //showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;
    
    /* set rtcm file and initialize rtcm struct */
    rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;
    
    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
            init_rtcm(&rtcm);
            break;
        }
    }
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    int i;
    
    trace(3,"freepreceph:\n");
    
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    free(lex->msgs); lex->msgs=NULL; lex->n =lex->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;
    
#ifdef EXTSTEC
    stec_free(nav);
#endif
    
    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;
    
    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);
    
    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;
    
    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;
        
        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n; 
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"no obs data\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"no nav data\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);
    
    /* delete duplicated ephemeris */
    uniqnav(nav);
    
    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            //settspan(ts,te);
        }
    }
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");
    
    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];
    
    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);
    
    for (i=0;i<3;i++) ra[i]=0.0;
    
    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {
        
        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */
        
        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) continue;
        
        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, char *name, double *r)
{
    FILE *fp;
    char buff[256],sname[256],*p,*q;
    double pos[3];
    
    trace(3,"getstapos: file=%s name=%s\n",file,name);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,'%'))) *p='\0';
        
        if (sscanf(buff,"%lf %lf %lf %s",pos,pos+1,pos+2,sname)<4) continue;
        
        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            pos[0]*=D2R;
            pos[1]*=D2R;
            pos2ecef(pos,r);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;
    
    trace(3,"antpos  : rcvno=%d\n",rcvno);
    
    if (postype==1) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            //showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==2) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
        if (!getstapos(posfile,name,rr)) {
            //showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
    }
    else if (postype==3) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) 
		{
            //showmsg("error : no position in rinex header");
            trace(1,"no position position in rinex header\n");
            return 0;
        }
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    const char *ext;
    
    trace(3,"openses :\n");
    
    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
        //showmsg("error : no sat ant pcv in %s",fopt->satantp);
        trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
        //showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    /* read dcb parameters */
    if (*fopt->dcb) {
        readdcb(fopt->dcb,nav);
    }
    /* read ionosphere data file */
    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I')) {
            readtec(fopt->iono,nav,0);
        }
#ifdef EXTSTEC
        else if (!strcmp(ext,".stec")||!strcmp(ext,".STEC")) {
            stec_read(fopt->iono,nav);
        }
#endif
    }
    /* open geoid data */
    if (sopt->geoid>0&&*fopt->geoid) {
        if (!opengeoid(sopt->geoid,fopt->geoid)) {
            //showmsg("error : no geoid data %s",fopt->geoid);
            trace(2,"no geoid data %s\n",fopt->geoid);
        }
    }
    /* read erp data */
    if (*fopt->eop) {
        if (!readerp(fopt->eop,&nav->erp)) {
            //showmsg("error : no erp data %s",fopt->eop);
            trace(2,"no erp data %s\n",fopt->eop);
        }
    }
    return 1;
}
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");
    
    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;
    
    /* close geoid data */
    closegeoid();
    
    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;
    
    /* close solution statistics and debug trace */
    rtkclosestat();
    traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv;
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];
    
    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(2,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    
    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt)
{
    FILE *fp=stdout;
    
    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);
    
    if (*outfile) {
        createdir(outfile);
        
        if (!(fp=fopen(outfile,"w"))) {
            //showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
    outheader(fp,infile,n,popt,sopt);
    
    if (*outfile) fclose(fp);
    
    return 1;
}

/* execute processing session ------------------------------------------------*/
/*/////////////////////////////////////////////////////////////////////////////
//函数：     static int execses（）
//功能：     事后数据处理
//参数：     gtime_t ts     %   解算开始时间
//           gtime_t te     %   解算结束时间
//           double  ti     %   采样间隔
//           prcopt_t *popt %   解算参数设置
//           solopt_t *sopt %   结果输出选项
//           filopt_t *fopt %   文件选项
//           int flag       %   1 
//           char **infile  %   输入文件路径
//           int *index     %   文件下标
//           int n          %   输入文件数
//           char *outfile  %   输出文件路径
return       O:(OK); 1(abort); 
// XGAO, 2014/9
/////////////////////////////////////////////////////////////////////////////*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                   const solopt_t *sopt, const filopt_t *fopt, int flag,
                   char **infile, const int *index, int n, char *outfile)
{
    FILE *fp;
    prcopt_t popt_=*popt;
    char tracefile[1024],statfile[1024];
    
    trace(3,"execses : n=%d outfile=%s\n",n,outfile);
    
    /* open debug trace */
	// 打开trace跟踪文件
	 
	//flag=1;  //Modified by CYJ 2014/10/16
 ////   if (flag&&sopt->trace>0) 
	//if(flag)
	//{
 //       if (*outfile) 
	//	{
 //           strcpy(tracefile,outfile);
 //           strcat(tracefile,".trace");
 //       }
 //       else 
	//	{
 //           strcpy(tracefile,fopt->trace);
 //       }
 //       traceclose();
 //       traceopen(tracefile); 
 // //      tracelevel(sopt->trace);
	//	tracelevel(7);
 //   }
	if (flag&&sopt->trace>0) 
		//if(flag)
	{
		if (*outfile) 
		{
			strcpy(tracefile,outfile);
			strcat(tracefile,".trace");
		}
		else 
		{
			strcpy(tracefile,fopt->trace);
		}
		traceclose();
		traceopen(tracefile); 
		tracelevel(sopt->trace);
		//tracelevel(4);
	}
	tracelevel(3);
    /* read obs and nav data */
	//读取观测数据与导航电文
    if (!readobsnav(ts,te,ti,infile,index,n,&popt_,&obss,&navs,stas)) return 0;
    
    /* set antenna paramters */
	//如果有精密星历或者SSR_COM数据，设置相应的天线改正参数
    if (popt_.sateph==EPHOPT_PREC||popt_.sateph==EPHOPT_SSRCOM) 
	{
        setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
               stas);
    }
    /* read ocean tide loading parameters */
	//读取海洋潮汐负荷参数
    if (popt_.mode>PMODE_SINGLE&&fopt->blq) {
        readotl(&popt_,fopt->blq,stas);
    }
    /* rover/reference fixed position */
	//Fixed 模式下进行基准站天线相位改正
    if (popt_.mode==PMODE_FIXED)
	{
        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos))
		{
            freeobsnav(&obss,&navs);
            return 0;
        }
    }
	//非PPP模式下进行流动站天线相位中心改正
    else if (PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC) 
	{
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos))
		{
            freeobsnav(&obss,&navs);
            return 0;
        }
    }

    /* open solution statistics */
	// 打开定位结果统计信息文件
    if (flag&&sopt->sstat>0) {
        strcpy(statfile,outfile);
        strcat(statfile,".stat");
        rtkclosestat();
        rtkopenstat(statfile,sopt->sstat);
    }

    /* write header to output file */
	// 写入结果文件文件头（解算参数等）
    if (flag&&!outhead(outfile,infile,n,&popt_,sopt)) {
        freeobsnav(&obss,&navs);
        return 0;
    }

    iobsu=iobsr=isbs=ilex=revs=aborts=0;
    
	// 如果是单点定位或者顺序解算
    if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) 
	{
        if ((fp=openfile(outfile)))
		{
			//定位计算
            procpos(fp,&popt_,sopt,0); /* forward */
            fclose(fp);
        }
    }
	//倒序解算
    else if (popt_.soltype==1) {
        if ((fp=openfile(outfile))) {
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
			//定位计算
            procpos(fp,&popt_,sopt,0); /* backward */
            fclose(fp);
        }
    }
	//顺序*倒序结合
    else { /* combined */
        solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        rbf=(double *)malloc(sizeof(double)*nepoch*3);
        rbb=(double *)malloc(sizeof(double)*nepoch*3);
        
        if (solf&&solb) {
            isolf=isolb=0;
			//定位计算
            procpos(NULL,&popt_,sopt,1); /* forward */
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
			//定位计算
            procpos(NULL,&popt_,sopt,1); /* backward */
            
            /* combine forward/backward solutions */
            if (!aborts&&(fp=openfile(outfile))) {
                combres(fp,&popt_,sopt);
                fclose(fp);
            }
        }
   //     else 
			//showmsg("error : memory allocation");
        free(solf);
        free(solb);
        free(rbf);
        free(rbb);
    }
    /* free obs and nav data */
    freeobsnav(&obss,&navs);
    
    return aborts?1:0;
}
/* execute processing session for each rover ---------------------------------*/
//流动站
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*rov_,*p,*q,s[64]="";
    
    trace(3,"execses_r: n=%d outfile=%s\n",n,outfile);
    
    for (i=0;i<n;i++) if (strstr(infile[i],"%r")) break;
    
    if (i<n) { /* include rover keywords */
        if (!(rov_=(char *)malloc(strlen(rov)+1))) return 0;
        strcpy(rov_,rov);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(rov_); for (;i>=0;i--) free(ifile[i]);
                return 0;
            }
        }
        for (p=rov_;;p=q+1) { /* for each rover */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_rov,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,p,"");
                reppath(outfile,ofile,t0,p,"");
                
                /* execute processing session */
                stat=execses(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile);
            }
            if (stat==1||!q) break;
        }
        free(rov_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        /* execute processing session *///
		//事后处理代码
        stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    }
    return stat;
}
/* execute processing session for each base station --------------------------*/
//基准站
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov, const char *base)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*base_,*p,*q,s[64];
    
    trace(3,"execses_b: n=%d outfile=%s\n",n,outfile);
    
    /* read prec ephemeris and sbas data */
	//读取精密星历和SBAS数据
    readpreceph(infile,n,popt,&navs,&sbss,&lexs);
    
    for (i=0;i<n;i++) 
		if (strstr(infile[i],"%b")) 
			break;
    
    if (i<n) 
	{ /* include base station keywords */
        if (!(base_=(char *)malloc(strlen(base)+1))) {
            freepreceph(&navs,&sbss,&lexs);
            return 0;
        }
        strcpy(base_,base);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(base_); 
				for (;i>=0;i--) 
					free(ifile[i]);
                freepreceph(&navs,&sbss,&lexs);
                return 0;
            }
        }
        for (p=base_;;p=q+1) { /* for each base station */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_base,p);
                if (ts.time) time2str(ts,s,0); 
				else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,"",p);
                reppath(outfile,ofile,t0,"",p);
                
                stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile,rov);
            }
            if (stat==1||!q) break;
        }
        free(base_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
		//流动站
        stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);
    }
    /* free prec ephemeris and sbas data */
    freepreceph(&navs,&sbss,&lexs);
    
    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise 
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first 
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
/*//////////////////////////////////////////////////////////////////////////////
**函数名：  extern int postpos();
**功能：    GNSS数据后处理算法
**输入参数; gtime_t ts;           % 设定数据解算开始时间
**          gtime_t te;           % 设定数据解算结束时间
**          double ti;            % 设定数据解算时间间隔（秒）
**          double tu;            % 设定数据解算单元块时间（秒）
**          const prcopt_t *popt; % 数据解算参数
**          const solopt_t *sopt; % 结果输出参数
**          const filopt_t *fopt; % 输入输出文件参数（天线改正、电离层、DCB等）
**          char **infile;        % 存储输入文件路径
**          int n;                % 输入文件数量
**          char *outfile;        % 存储结果输出文件路径
**          const char *rov;      % 流动站ID
**          const char *base;     % 基准站ID
**输出参数：定位结果：            % 0（成功）；>0(失败)；1（暂停）
**XGAO, 2014/9/22
/////////////////////////////////////////////////////////////////////////////*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile,
                   const char *rov, const char *base)
{
    gtime_t tts,tte,ttte;
    double tunit,tss;
    int i,j,k,nf,stat=0,week,flag=1,index[MAXINFILE]={0};
    char *ifile[MAXINFILE],ofile[1024],*ext;
    
    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);
    
    /* open processing session */
	//Step 1: 读取相关文件到指定的全局变量中
	//navs:  % 导航电文信息；
	//pcvss: % 卫星天线相位中心改正
	//pcrsr: % 接收机天线相位中心改正
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;
    
	//Step 2: 通过开始、结束时间以及解算时间块进行分段数据解算 
    if (ts.time!=0&&te.time!=0&&tu>=0.0) 
	{
		//如果设定的开始时间大于结束时间，释放全局变量，返回0；
        if (timediff(te,ts)<0.0) 
		{
            //showmsg("error : no period");
            closeses(&navs,&pcvss,&pcvsr);
            return 0;
        }

		//申请保存输入文件路径的空间，不成功，则返回-1；
        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024)))
			{
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }

		//如果单位时间块大于设定阈值（100天），则强制设定为该阈值
        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
        //settspan(ts,te);

		//子单元时间块tunit（当tu小于1天时则为tu,否则设定为1天）
        tunit=tu<86400.0?tu:86400.0;
		//tss：单位时间块开始时间（GPS周内秒）
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);
        

		//解算每个时间块
        for (i=0;;i++) 
		{ /* for each periods */
			//计算每个子单元时间块开始时间与结束时间
            tts=gpst2time(week,tss+i*tu);
            tte=timeadd(tts,tu-DTTOL);
            if (timediff(tts,te)>0.0) break;
            if (timediff(tts,ts)<0.0) tts=ts;
            if (timediff(tte,te)>0.0) tte=te;
            
            strcpy(proc_rov ,"");
            strcpy(proc_base,"");

            if (checkbrk("reading    : %s",time_str(tts,0))) 
			{
                stat=1;
                break;
            }
            for (j=k=nf=0;j<n;j++) 
			{
                
				//查找文件扩展名
                ext=strrchr(infile[j],'.');
                
				//是否为RTCM2格式
                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) 
				{
                    strcpy(ifile[nf++],infile[j]);
                }

				//文件非RTCM格式
                else 
				{
                    /* include next day precise ephemeris or rinex brdc nav */
                    ttte=tte;
					//如果有精密星历,星历文件时间扩展1小时
                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                        ttte=timeadd(ttte,3600.0);
                    }
					//如果为广播星历，星历文件时间扩展2小时
                    else if (strstr(infile[j],"brdc")) {
                        ttte=timeadd(ttte,7200.0);
                    }

					//根据单位时间块更新输入文件文件名，文件数自动增加变更文件数
                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
                }
                while (k<nf) index[k++]=j;
                
				//如果文件数量超限，则退出文件名更新循环
                if (nf>=MAXINFILE) {
                    trace(2,"too many input files. trancated\n");
                    break;
                }
            }

			//更新输出文件
            if (!reppath(outfile,ofile,tts,"","")&&i>0) flag=0;
            
            /* execute processing session */
			//事后定位计算，解算成功返回0，失败返回1；
            stat=execses_b(tts,tte,ti,popt,sopt,fopt,flag,ifile,index,nf,ofile,
                           rov,base);
            
            if (stat==1) break;
        }
		//解算完毕，释放动态申请的输入文件空间
        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }

	//Step 3: 开始时间进行了设定
    else if (ts.time!=0)
	{
        for (i=0;i<n&&i<MAXINFILE;i++) 
		{
            if (!(ifile[i]=(char *)malloc(1024))) 
			{
                for (;i>=0;i--) 
					free(ifile[i]);
				return -1;
            }

			//根据设定的开始时间更新输入文件名及其信息
            reppath(infile[i],ifile[i],ts,"","");
            index[i]=i;
        }
        reppath(outfile,ofile,ts,"","");
        
        /* execute processing session */
		//定位计算
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,ifile,index,n,ofile,rov,
                       base);
        
        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
    }

	//Step 4: 如果开始时间与结束时间均未人为设定，即直接解算读入文件
    else {
        for (i=0;i<n;i++) index[i]=i;
        
        /* execute processing session */
		//定位计算
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,
                       base);
    }
    /* close processing session */
	//释放全局变量
    closeses(&navs,&pcvss,&pcvsr);
    
    return stat;
}
