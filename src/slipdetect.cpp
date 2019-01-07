/* 历元间差分周跳探测与修复模板（多模）---------------------------------------------------*/
#include"../stdafx.h"
#include <stdarg.h>
#include "rtklib.h"
#include "Adjust\modern.h"
#include "Adjust\public.h"
#include "Probability\probability.h"

/* constants/macros ----------------------------------------------------------*/

#define MAXX        210                    //最多的参数个数（法方程迭代计算）

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define QFH(x)      ((x)<=0.0?(int)(-1.0):(int)(1.0))   //取符号函数，如果x<0.0，返回-1，否则+1
#define SSWR(x)     ((x)-(int)(x)>=0.5?(int)(x)+1:(int)(x))  //正数四舍五入

#define MAXSLIP     100      /* 周跳数预留数，即发生周跳时新增模糊度参数最大个数(20131111) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NC(opt)      0
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*(opt)->nf)

#define NB2(opt)    MAXSLIP    /* 预留参数：当周跳发生时，新增的模糊度参数个数 */


#define NR(opt)     (NP(opt)+NC(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt)+NB2(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */
#define IB2(index,opt) (NR(opt)+NB(opt)+index)   /* phase bias(slip satellite) (index:表示在slipflag中的下标(从0开始)) */


/* global variables ----------------------------------------------------------*/
/* 定义静态变量slipflag[]数组用于记录发生周跳的模糊度信息 */
static int nslip=0;                  /* slipflag数组非零数，即周跳发生数 */
static int slips[MAXX];           /* 存储周跳数 */
static int now_slipflg[MAXX];     /* 当前历元周跳发生情况标记 */ 
static int slips_doppler[MAXX];   /* 多普勒积分法探测的周跳近似值*/

extern void get_now_slipflg(int *slipflg)
{
	maticpy(slipflg,now_slipflg,MAXSLIP,1);
}

/* double-differenced measurement error covariance ---------------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;
    
    trace(3,"ddcov   : n=%d\n",n);
    
    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,6);
}

/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds) ----------------------*/
static int test_sys(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==0;
        case SYS_QZS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
    }
    return 0;
}

/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}
/* single-differenced measurement error variance -----------------------------*/
static double varerr(int sys, double el, double bl, double dt, int f,
                     const prcopt_t *opt)
{
    double a=opt->err[1],b=opt->err[2]/sin(el),c=opt->err[3]*bl/1E4;
    double d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    if (f>=opt->nf) fact=opt->eratio[f-opt->nf];
    if (fact<=0.0)  fact=opt->eratio[0];
    
    /* system factor affects to both of code and phase error (v.2.3.0) */
#if 1
    fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
#else
    if (f>=opt->nf) {
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    }
#endif
    return 2.0*(a*a*fact*fact+b*b*fact*fact+c*c)+d*d;
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
	/* if no phase observable, psudorange is also unusable */
	return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
		(f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* 单差观测量周跳修复 */
extern void slip_repaire(rtk_t *rtk, const nav_t *nav,const int *sat, const int ns,
						  const int *iu, double *y)
{
	prcopt_t *opt=&rtk->opt;
	int f,j,slip,k;
	double lamj;
	trace(2,"单差观测量周跳修复\n");
	for(f=0;f<opt->nf;f++)	
		for(j=0;j<ns;j++) {
			slip=slips[IB(sat[j],f,&rtk->opt)];
			if(slip==0) continue;
			k=f+iu[j]*2*opt->nf;
			lamj=nav->lam[sat[j]-1][f];
			y[f+iu[j]*2*opt->nf]-=slip*lamj;
			trace(2,"修复周跳：L%2d,sat:%3d,slip:%10d\n",f+1,sat[j],slip);
		}
}
/* detect cycle slip by doppler and phase difference -------------------------*/
static void detslp_dop_one(rtk_t *rtk, const obsd_t *obs, int i, int rcv,
	const nav_t *nav)
{
	/* detection with doppler disabled because of clock-jump issue (v.2.3.0) */
	int f,sat=obs[i].sat;
	double tt,dph,dpt,lam,thres,slp;

	/*FILE *fp;
	char *dop_file="E:\\论文4：单频周跳探测与修复\\detectAndRepairSlpByDt\\RESULT\\detslp_dop.data";

	trace(3,"detslp_dop: i=%d rcv=%d\n",i,rcv);

	fp=fopen(dop_file,"a+");*/
	/*if(fp==NULL) {trace(1,"%s文件打开失败！",dop_file);getchar();return;}*/

	for (f=0;f<rtk->opt.nf;f++) 
	{
		if (obs[i].L[f]==0.0||obs[i].D[f]==0.0||rtk->ssat[sat-1].ph[rcv-1][f]==0.0||
			rtk->ssat[sat-1].pd[rcv-1][f]==0.0) {
				continue;
		}
		if (fabs(tt=timediff(obs[i].time,rtk->ssat[sat-1].pt[rcv-1][f]))<DTTOL) continue;
		if ((lam=nav->lam[sat-1][f])<=0.0) continue;

		/* cycle slip threshold (cycle) */
		//        thres=MAXACC*tt*tt/2.0/lam+rtk->opt.err[4]*fabs(tt)*4.0;

		thres=0.03*fabs(tt)/lam*3;

		/* phase difference and doppler x time (cycle) */
		dph=obs[i].L[f]-rtk->ssat[sat-1].ph[rcv-1][f];
		dpt=-1*(obs[i].D[f]+rtk->ssat[sat-1].pd[rcv-1][f])*tt/2;

		//		trace(2,"obs[%d].D[%d]:%10.3lf\n",i,f,obs[i].D[f]);
		//		trace(2,"dph:%10.3lf,dpt:%10.3lf\n",dph,dpt);
		//		trace(2,"dph-dpt=%10.3lf,thres=%10.3lf\n",dph-dpt,thres);
		trace(2,"slip detected (sat=%2d rcv=%d L%d=%.3f %.3f dph-dpt=%10.3lf thres=%.3f)\n",
			sat,rcv,f+1,dph,dpt,dph-dpt,thres);

		//if(f==0) {  //接收机，时间，频率，卫星，d_phase-d_doppler,周跳判定阀值
		//	fprintf(fp,"%d %s %d  %3d %10.3lf %10.3lf\n",
		//		rcv,time_str(rtk->sol.time,3),f+1,sat,dph-dpt,thres );

		if (fabs(dph-dpt)<=thres) continue;

		slp=dph-dpt;

		now_slipflg[IB(sat,f,&rtk->opt)]=1;    //周跳记录

		rtk->ssat[sat-1].slip[f]|=1;

		if(rcv==2) slp*=-1;                       //以流动站为基准，后面双差为：流动站-基站 
		//多普勒探测到的可能周跳近似值
		if(slips_doppler[IB(sat,f,&rtk->opt)]<=QFH(slp)*SSWR(slp))
			slips_doppler[IB(sat,f,&rtk->opt)]=QFH(slp)*SSWR(slp);  

		errmsg(rtk,"detect slip by doppler (sat=%2d rcv=%d L%d=%.3f %.3f thres=%.3f)\n",
			sat,rcv,f+1,dph,dpt,thres);

	}/*fclose(fp);*/
}


/* 多普勒探测周跳：采用星间单差20140921-------------------------------------------------- */
static void detslp_doppler_dif_sat(rtk_t *rtk,const obsd_t *obs, const int *sat,
                   const int *iu, int ns, int rcv, const nav_t *nav,const double *azel) 
{
    int i,j,m,f,j1,j2,sysi,sysj;
	double lam,slp,dphi,dphj,dpti,dptj,thres,tt;

	int sati,nslp=0;

	/*FILE *fp;
	char *dop_file="E:\\Fixed\\DD模型\\detslp_dop.data";

	trace(3,"detslp_doppler_dif_sat:rcv=%d\n",rcv);
   
	fp=fopen(dop_file,"a+");
	if(fp==NULL) {trace(1,"%s文件打开失败！",dop_file);getchar();return;}*/
        
    for (i=0;i<ns;i++) 
	{
        
        /* detect cycle slip by LLI */
        for (f=0;f<rtk->opt.nf;f++) rtk->ssat[sat[i]-1].slip[f]&=0xFC;
		
	}

	for (m=0;m<2;m++) /* for each system (0:gps/qzss/sbas,1:glonass) */
		
	for (f=0;f<rtk->opt.nf;f++) {  //载波相位定位
		
		/* search reference satellite with highest elevation */
		for (i=-1,j=0;j<ns;j++) {
			sysi=rtk->ssat[sat[j]-1].sys;
			if ((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
			j1=iu[j];
			if (obs[j1].L[f]==0.0||obs[j1].D[f]==0.0||rtk->ssat[sat[j]-1].ph[rcv-1][f]==0.0||
				rtk->ssat[sat[j]-1].pd[rcv-1][f]==0.0) {
				continue;
			}
			if (fabs(tt=timediff(obs[j1].time,rtk->ssat[sat[j]-1].pt[rcv-1][f]))<DTTOL) continue;
			if ((lam=nav->lam[sat[j]-1][f])<=0.0) continue;
			if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
		}
		
		if (i<0) continue;
		
		trace(2,"参考卫星编号为：%d\n",sat[i]);
		
		sysi=rtk->ssat[sat[i]-1].sys;
		
		/* make double difference */
		for (j=0;j<ns;j++) {
			if (sat[i]==sat[j]) continue;
			sysi=rtk->ssat[sat[i]-1].sys;
			sysj=rtk->ssat[sat[j]-1].sys;
			if ((m==0&&sysj==SYS_GLO)||(m==1&&sysj!=SYS_GLO)) continue;
			j2=iu[j];
			if (obs[j2].L[f]==0.0||obs[j2].D[f]==0.0||rtk->ssat[sat[j]-1].ph[rcv-1][f]==0.0||
				rtk->ssat[sat[j]-1].pd[rcv-1][f]==0.0) {
				continue;
			}
			if (fabs(tt=timediff(obs[j2].time,rtk->ssat[sat[j]-1].pt[rcv-1][f]))<DTTOL) continue;
			if ((lam=nav->lam[sat[j]-1][f])<=0.0) continue;			
	
			/* phase difference and doppler x time (cycle) */			
			dphi=obs[iu[i]].L[f]-rtk->ssat[sat[i]-1].ph[rcv-1][f];
			dphj=obs[iu[j]].L[f]-rtk->ssat[sat[j]-1].ph[rcv-1][f];
			dpti=-1*(obs[iu[i]].D[f]+rtk->ssat[sat[i]-1].pd[rcv-1][f])*tt/2;
			dptj=-1*(obs[iu[j]].D[f]+rtk->ssat[sat[j]-1].pd[rcv-1][f])*tt/2;
			slp=(dphi-dphj)-(dpti-dptj);
			
			thres=0.03*fabs(tt)/lam*3;
			
			trace(2,"sat=%2d-%2d rcv=%d L%d=%10.3f %10.3f %10.3f %10.3f slp=%10.3lf thres=%10.3f\n",
				sat[i],sat[j],rcv,f+1,dphi,dpti,dphj,dptj,slp,thres);
/*			trace(2,"%10.3lf  %10.3lf  %10.3lf  %10.3lf\n",obs[iu[i]].L[f],
				rtk->ssat[sat[i]-1].ph[rcv-1][f],obs[iu[j]].L[f],rtk->ssat[sat[j]-1].ph[rcv-1][f]);
			trace(2,"%10.3lf  %10.3lf  %10.3lf  %10.3lf\n",obs[iu[i]].D[f],
				rtk->ssat[sat[i]-1].pd[rcv-1][f],obs[iu[j]].D[f],rtk->ssat[sat[j]-1].pd[rcv-1][f]);*/
		
			//if(f==0&&sat[j]==4&&rcv==1) {  //接收机，时间，频率，卫星，d_phase-d_doppler,周跳判定阀值
			//	fprintf(fp,"%d %s %d %d %3d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n",
			//		rcv,time_str(rtk->sol.time,3),f+1,sat[i],sat[j],dphi,dpti,dphj,dptj,slp,thres );
			//}
			
			if (fabs(slp)<=thres) continue;

			nslp++;  //记录周跳数

			now_slipflg[IB(sat[j],f,&rtk->opt)]=1;    //周跳记录
			
			rtk->ssat[sat[j]-1].slip[f]|=1;

			slp=-1*slp;    //  星间差分为   参考星-流动星，周跳符号要改正

			if(rcv==2) slp*=-1;                       //以流动站为基准，后面双差为：流动站-基站 
			//多普勒探测到的可能周跳近似值
			if(slips_doppler[IB(sat[j],f,&rtk->opt)]<=QFH(slp)*SSWR(slp)) //流动站和基站周跳标记取最大的
				slips_doppler[IB(sat[j],f,&rtk->opt)]=QFH(slp)*SSWR(slp);  
			
			errmsg(rtk,"slip detected (sat=%2d-%2d rcv=%d L%d=%10.3f %10.3f %10.3f %10.3f slp=%10.3lf thres=%.3f)\n",
				sat[i],sat[j],rcv,f+1,dphi,dphj,dpti,dptj,slp,thres);
		}
	}
	/*fclose(fp);*/
}
   

/* detect cycle slip by doppler------------------------------------------*/
extern void detslp_doppler(rtk_t *rtk, const obsd_t *obs, const int *sat,
                   const int *iu, const int *ir, int ns, const nav_t *nav,const double *azel)
{
    int i,f;
	
    trace(3,"detslp_doppler:\n");

	for(i=0;i<MAXX;i++) now_slipflg[i]=0;   //当前历元周跳标记数组清0
        
    for (i=0;i<ns;i++) {
        
        /* detect cycle slip by LLI */
        for (f=0;f<rtk->opt.nf;f++) {
			rtk->ssat[sat[i]-1].slip[f]&=0xFC;
			
			slips_doppler[IB(sat[i],f,&rtk->opt)]=0;
		}
		
		detslp_dop_one(rtk,obs,iu[i], 1, nav);
		detslp_dop_one(rtk,obs,ir[i], 2, nav);            
	}
//	detslp_doppler_dif_sat(rtk, obs,sat,iu,ns,1,nav,azel); 
//	detslp_doppler_dif_sat(rtk, obs,sat,ir,ns,2,nav,azel); 
}

/* 历元间差分观测方程方差阵 ---------------------------*/
static void dtcov(int nv, const double *Ri, const double *Rj, double *R)
{
    int i;
    
    trace(3,"dtcov   : nv=%d\n",nv);
    
    for (i=0;i<nv*nv;i++) R[i]=0.0;
	
    for (i=0;i<nv;i++) 	R[i+i*nv]=Ri[i]+Rj[i];
	
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,6);
}


static int ztss_sd(const int nl, const int nx, const int nc,const double *A,double *Lv,
				 const double *P,double *cycle2,int *W,int *is, double *s, double *vv)
{
	int i,j,index[4]={0};
	double *atp,*apa,*apl,*L,*x,*v,*vtp;
	double ratio,vpv[4]={10000,10000,10000,10000},vtv[4]={10000,10000,10000,10000},vpv_buf,vv_buf;
	double *slip;
	int num,number,number1,number2;

	trace(3,"周跳修复搜索：\n");

	atp=mat(nx,nl);apa=mat(nx,nx);apl=mat(nx,1);
	
	matmul("NN",nx,nl,nl,1.0,  A,P,0.0,atp);
	matmul("NT",nx,nx,nl,1.0,atp,A,0.0,apa);	
	if(matinv(apa,nx)){
		trace(2,"法方程矩阵求逆失败\n");
		free(atp);free(apa);free(apl);
		return 0;		
	}	
	x=mat(nx,1);v=mat(nl,1);vtp=mat(1,nl);slip=mat(nl-1,2);
	for(i=0;i<nc;i++){
		for(j=0,number=0;j<nl;j++){
				if(cycle2[j+i*nl]==0.0) number++;
		}
//		if(number<4) continue;    //未发生周跳的卫星数

		L=Lv+i*nl;
		//求解法方程
		matmul("NN",nx, 1,nl,1.0,atp,L,0.0,apl);
		matmul("NN",nx, 1,nx, 1.0,apa,apl,0.0,x);
		
		//计算v'pv
		matcpy(v,L,nl,1);
		matmul("TN",nl, 1,nx,-1.0,  A,  x,1.0,v);
		matmul("NN",1,nl,nl,1.0,  v,P,0.0,vtp);
		matmul("NN",1, 1,nl,1.0,vtp,v,0.0,&vpv_buf);
		//计算vv
		matmul("NN",1, 1,nl,1.0,  v,v,0.0,&vv_buf);

		//最优组合和次优组合(顶多只有三组单差模糊度周跳计算的vpv相同（最优），所以保留4组)
		if(vpv_buf<vpv[0]) {
			vpv[3]=vpv[2];
			vpv[2]=vpv[1];
			vpv[1]=vpv[0];
			vpv[0]=vpv_buf;
			index[3]=index[2];
			index[2]=index[1];
			index[1]=index[0];
			index[0]=i;
			vtv[3]=vtv[2];
			vtv[2]=vtv[1];
			vtv[1]=vtv[0];
			vtv[0]=vv_buf;
		}
		else if(vpv_buf<vpv[1]) {
			vpv[3]=vpv[2];
			vpv[2]=vpv[1];
			vpv[1]=vpv_buf;
			index[3]=index[2];
			index[2]=index[1];
			index[1]=i;
			vtv[3]=vtv[2];
			vtv[2]=vtv[1];
			vtv[1]=vv_buf;
		}
		else if(vpv_buf<vpv[2]){
			vpv[3]=vpv[2] ;
			vpv[2]=vpv_buf;
			index[3]=index[2];
			index[2]=i;
			vtv[3]=vtv[2];
			vtv[2]=vv_buf;
		}		
	}
	trace(2,"vpv对应最小的前4组周跳：\n");
	tracemat(2,cycle2+index[0]*nl,1,nl,10,3);
	tracemat(2,cycle2+index[1]*nl,1,nl,10,3);
	tracemat(2,cycle2+index[2]*nl,1,nl,10,3);
	tracemat(2,cycle2+index[3]*nl,1,nl,10,3);
	trace(2,"vpv依次为：\n");
	tracemat(2,vpv,1,4,10,3);
	num=0;
	while(num<2){  //最多比较两次
		if(fabs(vpv[1]-vpv[0])/vpv[0]<0.001){ //检验双差模糊度是否相同
			trace(2,"vpv[0]:%.3lf,vpv[1]:%.3lf\n",vpv[0],vpv[1]);
			tracemat(2,cycle2+index[0]*nl,1,nl,10,3);tracemat(2,cycle2+index[1]*nl,1,nl,10,3);
			for(j=0;j<nl-1;j++) slip[j]     =cycle2[j+1+index[0]*nl]-cycle2[0+index[0]*nl]; 
			for(j=0;j<nl-1;j++) slip[j+nl-1]=cycle2[j+1+index[1]*nl]-cycle2[0+index[1]*nl];
			trace(2,"slip:\n");tracemat(2,slip,1,nl-1,10,3);tracemat(2,slip+nl-1,1,nl-1,10,3);
			for(j=0;j<nl-1;j++) slip[j]-=slip[j+nl-1];  //做差比较
			trace(2,"slip:\n");tracemat(2,slip,1,nl-1,10,3);
			if(norm(slip,nl-1)<0.000001) { //如果双差模糊度相同，则替换掉次优模糊度备选组
				for(i=0,number1=0,number2=0;i<nl;i++) {
					if((cycle2[i+index[0]*nl]==0.0)&&(W[i]==1)) number1++;
					if((cycle2[i+index[1]*nl]==0.0)&&(W[i]==1)) number2++;
				}
				if(number1>=number2) {  //没有发生周跳的卫星数应该大于3				
					vpv[1]=vpv[2];
					vpv[2]=vpv[3];
					index[1]=index[2];
					index[2]=index[3];
					vtv[1]=vtv[2];
					vtv[2]=vtv[3];
				}
				else {
					vpv[0]=vpv[1];
					vpv[1]=vpv[2];
					vpv[2]=vpv[3];
					index[0]=index[1];
					index[1]=index[2];
					index[2]=index[3];
					vtv[0]=vtv[1];
					vtv[1]=vtv[2];
					vtv[2]=vtv[3];					
				}
			}

		}
		num++;	
	}

	matcpy(s,vpv,2,1);matcpy(vv,vtv,2,1);
	for(i=0;i<2;i++) {is[i]=index[i];}

	ratio=s[1]/s[0];

	trace(2,"模糊度备选组合搜索ratio=%6.3lf\n",ratio);
	trace(2,"最优组合与次优组合的下标依次为：%d,%d\n",index[0],index[1]);

	free(atp);free(apa);free(apl);free(x);
	free(v);free(vtp);free(slip);
	
	return 1;
}

static int ztss_dd(const int nl, const int nx, const int nc,const double *A,double *Lv,
				 const double *P,int *is, double *s, double *vv)
{
	int i,index[2]={0};
	double *atp,*apa,*apl,*L,*x,*v,*vtp;
	double ratio,vpv[2]={10000},vtv[2]={10000},vpv_buf,vv_buf;
	double *slip;
	
	trace(3,"周跳修复搜索：\n");

	atp=mat(nx,nl);apa=mat(nx,nx);apl=mat(nx,1);slip=mat(nl-1,2);
	
	matmul("NN",nx,nl,nl,1.0,  A,P,0.0,atp);
	matmul("NT",nx,nx,nl,1.0,atp,A,0.0,apa);	
	if(matinv(apa,nx)){
		trace(2,"法方程矩阵求逆失败\n");
		free(atp);free(apa);free(apl);
		return 0;		
	}	
	x=mat(nx,1);v=mat(nl,1);vtp=mat(1,nl);
	for(i=0;i<nc;i++){
		L=Lv+i*nl;
		//求解法方程
		matmul("NN",nx, 1,nl,1.0,atp,L,0.0,apl);
		matmul("NN",nx, 1,nx, 1.0,apa,apl,0.0,x);
		
		//计算v'pv
		matcpy(v,L,nl,1);
		matmul("TN",nl, 1,nx,-1.0,  A,  x,1.0,v);
		matmul("NN",1,nl,nl,1.0,  v,P,0.0,vtp);
		matmul("NN",1, 1,nl,1.0,vtp,v,0.0,&vpv_buf);
		//计算vv
		matmul("NN",1, 1,nl,1.0,  v,v,0.0,&vv_buf);
		//最优组合和次优组合
		if(vpv_buf<vpv[0]) {
			vpv[1]=vpv[0];
			vpv[0]=vpv_buf;
			index[1]=index[0];
			index[0]=i;
			vtv[1]=vtv[0];
			vtv[0]=vv_buf;
		}
		else if(vpv_buf<vpv[1]) {
			vpv[1]=vpv_buf;
			index[1]=i;
			vtv[1]=vv_buf;
		}
	}

	matcpy(s,vpv,2,1);matcpy(vv,vtv,2,1);
	for(i=0;i<2;i++) {is[i]=index[i];}

	ratio=s[1]/s[0];

	trace(2,"模糊度备选组合搜索ratio=%6.3lf\n",ratio);
	trace(2,"最优组合与次优组合的下标依次为：%d,%d\n",index[0],index[1]);

	free(atp);free(apa);free(apl);free(x);
	free(v);free(vtp);
	
	return 1;
}
//周跳搜索
extern void seach_slip(int nv,int nx,int *W,int *vflg,
	double *H,double *R,double *v,double *V,const nav_t *nav,double *ratio,
	double *rms_v,double *vpv,int *slp,double *w_ratio)
{
	int k,i,j,m,slip;
	int *index,sati,type,freq;
	char *stype;
	double vv,vtv[2],S[2],*cycle1,*cycle2,*Lv,*cycle3,ep=0.02,*E;
	double lami;
	int nc=0,t,is[2],flg=0;

	double *BP,*N,*NB,*BNB,*Qv,*Rv,*Qc,*dc,*dcQc,*Ql;
	double covd=0.0,Wd=0.0;

	trace(2,"seach_slip:\n");

	flg=0;k=0;nc=1;                     //记录周跳个数，当周跳个数为0时，结束解算
	cycle1=zeros(nv,1);                 //记录周跳发生情况（4舍5入值）
	index=imat(nv,1);
	for(i=0;i<nv;i++) {                 /* 记录周跳发生的状况 */
		if(W[i]==0) {
			index[k]=i;               //发生周跳的卫星在参数中的位置
			k++;                      //周跳的个数
			nc*=3;                    //发生周跳，搜索上下一周，共3个
			sati=(vflg[i]>> 8)&0xFF;  //卫星
			type=(vflg[i]>> 4)&0xF;   //观测量类型
			freq=vflg[i]&0xF;         //频率
			stype=type==0?"L":(type==1?"L":"C");
			lami=nav->lam[sati-1][freq];
			slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

			trace(1,"发生周跳卫星sat is :%d %s%d  V:%6.3lf  Slips:%d\n",
				sati,stype,freq+1,V[i],slip);

			cycle1[i]=slip*lami;    //转换为以米为单位
			vv=V[i];
			vv-=slip*lami;
			if(fabs(vv)>0.05)  {
				flg++;
				trace(1,"残差超限，周跳修复失败：周跳修复后残差V is :%6.3lf\n",vv);
				continue;
			}
		}
	}
	/*	if(flg>0){
	k=nv;
	for(i=0;i<nv;i++) {					
	index[i]=i;
	nc=(int)exn(3,k);
	}
	}*/
	if(k>0) {

		//周跳备选组合 、 周跳修正后的观测值残差
		cycle2=zeros(nv,nc);Lv=zeros(nv,nc);cycle3=zeros(nc,nv);E=eye(nv);

		for(i=0;i<nc;i++) {
			for(j=0;j<nv;j++) {
				cycle2[j+i*nv]=cycle1[j];
				Lv[j+i*nv]=v[j];
			}
		}

		for(i=0;i<k;i++) {   //全部的周跳备选组合(周跳的上下一周都为备选)
			sati=(vflg[i]>> 8)&0xFF;  //卫星
			type=(vflg[i]>> 4)&0xF;   //观测量类型
			freq=vflg[i]&0xF;         //频率
			lami=nav->lam[sati-1][freq];
			for(t=0,j=0,m=0;t<exn(3,i);t++) {
				m+=j;
				j=0;
				while(j<exn(3,k-i-1))   {cycle2[index[i]+(j+m)*nv]+=-1*lami;j++;}
				while(j<2*exn(3,k-i-1)) {cycle2[index[i]+(j+m)*nv]+= 0*lami;j++;}
				while(j<3*exn(3,k-i-1)) {cycle2[index[i]+(j+m)*nv]+= 1*lami;j++;}
			}
		}
		matmul("TN",nc,nv,nv,1.0,cycle2,E,0.0,cycle3);  //cycle3=cycle2';
		//		trace(2,"周跳备选组合cycle2：%d\n",nc);tracemat(2,cycle3,nc,nv,7,3);
		//根据备选组合，计算残差
		for(i=0;i<nc;i++) {
			for(j=0;j<nv;j++) {
				Lv[j+i*nv]-=cycle2[j+i*nv];
			}
		}	

		//搜索最优和次优解	
		ztss_sd(nv,nx,nc,H,Lv,R,cycle2,W,is,S,vtv);

		//输出最优解组合和次优解组合
		trace(2,"周跳最优解slip[0],vpv=%.3lf,RMS=%.3lf:\n",S[0],sqrt(vtv[0]/nv));
		tracemat(2,cycle2+is[0]*nv,1,nv,10,3);
		trace(2,"周跳次优解slip[1],vpv=%.3lf,RMS=%.3lf:\n",S[1],sqrt(vtv[1]/nv));
		tracemat(2,cycle2+is[1]*nv,1,nv,10,3);


		//计算周跳修复检查量《基于三差观测量的实时动态GPS周跳修复方法研究》20141006
		BP=mat(nx,nv);N=mat(nx,nx);NB=mat(nx,nv);BNB=mat(nv,nv);Ql=mat(nv,nv);
		Qv=mat(nv,nv);;Rv=mat(nv,nv);Qc=mat(nv,nv);dc=mat(nv,1);dcQc=mat(1,nv);

		matmul("NN",nx,nv,nv,1.0, H,R,0.0,BP);
		matmul("NT",nx,nx,nv,1.0,BP,H,0.0,N );
		matinv(N,nx);
		matmul("NN",nx,nv,nx,1.0,N, H,0.0, NB);
		matmul("TN",nv,nv,nx,1.0,H,NB,0.0,BNB);
		matcpy(Qv,BNB,nv,nv);
		matcpy(Ql,  R,nv,nv);   //Ql=inv(P);
		matinv(Ql,nv);

		matmul("NN",nv,nv,nv,1.0,E,Ql,-1.0,Qv);   //残差v的协因数阵Qv=Q-B*inv(N)*B'
		matmul("NN",nv,nv,nv,1.0,Qv,R,0.0,Rv);   //误差转残差转换矩阵Rv=Qv*P

		matmul("NN",nv,nv,nv,1.0, R,Rv,0.0,Qc);  //inv(Qc)=P-P*B*inv(N)*B'*P=P*Rv

		//		trace(2,"Qv:\n");tracemat(2,Qv,nv,nv,10,3);
		trace(2,"Rv:\n");tracemat(2,Rv,nv,nv,10,3);
		trace(2,"Qc:\n");tracemat(2,Qc,nv,nv,10,3);

		//dc=slip1-slip0 次优解-最优解
		matcpy(dc,cycle2+is[1]*nv,nv,1);
		matmul("NN",nv,1,nv,1.0,E,cycle2+is[0]*nv,-1.0,dc);

		trace(2,"dc:\n");tracemat(2,dc,1,nv,10,3);

		matmul("NN",1,nv,nv,1.0,  dc,Qc,0.0,dcQc);
		matmul("NN",1, 1,nv,4.0,dcQc,dc,0.0,&covd);
		Wd=(S[1]-S[0])/(sqrt(covd));                //周跳修复检查量

		trace(2,"vpv0:%10.3lf,vpv1=%10.3lf,covd=%10.3lf",S[0],S[1],covd);
		trace(2,"周跳修复检查量Wd:%10.3lf\n",Wd);

		*w_ratio=Wd;

		*ratio=S[1]/S[0];*vpv=S[0];*rms_v=vtv[0]/nv;
		for(i=0;i<nv;i++){
			sati=(vflg[i]>> 8)&0xFF;  //卫星
			type=(vflg[i]>> 4)&0xF;   //观测量类型
			freq=vflg[i]&0xF;         //频率
			stype=type==0?"L":(type==1?"L":"C");
			lami=nav->lam[sati-1][freq];
			vv=cycle2[i+is[0]*nv];  //周跳值，单位为m
			slip=QFH(vv)*SSWR(fabs(vv)/lami);  //4舍5入取整
			slp[i]=slip;
			if(slip!=0) {
				trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d\n",
					sati,stype,freq+1,V[i],slip);
			}								
		}
		free(cycle2);free(Lv);free(cycle3);free(E);free(index);	free(cycle1);

		free(BP);free(N);free(NB);free(BNB);free(Ql);
		free(Qv);free(Rv);free(Qc);free(dc);free(dcQc);
	}
}


/*- 三差法：双差+历元间差分  20140910--------------------- */
static int ddt2dx_(rtk_t *rtk, const nav_t *nav, double dt,
				  const double *x,const int *sat,const int ns, double *y,int n,double *e,double *azel,
				  const int *iu, const int *ir, double *dx, double *px)
{
	prcopt_t *opt=&rtk->opt;
	int i,j,k,m,f,sysi,sysj,nsat[2]={0};
	double lami,lamj,bl,dr[3];
	double *H,*Hi,*v,*R,*Ri,*Rj,*Q,*w,*ATP,*dxt;
	int nx=3,nv=0;
	double m2,*vtp;           //验后单位权中误差
	
	double m0=1,eps=2;        //粗差探测中使用
	double *P,*N,*V;
	int slip;              //slip:周跳数
	int *W;
	int vflg[MAXOBS*NFREQ*2+1],sati,freq,ff,type;   //记录站际历元间二次差分观测方程具体信息
	char *stype;
	double A[3*MAXOBS*2],L[MAXOBS*2],*Rl,*Rlk,*Rlj,*E;
	int satk,satj,t,nl=0;
	int nb[NFREQ*4]={0},b=0,lflg[MAXOBS*NFREQ*2+1];

	double vv,vtv[2],S[2],ratio=1.0,*cycle1,*cycle2,*Lv,*cycle3,ep=0.02;
	int nc=0,is[2],*index,flg=0;
	
	//	double m02=varerr(1,0.5,0,0,0,opt);	
	
	FILE *fp=fopen("E:\\Fixed\\DD模型\\dt_result.txt","a+");
	
	static double y0[MAXOBS*2],dx0[3],azel0[2*MAXSAT],bl0,dt0;
	static int   sat0[MAXSAT],ns0,n0,iu0[MAXSAT],ir0[MAXSAT];
	
	trace(3,"dt2dx_:%s\n",time_str(rtk->sol.time,3));

	for(i=0;i<MAXX;i++) now_slipflg[i]=0;   //当前历元周跳标记数组清0
	
	bl=baseline(x,rtk->rb,dr);
	
	if(bl0==0){  //初始历元
		matcpy(y0,y,2*opt->nf,n);        //站星距离
		matcpy(azel0,azel,2,n);          //卫星高度角
		memcpy(sat0,sat,MAXSAT);         //共视星序列
		memcpy(iu0,iu,MAXSAT);           //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);           //基准共视星序列
		bl0=bl;                          //基线长
		dt0=dt;
		ns0=ns;                          //共视卫星数
		n0=n;
		return 0;
	}
	
	trace(2,"y0:\n");   tracemat(2,   y0,2*opt->nf,n0,12,3);
	for(i=0;i<ns0;i++) {
		trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat0[i],iu0[i],ir0[i],azel0[1+iu0[i]*2]*R2D);
	}
	
	trace(2,"y:\n");   tracemat(2,   y,2*opt->nf,n,12,3);
	for(i=0;i<ns;i++) {
		trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat[i],iu[i],ir[i],azel[1+iu[i]*2]*R2D);
	}
	
	H=zeros(nx,opt->nf*ns);v=mat(opt->nf*ns,1);
	R=mat(opt->nf*ns,opt->nf*ns);
	Ri=mat(opt->nf*ns,1);Rj=mat(opt->nf*ns0,1);
	
	//站间历元二次差分
	for(m=0;m<2;m++) /* for each system (0:gps/qzss/sbas,1:glonass) */
		for(f=0;f<opt->nf;f++) {  //只是用了载波相位观测量
			for(i=0;i<ns;i++) { 
				sysi=rtk->ssat[sat[i]-1].sys;
				if((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
				if(!validobs(iu[i],ir[i],f,opt->nf,y))  continue;
				
				for(j=0;j<ns0;j++) {
					sysj=rtk->ssat[sat0[j]-1].sys;
					if((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
					if(sat[i]!=sat0[j])  continue;  //历元间差分，卫星匹配
					if(!validobs(iu0[j],ir0[j],f,opt->nf,y0))   continue;
										
					lami=nav->lam[sat[i]-1][f];
					lamj=nav->lam[sat0[j]-1][f];
					if(lami<0.0||lamj<0.0) continue;
					
					Hi=H+nv*nx;
					
					v[nv]=(y[f+iu[i]*2*opt->nf]-y[f+ir[i]*2*opt->nf])-(y0[f+iu0[j]*2*opt->nf]-y0[f+ir0[j]*2*opt->nf]);					
					
					for(k=0;k<3;k++) {      //每个频率上设一个相对钟差参数
						Hi[k]=-e[k+iu[i]*3];
					}  
					
					if(sysi==SYS_GLO) {nsat[1]++;}
					else              {nsat[0]++;}
					
					Ri[nv]=varerr(sysi,azel[1+iu[i]*2],bl,dt,f,opt);
                    Rj[nv]=varerr(sysj,azel0[1+iu0[j]*2],bl0,dt0,f,opt);				
					
					trace(4,"sat=%2d  v=%7.3f\n",sat[i],v[nv]);
					
					//vflg数据记录站际历元间二次差分观测方程具体信息
					//其结构依次为：卫星，观测量类型，频率
					vflg[nv++]=((sat[i]<<8)|((f<opt->nf?0:1)<<4)|(f%opt->nf));
					break;
				}
			}
		}

		/* 站际历元二次差分观测方程权阵 */
		dtcov(nv,Ri,Rj,R);

		trace(2,"站际历元二次差分：\n");
		trace(2,"H= %d*%d\n",nx,nv); tracemat(2,H,nx,nv,10,6);
		trace(2,"v=\n"); tracemat(2,v,1,nv,10,4);
		trace(2,"R= %d阶\n",nv); tracemat(2,R,nv,nv,10,6);

		Rlk=mat(nv,1);Rlj=mat(nv,1);

		//星间差分
		for (m=0;m<2;m++) /* for each system (0:gps/qzss/sbas,1:glonass) */	

		for (f=0;f<opt->nf;f++) {
			
    		/* search reference satellite with highest elevation */
			for (i=-1,j=0;j<ns;j++) {
				sysi=rtk->ssat[sat[j]-1].sys;
				if ((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
				if (!validobs(iu[j],ir[j],f,opt->nf,y)) continue;
				if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
			}
			if (i<0) continue;

			for(k=0;k<nv;k++) {    //找出参考星对应的残差
				satk=(vflg[k]>> 8)&0xFF;  //卫星
				type=(vflg[i]>> 4)&0xF;   //观测量类型
				freq=vflg[i]&0xF;         //频率
				if((sat[i]==satk)&&(f==type)&&(f%opt->nf==freq)) break; //匹配
			}
			
			for (j=0;j<nv;j++) {          //星间差分
				satj=(vflg[j]>> 8)&0xFF;  //卫星
				type=(vflg[j]>> 4)&0xF;   //观测量类型
				freq=vflg[j]&0xF;         //频率
				lamj=nav->lam[satj-1][freq];
				if (satk==satj) continue;
				sysi=rtk->ssat[satk-1].sys;
				sysj=rtk->ssat[satj-1].sys;
				if ((m==0&&sysj==SYS_GLO)||(m==1&&sysj!=SYS_GLO)) continue;
				
				ff=f%opt->nf;
				lami=nav->lam[satk-1][ff];
				lamj=nav->lam[satj-1][ff];
				if (lami<=0.0||lamj<=0.0) continue;
				
				/* 星间差分 */
				L[nl]=v[j]-v[k];    //非参考星-参考星(有利于后面周跳处理改正在非参考星上，
									//因为符号相同)
				
				/* partial derivatives by rover position */
				for (t=0;t<3;t++) {
					A[t+3*nl]=H[t+j*nx]-H[t+k*nx];
				}
				
	        	Rlk[nl]=Ri[k]+Rj[k];
				Rlj[nl]=Ri[j]+Rj[j];

				nb[b]++;
				//lflg数据记录三差观测方程具体信息
				//其结构依次为：卫星，观测量类型，频率
				lflg[nl++]=((satk<<16)|(satj<<8)|((f<opt->nf?0:1)<<4)|(f%opt->nf));
			}
			b++;
		}
		trace(2,"继续星间差分构成三差：\n");
		trace(2,"A'=%d*%d\n",3,nl);  tracemat(2,A,3,nl,10,6);
		trace(2,"l=\n"); tracemat(2,L,1,nl,10,4);
		
		if(nl<3) {
			trace(2,"历元间差分方程个数nv=%d,参数个数nx=%d,观测方程不够，无法解算！\n",nl,3);
			free(H);free(v);free(R);free(Ri);free(Rj);free(Rlk);free(Rlj);

			matcpy(y0,y,2*opt->nf,n);    //站星距离
			matcpy(azel0,azel,2,n);      //卫星高度角
			memcpy(sat0,sat,MAXSAT);     //共视星序列
			memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
			memcpy(ir0,ir,MAXSAT);       //基准共视星序列
			bl0=bl;                      //基线长
			dt0=dt;
			ns0=ns;                      //共视卫星数
			n0=n;
			return 0;
		}
		
		dxt=mat(nx,1);Q=zeros(nx,nx);w=mat(nv,1);ATP=zeros(nx,nv);vtp=mat(1,nv);
		Rl=zeros(nl,nl);
		/* 三差权阵 */
		ddcov(nb,b,Rlk,Rlj,nl,Rl);//k为参考星对应在nv中的位置，即在R对角上的位置

		/* 采用等权的方式定权20140729 */
//		R=eye(nl);
		
		trace(2,"Rl= %d阶\n",nl); tracemat(2,Rl,nl,nl,10,6);
		
		P=mat(nl,nl);N=mat(3*(3+1)/2,1);V=mat(nl,1);
		W=imat(nl,1);		
			
		if(!matinv(Rl,nl)) {

			 //先验方差给的过大，后面方差一致性检验无法通过，给予7倍权尝试，可以通过		 
			E=eye(nl);
			matmul("NN",nl,nl,nl,m0*m0,E,Rl,0.0,P);   //P=m0^2/D
			free(E);
			
			//粗差探测
			while(1) {

			DataSnooping_1(nl,3,A,L,P,eps,m0,dxt,N,V,W,&m2);	
		
			for(i=0;i<nl;i++) {           //当前历元周跳探测情况，周跳标记为1
				sati=(lflg[i]>> 8)&0xFF;  //卫星
				freq=lflg[i]&0xF;         //频率					
				f=freq;
				if(W[i]==0)	now_slipflg[IB(sati,f,&rtk->opt)]=1;  //当前历元发生周跳，标记为1
				else        now_slipflg[IB(sati,f,&rtk->opt)]=2;  //未发生周跳，标记为2	
				trace(2,"sat=%3d,now_slipflg=%3d\n",sati,now_slipflg[IB(sati,f,&rtk->opt)]);
			}
			/* now_slipflg数据标记：-------------------------------------------
			   0：  未参与站际历元二次差分模型的卫星，是否发生周跳未知；
			        或者站际历元二次差分模型计算失败（有效卫星数不够），
					周跳情况未知；
			   1：  表示站际历元二次差分模型探测出的周跳；
			   2：  通过站际历元二次差分模型探测未发生周跳或者发生周跳已经修复。
			-----------------------------------20140731-------------------------*/
			
			if(m2<0) {
				printf("计算失败\n");getchar();
				matcpy(y0,y,2*opt->nf,n);           //站星距离
				matcpy(azel0,azel,2,n);             //卫星高度角
				memcpy(sat0,sat,MAXSAT);            //共视星序列
				memcpy(iu0,iu,MAXSAT);              //流动站共视星序列
				memcpy(ir0,ir,MAXSAT);              //基准共视星序列
				bl0=bl;                             //基线长
				dt0=dt;
				ns0=ns;                             //共视卫星数
				n0=n;
				free(H);free(v);free(R);free(Ri);free(Rj);
				free(Q);free(w);free(ATP);free(dxt);
				free(P);free(V);free(W);free(N);
				free(Rl);free(Rlk);free(Rlj);
				return 0;
			} 
			flg=0;k=0;nc=1;                    //记录周跳个数，当周跳个数为0时，结束解算
			cycle1=zeros(nl,1);                 //记录周跳发生情况（4舍5入值）
			index=imat(nl,1);
			for(i=0;i<nl;i++) {                 /* 记录周跳发生的状况 */
				if(W[i]==0) {
					index[k]=i;               //发生周跳的卫星在参数中国的位置
					k++;                      //周跳的个数
					nc*=3;                    //发生周跳，搜索上下一周，共3个
					sati=(lflg[i]>> 8)&0xFF;  //卫星
					type=(lflg[i]>> 4)&0xF;   //观测量类型
					freq=lflg[i]&0xF;         //频率
					stype=type==0?"L":(type==1?"L":"C");
					lami=nav->lam[sati-1][freq];
					slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

					trace(1,"Slip detected:time:%s\n",time_str(rtk->sol.time,3));
					trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d\n",
						sati,stype,freq+1,V[i],slip);

					cycle1[i]=slip*lami;    //转换为以米为单位
					vv=V[i];
					vv-=slip*lami;
					if(fabs(vv)>0.02)  {
						flg++;
						trace(1,"残差超限，周跳修复失败：周跳修复后残差V:%6.3lf\n",vv);
						continue;
					}
				}
			}
			if(flg>0){
				k=nl;
				for(i=0;i<nl;i++) {					
					index[i]=i;
					nc=(int)exn(3,k);
				}
			}
			if(k>0) {
				
				//周跳备选组合 、 周跳修正后的观测值残差
				cycle2=zeros(nl,nc);Lv=zeros(nl,nc);cycle3=zeros(nc,nl);E=eye(nl);
				
				for(i=0;i<nc;i++) {
					for(j=0;j<nl;j++) {
						cycle2[j+i*nl]=cycle1[j];
						Lv[j+i*nl]=L[j];
					}
				}
				
				for(i=0;i<k;i++) {   //全部的周跳备选组合(周跳的上下一周都为备选)
					sati=(lflg[i]>> 8)&0xFF;  //卫星
					type=(lflg[i]>> 4)&0xF;   //观测量类型
					freq=lflg[i]&0xF;         //频率
					lami=nav->lam[sati-1][freq];
					for(t=0,j=0,m=0;t<exn(3,i);t++) {
						m+=j;
						j=0;
						while(j<exn(3,k-i-1))   {cycle2[index[i]+(j+m)*nl]+=-1*lami;j++;}
						while(j<2*exn(3,k-i-1)) {cycle2[index[i]+(j+m)*nl]+= 0*lami;j++;}
						while(j<3*exn(3,k-i-1)) {cycle2[index[i]+(j+m)*nl]+= 1*lami;j++;}
					}
				}
				matmul("TN",nc,nl,nl,1.0,cycle2,E,0.0,cycle3);  //cycle3=cycle2';
				trace(2,"周跳备选组合cycle2：%d\n",nc);tracemat(2,cycle3,nc,nl,7,3);
				//根据备选组合，计算残差
				for(i=0;i<nc;i++) {
					for(j=0;j<nl;j++) {
						Lv[j+i*nl]-=cycle2[j+i*nl];
					}
				}	
				
				//搜索最优和次优解
				ztss_dd(nl,nx,nc,A,Lv,P,is,S,vtv);
				
				//输出最优解组合和次优解组合
				trace(2,"周跳最优解slip[0],vpv=%.3lf,RMS=%.3lf:\n",S[0],sqrt(vtv[0]/nl));
				tracemat(2,cycle2+is[0]*nl,1,nl,10,3);
				trace(2,"周跳次优解slip[1],vpv=%.3lf,RMS=%.3lf:\n",S[1],sqrt(vtv[1]/nl));
				tracemat(2,cycle2+is[1]*nl,1,nl,10,3);
				
				ratio=S[1]/S[0];
				if(ratio>5&&sqrt(vtv[0]/nl)<ep) {
					for(i=0;i<nl;i++){
						sati=(lflg[i]>> 8)&0xFF;  //卫星
						type=(lflg[i]>> 4)&0xF;   //观测量类型
						freq=lflg[i]&0xF;         //频率
						stype=type==0?"L":(type==1?"L":"C");
						lami=nav->lam[sati-1][freq];
						vv=cycle2[i+is[0]*nl];  //周跳值，单位为m
						slip=QFH(vv)*SSWR(fabs(vv)/lami);  //4舍5入取整
						if(slip==0) k--;
						L[i]-=slip*lami;      //用于再进行一次平差求得更为精确的坐标差信息
						if(slip!=0) {
							trace(1,"Slip detected:time:%s\n",time_str(rtk->sol.time,3));
							trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d\n",
								sati,stype,freq+1,V[i],slip);
						}
						
						//将周跳数记录到slips数组中
						f=freq;
						slips[IB(sati,f,&rtk->opt)]+=slip;	
						
						//当前历元周跳修复，修复标记为2
						now_slipflg[IB(sati,f,&rtk->opt)]=2;
						
						W[i]=1;					
					}
				}
				free(cycle2);free(Lv);free(cycle3);free(E);free(index);			
			}

			//对残差进行
			for(i=0;i<nl;i++) {                /* 记录周跳发生的状况 */
				if(W[i]==0) {
					k++;
					sati=(lflg[i]>> 8)&0xFF;  //卫星
					type=(lflg[i]>> 4)&0xF;   //观测量类型
					freq=lflg[i]&0xF;         //频率
					stype=type==0?"L":(type==1?"L":"C");
					lami=nav->lam[sati-1][freq];
					slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整
		     		if(slip==0) k--;
					L[i]-=slip*lami;      //用于再进行一次平差求得更为精确的坐标差信息
					trace(1,"Slip detected:time:%s\n",time_str(rtk->sol.time,3));
					trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d\n",
						sati,stype,freq+1,V[i],slip);
					
					//将周跳数记录到slips数组中
					f=freq;
					slips[IB(sati,f,&rtk->opt)]+=slip;	

				    //当前历元周跳修复，修复标记为2
					now_slipflg[IB(sati,f,&rtk->opt)]=2;
				}
			}
			fprintf(fp,"%s %.4lf\n",time_str(rtk->sol.time,3),m2);
			fclose(fp);
	//		break;           //20140708 存在周跳误判的情况，有待进一步完善
			if(k==0) break;  //当不存在周跳后，结束计算
			}

			for(i=0;i<3;i++) dx[i]=dxt[i];
        
/*			m2=m2*m2;         //使用验后单位权中误差		
			px[0]=m2*N[ij(0,0)];px[1]=m2*N[ij(1,1)];px[2]=m2*N[ij(2,2)];
			px[3]=m2*N[ij(1,0)];px[4]=m2*N[ij(2,0)];px[5]=m2*N[ij(2,1)];*/

			m0=m0*m0;   //使用验前单位权中误差		
			px[0]=m0*N[ij(0,0)];px[1]=m0*N[ij(1,1)];px[2]=m0*N[ij(2,2)];
			px[3]=m0*N[ij(1,0)];px[4]=m0*N[ij(2,0)];px[5]=m0*N[ij(2,1)];  
						
			trace(2,"Dx:\n");    tracemat(2,px,1,6,9,7);
		}
		else {
			matcpy(y0,y,2*opt->nf,n);           //站星距离
			matcpy(azel0,azel,2,n);     //卫星高度角
			memcpy(sat0,sat,MAXSAT);     //共视星序列
			memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
			memcpy(ir0,ir,MAXSAT);       //基准共视星序列
			bl0=bl;                      //基线长
			dt0=dt;
			ns0=ns;                       //共视卫星数
			n0=n;  
			free(H);free(v);free(R);free(Ri);free(Rj);
			free(Q);free(w);free(ATP);free(dxt);
			free(P);free(V);free(W);free(N);
			free(Rl);free(Rlk);free(Rlj);
			return 0;
		}
		
		/* 历元间差分计算成功，记录当前历元信息 */
		matcpy(y0,y,2*opt->nf,n);           //站星距离
		matcpy(azel0,azel,2,n);     //卫星高度角
		memcpy(sat0,sat,MAXSAT);     //共视星序列
		memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);       //基准共视星序列
		bl0=bl;                      //基线长
		dt0=dt;
		ns0=ns;                      //共视卫星数
		n0=n;
		
		free(H);free(v);free(R);free(Ri);free(Rj);
		free(Q);free(w);free(ATP);free(dxt);
		free(P);free(V);free(W);free(N);
		free(Rl);free(Rlk);free(Rlj);
		
		return 1;
}


/* 站间单差，历元间差分计算坐标增量(单历元—法方程叠加) 2013.12.09 */
static int dt2dx_(rtk_t *rtk, const nav_t *nav, double dt,
	const double *x,const int *sat,const int ns, double *y,int n,double *e,double *azel,
	const int *iu, const int *ir, double *dx, double *px)
{
	prcopt_t *opt=&rtk->opt;
	int i,j,k,m,f,sysi,sysj,nsat[2]={0};
	double lami,lamj,bl,dr[3];
	double *H,*Hi,*v,*R,*Ri,*Rj,*Q,*w,*ATP,*dxt;
	int nx=5,nv=0;
	double m2,*vtp;           //验后单位权中误差

	double m0=1,eps=2;    //粗差探测中使用
	double *P,*N,*V;
	int *W;
	int slip;    //slip:周跳数
	int vflg[MAXOBS*NFREQ*2+1],sati,freq,type;   //记录站际历元间二次差分观测方程具体信息
	char *stype;
	double *Ws;   //Co_bust中使用

	double vv,ratio=1.0,ep=0.02;
	int  flg=0;
	int info=0; //粗差探测是否成功标志
	int suc_info=1;
	double rms_v=0.0,vpv=0.0;
	int *slp;
	int info_LEGE=1;
	double *vc;
	double w_ratio=1.0;

	static double y0[MAXOBS*2],dx0[3],azel0[2*MAXSAT],bl0,dt0;
	static int   sat0[MAXSAT],ns0,n0,iu0[MAXSAT],ir0[MAXSAT];
	static double e0[MAXOBS*3];    //用于随机模型补偿二次差分模型舍去项误差


	//抗差估计
	//fact:权因子函数选择：0---IGG1;1---IGG2;2---Huber;
	//k0:等价权函数保权界值；k1:等价权函数零权界值
	//eps: 限差；
	double fact=0,eps_x=0.001,k0=1.5,k1=3.0,*R1; 

	trace(3,"dt2dx_:%s\n",time_str(rtk->sol.time,3));

	for(i=0;i<MAXSLIP;i++) now_slipflg[i]=0;   //当前历元周跳标记数组清0

	bl=baseline(x,rtk->rb,dr);

	if(bl0==0){  //初始历元
		matcpy(y0,y,2*opt->nf,n);        //站星距离
		matcpy(azel0,azel,2,n);          //卫星高度角
		memcpy(sat0,sat,MAXSAT);         //共视星序列
		memcpy(iu0,iu,MAXSAT);           //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);           //基准共视星序列
		bl0=bl;                          //基线长
		dt0=dt;
		ns0=ns;                          //共视卫星数
		n0=n;
		return 0;
	}

	trace(2,"y0:\n");   tracemat(2,   y0,2*opt->nf,n0,12,3);
	for(i=0;i<ns0;i++) {
		trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat0[i],iu0[i],ir0[i],azel0[1+iu0[i]*2]*R2D);
	}

	trace(2,"y:\n");   tracemat(2,   y,2*opt->nf,n,12,3);
	for(i=0;i<ns;i++) {
		trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat[i],iu[i],ir[i],azel[1+iu[i]*2]*R2D);
	}

	H=zeros(nx,opt->nf*ns);v=mat(opt->nf*ns,1);
	R=mat(opt->nf*ns,opt->nf*ns);
	Ri=mat(opt->nf*ns,1);Rj=mat(opt->nf*ns0,1);

	for(m=0;m<2;m++) for(f=0;f<opt->nf;f++) {  //只是用了载波相位观测量
		for(i=0;i<ns;i++) 
		{ 
			sysi=rtk->ssat[sat[i]-1].sys;
			if((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
			if(!validobs(iu[i],ir[i],f,opt->nf,y))  continue;

			for(j=0;j<ns0;j++) {
				sysj=rtk->ssat[sat0[j]-1].sys;
				if((m==0&&sysi==SYS_GLO)||(m==1&&sysi!=SYS_GLO)) continue;
				if(sat[i]!=sat0[j])  continue;  //历元间差分，卫星匹配
				if(!validobs(iu0[j],ir0[j],f,opt->nf,y0))   continue;

				lami=nav->lam[sat[i]-1][f];
				lamj=nav->lam[sat0[j]-1][f];
				if(lami<0.0||lamj<0.0) continue;

				Hi=H+nv*nx;

				v[nv]=(y[f+iu[i]*2*opt->nf]-y[f+ir[i]*2*opt->nf])-(y0[f+iu0[j]*2*opt->nf]-y0[f+ir0[j]*2*opt->nf]);

				for(k=0;k<3;k++) {   //每个频率上设一个相对钟差参数
					Hi[k]=-e[k+iu[i]*3];
				}  


				Hi[3]=1.0;   //每个频率上设一个相对钟差参数
				if(sysi==SYS_GLO) {Hi[4]=1.0;nsat[1]++;}
				else              {Hi[4]=0.0;nsat[0]++;}

				Ri[nv]=varerr(sysi,azel[1+iu[i]*2],bl,dt,f,opt);
				Rj[nv]=varerr(sysj,azel0[1+iu0[j]*2],bl0,dt0,f,opt);				

				trace(4,"sat=%2d  v=%7.3f\n",sat[i],v[nv]);

				//vflg数据记录站际历元间二次差分观测方程具体信息
				//其结构依次为：卫星，观测量类型，频率
				vflg[nv++]=((sat[i]<<8)|((f<opt->nf?0:1)<<4)|(f%opt->nf));
				break;
			}
		}
	}

	/* shrink design matrix:nx=5->4 */
	if(!nsat[0]||!nsat[1]) {
		nx=4;
		for(i=0;i<nv;i++) {
			for(j=0;j<4;j++) {
				H[j+i*4]=H[j+i*5];				
			}
		}
	} 

	trace(2,"H= %d*%d\n",nx,nv); tracemat(2,H,nx,nv,10,6);
	trace(2,"v=\n"); tracemat(2,v,nv,1,10,4);

	if(nv<nx) {
		trace(2,"历元间差分方程个数nv=%d,参数个数nx=%d,观测方程不够，无法解算！\n",nv,nx);
		free(H);free(v);free(R);free(Ri);free(Rj);

		matcpy(y0,y,2*opt->nf,n);           //站星距离
		matcpy(azel0,azel,2,n);     //卫星高度角
		memcpy(sat0,sat,MAXSAT);     //共视星序列
		memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);       //基准共视星序列
		bl0=bl;                      //基线长
		dt0=dt;
		ns0=ns;                       //共视卫星数
		n0=n;
		return 0;
	}

	dxt=mat(nx,1);Q=zeros(nx,nx);w=mat(nv,1);ATP=zeros(nx,nv);vtp=mat(1,nv);

	/* 历元间差分观测方程权阵 */
	dtcov(nv,Ri,Rj,R);

	trace(2,"R= %d阶\n",nv); tracemat(2,R,nv,nv,10,6);

	P=mat(nv,1);N=mat(nx*(nx+1)/2,1);V=mat(nv,1);
	W=imat(nv,1);Ws=mat(nv,1);vc=mat(nv,1);

	//多普勒积分法探测到的周跳标记，周跳数>1才标记(有时可能出现虚报的情况)
	for(i=0;i<nv;i++) {           
		sati=(vflg[i]>> 8)&0xFF;  //卫星
		freq=vflg[i]&0xF;         //频率					
		f=freq;
		W[i]=1;
		//		if((now_slipflg[IB(sati,f,&rtk->opt)]==1)&&(slips_doppler[IB(sati,f,&rtk->opt)])>1)
		/*			if(now_slipflg[IB(sati,f,&rtk->opt)]==1)
		W[i]=0; //由多普勒积分法探测到的周跳
		else W[i]=1;	*/			
		trace(2,"sat=%3d,now_slipflg=%3d,slips_doppler=%3d\n",
			sati,now_slipflg[IB(sati,f,&rtk->opt)],slips_doppler[IB(sati,f,&rtk->opt)]);
	}

	if(!matinv(R,nv)) {

		//先验方差给的过大，后面方差一致性检验无法通过，给予7倍权尝试，可以通过		
		for(i=0;i<nv;i++) P[i]=R[i+i*nv]*m0*m0;  

		for(i=0;i<nv;i++) R[i+i*nv]=R[i+i*nv]*m0*m0; 

		R1=mat((nv+1)*nv/2,1);
		for(i=0;i<nv;i++) 
			for(j=0;j<=i;j++) R1[ij(i,j)]=R[i+j*nv];   //取下三角矩阵


		//info=0为初次平差没有发现粗差；1,2,3：依次表示粗差探测第1，第2，第3种情况
		trace(2,"v:\n");tracemat(2,v,1,nv,10,3);
		info=DataSnooping(nv,nx,H,v,P,eps,m0,dxt,N,V,W,&m2);	//粗差探测
		if(info>0) {  //探测到粗差
			trace(1,"\nTime:%s\n",time_str(rtk->sol.time,3));
			trace(1,"Detect gloss by DataSnooping,Situation %d\n",info);
			trace(2,"粗差探测结果：W=");tracemati(2,W,1,nv,4);			
		}

		info_LEGE=1;  //粗差探测2默认设置为1
		if(info>1) {   //情况2,3
			flg=0;
			for(i=0;i<nv;i++) {                 /* 记录周跳发生的状况 */
				if(W[i]==0) {							
					sati=(vflg[i]>> 8)&0xFF;  //卫星
					type=(vflg[i]>> 4)&0xF;   //观测量类型
					freq=vflg[i]&0xF;         //频率
					stype=type==0?"L":(type==1?"L":"C");
					lami=nav->lam[sati-1][freq];
					slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

					trace(1,"Slip detected by DataSnooping:time:%s\n",time_str(rtk->sol.time,3));
					trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d,周跳修复后残差V:%6.3lf\n",
						sati,stype,freq+1,V[i],slip,V[i]-lami*slip);

					vv=V[i];
					vv-=slip*lami;
					if(fabs(vv)>0.02)  {
						flg++;
						trace(1,"残差超限，DataSnooping探测周跳定位可能错误，周跳修复后残差V:%6.3lf\n",vv);
						continue;
					}
				}
			}

			if(flg>0) { //探测粗差

				if(info==3&&flg>0) {             //粗差探测第3种情况
					for(i=0;i<nv;i++) W[i]=1;    //标记卫星清零
				}
			}
			//如果有效观测量只有5,采用搜索法确定粗差(这种情况可能已经处理)
			if(nv==(nx+1)) for(i=0;i<nv;i++) W[i]=1;   

			trace(1,"Detect gloss by DataSnooping_single:\n");
			DataSnooping_single(nv,nx,H,v,P,m0,vflg,nav,dxt,N,V,W,&m2);

		}

		//如果有效观测量只有5,采用搜索法确定粗差(Datasnooping可能无法探测粗差失败)
		if(nv==(nx+1)) {
			for(i=0;i<nv;i++) W[i]=1;   
			trace(1,"Detect gloss by DataSnooping_single:\n");
			DataSnooping_single(nv,nx,H,v,P,m0,vflg,nav,dxt,N,V,W,&m2); 		
		}

		if(nv>nx&&info_LEGE>0) {
			matcpy(vc,v,nv,1);
			for(i=0,flg=0;i<nv;i++) {         /* 记录周跳发生的状况 */
				if(W[i]==0) {							
					sati=(vflg[i]>> 8)&0xFF;  //卫星
					type=(vflg[i]>> 4)&0xF;   //观测量类型
					freq=vflg[i]&0xF;         //频率
					stype=type==0?"L":(type==1?"L":"C");
					lami=nav->lam[sati-1][freq];
					slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

					trace(2,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d,周跳修复后残差V:%6.3lf\n",
						sati,stype,freq+1,V[i],slip,V[i]-lami*slip);
					vc[i]-=slip*lami;
				}

			}
			m2=CoRobust(nv,nx,H,vc,R1,fact,k0,k1,m0,eps_x,dxt,N,V,Ws); //抗差估计
			trace(2,"抗差估计 by CoRobust,W：");tracemat(2,Ws,1,nv,10,3);
			for(i=0;i<nv;i++) {
				V[i]+=v[i]-vc[i];
				if(Ws[i]==0.0) W[i]=0;
				else W[i]=1;
				if(fabs(v[i]-vc[i])>0.01) W[i]=0;//含有周跳
			}
		}
		else{ 
			m2=-1;
		}
		for(i=0;i<3;i++) dx[i]=dxt[i];	

		for(i=0;i<nv;i++) {         /* 当前历元周跳探测情况，周跳标记为1 */
			sati=(vflg[i]>> 8)&0xFF;  //卫星
			freq=vflg[i]&0xF;         //频率					
			f=freq;
			if(W[i]==0)	now_slipflg[IB(sati,f,&rtk->opt)]=1;  //当前历元发生周跳，标记为1
			else        now_slipflg[IB(sati,f,&rtk->opt)]=2;  //未发生周跳，标记为2	
			trace(2,"sat=%3d,now_slipflg=%3d\n",sati,now_slipflg[IB(sati,f,&rtk->opt)]);
		}
		/* now_slipflg数据标记：-------------------------------------------
		0：  未参与站际历元二次差分模型的卫星，是否发生周跳未知；
		或者站际历元二次差分模型计算失败（有效卫星数不够），
		周跳情况未知；
		1：  表示站际历元二次差分模型探测出的周跳；
		2：  通过站际历元二次差分模型探测未发生周跳或者发生周跳已经修复。
		-----------------------------------20140731-------------------------*/			
		for(i=0,k=0;i<nv;i++) {
			if(W[i]==0) 	k++;   //可能发生周跳的个数
		}
		//Modified by CYJ 2014/10/15
		//如果有效观测为参数个数，即nv=nx时，Datasnooping会探测失败，本文周跳探测模型无效
		//周跳标记应该全部为0，即周跳是否发生无从而知，考虑到rtklib针对单频数据没有进行周跳
		//探测，建议将其标记为1.若以后补充了其他模型可以考虑将其标记为0(即是否周跳未知)
		if(nv==nx) {  
			for(i=0;i<nv;i++){
				now_slipflg[IB(sati,f,&rtk->opt)]=1;
				slips[IB(sati,f,&rtk->opt)]=0;	
				suc_info=0;   //周跳探测修复失败，需跳出循环
			}
			k=0;  
		}

		if(k>0) {			
			slp=imat(nv,1);for(i=0;i<nv;i++) slp[i]=0;
			//周跳搜索计算
			rms_v=10;  //默认设置
			ratio=1;
			if(info_LEGE>0) {  //粗差探测成功才搜索，否则没必要搜索计算
				seach_slip(nv,nx,W,vflg,H,R,v,V,nav,&ratio,&rms_v,&vpv,slp,&w_ratio);//周跳探测、修复
			}


			trace(2,"w_ratio:%10.3lf\n",w_ratio);
			trace(2,"slp:\n");tracemati(2,slp,1,nv,4);
			if((rms_v<ep&&ratio>20)||(w_ratio>2.01)) {

				for(i=0;i<nv;i++){
					sati=(vflg[i]>> 8)&0xFF;  //卫星
					type=(vflg[i]>> 4)&0xF;   //观测量类型
					freq=vflg[i]&0xF;         //频率
					stype=type==0?"L":(type==1?"L":"C");
					lami=nav->lam[sati-1][freq];
					slip=slp[i];  //4舍5入取整
					v[i]-=slip*lami;      //用于再进行一次平差求得更为精确的坐标差信息
					if(slip==0) k--;
					if(slip!=0) {
						trace(1,"Slip repaired success:time:%s\n",time_str(rtk->sol.time,3));
						trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d   修复后残差V:%6.3lf\n",
							sati,stype,freq+1,V[i],slip,V[i]-slip*lami);
					}

					//将周跳数记录到slips数组中
					f=freq;
					slips[IB(sati,f,&rtk->opt)]+=slip;	

					//当前历元周跳修复，修复标记为2
					now_slipflg[IB(sati,f,&rtk->opt)]=2;					
					W[i]=1;	

				}
			}
			else {
				trace(1,"Slip repaired failed:time:%s\n",time_str(rtk->sol.time,3));

				suc_info=0;   //周跳探测修复失败，需跳出循环
				for(i=0;i<nv;i++){
					sati=(vflg[i]>> 8)&0xFF;  //卫星
					type=(vflg[i]>> 4)&0xF;   //观测量类型
					freq=vflg[i]&0xF;         //频率
					f=freq;
					if(info==1) { //对于粗差第一种情况，只对周跳卫星重新初始化
						if(!W[i])//周跳探测失败，存在多处周跳
						{
							now_slipflg[IB(sati,f,&rtk->opt)]=1;
							slips[IB(sati,f,&rtk->opt)]=0; 
							rtk->ssat[sati-1].slipN[f]++;

						}						
					}
					else {     //对于粗差第二和第三种情况，对全部卫星重新初始化
						now_slipflg[IB(sati,f,&rtk->opt)]=1;
						slips[IB(sati,f,&rtk->opt)]=0; 
						//rtk->ssat[sati-1].slipN[f]++;
					}
				}
			}
			free(slp);
		}
		free(R1);

//		if(k==0||suc_info==0) {free(R1);}  //当不存在周跳后或者修复失败，结束计算

		if(m2<0) {  
			printf("计算失败\n");getchar();
			matcpy(y0,y,2*opt->nf,n);           //站星距离
			matcpy(azel0,azel,2,n);             //卫星高度角
			memcpy(sat0,sat,MAXSAT);            //共视星序列
			memcpy(iu0,iu,MAXSAT);              //流动站共视星序列
			memcpy(ir0,ir,MAXSAT);              //基准共视星序列
			//matcpy( e0,e,3,MAXOBS);             //站星方向矢量
			bl0=bl;                             //基线长
			dt0=dt;
			ns0=ns;                             //共视卫星数
			n0=n;
			free(H);free(v);free(R);free(Ri);free(Rj);
			free(Q);free(w);free(ATP);free(dxt);free(vtp);
			free(P);free(V);free(W);free(N);free(Ws);free(vc);
			return 0;
		}

		for(i=0;i<3;i++) dx[i]=dxt[i];		

		m0=m0*m0;   //使用验前单位权中误差		
		px[0]=m0*N[ij(0,0)];px[1]=m0*N[ij(1,1)];px[2]=m0*N[ij(2,2)];
		px[3]=m0*N[ij(1,0)];px[4]=m0*N[ij(2,0)];px[5]=m0*N[ij(2,1)];  

		trace(2,"Dx:\n");    tracemat(2,px,1,6,9,7);
		trace(2,"钟差参数dt:%10.3lf\n",dxt[3]);    
	}
	else {
		matcpy(y0,y,2*opt->nf,n);           //站星距离
		matcpy(azel0,azel,2,n);     //卫星高度角
		memcpy(sat0,sat,MAXSAT);     //共视星序列
		memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);       //基准共视星序列
		//			matcpy( e0,e,3,MAXOBS);      //站星方向矢量
		bl0=bl;                      //基线长
		dt0=dt;
		ns0=ns;                       //共视卫星数
		n0=n;  
		free(H);free(v);free(R);free(Ri);free(Rj);
		free(Q);free(w);free(ATP);free(dxt);free(vtp);
		free(P);free(V);free(W);free(N);free(Ws);free(vc);
		return 0;
	}


	/* 历元间差分计算成功，记录当前历元信息 */
	matcpy(y0,y,2*opt->nf,n);           //站星距离
	matcpy(azel0,azel,2,n);     //卫星高度角
	memcpy(sat0,sat,MAXSAT);     //共视星序列
	memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
	memcpy(ir0,ir,MAXSAT);       //基准共视星序列
	bl0=bl;                      //基线长
	dt0=dt;
	ns0=ns;                       //共视卫星数
	n0=n;

	free(H);free(v);free(R);free(Ri);free(Rj);
	free(Q);free(w);free(ATP);free(dxt);free(vtp);
	free(P);free(V);free(W);free(N);free(Ws);free(vc);

	return 1;
}



/* 站间单差，历元间差分计算坐标增量(单历元—法方程叠加) 2013.12.09 */
static int dt2dx_Mod(rtk_t *rtk, const nav_t *nav, double dt,
	const double *x,const int *sat,const int ns, double *y,int n,double *e,double *azel,
	const int *iu, const int *ir, double *dx, double *px)
{
	prcopt_t *opt=&rtk->opt;
	int i,j,k,m,f,sysi,sysj,nsat[2]={0};
	double lami,lamj,bl,dr[3];
	double *H,*Hi,*v,*R,*Ri,*Rj,*Q,*w,*ATP,*dxt;
	int nx=5,nv=0;
	double m2,*vtp;           //验后单位权中误差
	int mask[4]={0};
	double m0=1,eps=2;    //粗差探测中使用
	double *P,*N,*V;
	int *W;
	int slip;    //slip:周跳数
	int vflg[MAXOBS*NFREQ*2+1],sati,freq,type;   //记录站际历元间二次差分观测方程具体信息
	char *stype;
	double *Ws;   //Co_bust中使用
	int NX=4+3;
	double vv,ratio=1.0,ep=0.02;
	int  flg=0;
	int info=0; //粗差探测是否成功标志
	int suc_info=1;
	double rms_v=0.0,vpv=0.0;
	int *slp;
	int info_LEGE=1;
	double *vc;
	double w_ratio=1.0;

	static double y0[MAXOBS*2],dx0[3],azel0[2*MAXSAT],bl0,dt0;
	static int   sat0[MAXSAT],ns0,n0,iu0[MAXSAT],ir0[MAXSAT];

	//抗差估计
	//fact:权因子函数选择：0---IGG1;1---IGG2;2---Huber;
	//k0:等价权函数保权界值；k1:等价权函数零权界值
	//eps: 限差；
	double fact=0,eps_x=0.001,k0=3,k1=8,*R1; 

	trace(3,"dt2dx_:%s\n",time_str(rtk->sol.time,3));

	for(i=0;i<MAXSLIP;i++) now_slipflg[i]=0;   //当前历元周跳标记数组清0

	bl=baseline(x,rtk->rb,dr);

	if(bl0==0){  //初始历元
		matcpy(y0,y,2*opt->nf,n);        //站星距离
		matcpy(azel0,azel,2,n);          //卫星高度角
		memcpy(sat0,sat,MAXSAT);         //共视星序列
		memcpy(iu0,iu,MAXSAT);           //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);           //基准共视星序列
		bl0=bl;                          //基线长
		dt0=dt;
		ns0=ns;                          //共视卫星数
		n0=n;
		return 0;
	}

	//trace(2,"y0:\n");   tracemat(2,   y0,2*opt->nf,n0,12,3);
	//for(i=0;i<ns0;i++) 
	//{
	//	trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat0[i],iu0[i],ir0[i],azel0[1+iu0[i]*2]*R2D);
	//}

	//trace(2,"y:\n");   tracemat(2,   y,2*opt->nf,n,12,3);
	//for(i=0;i<ns;i++) {
	//	trace(2,"sat=%2d %3d %3d azel=%5.1f\n ",sat[i],iu[i],ir[i],azel[1+iu[i]*2]*R2D);
	//}

	H=zeros(nx,opt->nf*ns);v=mat(opt->nf*ns,1);
	R=mat(opt->nf*ns,opt->nf*ns);
	Ri=mat(opt->nf*ns,1);Rj=mat(opt->nf*ns0,1);
	/* for each system (0:gps/qzss/sbas,1:glonass) */
	for(m=0;m<3;m++) for(f=0;f<opt->nf;f++) 
	{  //只是用了载波相位观测量
		for(i=0;i<ns;i++) 
		{ 
			sysi=rtk->ssat[sat[i]-1].sys;
			if (!test_sys(sysi,m)) continue;
			if(!validobs(iu[i],ir[i],f,opt->nf,y))  continue;

			for(j=0;j<ns0;j++) {
				sysj=rtk->ssat[sat0[j]-1].sys;
				if (!test_sys(sysi,m)) continue;
				if(sat[i]!=sat0[j])  continue;  //历元间差分，卫星匹配
				if(!validobs(iu0[j],ir0[j],f,opt->nf,y0))   continue;

				lami=nav->lam[sat[i]-1][f];
				lamj=nav->lam[sat0[j]-1][f];
				if(lami<0.0||lamj<0.0) continue;

				Hi=H+nv*nx;

				v[nv]=(y[f+iu[i]*2*opt->nf]-y[f+ir[i]*2*opt->nf])-(y0[f+iu0[j]*2*opt->nf]-y0[f+ir0[j]*2*opt->nf]);


				//for(k=0;k<3;k++) {   //每个频率上设一个相对钟差参数
				//	Hi[k]=-e[k+iu[i]*3];
				//}  

				for (k=0;k<NX;k++) Hi[k]=k<3?-e[k+iu[i]*3]:(k==3?1.0:0.0);

				/* time system and receiver bias offset correction */
				//时间系统和接收机的偏移校正
				//Hi[3]:GPS; Hi[4]:GLO; Hi[5]:GAL; Hi[6]:CMP.
				if      (sysj==SYS_GLO) { Hi[4]=1.0; mask[1]=1;}
				else if (sysj==SYS_GAL) { Hi[5]=1.0; mask[2]=1;}
				else if (sysj==SYS_CMP) { Hi[6]=1.0; mask[3]=1;}
				else mask[0]=1;

				Ri[nv]=varerr(sysi,azel[1+iu[i]*2],bl,dt,f,opt);
				Rj[nv]=varerr(sysj,azel0[1+iu0[j]*2],bl0,dt0,f,opt);				

				trace(4,"sat=%2d  v=%7.3f\n",sat[i],v[nv]);

				//vflg数据记录站际历元间二次差分观测方程具体信息
				//其结构依次为：卫星，观测量类型，频率
				vflg[nv++]=((sat[i]<<8)|((f<opt->nf?0:1)<<4)|(f%opt->nf));
				break;
			}
		}
	}

	/* shrink design matrix:nx=5->4 */
	//if(!nsat[0]||!nsat[1]) {
	//	nx=4;
	//	for(i=0;i<nv;i++) {
	//		for(j=0;j<4;j++) {
	//			H[j+i*4]=H[j+i*5];				
	//		}
	//	}
	//} 

	//Modified by XGAO,2014/10/15
	int sysnum=0;
	for(i=0;i<4;i++) sysnum+=mask[i];
	//Single System



	trace(2,"H= %d*%d\n",nx,nv); tracemat(2,H,nx,nv,10,6);
	trace(2,"v=\n"); tracemat(2,v,nv,1,10,4);

	if(nv<nx) 
	{
		//trace(2,"历元间差分方程个数nv=%d,参数个数nx=%d,观测方程不够，无法解算！\n",nv,nx);
		free(H);free(v);free(R);free(Ri);free(Rj);

		matcpy(y0,y,2*opt->nf,n);           //站星距离
		matcpy(azel0,azel,2,n);     //卫星高度角
		memcpy(sat0,sat,MAXSAT);     //共视星序列
		memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);       //基准共视星序列
		bl0=bl;                      //基线长
		dt0=dt;
		ns0=ns;                       //共视卫星数
		n0=n;
		return 0;
	}

	dxt=mat(nx,1);Q=zeros(nx,nx);w=mat(nv,1);ATP=zeros(nx,nv);vtp=mat(1,nv);

	/* 历元间差分观测方程权阵 */
	dtcov(nv,Ri,Rj,R);

	//trace(2,"R= %d阶\n",nv); tracemat(2,R,nv,nv,10,6);

	P=mat(nv,1);N=mat(nx*(nx+1)/2,1);V=mat(nv,1);
	W=imat(nv,1);Ws=mat(nv,1);vc=mat(nv,1);

	//多普勒积分法探测到的周跳标记，周跳数>1才标记(有时可能出现虚报的情况)
	for(i=0;i<nv;i++) {           
		W[i]=1;
	}

	if(!matinv(R,nv)) 
	{

		//先验方差给的过大，后面方差一致性检验无法通过，给予7倍权尝试，可以通过		
		for(i=0;i<nv;i++) P[i]=R[i+i*nv]*m0*m0;  

		for(i=0;i<nv;i++) R[i+i*nv]=R[i+i*nv]*m0*m0; 

		R1=mat((nv+1)*nv/2,1);
		for(i=0;i<nv;i++) 
			for(j=0;j<=i;j++) R1[ij(i,j)]=R[i+j*nv];   //取下三角矩阵

		while(1) {

			//info=0为初次平差没有发现粗差；1,2,3：依次表示粗差探测第1，第2，第3种情况
			trace(2,"v:\n");tracemat(2,v,1,nv,10,3);
			//*m2=sqrt(pvv/(nn-t))
			info=DataSnooping(nv,nx,H,v,P,eps,m0,dxt,N,V,W,&m2);	//粗差探测
			if(info>0) {  //探测到粗差
				trace(1,"\nTime:%s\n",time_str(rtk->sol.time,3));
				trace(1,"Detect gloss by DataSnooping,Situation %d\n",info);
				trace(2,"粗差探测结果：W=");tracemati(2,W,1,nv,4);			
			}

			info_LEGE=1;  //粗差探测2默认设置为1
			if(info>1) {   //情况2,3
				flg=0;
				for(i=0;i<nv;i++) {                 /* 记录周跳发生的状况 */
					if(W[i]==0) {							
						sati=(vflg[i]>> 8)&0xFF;  //卫星
						type=(vflg[i]>> 4)&0xF;   //观测量类型
						freq=vflg[i]&0xF;         //频率
						stype=type==0?"L":(type==1?"L":"C");
						lami=nav->lam[sati-1][freq];
						slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

						trace(1,"Slip detected by DataSnooping:time:%s\n",time_str(rtk->sol.time,3));
						trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d,周跳修复后残差V:%6.3lf\n",
							sati,stype,freq+1,V[i],slip,V[i]-lami*slip);

						vv=V[i];
						vv-=slip*lami;
						if(fabs(vv)>0.02)  {
							flg++;
							trace(1,"残差超限，DataSnooping探测周跳定位可能错误，周跳修复后残差V:%6.3lf\n",vv);
							continue;
						}
					}
				}

				if(flg>0) { //探测粗差

					if(info==3&&flg>0) {             //粗差探测第3种情况
						for(i=0;i<nv;i++) W[i]=1;    //标记卫星清零
					}
				}
				//如果有效观测量只有5,采用搜索法确定粗差(这种情况可能已经处理)
				if(nv==(nx+1)) for(i=0;i<nv;i++) W[i]=1;   

				trace(1,"Detect gloss by DataSnooping_single:\n");
				DataSnooping_single(nv,nx,H,v,P,m0,vflg,nav,dxt,N,V,W,&m2); 

				if(info_LEGE>0) {

					matcpy(vc,v,nv,1);
					for(i=0,flg=0;i<nv;i++) {         /* 记录周跳发生的状况 */
						if(W[i]==0) {							
							sati=(vflg[i]>> 8)&0xFF;  //卫星
							type=(vflg[i]>> 4)&0xF;   //观测量类型
							freq=vflg[i]&0xF;         //频率
							stype=type==0?"L":(type==1?"L":"C");
							lami=nav->lam[sati-1][freq];
							slip=QFH(V[i])*SSWR(fabs(V[i])/lami);  //4舍5入取整

							trace(2,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d,周跳修复后残差V:%6.3lf\n",
								sati,stype,freq+1,V[i],slip,V[i]-lami*slip);
							vc[i]-=slip*lami;
						}
					}
					m2=CoRobust(nv,nx,H,vc,R1,fact,k0,k1,m0,eps_x,dxt,N,V,Ws); //抗差估计
					trace(2,"抗差估计 by CoRobust,W：");tracemat(2,Ws,1,nv,10,3);
					for(i=0;i<nv;i++) {
						V[i]+=v[i]-vc[i];
						if(Ws[i]==0.0) W[i]=0;
						else W[i]=1;
						if(fabs(v[i]-vc[i])>0.01) W[i]=0;//含有周跳
					}

				}
				else{ 
					m2=-1;
				}
				for(i=0;i<3;i++) dx[i]=dxt[i];	
			}

			for(i=0;i<nv;i++) {         /* 当前历元周跳探测情况，周跳标记为1 */
				sati=(vflg[i]>> 8)&0xFF;  //卫星
				freq=vflg[i]&0xF;         //频率					
				f=freq;
				if(W[i]==0)	now_slipflg[IB(sati,f,&rtk->opt)]=1;  //当前历元发生周跳，标记为1
				else        now_slipflg[IB(sati,f,&rtk->opt)]=2;  //未发生周跳，标记为2	
				trace(2,"sat=%3d,now_slipflg=%3d\n",sati,now_slipflg[IB(sati,f,&rtk->opt)]);
			}
			/* now_slipflg数据标记：-------------------------------------------
			0：  未参与站际历元二次差分模型的卫星，是否发生周跳未知；
			或者站际历元二次差分模型计算失败（有效卫星数不够），
			周跳情况未知；
			1：  表示站际历元二次差分模型探测出的周跳；
			2：  通过站际历元二次差分模型探测未发生周跳或者发生周跳已经修复。
			-----------------------------------20140731-------------------------*/			
			for(i=0,k=0;i<nv;i++) {
				if(W[i]==0) 	k++;   //可能发生周跳的个数
			}
			//Modified by CYJ 2014/10/15
			//如果有效观测为参数个数，即nv=nx时，Datasnooping会探测失败，本文周跳探测模型无效
			//周跳标记应该全部为0，即周跳是否发生无从而知，考虑到rtklib针对单频数据没有进行周跳
			//探测，建议将其标记为1.若以后补充了其他模型可以考虑将其标记为0(即是否周跳未知)
			if(nv==nx) {  
				for(i=0;i<nv;i++){
					now_slipflg[IB(sati,f,&rtk->opt)]=1;
					slips[IB(sati,f,&rtk->opt)]=0;				
				}				
			}

			if(k>0) {			
				slp=imat(nv,1);for(i=0;i<nv;i++) slp[i]=0;
				//周跳搜索计算
				rms_v=10;  //默认设置
				ratio=1;
				if(info_LEGE>0) {  //粗差探测成功才搜索，否则没必要搜索计算
					seach_slip(nv,nx,W,vflg,H,R,v,V,nav,&ratio,&rms_v,&vpv,slp,&w_ratio);
				}

				trace(2,"w_ratio:%10.3lf\n",w_ratio);
				trace(2,"slp:\n");tracemati(2,slp,1,nv,4);
				if((rms_v<ep&&ratio>30)||(w_ratio>3.0)) 
				{

					for(i=0;i<nv;i++){
						sati=(vflg[i]>> 8)&0xFF;  //卫星
						type=(vflg[i]>> 4)&0xF;   //观测量类型
						freq=vflg[i]&0xF;         //频率
						stype=type==0?"L":(type==1?"L":"C");
						lami=nav->lam[sati-1][freq];
						slip=slp[i];  //4舍5入取整
						v[i]-=slip*lami;      //用于再进行一次平差求得更为精确的坐标差信息
						if(slip==0) k--;
						if(slip!=0) {
							trace(1,"Slip repaired success:time:%s\n",time_str(rtk->sol.time,3));
							trace(1,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d   修复后残差V:%6.3lf\n",
								sati,stype,freq+1,V[i],slip,V[i]-slip*lami);

							//fprintf(fp_slp,"Slip repaired success:time:%s\n",time_str(rtk->sol.time,3));
							//fprintf(fp_slp,"发生周跳卫星sat:%d %s%d  V:%6.3lf  Slips:%d   修复后残差V:%6.3lf\n",
							//	sati,stype,freq+1,V[i],slip,V[i]-slip*lami);
						}

						//将周跳数记录到slips数组中
						f=freq;
						slips[IB(sati,f,&rtk->opt)]+=slip;	

						//当前历元周跳修复，修复标记为2
						now_slipflg[IB(sati,f,&rtk->opt)]=2;					
						W[i]=1;	

					}
				}
				else {
					trace(1,"Slip repaired failed:time:%s\n",time_str(rtk->sol.time,3));
					//fprintf(fp_slp,"Slip repaired failed:time:%s\n",time_str(rtk->sol.time,3));
					suc_info=0;   //周跳探测修复失败，需跳出循环
					for(i=0;i<nv;i++){
						sati=(vflg[i]>> 8)&0xFF;  //卫星
						type=(vflg[i]>> 4)&0xF;   //观测量类型
						freq=vflg[i]&0xF;         //频率
						f=freq;
						if(info==1) { //对于粗差第一种情况，只对周跳卫星重新初始化
							if(!W[i])//周跳探测失败，存在多处周跳
							{
								now_slipflg[IB(sati,f,&rtk->opt)]=1;
								slips[IB(sati,f,&rtk->opt)]=0; 
								//fprintf(fp_slp,"需要重新初始化卫星模糊度参数：%3d\n",sati);
							}						
						}
						else {     //对于粗差第二和第三种情况，对全部卫星重新初始化
							now_slipflg[IB(sati,f,&rtk->opt)]=1;
							slips[IB(sati,f,&rtk->opt)]=0; 
							//fprintf(fp_slp,"需要重新初始化卫星模糊度参数：%3d\n",sati);
						}
					}
				}
				free(slp);
			}

			if(k==0||suc_info==0) {free(R1);break;}  //当不存在周跳后或者修复失败，结束计算
		}

		if(m2<0) 
		{  
			printf("计算失败\n");getchar();
			matcpy(y0,y,2*opt->nf,n);           //站星距离
			matcpy(azel0,azel,2,n);             //卫星高度角
			memcpy(sat0,sat,MAXSAT);            //共视星序列
			memcpy(iu0,iu,MAXSAT);              //流动站共视星序列
			memcpy(ir0,ir,MAXSAT);              //基准共视星序列
			bl0=bl;                             //基线长
			dt0=dt;
			ns0=ns;                             //共视卫星数
			n0=n;
			free(H);free(v);free(R);free(Ri);free(Rj);
			free(Q);free(w);free(ATP);free(dxt);free(vtp);
			free(P);free(V);free(W);free(N);free(Ws);free(vc);
			return 0;
		}

		for(i=0;i<3;i++) dx[i]=dxt[i];		

		m0=m0*m0;   //使用验前单位权中误差		
		px[0]=m0*N[ij(0,0)];px[1]=m0*N[ij(1,1)];px[2]=m0*N[ij(2,2)];
		px[3]=m0*N[ij(1,0)];px[4]=m0*N[ij(2,0)];px[5]=m0*N[ij(2,1)];  

		trace(2,"Dx:\n");    tracemat(2,px,1,6,9,7);
		trace(2,"钟差参数dt:%10.3lf\n",dxt[3]);    
	}
	else {
		matcpy(y0,y,2*opt->nf,n);           //站星距离
		matcpy(azel0,azel,2,n);     //卫星高度角
		memcpy(sat0,sat,MAXSAT);     //共视星序列
		memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
		memcpy(ir0,ir,MAXSAT);       //基准共视星序列
		bl0=bl;                      //基线长
		dt0=dt;
		ns0=ns;                       //共视卫星数
		n0=n;  
		free(H);free(v);free(R);free(Ri);free(Rj);
		free(Q);free(w);free(ATP);free(dxt);free(vtp);
		free(P);free(V);free(W);free(N);free(Ws);free(vc);
		return 0;
	}


	/* 历元间差分计算成功，记录当前历元信息 */
	matcpy(y0,y,2*opt->nf,n);           //站星距离
	matcpy(azel0,azel,2,n);     //卫星高度角
	memcpy(sat0,sat,MAXSAT);     //共视星序列
	memcpy(iu0,iu,MAXSAT);       //流动站共视星序列
	memcpy(ir0,ir,MAXSAT);       //基准共视星序列
	bl0=bl;                      //基线长
	dt0=dt;
	ns0=ns;                       //共视卫星数
	n0=n;

	free(H);free(v);free(R);free(Ri);free(Rj);
	free(Q);free(w);free(ATP);free(dxt);free(vtp);
	free(P);free(V);free(W);free(N);free(Ws);free(vc);

	return 1;
}

/* 站间单差，历元间差分计算坐标增量（法方程叠加） 2013.12.09 */
extern  int dt2dx_one(rtk_t *rtk, const nav_t *nav, double dt,
					 const double *x,const int *sat, const int ns,double *y, int n, double *e,
					 double *azel, const int *iu, const int *ir, double *dx0,double *Qx)
{
	int i;
	int info;
	double dx[3],px[6],dxt[3];
	static double sx[3];
	static double X0[3];
    /*FILE *fp;
	char *file="E:\\论文4：单频周跳探测与修复\\detectAndRepairSlpByDt\\RESULT\\sx.txt";
    fp=fopen(file,"a+");*/
	trace(3,"dt2dx_one:time=%s\n",time_str(rtk->sol.time,3));
	
	if(norm(X0,3)<=0) {    //X0 记录上一历元坐标初值
		for(i=0;i<3;i++) X0[i]=x[i];
	}
	//x用于计算基线长，dt用于计算观测方差
	info=dt2dx_(rtk,nav,dt,x,sat,ns,y,n,e,azel,iu,ir,dx,px);
	
	if(!info) {   //历元间差分计算失败
		// 更新X0到当前历元
		for(i=0;i<3;i++)  X0[i]=x[i];
		return 0;
	}  
	
	//更新dx0
	for(i=0;i<3;i++) dx0[i]=dx[i];
	
	//（静态模拟，反映历元间差分计算精度情况）
	for(i=0;i<3;i++)  dxt[i]=x[i]-X0[i]+dx[i];    //计算相邻历元坐标差值
	for(i=0;i<3;i++)  sx[i]+=dxt[i];              //计算坐标差值的和
	// 更新X0到当前历元
	for(i=0;i<3;i++)  X0[i]=x[i];
	
	//fprintf(fp,"%s %.4lf  %.4lf %.4lf  %.4lf  %.4lf  %.4lf\n",
	//	time_str(rtk->sol.time,3),dxt[0],dxt[1],dxt[2],sx[0],sx[1],sx[2]);
	//fclose(fp);
	
	matcpy(Qx,px,6,1);
	
	trace(2,"dx:");tracemat(2,dx,1,3,10,4);
	trace(2,"Qx:");tracemat(2,Qx,1,3,10,4);
	
	return 1;
}
