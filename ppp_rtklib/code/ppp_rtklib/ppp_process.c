#include "rtklib.h"

/*存储观测的数据*/
static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */
static obs_t obss = { 0 };          /* observation data */
static nav_t navs = { 0 };          /* navigation data */
static sta_t stas[MAXRCV];      /* station infomation */
static int iobsu = 0;            /* current rover observation data index */
static int iobsr = 0;            /* current reference observation data index */
static int revs = 0;            /* analysis direction (0:forward,1:backward) */
static int aborts = 0;            /* abort status */
static int nepoch = 0;            /* number of observation epochs */


prcopt_t popt_ = { 0 };

const prcopt_t prcopt_default = { /* defaults processing options */
	PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
	15.0*D2R,{ { 0,0 } },           /* elmin,snrmask */
	0,1,1,5,0,10,               /* sateph,modear,glomodear,maxout,minlock,minfix */
	0,0,0,0,                    /* estion,esttrop,dynamics,tidecorr */
	1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
	0,0,                        /* rovpos,refpos */
	{ 100.0,100.0 },              /* eratio[] */
	{ 100.0,0.003,0.003,0.0,1.0 }, /* err[] */
	{ 30.0,0.03,0.3 },            /* std[] */
	{ 1E-4,1E-3,1E-4,1E-1,1E-2 }, /* prn[] */
	5E-12,                      /* sclkstab */
	{ 3.0,0.9999,0.20 },          /* thresar */
	0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
	30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
	{ 0 },{ 0 },{ 0 },                /* baseline,ru,rb */
	{ "","" },                    /* anttype */
	{ { 0 } },{ { 0 } },{ 0 }             /* antdel,pcv,exsats */
};
const solopt_t solopt_default = { /* defaults solution output options */
	SOLF_LLH,TIMES_GPST,1,3,    /* posf,times,timef,timeu */
	0,1,0,0,0,0,                /* degf,outhead,outopt,datum,height,geoid */
	0,0,0,                      /* solstatic,sstat,trace */
	{ 0.0,0.0 },                  /* nmeaintv */
	" ",""                      /* separator/program name */
};

/*pntpos里面用到了*/
const double chisqr[100] = {      /* chi-sqr(n) (alpha=0.001) */
	10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
	31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
	46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
	61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
	74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
	88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
	101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
	113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
	126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
	138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};

/*功能函数--------------------------------------------------------------------*/
/*输出文件头函数*/
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
	double pos[3], dms1[3], dms2[3];
	const char *sep = opt->sep;

	trace(3, "outrpos :\n");

	if (opt->posf == SOLF_LLH || opt->posf == SOLF_ENU) {
		ecef2pos(r, pos);
		if (opt->degf) {
			deg2dms(pos[0] * R2D, dms1);
			deg2dms(pos[1] * R2D, dms2);
			fprintf(fp, "%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
				dms1[0], sep, dms1[1], sep, dms1[2], sep, dms2[0], sep, dms2[1],
				sep, dms2[2], sep, pos[2]);
		}
		else {
			fprintf(fp, "%13.9f%s%14.9f%s%10.4f", pos[0] * R2D, sep, pos[1] * R2D,
				sep, pos[2]);
		}
	}
	else if (opt->posf == SOLF_XYZ) {
		fprintf(fp, "%14.4f%s%14.4f%s%14.4f", r[0], sep, r[1], sep, r[2]);
	}
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
	const solopt_t *sopt)
{
	const char *s1[] = { "GPST","UTC","JST" };
	gtime_t ts, te;
	double t1, t2;
	int i, j, w1, w2;
	char s2[32], s3[32];

	trace(3, "outheader: n=%d\n", n);

	if (sopt->posf == SOLF_NMEA) return;

	if (sopt->outhead) {
		if (!*sopt->prog) {
			fprintf(fp, "%s program   : RTKLIB ver.%s\n", COMMENTH, VER_RTKLIB);
		}
		else {
			fprintf(fp, "%s program   : %s\n", COMMENTH, sopt->prog);
		}
		for (i = 0; i<n; i++) {
			fprintf(fp, "%s inp file  : %s\n", COMMENTH, file[i]);
		}
		for (i = 0; i<obss.n; i++)    if (obss.data[i].rcv == 1) break;
		for (j = obss.n - 1; j >= 0; j--) if (obss.data[j].rcv == 1) break;
		if (j<i) { fprintf(fp, "\n%s no rover obs data\n", COMMENTH); return; }
		ts = obss.data[i].time;
		te = obss.data[j].time;
		t1 = time2gpst(ts, &w1);
		t2 = time2gpst(te, &w2);
		if (sopt->times >= 1) ts = gpst2utc(ts);
		if (sopt->times >= 1) te = gpst2utc(te);
		if (sopt->times == 2) ts = timeadd(ts, 9 * 3600.0);
		if (sopt->times == 2) te = timeadd(te, 9 * 3600.0);
		time2str(ts, s2, 1);
		time2str(te, s3, 1);
		fprintf(fp, "%s obs start : %s %s (week%04d %8.1fs)\n", COMMENTH, s2, s1[sopt->times], w1, t1);
		fprintf(fp, "%s obs end   : %s %s (week%04d %8.1fs)\n", COMMENTH, s3, s1[sopt->times], w2, t2);
	}
	if (sopt->outopt) {
		outprcopt(fp, popt);
	}
	if (PMODE_DGPS <= popt->mode&&popt->mode <= PMODE_FIXED&&popt->mode != PMODE_MOVEB) {
		fprintf(fp, "%s ref pos   :", COMMENTH);
		outrpos(fp, popt->rb, sopt);
		fprintf(fp, "\n");
	}
	if (sopt->outhead || sopt->outopt) fprintf(fp, "%s\n", COMMENTH);

	outsolhead(fp, sopt);
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, char **infile, int n,
	const prcopt_t *popt, const solopt_t *sopt)
{
	FILE *fp = stdout;

	trace(3, "outhead: outfile=%s n=%d\n", outfile, n);

	if (*outfile) {
		createdir(outfile);

		if (!(fp = fopen(outfile, "w"))) {
			showmsg("error : open output file %s", outfile);
			return 0;
		}
	}
	/* output header */
	outheader(fp, infile, n, popt, sopt);

	if (*outfile) fclose(fp);

	return 1;
}
/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
	va_list arg;
	char buff[1024], *p = buff;
	if (!*format) return showmsg("");
	va_start(arg, format);
	p += vsprintf(p, format, arg);
	va_end(arg);
	return showmsg(buff);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
	double tt;
	int n;

	for (; *i<obs->n; (*i)++) if (obs->data[*i].rcv == rcv) break;
	for (n = 0; *i + n<obs->n; n++) {
		tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
		if (obs->data[*i + n].rcv != rcv || tt>DTTOL) break;
	}
	return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
	double tt;
	int n;

	for (; *i >= 0; (*i)--) if (obs->data[*i].rcv == rcv) break;
	for (n = 0; *i - n >= 0; n++) {
		tt = timediff(obs->data[*i - n].time, obs->data[*i].time);
		if (obs->data[*i - n].rcv != rcv || tt<-DTTOL) break;
	}
	return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
	gtime_t time = { 0 };
	char path[1024];
	int i, nu, nr, n = 0;

	trace(3, "infunc  : revs=%d iobsu=%d iobsr=%d \n", revs, iobsu, iobsr);

	if (0 <= iobsu&&iobsu<obss.n) {
		settime((time = obss.data[iobsu].time));
		if (checkbrk("processing : %s Q=%d", time_str(time, 0), solq)) {
			aborts = 1; showmsg("aborted"); return -1;
		}
	}
	if (!revs) { /* input forward data */
		if ((nu = nextobsf(&obss, &iobsu, 1)) <= 0) return -1;
		if (popt->intpref) {
			for (; (nr = nextobsf(&obss, &iobsr, 2))>0; iobsr += nr)
				if (timediff(obss.data[iobsr].time, obss.data[iobsu].time)>-DTTOL) break;
		}
		else {
			for (i = iobsr; (nr = nextobsf(&obss, &i, 2))>0; iobsr = i, i += nr)
				if (timediff(obss.data[i].time, obss.data[iobsu].time)>DTTOL) break;
		}
		nr = nextobsf(&obss, &iobsr, 2);
		for (i = 0; i<nu&&n<MAXOBS; i++) obs[n++] = obss.data[iobsu + i];
		for (i = 0; i<nr&&n<MAXOBS; i++) obs[n++] = obss.data[iobsr + i];
		iobsu += nu;

		
	}
	return n;
}
/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
	int mode)
{
	gtime_t time = { 0 };
	sol_t sol = { { 0 } };
	rtk_t rtk;
	obsd_t obs[MAXOBS];
	double rb[3] = { 0 };
	int i, nobs, n, solstatic, pri[] = { 0,1,2,3,4,5,1,6 };

	trace(3, "procpos : mode=%d\n", mode);

	solstatic = sopt->solstatic &&
		(popt->mode == PMODE_STATIC || popt->mode == PMODE_PPP_STATIC);

	rtkinit(&rtk, popt);

	while ((nobs = inputobs(obs, rtk.sol.stat, popt)) >= 0) {

		/* exclude satellites */
		for (i = n = 0; i<nobs; i++) {
			if ((satsys(obs[i].sat, NULL)&popt->navsys) &&
				popt->exsats[obs[i].sat - 1] != 1) obs[n++] = obs[i];
		}
		if (n <= 0) continue;

		if (!rtkpos(&rtk, obs, n, &navs)) continue;

		if (mode == 0) { /* forward/backward */
			if (!solstatic) {
				outsol(fp, &rtk.sol, rtk.rb, sopt);
			}
			else if (time.time == 0 || pri[rtk.sol.stat] <= pri[sol.stat]) {
				sol = rtk.sol;
				for (i = 0; i<3; i++) rb[i] = rtk.rb[i];
				if (time.time == 0 || timediff(rtk.sol.time, time)<0.0) {
					time = rtk.sol.time;
				}
			}
		}
	}
	if (mode == 0 && solstatic&&time.time != 0.0) {
		sol.time = time;
		outsol(fp, &sol, rb, sopt);
	}
	rtkfree(&rtk);
}


/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
	trace(3, "openfile: outfile=%s\n", outfile);

	return !*outfile ? stdout : fopen(outfile, "a");
}
/*-----------------------------------------------------------------------------*/


/*模块一：输入数据*/
int importData(const filopt_t *fopt, const prcopt_t *popt,pcvs_t *pcvs, pcvs_t *pcvr, 
	obs_t *obs, nav_t *nav,sta_t *sta)
{
	gtime_t ts = { 0 }, te = { 0 };
	int ret;
	popt_ = *popt;

	trace(3, "Import Data");

	/* read satellite antenna parameters */
	if (*fopt->satantp&&!(readpcv(fopt->satantp, pcvs))) {
		showmsg("error : no sat ant pcv in %s", fopt->satantp);
		trace(1, "sat antenna pcv read error: %s\n", fopt->satantp);
		return 0;
	}
	/* read receiver antenna parameters */
	if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp, pcvr))) {
		showmsg("error : no rec ant pcv in %s", fopt->rcvantp);
		trace(1, "rec antenna pcv read error: %s\n", fopt->rcvantp);
		return 0;
	}

	/*读取sp3精密星历*/
	readsp3(fopt->sp3, nav, 0);

	/*读取钟差文件*/
	readrnxc(fopt->clk, nav);

	/*读取观测值文件*/
	obs->data = NULL; obs->n = obs->nmax = 0;
	ret = readrnxt(fopt->obs, 1, ts, te, 0.0, (&popt_)->rnxopt[0], obs,nav,sta);

	/*读取导航文件*/
	nav->eph = NULL; nav->n = nav->nmax = 0;
	nav->geph = NULL; nav->ng = nav->ngmax = 0;
	nav->seph = NULL; nav->ns = nav->nsmax = 0;
	nepoch = 0;
	ret = readrnxt(fopt->nav, 2, ts, te, 0.0, (&popt_)->rnxopt[1], obs, nav, sta+1);
	if (obs->n <= 0) {
		checkbrk("error : no obs data");
		trace(1, "no obs data\n");
		return 0;
	}
	if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
		checkbrk("error : no nav data");
		trace(1, "no nav data\n");
		return 0;
	}
	/* sort observation data */
	nepoch = sortobs(obs);

	/* delete duplicated ephemeris */
	uniqnav(nav);

}
/*模块二：预处理及处理入口procpos*/
int process(const filopt_t *fopt, const prcopt_t *popt, pcvs_t *pcvs, pcvs_t *pcvr,
	obs_t *obs, nav_t *nav, sta_t *sta, char *outfile, const solopt_t *sopt)
{
	FILE *fp;

	/* set antenna paramters */
	if (popt_.sateph == EPHOPT_PREC || popt_.sateph == EPHOPT_SSRCOM) {
		setpcv(obss.n>0 ? obss.data[0].time : timeget(), &popt_, &navs, &pcvss, &pcvsr,
			stas);
	}
	char* infile[] = { fopt->obs ,fopt->nav,fopt->sp3,fopt->clk };

	/* write header to output file */
	outhead(outfile, infile, 4, &popt_, sopt);


	iobsu = iobsr = revs = aborts = 0;
	if (popt_.mode == PMODE_SINGLE || popt_.soltype == 0) {
		if ((fp = openfile(outfile))) {
			procpos(fp, &popt_, sopt, 0); /* forward */
			fclose(fp);
		}
	}

}
/*ppp处理过程：ppp_process*/
int ppp_process(const prcopt_t *popt, const solopt_t *sopt,const filopt_t *fopt,char *outfile)
{
	popt_ = *popt;

	if (!importData(fopt,popt, &pcvss, &pcvsr, &obss, &navs, &stas))
	{
		showmsg("导入数据错误!\n");
	}
	printf("---------------input file success!---------------\n");

	if (!process(fopt, popt, &pcvss, &pcvsr, &obss, &navs, &stas,outfile,sopt))
	{
		showmsg("预处理错误!\n");
	}
	printf("\n---------------数据处理成功!---------------\n");

	/*释放数据占用的空间*/
	freeData(&pcvss, &pcvsr, &navs, &obss);
}