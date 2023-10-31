#include "rtklib.h"

/* global variables ----------------------------------------------------------*/
static int statlevel = 0;          /* rtk status output level (0:off) */
static FILE *fp_stat = NULL;       /* rtk status file pointer */


/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */



/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
	char buff[256], tstr[32];
	int n;
	va_list ap;
	time2str(rtk->sol.time, tstr, 2);
	n = sprintf(buff, "%s: ", tstr + 11);
	va_start(ap, format);
	n += vsprintf(buff + n, format, ap);
	va_end(ap);
	n = n<MAXERRMSG - rtk->neb ? n : MAXERRMSG - rtk->neb;
	memcpy(rtk->errbuf + rtk->neb, buff, n);
	rtk->neb += n;
	trace(2, "%s", buff);
}
								   /* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t *rtk)
{
	ssat_t *ssat;
	double tow, pos[3], vel[3], acc[3], vela[3] = { 0 }, acca[3] = { 0 }, xa[3];
	int i, j, week, est, nfreq, nf = NF(&rtk->opt);
	char id[32];

	if (statlevel <= 0 || !fp_stat) return;

	trace(3, "outsolstat:\n");

	est = rtk->opt.mode >= PMODE_DGPS;
	nfreq = est ? nf : 1;
	tow = time2gpst(rtk->sol.time, &week);

	/* receiver position */
	if (est) {
		for (i = 0; i<3; i++) xa[i] = i<rtk->na ? rtk->xa[i] : 0.0;
		fprintf(fp_stat, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
			rtk->sol.stat, rtk->x[0], rtk->x[1], rtk->x[2], xa[0], xa[1], xa[2]);
	}
	else {
		fprintf(fp_stat, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
			rtk->sol.stat, rtk->sol.rr[0], rtk->sol.rr[1], rtk->sol.rr[2],
			0.0, 0.0, 0.0);
	}
	/* receiver velocity and acceleration */
	if (est&&rtk->opt.dynamics) {
		ecef2pos(rtk->sol.rr, pos);
		ecef2enu(pos, rtk->x + 3, vel);
		ecef2enu(pos, rtk->x + 6, acc);
		if (rtk->na >= 6) ecef2enu(pos, rtk->xa + 3, vela);
		if (rtk->na >= 9) ecef2enu(pos, rtk->xa + 6, acca);
		fprintf(fp_stat, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
			week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], acc[0], acc[1], acc[2],
			vela[0], vela[1], vela[2], acca[0], acca[1], acca[2]);
	}
	else {
		ecef2pos(rtk->sol.rr, pos);
		ecef2enu(pos, rtk->sol.rr + 3, vel);
		fprintf(fp_stat, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
			week, tow, rtk->sol.stat, vel[0], vel[1], vel[2],
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}
	/* receiver clocks */
	fprintf(fp_stat, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
		week, tow, rtk->sol.stat, 1, rtk->sol.dtr[0] * 1E9,
		(rtk->sol.dtr[0] + rtk->sol.dtr[1])*1E9, 0.0, 0.0);

	/* ionospheric parameters */
	if (est&&rtk->opt.ionoopt == IONOOPT_EST) {
		for (i = 0; i<MAXSAT; i++) {
			ssat = rtk->ssat + i;
			if (!ssat->vs) continue;
			satno2id(i + 1, id);
			j = II(i + 1, &rtk->opt);
			xa[0] = j<rtk->na ? rtk->xa[j] : 0.0;
			fprintf(fp_stat, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow, rtk->sol.stat,
				id, ssat->azel[0] * R2D, ssat->azel[1] * R2D, rtk->x[j], xa[0]);
		}
	}
	/* tropospheric parameters */
	if (est && (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG)) {
		for (i = 0; i<2; i++) {
			j = IT(i, &rtk->opt);
			xa[0] = j<rtk->na ? rtk->xa[j] : 0.0;
			fprintf(fp_stat, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
				i + 1, rtk->x[j], xa[0]);
		}
	}
	/* receiver h/w bias */
	if (est&&rtk->opt.glomodear == 2) {
		for (i = 0; i<nfreq; i++) {
			j = IL(i, &rtk->opt);
			xa[0] = j<rtk->na ? rtk->xa[j] : 0.0;
			fprintf(fp_stat, "$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
				i + 1, rtk->x[j], xa[0]);
		}
	}
	if (rtk->sol.stat == SOLQ_NONE || statlevel <= 1) return;

	/* residuals and status */
	for (i = 0; i<MAXSAT; i++) {
		ssat = rtk->ssat + i;
		if (!ssat->vs) continue;
		satno2id(i + 1, id);
		for (j = 0; j<nfreq; j++) {
			fprintf(fp_stat, "$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
				week, tow, id, j + 1, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
				ssat->resp[j], ssat->resc[j], ssat->vsat[j], ssat->snr[j] * 0.25,
				ssat->fix[j], ssat->slip[j] & 3, ssat->lock[j], ssat->outc[j],
				ssat->slipc[j], ssat->rejc[j]);
		}
	}
}

/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
	sol_t sol0 = { { 0 } };
	ambc_t ambc0 = { { { 0 } } };
	ssat_t ssat0 = { 0 };
	int i;

	trace(3, "rtkinit :\n");

	rtk->sol = sol0;
	for (i = 0; i<6; i++) rtk->rb[i] = 0.0;
	rtk->nx = opt->mode <= PMODE_FIXED ? NX(opt) : pppnx(opt);
	rtk->na = opt->mode <= PMODE_FIXED ? NR(opt) : 0;
	rtk->tt = 0.0;
	rtk->x = zeros(rtk->nx, 1);
	rtk->P = zeros(rtk->nx, rtk->nx);
	rtk->xa = zeros(rtk->na, 1);
	rtk->Pa = zeros(rtk->na, rtk->na);
	rtk->nfix = rtk->neb = 0;
	for (i = 0; i<MAXSAT; i++) {
		rtk->ambc[i] = ambc0;
		rtk->ssat[i] = ssat0;
	}
	for (i = 0; i<MAXERRMSG; i++) rtk->errbuf[i] = 0;
	rtk->opt = *opt;
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
{
	trace(3, "rtkfree :\n");

	rtk->nx = rtk->na = 0;
	free(rtk->x); rtk->x = NULL;
	free(rtk->P); rtk->P = NULL;
	free(rtk->xa); rtk->xa = NULL;
	free(rtk->Pa); rtk->Pa = NULL;
}
/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by
* precise positioning
* args   : rtk_t *rtk       IO  rtk control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .dtr[0]    O   receiver clock bias (s)
*                .dtr[1]    O   receiver glonass-gps time offset (s)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellite
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  sat(s+1) status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
extern int rtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt = &rtk->opt;
	sol_t solb = { { 0 } };
	gtime_t time;
	int i, nu, nr;
	char msg[128] = "";

	trace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);
	trace(4, "obs=\n"); traceobs(4, obs, n);
	/*trace(5,"nav=\n"); tracenav(5,nav);*/

	/* set base staion position */
	if (opt->refpos <= 3 && opt->mode != PMODE_SINGLE&&opt->mode != PMODE_MOVEB) {
		for (i = 0; i<6; i++) rtk->rb[i] = i<3 ? opt->rb[i] : 0.0;
	}
	/* count rover/base station observations */
	for (nu = 0; nu <n&&obs[nu].rcv == 1; nu++);
	for (nr = 0; nu + nr<n&&obs[nu + nr].rcv == 2; nr++);

	time = rtk->sol.time; /* previous epoch */

						  /* rover position by single point positioning */
	if (!pntpos(obs, nu, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {
		errmsg(rtk, "point pos error (%s)\n", msg);

		if (!rtk->opt.dynamics) {
			outsolstat(rtk);
			return 0;
		}
	}
	if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);

	/* single point positioning */
	if (opt->mode == PMODE_SINGLE) {
		outsolstat(rtk);
		return 1;
	}
	/* precise point positioning */
	if (opt->mode >= PMODE_PPP_KINEMA) {
		pppos(rtk, obs, nu, nav);
		pppoutsolstat(rtk, statlevel, fp_stat);
		return 1;
	}
}