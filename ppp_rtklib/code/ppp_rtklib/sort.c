#include "rtklib.h"

/*Input/output?functions-------------------------------------------------------*/

/* compare ephemeris ---------------------------------------------------------*/
static int cmpeph(const void *p1, const void *p2)
{
	eph_t *q1 = (eph_t *)p1, *q2 = (eph_t *)p2;
	return q1->ttr.time != q2->ttr.time ? (int)(q1->ttr.time - q2->ttr.time) :
		(q1->toe.time != q2->toe.time ? (int)(q1->toe.time - q2->toe.time) :
			q1->sat - q2->sat);
}
/* sort and unique ephemeris -------------------------------------------------*/
static void uniqeph(nav_t *nav)
{
	eph_t *nav_eph;
	int i, j;

	trace(3, "uniqeph: n=%d\n", nav->n);

	if (nav->n <= 0) return;

	qsort(nav->eph, nav->n, sizeof(eph_t), cmpeph);

	for (i = j = 0; i<nav->n; i++) {
		if (nav->eph[i].sat != nav->eph[j].sat ||
			nav->eph[i].iode != nav->eph[j].iode) {
			nav->eph[++j] = nav->eph[i];
		}
	}
	nav->n = j + 1;

	if (!(nav_eph = (eph_t *)realloc(nav->eph, sizeof(eph_t)*nav->n))) {
		trace(1, "uniqeph malloc error n=%d\n", nav->n);
		free(nav->eph); nav->eph = NULL; nav->n = nav->nmax = 0;
		return;
	}
	nav->eph = nav_eph;
	nav->nmax = nav->n;

	trace(4, "uniqeph: n=%d\n", nav->n);
}
/* compare glonass ephemeris -------------------------------------------------*/
static int cmpgeph(const void *p1, const void *p2)
{
	geph_t *q1 = (geph_t *)p1, *q2 = (geph_t *)p2;
	return q1->tof.time != q2->tof.time ? (int)(q1->tof.time - q2->tof.time) :
		(q1->toe.time != q2->toe.time ? (int)(q1->toe.time - q2->toe.time) :
			q1->sat - q2->sat);
}
/* sort and unique glonass ephemeris -----------------------------------------*/
static void uniqgeph(nav_t *nav)
{
	geph_t *nav_geph;
	int i, j;

	trace(3, "uniqgeph: ng=%d\n", nav->ng);

	if (nav->ng <= 0) return;

	qsort(nav->geph, nav->ng, sizeof(geph_t), cmpgeph);

	for (i = j = 0; i<nav->ng; i++) {
		if (nav->geph[i].sat != nav->geph[j].sat ||
			nav->geph[i].toe.time != nav->geph[j].toe.time ||
			nav->geph[i].svh != nav->geph[j].svh) {
			nav->geph[++j] = nav->geph[i];
		}
	}
	nav->ng = j + 1;

	if (!(nav_geph = (geph_t *)realloc(nav->geph, sizeof(geph_t)*nav->ng))) {
		trace(1, "uniqgeph malloc error ng=%d\n", nav->ng);
		free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
		return;
	}
	nav->geph = nav_geph;
	nav->ngmax = nav->ng;

	trace(4, "uniqgeph: ng=%d\n", nav->ng);
}
/* compare sbas ephemeris ----------------------------------------------------*/
static int cmpseph(const void *p1, const void *p2)
{
	seph_t *q1 = (seph_t *)p1, *q2 = (seph_t *)p2;
	return q1->tof.time != q2->tof.time ? (int)(q1->tof.time - q2->tof.time) :
		(q1->t0.time != q2->t0.time ? (int)(q1->t0.time - q2->t0.time) :
			q1->sat - q2->sat);
}
/* sort and unique sbas ephemeris --------------------------------------------*/
static void uniqseph(nav_t *nav)
{
	seph_t *nav_seph;
	int i, j;

	trace(3, "uniqseph: ns=%d\n", nav->ns);

	if (nav->ns <= 0) return;

	qsort(nav->seph, nav->ns, sizeof(seph_t), cmpseph);

	for (i = j = 0; i<nav->ns; i++) {
		if (nav->seph[i].sat != nav->seph[j].sat ||
			nav->seph[i].t0.time != nav->seph[j].t0.time) {
			nav->seph[++j] = nav->seph[i];
		}
	}
	nav->ns = j + 1;

	if (!(nav_seph = (seph_t *)realloc(nav->seph, sizeof(seph_t)*nav->ns))) {
		trace(1, "uniqseph malloc error ns=%d\n", nav->ns);
		free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
		return;
	}
	nav->seph = nav_seph;
	nav->nsmax = nav->ns;

	trace(4, "uniqseph: ns=%d\n", nav->ns);
}
/* satellite carrier wave length -----------------------------------------------
* get satellite carrier wave lengths
* args   : int    sat       I   satellite number
*          int    frq       I   frequency index (0:L1,1:L2,2:L5/3,...)
*          nav_t  *nav      I   navigation messages
* return : carrier wave length (m) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double satwavelen(int sat, int frq, const nav_t *nav)
{
	const double freq_glo[] = { FREQ1_GLO,FREQ2_GLO,FREQ3_GLO };
	const double dfrq_glo[] = { DFRQ1_GLO,DFRQ2_GLO,0.0 };
	int i, sys = satsys(sat, NULL);

	if (sys == SYS_GLO) {
		if (0 <= frq&&frq <= 2) {
			for (i = 0; i<nav->ng; i++) {
				if (nav->geph[i].sat != sat) continue;
				return CLIGHT / (freq_glo[frq] + dfrq_glo[frq] * nav->geph[i].frq);
			}
		}
	}
	else if (sys == SYS_CMP) {
		if (frq == 1) return CLIGHT / FREQ2_CMP; /* B1 */
		else if (frq == 3) return CLIGHT / FREQ6_CMP; /* B3 */
		else if (frq == 4) return CLIGHT / FREQ7_CMP; /* B2 */
	}
	else {
		if (frq == 0) return CLIGHT / FREQ1; /* L1/E1 */
		else if (frq == 1) return CLIGHT / FREQ2; /* L2 */
		else if (frq == 2) return CLIGHT / FREQ5; /* L5/E5a */
		else if (frq == 3) return CLIGHT / FREQ6; /* L6/LEX */
		else if (frq == 4) return CLIGHT / FREQ7; /* E5b */
		else if (frq == 5) return CLIGHT / FREQ8; /* E5a+b */
	}
	return 0.0;
}
/* unique ephemerides ----------------------------------------------------------
* unique ephemerides in navigation data and update carrier wave length
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
	int i, j;

	trace(3, "uniqnav: neph=%d ngeph=%d nseph=%d\n", nav->n, nav->ng, nav->ns);

	/* unique ephemeris */
	uniqeph(nav);
	uniqgeph(nav);
	uniqseph(nav);

	/* update carrier wave length */
	for (i = 0; i<MAXSAT; i++) for (j = 0; j<NFREQ; j++) {
		nav->lam[i][j] = satwavelen(i + 1, j, nav);
	}
}
/* compare observation data -------------------------------------------------*/
static int cmpobs(const void *p1, const void *p2)
{
	obsd_t *q1 = (obsd_t *)p1, *q2 = (obsd_t *)p2;
	double tt = timediff(q1->time, q2->time);
	if (fabs(tt)>DTTOL) return tt<0 ? -1 : 1;
	if (q1->rcv != q2->rcv) return (int)q1->rcv - (int)q2->rcv;
	return (int)q1->sat - (int)q2->sat;
}
/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t *obs    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
	int i, j, n;

	trace(3, "sortobs: nobs=%d\n", obs->n);

	if (obs->n <= 0) return 0;

	qsort(obs->data, obs->n, sizeof(obsd_t), cmpobs);

	/* delete duplicated data */
	for (i = j = 0; i<obs->n; i++) {
		if (obs->data[i].sat != obs->data[j].sat ||
			obs->data[i].rcv != obs->data[j].rcv ||
			timediff(obs->data[i].time, obs->data[j].time) != 0.0) {
			obs->data[++j] = obs->data[i];
		}
	}
	obs->n = j + 1;

	for (i = n = 0; i<obs->n; i = j, n++) {
		for (j = i + 1; j<obs->n; j++) {
			if (timediff(obs->data[j].time, obs->data[i].time)>DTTOL) break;
		}
	}
	return n;
}
/* screen by time --------------------------------------------------------------
* screening by time start, time end, and time interval
* args   : gtime_t time  I      time
*          gtime_t ts    I      time start (ts.time==0:no screening by ts)
*          gtime_t te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-----------------------------------------------------------------------------*/
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint)
{
	return (tint <= 0.0 || fmod(time2gpst(time, NULL) + DTTOL, tint) <= DTTOL*2.0) &&
		(ts.time == 0 || timediff(time, ts) >= -DTTOL) &&
		(te.time == 0 || timediff(time, te)<  DTTOL);
}