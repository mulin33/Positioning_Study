#include "rtklib.h"

#define NUMSYS      6                   /* number of systems */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */


/* read rinex clock ----------------------------------------------------------*/
static int readrnxclk(FILE *fp, const char *opt, int index, nav_t *nav)
{
	pclk_t *nav_pclk;
	gtime_t time;
	double data[2];
	int i, j, sat;
	char buff[MAXRNXLEN], satid[8] = "";

	trace(3, "readrnxclk: index=%d\n", index);

	if (!nav) return 0;

	while (fgets(buff, sizeof(buff), fp)) {

		if (str2time(buff, 8, 26, &time)) {
			trace(2, "rinex clk invalid epoch: %34.34s\n", buff);
			continue;
		}
		strncpy(satid, buff + 3, 4);

		/* only read AS (satellite clock) record */
		if (strncmp(buff, "AS", 2) || !(sat = satid2no(satid))) continue;

		for (i = 0, j = 40; i<2; i++, j += 20) data[i] = str2num(buff, j, 19);

		if (nav->nc >= nav->ncmax) {
			nav->ncmax += 1024;
			if (!(nav_pclk = (pclk_t *)realloc(nav->pclk, sizeof(pclk_t)*(nav->ncmax)))) {
				trace(1, "readrnxclk malloc error: nmax=%d\n", nav->ncmax);
				free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
				return -1;
			}
			nav->pclk = nav_pclk;
		}
		if (nav->nc <= 0 || fabs(timediff(time, nav->pclk[nav->nc - 1].time))>1E-9) {
			nav->nc++;
			nav->pclk[nav->nc - 1].time = time;
			nav->pclk[nav->nc - 1].index = index;
			for (i = 0; i<MAXSAT; i++) {
				nav->pclk[nav->nc - 1].clk[i][0] = 0.0;
				nav->pclk[nav->nc - 1].std[i][0] = 0.0f;
			}
		}
		nav->pclk[nav->nc - 1].clk[sat - 1][0] = data[0];
		nav->pclk[nav->nc - 1].std[sat - 1][0] = (float)data[1];
	}
	return nav->nc>0;
}


/* uncompress and read rinex file --------------------------------------------*/
static int readrnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
	const char *opt, int flag, int index, char *type,
	obs_t *obs, nav_t *nav, sta_t *sta)
{
	FILE *fp;
	int cstat, stat;
	char tmpfile[1024];

	trace(3, "readrnxfile: file=%s flag=%d index=%d\n", file, flag, index);


	/* uncompress file */
	if ((cstat = uncompress(file, tmpfile))<0) {
		trace(2, "rinex file uncompact error: %s\n", file);
		return 0;
	}
	if (!(fp = fopen(cstat ? tmpfile : file, "r"))) {
		trace(2, "rinex file open error: %s\n", cstat ? tmpfile : file);
		return 0;
	}
	/* read rinex file */
	stat = readrnxclk(fp, opt, index, nav);

	fclose(fp);

	/* delete temporary file */
	if (cstat) remove(tmpfile);

	return stat;
}


/* compare precise clock -----------------------------------------------------*/
static int cmppclk(const void *p1, const void *p2)
{
	pclk_t *q1 = (pclk_t *)p1, *q2 = (pclk_t *)p2;
	double tt = timediff(q1->time, q2->time);
	return tt<-1E-9 ? -1 : (tt>1E-9 ? 1 : q1->index - q2->index);
}
/* combine precise clock -----------------------------------------------------*/
static void combpclk(nav_t *nav)
{
	pclk_t *nav_pclk;
	int i, j, k;

	trace(3, "combpclk: nc=%d\n", nav->nc);

	if (nav->nc <= 0) return;

	qsort(nav->pclk, nav->nc, sizeof(pclk_t), cmppclk);

	for (i = 0, j = 1; j<nav->nc; j++) {
		if (fabs(timediff(nav->pclk[i].time, nav->pclk[j].time))<1E-9) {
			for (k = 0; k<MAXSAT; k++) {
				if (nav->pclk[j].clk[k][0] == 0.0) continue;
				nav->pclk[i].clk[k][0] = nav->pclk[j].clk[k][0];
				nav->pclk[i].std[k][0] = nav->pclk[j].std[k][0];
			}
		}
		else if (++i<j) nav->pclk[i] = nav->pclk[j];
	}
	nav->nc = i + 1;

	if (!(nav_pclk = (pclk_t *)realloc(nav->pclk, sizeof(pclk_t)*nav->nc))) {
		free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
		trace(1, "combpclk malloc error nc=%d\n", nav->nc);
		return;
	}
	nav->pclk = nav_pclk;
	nav->ncmax = nav->nc;

	trace(4, "combpclk: nc=%d\n", nav->nc);
}
/* read rinex clock files ------------------------------------------------------
* read rinex clock files
* args   : char *file    I      file (wild-card * expanded)
*          nav_t *nav    IO     navigation data    (NULL: no input)
* return : number of precise clock
*-----------------------------------------------------------------------------*/
extern int readrnxc(const char *file, nav_t *nav)
{
	gtime_t t = { 0 };
	int i, n, index = 0, stat = 1;
	char *files[MAXEXFILE] = { 0 }, type;

	trace(3, "readrnxc: file=%s\n", file);

	for (i = 0; i<MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]); return 0;
		}
	}
	/* expand wild-card */
	n = expath(file, files, MAXEXFILE);

	/* read rinex clock files */
	for (i = 0; i<n; i++) {
		if (readrnxfile(files[i], t, t, 0.0, NULL, 1, index++, &type, NULL, nav, NULL)) {
			continue;
		}
		stat = 0;
		break;
	}
	for (i = 0; i<MAXEXFILE; i++) free(files[i]);

	if (!stat) return 0;

	/* unique and combine ephemeris and precise clock */
	combpclk(nav);

	return nav->nc;
}
