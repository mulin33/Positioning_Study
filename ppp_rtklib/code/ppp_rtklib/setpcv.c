#include "rtklib.h"

static sta_t stas[MAXRCV];      /* station infomation */
/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args   : int    sat         I   satellite number (0: receiver antenna)
*          char   *type       I   antenna type for receiver antenna
*          gtime_t time       I   time to search parameters
*          pcvs_t *pcvs       IO  antenna parameters
* return : antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time,
	const pcvs_t *pcvs)
{
	pcv_t *pcv;
	char buff[MAXANT], *types[2], *p;
	int i, j, n = 0;

	trace(3, "searchpcv: sat=%2d type=%s\n", sat, type);

	if (sat) { /* search satellite antenna */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			if (pcv->sat != sat) continue;
			if (pcv->ts.time != 0 && timediff(pcv->ts, time)>0.0) continue;
			if (pcv->te.time != 0 && timediff(pcv->te, time)<0.0) continue;
			return pcv;
		}
	}
	else {
		strcpy(buff, type);
		for (p = strtok(buff, " "); p&&n<2; p = strtok(NULL, " ")) types[n++] = p;
		if (n <= 0) return NULL;

		/* search receiver antenna with radome at first */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			for (j = 0; j<n; j++) if (!strstr(pcv->type, types[j])) break;
			if (j >= n) return pcv;
		}
		/* search receiver antenna without radome */
		for (i = 0; i<pcvs->n; i++) {
			pcv = pcvs->pcv + i;
			if (strstr(pcv->type, types[0]) != pcv->type) continue;

			trace(2, "pcv without radome is used type=%s\n", type);
			return pcv;
		}
	}
	return NULL;
}

/* set antenna parameters ----------------------------------------------------*/
 void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
	const pcvs_t *pcvr, const sta_t *sta)
{
	pcv_t *pcv;
	double pos[3], del[3];
	int i, j, mode = PMODE_DGPS <= popt->mode&&popt->mode <= PMODE_FIXED;
	char id[64];

	/* set satellite antenna parameters */
	for (i = 0; i<MAXSAT; i++) {
		if (!(satsys(i + 1, NULL)&popt->navsys)) continue;
		if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
			satno2id(i + 1, id);
			trace(2, "no satellite antenna pcv: %s\n", id);
			continue;
		}
		nav->pcvs[i] = *pcv;
	}
	for (i = 0; i<(mode ? 2 : 1); i++) {
		if (!strcmp(popt->anttype[i], "*")) { /* set by station parameters */
			strcpy(popt->anttype[i], sta[i].antdes);
			if (sta[i].deltype == 1) { /* xyz */
				if (norm(sta[i].pos, 3)>0.0) {
					ecef2pos(sta[i].pos, pos);
					ecef2enu(pos, sta[i].del, del);
					for (j = 0; j<3; j++) popt->antdel[i][j] = del[j];
				}
			}
			else { /* enu */
				for (j = 0; j<3; j++) popt->antdel[i][j] = stas[i].del[j];
			}
		}
		if (!(pcv = searchpcv(0, popt->anttype[i], time, pcvr))) {
			trace(2, "no receiver antenna pcv: %s\n", popt->anttype[i]);
			*popt->anttype[i] = '\0';
			continue;
		}
		strcpy(popt->anttype[i], pcv->type);
		popt->pcvr[i] = *pcv;
	}
}