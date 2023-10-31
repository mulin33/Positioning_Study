#include "rtklib.h"

int freeData(pcvs_t *pcvs, pcvs_t *pcvr,nav_t *nav,obs_t *obs)
{
	/* free antenna parameters */
	free(pcvs->pcv); pcvs->pcv = NULL; pcvs->n = pcvs->nmax = 0;
	free(pcvr->pcv); pcvr->pcv = NULL; pcvr->n = pcvr->nmax = 0;

	/*释放nav结构体*/
	int i = 0;
	free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
	free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
	free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
	for (i = 0; i<nav->nt; i++) {
		free(nav->tec[i].data);
		free(nav->tec[i].rms);
	}
	free(nav->tec); nav->tec = NULL; nav->nt = nav->ntmax = 0;

	/*free obs and nav data*/
	free(obs->data); obs->data = NULL; obs->n = obs->nmax = 0;
	free(nav->eph); nav->eph = NULL; nav->n = nav->nmax = 0;
	free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
	free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
}

