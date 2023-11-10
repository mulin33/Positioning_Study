#include "rtklib.h"

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
	const solopt_t *sopt,obs_t obss)
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
extern int outhead(const char *outfile, char **infile, int n,
	const prcopt_t *popt, const solopt_t *sopt,obs_t obss)
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
	outheader(fp, infile, n, popt, sopt, obss);

	if (*outfile) fclose(fp);

	return 1;
}
