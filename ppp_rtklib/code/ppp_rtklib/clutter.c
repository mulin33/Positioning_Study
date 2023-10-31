#include "rtklib.h"

extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* expand file path ------------------------------------------------------------
* expand file path with wild-card (*) in file
* args   : char   *path     I   file path to expand (captal insensitive)
*          char   *paths    O   expanded file paths
*          int    nmax      I   max number of expanded file paths
* return : number of expanded file paths
* notes  : the order of expanded files is alphabetical order
*-----------------------------------------------------------------------------*/
extern int expath(const char *path, char *paths[], int nmax)
{
	int i, j, n = 0;
	char tmp[1024];
#ifdef WIN32
	WIN32_FIND_DATA file;
	HANDLE h;
	char dir[1024] = "", *p;

	trace(3, "expath  : path=%s nmax=%d\n", path, nmax);

	if ((p = strrchr(path, '\\'))) {
		strncpy(dir, path, p - path + 1); dir[p - path + 1] = '\0';
	}
	if ((h = FindFirstFile((LPCTSTR)path, &file)) == INVALID_HANDLE_VALUE) {
		strcpy(paths[0], path);
		return 1;
	}
	sprintf(paths[n++], "%s%s", dir, file.cFileName);
	while (FindNextFile(h, &file) && n<nmax) {
		if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
		sprintf(paths[n++], "%s%s", dir, file.cFileName);
	}
	FindClose(h);
#else
	struct dirent *d;
	DIR *dp;
	const char *file = path;
	char dir[1024] = "", s1[1024], s2[1024], *p, *q;

	trace(3, "expath  : path=%s nmax=%d\n", path, nmax);

	if ((p = strrchr(path, '/')) || (p = strrchr(path, '\\'))) {
		file = p + 1; strncpy(dir, path, p - path + 1); dir[p - path + 1] = '\0';
	}
	if (!(dp = opendir(*dir ? dir : "."))) return 0;
	while ((d = readdir(dp))) {
		if (*(d->d_name) == '.') continue;
		sprintf(s1, "^%s$", d->d_name);
		sprintf(s2, "^%s$", file);
		for (p = s1; *p; p++) *p = (char)tolower((int)*p);
		for (p = s2; *p; p++) *p = (char)tolower((int)*p);
		for (p = s1, q = strtok(s2, "*"); q; q = strtok(NULL, "*")) {
			if ((p = strstr(p, q))) p += strlen(q); else break;
		}
		if (p&&n<nmax) sprintf(paths[n++], "%s%s", dir, d->d_name);
	}
	closedir(dp);
#endif
	/* sort paths in alphabetical order */
	for (i = 0; i<n - 1; i++) {
		for (j = i + 1; j<n; j++) {
			if (strcmp(paths[i], paths[j])>0) {
				strcpy(tmp, paths[i]);
				strcpy(paths[i], paths[j]);
				strcpy(paths[j], tmp);
			}
		}
	}
	for (i = 0; i<n; i++) trace(3, "expath  : file=%s\n", paths[i]);

	return n;
}
/* create directory ------------------------------------------------------------
* create directory if not exist
* args   : char   *path     I   file path to be saved
* return : none
* notes  : not recursive. only one level
*-----------------------------------------------------------------------------*/
extern void createdir(const char *path)
{
	char buff[1024], *p;

	tracet(3, "createdir: path=%s\n", path);

	strcpy(buff, path);
	if (!(p = strrchr(buff, FILEPATHSEP))) return;
	*p = '\0';

#ifdef WIN32
	CreateDirectory(buff, NULL);
#else
	mkdir(buff, 0777);
#endif
}

/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	if (norm(rs, 3)<RE_WGS84) return -1.0;
	for (i = 0; i<3; i++) e[i] = rs[i] - rr[i];
	r = norm(e, 3);
	for (i = 0; i<3; i++) e[i] /= r;
	return r + OMGE*(rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT;
}
/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern double satazel(const double *pos, const double *e, double *azel)
{
	double az = 0.0, el = PI / 2.0, enu[3];

	if (pos[2]>-RE_WGS84) {
		ecef2enu(pos, e, enu);
		az = dot(enu, enu, 2)<1E-12 ? 0.0 : atan2(enu[0], enu[1]);
		if (az<0.0) az += 2 * PI;
		el = asin(enu[2]);
	}
	if (azel) { azel[0] = az; azel[1] = el; }
	return el;
}
/* compute dops ----------------------------------------------------------------
* compute DOP (dilution of precision)
* args   : int    ns        I   number of satellites
*          double *azel     I   satellite azimuth/elevation angle (rad)
*          double elmin     I   elevation cutoff angle (rad)
*          double *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
* return : none
* notes  : dop[0]-[3] return 0 in case of dop computation error
*-----------------------------------------------------------------------------*/
#define SQRT(x)     ((x)<0.0?0.0:sqrt(x))

extern void dops(int ns, const double *azel, double elmin, double *dop)
{
	double H[4 * MAXSAT], Q[16], cosel, sinel;
	int i, n;

	for (i = 0; i<4; i++) dop[i] = 0.0;
	for (i = n = 0; i<ns&&i<MAXSAT; i++) {
		if (azel[1 + i * 2]<elmin || azel[1 + i * 2] <= 0.0) continue;
		cosel = cos(azel[1 + i * 2]);
		sinel = sin(azel[1 + i * 2]);
		H[4 * n] = cosel*sin(azel[i * 2]);
		H[1 + 4 * n] = cosel*cos(azel[i * 2]);
		H[2 + 4 * n] = sinel;
		H[3 + 4 * n++] = 1.0;
	}
	if (n<4) return;

	matmul("NT", 4, 4, n, 1.0, H, H, 0.0, Q);
	if (!matinv(Q, 4)) {
		dop[0] = SQRT(Q[0] + Q[5] + Q[10] + Q[15]); /* GDOP */
		dop[1] = SQRT(Q[0] + Q[5] + Q[10]);       /* PDOP */
		dop[2] = SQRT(Q[0] + Q[5]);             /* HDOP */
		dop[3] = SQRT(Q[10]);                 /* VDOP */
	}
}