#include "rtklib.h"

/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : gtime_t time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(gtime_t time, const double *pos, const double *azel,
	double humi)
{
	const double temp0 = 15.0; /* temparature at sea level */
	double hgt, pres, temp, e, z, trph, trpw;

	if (pos[2]<-100.0 || 1E4<pos[2] || azel[1] <= 0) return 0.0;

	/* standard atmosphere */
	hgt = pos[2]<0.0 ? 0.0 : pos[2];

	pres = 1013.25*pow(1.0 - 2.2557E-5*hgt, 5.2568);
	temp = temp0 - 6.5E-3*hgt + 273.16;
	e = 6.108*humi*exp((17.15*temp - 4684.0) / (temp - 38.45));

	/* saastamoninen model */
	z = PI / 2.0 - azel[1];
	trph = 0.0022768*pres / (1.0 - 0.00266*cos(2.0*pos[0]) - 0.00028*hgt / 1E3) / cos(z);
	trpw = 0.002277*(1255.0 / temp + 0.05)*e / cos(z);
	return trph + trpw;
}

static double interpc(const double coef[], double lat)
{
	int i = (int)(lat / 15.0);
	if (i<1) return coef[0]; else if (i>4) return coef[4];
	return coef[i - 1] * (1.0 - lat / 15.0 + i) + coef[i] * (lat / 15.0 - i);
}
static double mapf(double el, double a, double b, double c)
{
	double sinel = sin(el);
	return (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
}
static double nmf(gtime_t time, const double pos[], const double azel[],
	double *mapfw)
{
	/* ref [5] table 3 */
	/* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
	const double coef[][5] = {
		{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3 },
		{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3 },
		{ 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3 },

		{ 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5 },
		{ 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5 },
		{ 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5 },

		{ 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4 },
		{ 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3 },
		{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2 }
	};
	const double aht[] = { 2.53E-5, 5.49E-3, 1.14E-3 }; /* height correction */

	double y, cosy, ah[3], aw[3], dm, el = azel[1], lat = pos[0] * R2D, hgt = pos[2];
	int i;

	if (el <= 0.0) {
		if (mapfw) *mapfw = 0.0;
		return 0.0;
	}
	/* year from doy 28, added half a year for southern latitudes */
	y = (time2doy(time) - 28.0) / 365.25 + (lat<0.0 ? 0.5 : 0.0);

	cosy = cos(2.0*PI*y);
	lat = fabs(lat);

	for (i = 0; i<3; i++) {
		ah[i] = interpc(coef[i], lat) - interpc(coef[i + 3], lat)*cosy;
		aw[i] = interpc(coef[i + 6], lat);
	}
	/* ellipsoidal height is used instead of height above sea level */
	dm = (1.0 / sin(el) - mapf(el, aht[0], aht[1], aht[2]))*hgt / 1E3;

	if (mapfw) *mapfw = mapf(el, aw[0], aw[1], aw[2]);

	return mapf(el, ah[0], ah[1], ah[2]) + dm;
}

/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref [5] (NMF) and [9] (GMF)
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
extern double tropmapf(gtime_t time, const double pos[], const double azel[],
	double *mapfw)
{
#ifdef IERS_MODEL
	const double ep[] = { 2000,1,1,12,0,0 };
	double mjd, lat, lon, hgt, zd, gmfh, gmfw;
#endif
	trace(4, "tropmapf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",
		pos[0] * R2D, pos[1] * R2D, pos[2], azel[0] * R2D, azel[1] * R2D);

	if (pos[2]<-1000.0 || pos[2]>20000.0) {
		if (mapfw) *mapfw = 0.0;
		return 0.0;
	}
#ifdef IERS_MODEL
	mjd = 51544.5 + (timediff(time, epoch2time(ep))) / 86400.0;
	lat = pos[0];
	lon = pos[1];
	hgt = pos[2] - geoidh(pos); /* height in m (mean sea level) */
	zd = PI / 2.0 - azel[1];

	/* call GMF */
	gmf_(&mjd, &lat, &lon, &hgt, &zd, &gmfh, &gmfw);

	if (mapfw) *mapfw = gmfw;
	return gmfh;
#else
	return nmf(time, pos, azel, mapfw); /* NMF */
#endif
}


/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
	const double *azel)
{
	const double ion_default[] = { /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	double tt, f, psi, phi, lam, amp, per, x;
	int week;

	if (pos[2]<-1E3 || azel[1] <= 0) return 0.0;
	if (norm(ion, 8) <= 0.0) ion = ion_default;

	/* earth centered angle (semi-circle) */
	psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi = pos[0] / PI + psi*cos(azel[0]);
	if (phi> 0.416) phi = 0.416;
	else if (phi<-0.416) phi = -0.416;
	lam = pos[1] / PI + psi*sin(azel[0]) / cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
	phi += 0.064*cos((lam - 1.617)*PI);

	/* local time (s) */
	tt = 43200.0*lam + time2gpst(t, &week);
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

									   /* slant factor */
	f = 1.0 + 16.0*pow(0.53 - azel[1] / PI, 3.0);

	/* ionospheric delay */
	amp = ion[0] + phi*(ion[1] + phi*(ion[2] + phi*ion[3]));
	per = ion[4] + phi*(ion[5] + phi*(ion[6] + phi*ion[7]));
	amp = amp<    0.0 ? 0.0 : amp;
	per = per<72000.0 ? 72000.0 : per;
	x = 2.0*PI*(tt - 50400.0) / per;

	return CLIGHT*f*(fabs(x)<1.57 ? 5E-9 + amp*(1.0 + x*x*(-0.5 + x*x / 24.0)) : 5E-9);
}

/* phase windup correction -----------------------------------------------------
* phase windup correction (ref [7] 5.1.2)
* args   : gtime_t time     I   time (GPST)
*          double  *rs      I   satellite position (ecef) {x,y,z} (m)
*          double  *rr      I   receiver  position (ecef) {x,y,z} (m)
*          double  *phw     IO  phase windup correction (cycle)
* return : none
* notes  : the previous value of phase windup correction should be set to *phw
*          as an input. the function assumes windup correction has no jump more
*          than 0.5 cycle.
*-----------------------------------------------------------------------------*/
extern void windupcorr(gtime_t time, const double *rs, const double *rr,
	double *phw)
{
	double ek[3], exs[3], eys[3], ezs[3], ess[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], rsun[3], cosp, ph, erpv[5] = { 0 };
	int i;

	trace(4, "windupcorr: time=%s\n", time_str(time, 0));

	/* sun position in ecef */
	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

	/* unit vector satellite to receiver */
	for (i = 0; i<3; i++) r[i] = rr[i] - rs[i];
	if (!normv3(r, ek)) return;

	/* unit vectors of satellite antenna */
	for (i = 0; i<3; i++) r[i] = -rs[i];
	if (!normv3(r, ezs)) return;
	for (i = 0; i<3; i++) r[i] = rsun[i] - rs[i];
	if (!normv3(r, ess)) return;
	cross3(ezs, ess, r);
	if (!normv3(r, eys)) return;
	cross3(eys, ezs, exs);

	/* unit vectors of receiver antenna */
	ecef2pos(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

													/* phase windup effect */
	cross3(ek, eys, eks);
	cross3(ek, eyr, ekr);
	for (i = 0; i<3; i++) {
		ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];
		dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];
	}
	cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
	if (cosp<-1.0) cosp = -1.0;
	else if (cosp> 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	cross3(ds, dr, drs);
	if (dot(ek, drs, 3)<0.0) ph = -ph;

	*phw = ph + floor(*phw - ph + 0.5); /* in cycle */
}

/* interpolate antenna phase center variation --------------------------------*/
static double interpvar(double ang, const double *var)
{
	double a = ang / 5.0; /* ang=0-90 */
	int i = (int)a;
	if (i<0) return var[0]; else if (i >= 18) return var[18];
	return var[i] * (1.0 - a + i) + var[i + 1] * (a - i);
}
/* receiver antenna model ------------------------------------------------------
* compute antenna offset by antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double *azel     I   azimuth/elevation for receiver {az,el} (rad)
*          int     opt      I   option (0:only offset,1:offset+pcv)
*          double *dant     O   range offsets for each frequency (m)
* return : none
* notes  : current version does not support azimuth dependent terms
*-----------------------------------------------------------------------------*/
extern void antmodel(const pcv_t *pcv, const double *del, const double *azel,
	int opt, double *dant)
{
	double e[3], off[3], cosel = cos(azel[1]);
	int i, j;

	trace(4, "antmodel: azel=%6.1f %4.1f opt=%d\n", azel[0] * R2D, azel[1] * R2D, opt);

	e[0] = sin(azel[0])*cosel;
	e[1] = cos(azel[0])*cosel;
	e[2] = sin(azel[1]);

	for (i = 0; i<NFREQ; i++) {
		for (j = 0; j<3; j++) off[j] = pcv->off[i][j] + del[j];

		dant[i] = -dot(off, e, 3) + (opt ? interpvar(90.0 - azel[1] * R2D, pcv->var[i]) : 0.0);
	}
	trace(5, "antmodel: dant=%6.3f %6.3f\n", dant[0], dant[1]);
}
/* satellite antenna model ------------------------------------------------------
* compute satellite antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double nadir     I   nadir angle for satellite (rad)
*          double *dant     O   range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void antmodel_s(const pcv_t *pcv, double nadir, double *dant)
{
	int i;

	trace(4, "antmodel_s: nadir=%6.1f\n", nadir*R2D);

	for (i = 0; i<NFREQ; i++) {
		dant[i] = interpvar(nadir*R2D / 5.0, pcv->var[i]);
	}
	trace(5, "antmodel_s: dant=%6.3f %6.3f\n", dant[0], dant[1]);
}

/*ÓÃÓÚtidedisp*/
/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int geterp(const erp_t *erp, gtime_t time, double *erpv)
{
	const double ep[] = { 2000,1,1,12,0,0 };
	double mjd, day, a;
	int i = 0, j, k;

	trace(4, "geterp:\n");

	if (erp->n <= 0) return 0;

	mjd = 51544.5 + (timediff(gpst2utc(time), epoch2time(ep))) / 86400.0;

	if (mjd <= erp->data[0].mjd) {
		day = mjd - erp->data[0].mjd;
		erpv[0] = erp->data[0].xp + erp->data[0].xpr*day;
		erpv[1] = erp->data[0].yp + erp->data[0].xpr*day;
		erpv[2] = erp->data[0].ut1_utc - erp->data[0].lod*day;
		erpv[3] = erp->data[0].lod;
		return 1;
	}
	if (mjd >= erp->data[erp->n - 1].mjd) {
		day = mjd - erp->data[erp->n - 1].mjd;
		erpv[0] = erp->data[erp->n - 1].xp + erp->data[erp->n - 1].xpr*day;
		erpv[1] = erp->data[erp->n - 1].yp + erp->data[erp->n - 1].ypr*day;
		erpv[2] = erp->data[erp->n - 1].ut1_utc - erp->data[erp->n - 1].lod*day;
		erpv[3] = erp->data[erp->n - 1].lod;
		return 1;
	}
	for (j = 0, k = erp->n - 1; j <= k;) {
		i = (j + k) / 2;
		if (mjd<erp->data[i].mjd) k = i - 1; else j = i + 1;
	}
	if (erp->data[i].mjd == mjd - erp->data[i + 1].mjd) {
		a = 0.5;
	}
	else {
		a = (mjd - erp->data[i + 1].mjd) / (erp->data[i].mjd - mjd - erp->data[i + 1].mjd);
	}
	erpv[0] = (1.0 - a)*erp->data[i].xp + a*erp->data[i + 1].xp;
	erpv[1] = (1.0 - a)*erp->data[i].yp + a*erp->data[i + 1].yp;
	erpv[2] = (1.0 - a)*erp->data[i].ut1_utc + a*erp->data[i + 1].ut1_utc;
	erpv[3] = (1.0 - a)*erp->data[i].lod + a*erp->data[i + 1].lod;
	return 1;
}