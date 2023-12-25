# day06_最小二乘和单点定位

## 1 线性和非线性最小二乘公式

首先参考[最小二乘原理.ppt]()

### 1.1 Linear LSE(Least Square Estimation)

原式：
$$
y=Hx+v
$$
结果：
$$
x=(H^TWH)^{-1}H^TWy	\\
W=diag(\sigma^{-2}_1,\sigma^{-2}_2,...,\sigma^{-2}_m)
$$

### 1.2 Gauss-Newton iteration for non-linear LSE

原式：
$$
y=h(x)+v	\\
y=h(x_0)+H(x-x_0)+v	& ——泰勒展开\\
H=\frac{\partial h(x)}{\partial x}|_{x=x_0}
$$

结果：

**如果*x0*不是真值，需要通过*Gauss-Newton*迭代获得**
$$
x=x_0+(H^TWH)^{-1}H^TW(y-h(x_0))	\\
x_{i+1}=x_i+(H^TWH)^{-1}H^TW(y-h(x_0)) & iteration
$$

## 2 单点定位SPP算法

### 2.1 原式

根据[Rinex 文件格式简介]() Chapter3.Basic Definitions:

**上标(n)表示卫星标识**
$$
Pseudo Range(PR)=\rho^{(n)}+cdt_r-cdT^{(n)}+I^{(n)}+T^{(n)}+\varepsilon	\\
\rho^{(n)}=\sqrt{(x^{(n)}-x)^2+(y^{(n)}-y)^2+(z^{(n)}-z)^2}
$$

$$
y=h(x)+v	\\
\\
x=(x,y,z,cdt_r)^T	\\
y=(P^1_r,P^2_r,...,P^n_r)^T	\\

h(x)=\left[
\begin{matrix}
\rho^1+cdt_r-cdT^1+I^1+T^1	\\
\rho^2+cdt_r-cdT^2+I^2+T^2	\\
.							\\
.							\\
\rho^n+cdt_r-cdT^n+I^n+T^n
\end{matrix}
\right]
$$

### 2.2 线性化

**雅可比矩阵*H***只与各颗卫星相对于用户的几何位置有关，因此***H***通常被称为**几何矩阵**
$$
y=h(x_0)+H(x-x_0)+v	\\
\\
\frac{\partial \rho^{(n)}}{\partial x}=\frac{-(x^{(n)}-x)}{\sqrt{(x^{(n)}-x)^2+(y^{(n)}-y)^2+(z^{(n)}-z)^2}}=\frac{-(x^{(n)}-x)}{\rho^{(n)}}	\\
\therefore H\cdot(x-x_0)=\left[
\begin{matrix}
-\frac{x^1-x}{\rho^1} & -\frac{y^1-y}{\rho^1} & -\frac{z^1-z}{\rho^1} & 1	\\
-\frac{x^2-x}{\rho^2} & -\frac{y^2-y}{\rho^2} & -\frac{z^2-z}{\rho^2} & 1	\\
.																		\\
.																		\\
-\frac{x^n-x}{\rho^n} & -\frac{y^n-y}{\rho^n} & -\frac{z^n-z}{\rho^n} & 1
\end{matrix}
\right]

\cdot	%%点乘符号

\left[
\begin{matrix}
\Delta x \\
\Delta y \\
\Delta z \\
\Delta cdt_r
\end{matrix}
\right]
$$

## 3 权阵和误差

### 3.1 公式

在RTKLIB中，***weight matrix W***定义如下：

![image-20231123085054276](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231123085054276.png)

### 3.2 代码

1. **上面误差代码实现如下：**

```c
/* varerr()函数来自pntpos.c and 代码片段来自pntpos.c->rescode()------------------------*/
/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);	//这里GPS、Galileo、QZSS作者设定误差是一样的
    if (el<MIN_EL) el=MIN_EL;	//第一次迭代时，el=0,将其赋值MIN_EL(5°)
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free 无电离层组合时还需要再乘以SQR(3.0)*/
    return SQR(fact)*varr;
}

/* variance of pseudorange error */
var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
```

2. **各方差源获得代码如下链接：**

$F^s$​[卫星系统误差因子](#### satellite system error factor)

$R_r,a_\sigma,b_\sigma$[伪距载波ration和载波方差](####载波error ratio/factor)

$\sigma_{eph}$  [轨道和钟差方差](####pos和clk方差)

$\sigma_{bias}$  [硬件延迟]()

-----

#### $(F^s)$satellite system error factor

```c
/* 代码片段来自rtklib.h--------------------------------------------------------------*/
#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_IRN   1.5                 /* error factor: IRNSS */
#define EFACT_SBS   3.0                 /* error factor: SBAS */
```

#### $(R_r,a_\sigma,b_\sigma)$伪距载波error ratio/factor

伪距单点定位计算权重时，*varerr()*函数调用了*prcopt_t->err[5]*，该err[]值初始化在***main()***函数里的下面这行代码(也可以自己定义)：

```c
//prcopt_default 来自rtkcmn.c
prcopt_t prcopt = prcopt_default;
```

```c
/* 代码片段来自rtklib.h->prcopt_t结构体定义->err[5]-------------------------------------*/
double err[5];      /* measurement error factor */
                    /* [0]:reserved */
                    /* [1-3]:error factor a/b/c of phase (m) */
                    /* [4]:doppler frequency (hz) */
```

***err[5]***各参数代表含义：

| *prcopt_t->err[5]* |                        *descriptions*                        | *defalt value* |
| :----------------: | :----------------------------------------------------------: | :------------: |
|      *err[0]*      | *the ratio of std of pseudorange errors to carrier-phase for L1/L2/L5* |      100       |
|      *err[1]*      |       *base term of carrier-phase error std* ***(m)***       |     0.003      |
|      *err[2]*      | *the elevation dependent term of carrier-phase error std* ***(m/sin(el))*** |     0.003      |
|      *err[3]*      | *the baseline-length dependent term of carrier-phase error std* ***(m/10km)*** |       0        |
|      *err[4]*      | (Current version does not use the value) *the std of Doppler errors* ***(Hz)*** |       1        |

![image-20231122110342259](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231122110342259.png)

#### $(\sigma_{eph})$pos和clk方差


$$
\sigma_{eph}^2=vare+varc
$$

>  在SPP定位中，从n文件获得星历误差

(留坑待补)。。。。。

> 在PPP定位中，从sp3+clk文件(or 仅sp3文件)获得星历误差

1. **预定义代码片段**

```c
/* 预定义代码来自preceph.c头部---------------------------------------------------------*/
#define EXTERR_CLK  1E-3     /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7     /* extrapolation error for ephem (m/s^2) */
```

2. **$vare$获得：**

从sp3文件读取到某个历元某卫星的x、y、z标准差，正常情况下，$vare$赋值为**靠近插值时间的前一个历元卫星位置方差**,没有进行插值获得

```c
/* 代码片段来自preceph.c->pephpos()函数---------------------------------------------------------------*/
if (vare) {	//如果定义了vare浮点数指针则取值精密星历里的标准差，vare(variance of ephemeries)
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
        
        /* extrapolation error for orbit 推断轨道误差*/
		/* 当插值历元位于轨道弧的首或末端，因为其不处于轨道弧中间，其误差会增大，所以需额外加上推断误差
		 * vare赋值是靠近插值时间的历元的方差，没有对轨道方差进行插值，why?----------------------------*/
		/* question:EXTERR_EPH和EXTERR_CLK是怎么赋值的呢-----------------------------------------------*/
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
```

3. **$varc$获得(从sp3和clk)：**

从sp3和clk文件获得某历元某卫星钟差标准差，正常情况下，$varc$标准差赋值为**靠近插值时间的前一个历元卫星位置方差加上预定义系数*时间差 **

```c
 /* linear interpolation for clock 钟差线性插值*/
	/* dts存储钟差，std包括卫星位置和钟差的标准差*/
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
    if (t[0]<=0.0) {	//插值历元在轨道弧前面
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {	//插值历元在轨道弧后面
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
	/* question:下面这段代码应该有一处问题
	 * 1、计算std那里前面的钟差没有乘以光速CLIGHT，根据作者注释std单位为m*/
    else if (c[0]!=0.0&&c[1]!=0.0) {	//插值历元在轨道弧上
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);	//插值钟差
        i=t[0]<-t[1]?0:1;
		/* 将std存储单位修从->m，如下：----------------------------------*/
        //std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
		std = nav->peph[index + i].std[sat - 1][3] * CLIGHT + EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
```

#### $(\sigma_{bias})$硬件延迟

1. 如果是无电离层组合,$\sigma_{vare}=0$
2. 如果是单频，$\sigma_{vare}=SQR(ERR\_CBIAS);$# define ERR_CBIAS   0.3 

```c
/* 代码片段来自pntpos.c->prange()函数-----------------------------------------------------------*/
/* psendorange with code bias correction 改正硬件延迟的代码-------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    ...//省略，只保留关键代码
    *var=0.0;	//vmeas初始化为0
    
    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency 如果是无电离层组合*/
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
            gamma=SQR(FREQ1/FREQ2);
            return (P2-gamma*P1)/(1.0-gamma);	//返回无电离层组合观测值,vare=0
        }
        ...//此处省略其他卫星系统
    }
    else { /* single-freq (L1/E1/B1) 单频，卫星钟差包含码偏差，需要调用gettgd()进行纠正*/
        *var=SQR(ERR_CBIAS);	//#define ERR_CBIAS   0.3         /* code bias error Std (m) */
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;	
        }
         ...//此处省略其他卫星系统
    }
}
```



```c
/* 代码片段来自pntpos.c--------------------------------------------------------------*/
/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#if 0 /* enable GPS-QZS time offset estimation */
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
```

### 3.3 误差计算代码

```c

```

$$
W=diag(\sigma^{-2}_1,\sigma^{-2}_2,...,\sigma^{-2}_m)	\\
\sigma^{-2}=F^sR_r(a_\sigma^2+b_\sigma^2/sinEl^s_r)+\sigma^2_{eph}+\sigma^2_{ion}+\sigma^2_{trop}+\sigma^2_{bias}	\\


F^s ——\text{satellite system error factor}\text{(1:GPS,Galileo,QZSS and BeiDou; 1.5:GLONASS; 3:SBAS)}	\\
\leftline{R_r——\text{code/carrier-phase error ration}	}\\
a_\sigma,b_\sigma——\text{carrier-phase error factor a and b(m)}	\\
$$






## 参考文献

[RTKLIB manual_2.4.2]()

