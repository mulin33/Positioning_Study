# day02_文件读取

## 1 天线文件(.atx)

### 1.1 参考链接

> [igs 14 introduction]([files.igs.org/pub/data/format/antex14.txt](https://files.igs.org/pub/data/format/antex14.txt))
>
> [antenna calibration introduction](https://www.ngs.noaa.gov/ANTCAL/FAQ.xhtml#faq4)

### 1.2 前置知识

天线相位影响精度在cm级

**名词解释**

```
ARP：antenna reference point
PCO:phase center offset（和频率有关）
PCV：phase center variations（和频率有关）
```

**接收机天线相位改正部分(.atx)**

![image-20231104110823521](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231104110823521.png)

**卫星天线相位改正部分(.atx)**

![image-20231104111058094](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231104111058094.png)

### 1.4 问题

1. **question:难道pcvs->pcv定义的是指向结构体数组的指针?那么realloc函数在分配空间时难道有创建结构体数组的功能？**

```c
typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    double off[NFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
    double var[NFREQ][19]; /* phase center variation (m) */
                        /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

/* add antenna parameter -----------------------------------------------------*/
static void addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
    pcv_t *pcvs_pcv;
    
    if (pcvs->nmax<=pcvs->n) {
        pcvs->nmax+=256;
        if (!(pcvs_pcv=(pcv_t *)realloc(pcvs->pcv,sizeof(pcv_t)*pcvs->nmax))) {
            trace(1,"addpcv: memory allocation error\n");
            free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
            return;
        }
        pcvs->pcv=pcvs_pcv;
    }
    pcvs->pcv[pcvs->n++]=*pcv;	//question:难道pcvs->pcv定义的是指向结构体数组的指针
}
```

2. 时间转换函数还需要阅读

## 2 观测值文件(obs)

LLI失锁标志还得再研究研究

```c
/*这是一个用来存储obs头文件中（# / TYPES OF OBSERV）的三维数组----------*/
char tobs[NUMSYS][MAXOBSTYPE][4]={{""}}; //from readrnxfp.c

sigind_t index[NUMSYS]={{0}};
/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int idx[MAXOBSTYPE];                /* signal freq-index */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    uint8_t pri [MAXOBSTYPE];           /* signal priority (15-0) */
    uint8_t type[MAXOBSTYPE];           /* type (0:C,1:L,2:D,3:S) */
    uint8_t code[MAXOBSTYPE];           /* obs-code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

/* set_index函数*/
if (nsys>=1) set_index(ver,SYS_GPS,opt,tobs[0],index  );

static char *obscodes[]={       /* observation code strings */
    
    ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    "6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
    "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};
static char codepris[7][MAXFREQ][16]={  /* code priority for each freq-index */
   /*    0         1          2          3         4         5     */
    {"CPYWMNSL","PYWCMNDLSX","IQX"     ,""       ,""       ,""      ,""}, /* GPS */
    {"CPABX"   ,"PCABX"     ,"IQX"     ,""       ,""       ,""      ,""}, /* GLO */
    {"CABXZ"   ,"IQX"       ,"IQX"     ,"ABCXZ"  ,"IQX"    ,""      ,""}, /* GAL */
    {"CLSXZ"   ,"LSX"       ,"IQXDPZ"  ,"LSXEZ"  ,""       ,""      ,""}, /* QZS */
    {"C"       ,"IQX"       ,""        ,""       ,""       ,""      ,""}, /* SBS */
    {"IQXDPAN" ,"IQXDPZ"    ,"DPX"     ,"IQXA"   ,"DPX"    ,""      ,""}, /* BDS */
    {"ABCX"    ,"ABCX"      ,""        ,""       ,""       ,""      ,""}  /* IRN */
};
```







## 卫星系统标识[^卫星系统参考链接]

### 1 缩略词

**卫星标识**包括系统标识、SVN号、COSPAR-ID号、PRN四类：

1. **GNSS系统标识**：不同卫星导航系统的标识
   ——C:BDS
   ——G:GPS
   ——R:GLONASS
   ——E:GLONASS
2. **SVN**号：（Space Vehicle Number）空间飞行器编号，表示导航卫星的唯一编号
3. **COSPAR-ID**：（Committee on space research-ID）国际卫星标识符，用于命名、标识人造卫星，由两排数字与一排字母组成。第一排数字为该卫星的发射年，第二排数字为该卫星在其发射年的全球发射顺序，跟在第二排数字右侧的字母是在该次发射任务中分离出多个部分时用于标识每一部分;
4. **PRN号**：利用伪随机码标识导航卫星的编号

一般sNNN标识SVN号，sNN标识prn号

* 对于不同卫星系统，sNN一般从s01开始，常见的是我们说的prn号(不同卫星系统叫法不同)

——**GPS、BD**：the PRN number

——**GLONASS**：the Slot number

——**Galileo**：the SVID number

——**QZSS**：PRN number minus 192(.atx是这么定义的)

——**SBAS**：PRN number minus 100(.atx是这么定义的)

* sNNN一般表示SVN号

——GPS：SVN number

——GLONASS：GLONASS number

——Galileo：GSAT number

——QZSS：SVN number（183-202）

——SBAS：120-158

* 在RTKLIB中各卫星系统sNN及其number如下
* **注意**：下面number是在各卫星系统都预定义了的前提下，如果缺少某个卫星系统其number会减去对应卫星系统卫星数量

| 卫星系统 |        sNN        | number  |
| :------: | :---------------: | :-----: |
|   GPS    |  G01-G32**(32)**  |  1-32   |
| GLONASS  |  R01-R27**(27)**  |  33-59  |
| Galileo  |  E01-E36**(36)**  |  60-95  |
|   QZSS   | J193-J202**(10)** | 96-105  |
|   BDS    |  C01-C63**(63)**  | 106-168 |
|  IRNSS   |  I01-I14**(14)**  | 169-182 |
|   LEO    |  L01-L10**(10)**  | 183-192 |
|   SBAS   |  S20-S58**(39)**  | 193-231 |





[^卫星系统参考链接]:[System Introduction (csno-tarc.cn)](http://www.csno-tarc.cn/en/system/introduction)
