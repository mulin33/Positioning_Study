# day01 rtklib认识及编译

## 1.RTKLIB认识

目前有***2.4.2 p13***的**稳定**版本和***2.4.3***的**开发测试**版本，

 2.4.2 版本有PPP-AR部分，而2.4.3没有PPP-AR

下载下来的目录文件如下

![image-20231021211132637](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231021211132637.png)

## 2.RTKLIB在VS上编译

### 2.1 移植到VS项目

1. 创建**C++空项目**，可以勾选“解决方案和项目放在统一目录中”，记住创建的项目目录。

2. 把 RTKLIB 源码文件中**整个src文件夹**复制到创建的**项目**文件目录中。

3. 把 RTKLIB 源码文件中 **\app\consapp** 中的 **rnx2rtkp.c** 放到刚刚复制过去的 **src文件夹**，重命名为main.c。

4. 在解决源文件中添加名为 “**src**” 的筛选器，再在 **src** 筛选器下面添加名为 “**rcv**” 的筛选器 。

   右键添加现有项目把 **src/rcv文件夹** 中的所有文件加到 **src/rcv筛选器** 中，**src** 中所有代码文件加到 **src** 筛选器中。

### 2.2 修改bug

1. > **error C1083: 无法打开包括文件:“rtklib.h”: No such file ordirectory**
   >
   > **解决**：这是因为rcv里的文件找不到rtklib.h这个头文件。把在 **src/rcv文件夹几个的.c文件** 中的 **#include "rtklib.h"** 修改为 **#include "../rtklib.h”**

2. >  **error C1083：无法打开包括文件: "pthread.h": No such file or directory**
   >
   > **解决**：项目要在Windows下编译运行的，加**WIN32**，加了这一项后RTKLIB就不会用Linux下的<pthread.h>和<sys/select.h>

3. > **error C4703: 使用了可能未初始化的本地指针变量“sbs”**
   >
   > **解决**：对指针变量进行初始化，根据VS提示将ephemeris.c中`const sbssatp_t *sbs;`改为`const sbssatp_t *sbs = NULL;`

4. > **error C2466: 不能分配常量大小为 0 的数组**
   >
   > 根据报错发现没有预定义宏时数组大小为0，因此需加上`ENAGLO`,为了使代码能运行加上其他卫星系统的预定义

![image-20231021221656828](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231021221656828.png)

----------

***2-4bug修改如下***：项目属性-->C/C++-->预处理器-->预处理器定义-->粘贴下面

```
_CRT_SECURE_NO_WARNINGS
_WINSOCK_DEPRECATED_NO_WARNINGS             
WIN32
ENAGLO
ENAGAL
ENACMP
ENAQZS
ENAIRN
```

--------

5. > **error C4146: 一元负运算符应用于无符号类型，结果仍为无符号类型**
   >
   > **解决**：出现这个问题是因为勾选了VS的SDL检查[^SDL]，关闭SDL检查：项目属性-->C/C++-->SDL检查，选测否，就可以将其关闭了

6. > LNK2019 无法解析的外部符号：...，该符号在函数的...中被引用 
   >
   > **解决**：该函数引用了外部DLL文件，需要链接两个DLL文件,项目属性--> 链接器 --> 附加依赖项-->编辑-->粘贴如下
   >
   > ```
   > winmm.lib
   > ws2_32.lib
   > ```

###  2.3 生成项目

出现如下结果即生成成功

![生成项目](https://raw.githubusercontent.com/mulin33/ImageHost/main/blogImg/image-20231022104302770.png)

## 3.procpos时序图



```mermaid
%% ppp_procpos后面的时序图
sequenceDiagram
autonumber

participant p as postpos.c
participant r as rtkpos.c
participant spp as pntpos.c
participant e as ephemeris.c
participant sp3 as preceph.c
participant ppp as ppp.c


p->>+p:procpos(mode=0即froward)
p->>+r:rtkinit(rtk,popt)
r-->>-p:初始化rtk_t结构体

%% 开始单历元解算loop循环
rect rgb(240 248 255)
loop while(nobs=inputobs()>=0)

p->>p:inputobs输出一个历元观测数据obsd_t obs

p->>+r:rtkpos(rtk,obs,nobs,navs)<br/>ok
%%pntpos()单历元单点定位
r->>+spp:pntpos()

%%satposs计算一个历元的卫星位置
spp->>+e:satposs()<br/>计算卫星pos/v/clk/drift/var
rect rgb(255 255 224)
loop 循环一个历元观测卫星(for i=0,i<nobs,i++)
	note over e:获得信号发射钟面时间t1(接收机时间-伪距/光速)，包含tr和ts
	e->>e:ephclk()由钟面时间从广播星历获得卫星钟差dt(not 相对论和tgd改正)
	note over e:获得信号发射时间真值t2(钟面时-钟差dt)，包含tr
	e->>+e:satpos()由t2，switch(n/sp3)计算卫星位置和钟差
	
	rect rgb(255 228 225)
	%%广播星历n获得卫星位置
	alt prcopt.sateph=EPHOPT_BRDC
		e->>+e:ephpos()先调用seleph()匹配对应星历eph
		e->>-e:后调用eph2pos()由n计算卫星位置和钟差(相对论改正 not tgd)
		note over e:扰动法计算卫星速度和钟漂:t2+0.001s再次调用eph2pos获得新的位置和钟差
		e->>-e:ephpos()返回从广播星历获得卫星pos/v/bias(去除相对论包含tgd)/drift/var
	%%精密星历sp3获得卫星位置	
	else prcopt.sateph=EPHOPT_PREC
		e->>+sp3:peph2pos()从精密星历获得卫星pos/clk
		sp3->>sp3:pephpos()Neville插值10个历元获得pos
		sp3->>sp3:pephclk()获得clk
		sp3->>sp3:satantoff()从COM->Ant phase Center
		sp3->>sp3:去除卫星时钟相对论效应sat clock {bias,drift}
		sp3-->>-e:return peph2pos(1)
	end
	end
	e-->>-spp:return satposs(void)
end
end

%%estpos()解算接收机pos和dtr
spp->>+spp:estpos()lsq伪距单点定位
rect rgb(255 255 240)
loop 迭代lsq10次，for(i<MAXITR,i++)
%%estpos()调用rescode()
spp->>+spp:rescode()pseudorange residuals
rect rgb(255 228 225)
loop 循环每一个卫星
	spp->>spp:satexlude()剔除不符合条件的卫星
	spp->>spp:geodist()消除earth rotation得卫地距离和r->s单位向量(构建H矩阵)
	spp->>spp:satazel()计算方位角
	
	spp->>spp:prange()伪距硬件延迟改正(dcb and tgd)
end
end
spp-->>-spp:return:rescode()

spp-->>-spp:return:estpos()
end
end

spp->>spp:raim_fde()接收机自主完好性检验
spp-->>-r:return pntpos(1)
%%ppp定位
r->>+ppp:pppos(),if pntpos返回1，才会进行pppos()
ppp-->>-r:return pppos()


r-->>-p:return rtkpos
end
end
p->>-p: return end





```







-------

[^SDL]: Security  Development Lifecycle 安全开发生命周期。这个是为了检测开发者的代码安全，如果我们不去掉它，它会按照SDL的规范来编译代码，就会有一些以前常用的东西报错。比如**scanf**函数，可以将其变成 **scanf_s**函数 就可以使用了。