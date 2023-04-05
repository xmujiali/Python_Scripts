# 一、拓扑分析
通过以下命令
```
0-topo.py
```
可以获得各异构体的顶点、边、环类型及数目。

# 二、半经验和DFT的几何优化
## 1、选择需要的异构体并进行半经验/DFT计算

运行
```
1-filter.py
```
它首先会打印出所有可选的筛选方案（ALL全选/GT某字段大于阈值/LT某字段小于阈值）
- 如果选择ALL，则考虑所有异构体
- 如果选择GT或者LT，则程序首先打印csv的数据的列名，待用户选择某个列后，再打印该列的统计摘要并提示用户输入阈值进行筛选。
  
最后，将任务分解到数个bash文件分别运行，即可
```
misc-partition_jobs.py
bash run_partition_x.sh &
```


示例：
假定
- 用户已经使用半经验计算完成了对所有异构体的几何优化，chk文件在./semi_opt/C_28_异构体编号.chk
- 用户当前目录下有一个result.csv
- 希望能够提取所有PM6_beta大于2000的异构体先进行DFT水平的几何优化，再计算非线性光学
- 后续计算希望放到./dft_opt目录和./dft_nlo下)
预先在准备./samples/dft_opt.gjf和./samples/dft_nlo.gjf，内容如下：
```
%oldchk=../semi_opt/C_28_(isomer).chk
%chk=C_28_(isomer).chk
%mem=4000MB
%CPU=8
#p M062X/6-31G* scf=(xqc,tight,maxcycle=512) opt=tight geom=check 

M062X/6-31G* optimization for isomer no. (isomer)

0   1

```
```
%oldchk=../dft_opt/C_28_(isomer).chk
%chk=C_28_(isomer).chk
%mem=4000MB
%CPU=8
#p M062X/6-31G* scf=(xqc,tight,maxcycle=512) opt=tight geom=check guess=read polar

M062X/6-31G* nlo computation for isomer no. (isomer)

0   1
```
操作如下：
```
[jli@chem174 C28] $ 1-filter.py
Choose a filter from following:

--------------------
key =  all
choose all isomers.
--------------------
key =  GT
returns isomers with value larger than you give.
--------------------
key =  LT
returns isomers with value less than you give.
--------------------
Which one do you want to use: GT
Give men csv filename[Enter for ./result.csv]:
Index(['PM6_opt', 'PM6_beta'], dtype='object')
col_name = PM6_beta
count      7903.000000
mean        553.137543
std        2954.252357
min           0.000000
25%         221.646216
50%         363.318116
75%         587.595349
max      230853.643861
Name: PM6_beta, dtype: float64
Only isomers with value larger than you give will be considered.
value = 2000
Give me the template file for gaussian inputs:
samples/dft_opt.gjf
Input the directory name where computational jobs located.
dft_opt
Job done!

```

```

```

## 2、检查优化是否全部成功(半经验/DFT)
```
2-check_opt.py
```
如果全部成功会给出提示。如果有部分失败的任务也会记录到failed_jobs和redo_jobs.sh中，手工修改输入文件（提示：可能CAGE给出的结构本身是错误的！需要自行根据邻接表重新连接原子。），然后运行
```
bash redo_jobs.sh
```

## 3、统计能量(半经验/DFT)

```
3-sum_opt.py
```
将为你统计相关的半经验计算结果






## 统计NLO计算结果

```
4-sum_nlo.py
```


