---
title: RNAseq 数据预处理的常规流程
date: 2024-07-12 14:54:56
tags: 
    - 生物信息学
    - bulk seq
categories:
    - 生物信息学
cover: /img/RNAseq_数据预处理的常规流程/picPcaAfter.jpeg
mathjax: true
---

# RNAseq 数据预处理的常规流程

## 目录

0. <a href='#0'>环境设置</a>
1. <a href='#1'>数据展示</a>
2. <a href='#2'>数据质量检测</a>
3. <a href='#3'>过滤低表达基因</a>
4. <a href='#4'>数据标准化</a>
5. <a href='#5'>去批次</a>
6. <a href='#6'>处理后质量检测</a>

## <div id='0'>0. 环境设置</div>
读入本次实验需要用到的数据包。

```R
library(tidyverse) # 格式整理
library(data.table) # 读入文件
library(clusterProfiler) # Gene ID 转换为 Gene Symbol
library(org.Hs.eg.db) # 人类基因组注释文件
library(FactoMineR) # 计算 PCA
library(factoextra) # PCA 画图
library(pheatmap) # 热图
library(ggplot2) # 画图
library(reshape2) # 长宽表转换
library(edgeR) # 去除低表达基因 
```

## <div id='1'>1. 数据展示</div>

此次演示中使用来自 GSE222699 的 RNAseq 数据，数据从 GEO 官网下载。

```tree
GSE222699_RAW
├─ GSM6929162_HDMB03_CTRL_1.txt.gz
├─ GSM6929163_HDMB03_CTRL_2.txt.gz
├─ GSM6929164_HDMB03_CTRL_3.txt.gz
├─ GSM6929165_HDMB03_CTRL_4.txt.gz
├─ GSM6929166_HDMB03_CTRL_1.txt.gz
├─ GSM6929167_HDMB03_CTRL_2.txt.gz
├─ GSM6929168_HDMB03_CTRL_3.txt.gz
└─ GSM6929169_HDMB03_CTRL_4.txt.gz
```

### 1.1 读入数据

这组数据并没有将所有样本整合为一个表达矩阵，而是每个样本各有一个文件，因此我们首先整理表达矩阵。

```R
# 将所有文件列入一个列表
files <- list.files(
    path = "data/GSE222699_RAW/",
    full.names = TRUE
)

# 读入列表中的第一个文件作为模板
# 因为文件属于 .gz 压缩文件，直接用 fread 读入后再转换为 data.frame，这是一种更方便的操作
counts <- fread(files[[1]]) %>%
    as.data.frame()

# 从第二个文件开始，读入后与前几个文件合并
for (file in files[2:length(files)]) {
    tmp <- fread(file) %>%
        as.data.frame()
    counts <- merge(counts, tmp, by = "V1")
}

# 查看合并后的表达矩阵
View(counts)
```

|   |        V1          | HDMB03_shCTRL_1 | HDMB03_shCTRL_2 | HDMB03_shCTRL_3 | HDMB03_shCTRL_4 | HDMB03_shPPHLN1_1 | HDMB03_shPPHLN1_2 | HDMB03_shPPHLN1_3 | HDMB03_shPPHLN1_4 |
|---|--------------------|-----------------|-----------------|-----------------|-----------------|-------------------|-------------------|-------------------|--------|
| 1 | ENSG00000000003    |        1041     |        966      |       828       |     1053        |    1334           |   1563            |  1157             | 1258   |
| 2 | ENSG00000000005    |           0     |          0      |         0       |        0        |       0           |      0            |     0             |    0   |
| 3 | ENSG00000000419    |        1118     |       1382      |       999       |     1414        |    1140           |   1648            |  1067             | 1168   |
| 4 | ENSG00000000457    |         315     |        394      |       369       |      344        |     482           |    595            |   435             |  388   |
| 5 | ENSG00000000460    |         447     |        436      |       439       |      411        |     479           |    635            |   496             |  339   |
| 6 | ENSG00000000938    |           0     |          0      |         0       |        0        |       0           |      0            |     0             |    0   |

### 1.2 Gene ID 转换为 Gene Symbol

可以看到，合并后的矩阵以 Gene ID 表示基因，这不利于我们后续对基因的阅读，需要转换为方便阅读的 Gene Symbol。

```R
# 通过 clusterProfiler 中的 bitr() 函数将 Gene ID 与 Symbol 一一对应
# drop = FALSE 是指无法 mapping 到 Symbol 的 Gene ID 依然保留在表中，方便后续的合并
# 如果 drop = TRUE，那么无法 mapping 的 Gene ID 会被舍去，id2symbol 的顺序就会和 counts 不同，不利于合并
id2symbol <- bitr(
    geneID = counts$V1,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    drop = FALSE
)

# 消除对应到多个 Symbol 的 Gene ID，只保留 mapping 到的第一个 Symbol
# duplicated() 代表重复的元素，对其取反 !duplicated() 代表不重复的元素
id2symbol <- id2symbol[which(!duplicated(id2symbol$ENSEMBL)), ]

# 比对 id2symbol 与 counts 是否一一对应
table(id2symbol$ENSEMBL == counts$V1) # TRUE

# 将 Gene Symbol 合并入 counts
counts$Symbol <- id2symbol$SYMBOL

# 去除没有 mapping 到的 Gene ID
# na.omit() 去除包含 NA 的行
counts <- na.omit(counts)

# 查看 ID 转换后的表达矩阵
View(counts)
```

|   |        V1          | HDMB03_shCTRL_1 | HDMB03_shCTRL_2 | HDMB03_shCTRL_3 | HDMB03_shCTRL_4 | HDMB03_shPPHLN1_1 | HDMB03_shPPHLN1_2 | HDMB03_shPPHLN1_3 | HDMB03_shPPHLN1_4 | Symbol |
|---|--------------------|-----------------|-----------------|-----------------|-----------------|-------------------|-------------------|-------------------|--------|-------|
| 1 | ENSG00000000003    |        1041     |        966      |       828       |     1053        |    1334           |   1563            |  1157             | 1258   |TSPAN6 |
| 2 | ENSG00000000005    |           0     |          0      |         0       |        0        |       0           |      0            |     0             |    0   |  TNMD |
| 3 | ENSG00000000419    |        1118     |       1382      |       999       |     1414        |    1140           |   1648            |  1067             | 1168   |  DPM1 |
| 4 | ENSG00000000457    |         315     |        394      |       369       |      344        |     482           |    595            |   435             |  388   | SCYL3 |
| 5 | ENSG00000000460    |         447     |        436      |       439       |      411        |     479           |    635            |   496             |  339   | FIRRM |
| 6 | ENSG00000000938    |           0     |          0      |         0       |        0        |       0           |      0            |     0             |    0   |   FGR |

### 1.3 消除重复的基因

多个 Gene ID 常常对应到同一个 Gene Symbol，导致矩阵中出现重复的基因名，需要想办法去除。

```R
# 去除 Gene ID 这一列
# [, -1] 代表去除第一列
counts <- counts[, -1]

# 取出重复的基因与不重复的基因
# 查看有多少重复的基因
table(duplicated(counts$Symbol))

# 取出重复基因所在的行号
dupl <- duplicated(counts$Symbol)

# 取出重复的基因名
geneDupl <- counts$Symbol[dupl]

# 取出重复的基因所在的行
# 这边不直接使用 duplicated() 函数，是因为 duplicated() 函数不会取出第一项重复的元素
# 比如如果一个元素共重复了 5 次， duplicated() 只会取出后 4 次，而忽略第 1 次
countsDupl <- subset(
    x = counts,
    subset = Symbol %in% geneDupl
)

# 不在 geneDupl 中的基因则认为是 unique 的基因
countsUniq <- subset(
    x = counts,
    subset = !(Symbol %in% geneDupl)
)

# 合并重复的基因
# 合并的方法有很多种，包括：均数 mean，中位数 median，最大值 max
# 因为此处的数据类型为 Raw counts，为了尽量保护其整数性，此处选择使用 median
countsDupl <- aggregate(
    x = .~Symbol,
    data = countsDupl,
    FUN = median
)

# 将合并好的重复基因与非重复基因合并，得到完整的表达矩阵
counts <- rbind(
    countsDupl, countsUniq,
    deparse.level = 2
)

# 将 Gene Symbol 转换为行名
counts <- counts %>%
    remove_rownames() %>%
    column_to_rownames(var = "Symbol")

# Counts 取整数
# 即使是使用 median 处理重复基因，也难以避免 .5 的出现，因此依然需要四舍五入 round()
counts <- round(counts)

# 查看处理后的表达矩阵
View(counts)
```

|       | HDMB03_shCTRL_1 | HDMB03_shCTRL_2 | HDMB03_shCTRL_3 | HDMB03_shCTRL_4 | HDMB03_shPPHLN1_1 | HDMB03_shPPHLN1_2 | HDMB03_shPPHLN1_3 | HDMB03_shPPHLN1_4 |
|-------|-----------------|-----------------|-----------------|-----------------|-------------------|-------------------|-------------------|-------------------|
|ARL17A |             84  |            90   |           88    |          70     |           98      |         103       |       102         |       78          |
|BAZ2B  |            694  |           730   |          674    |         734     |          672      |         930       |       778         |      752          |
|BOLA2  |              2  |             0   |            2    |           6     |            0      |           0       |         0         |        0          |
|CCNYL6 |              1  |             2   |            0    |           1     |            0      |           0       |         1         |        3          |
|CKMT1B |            490  |           726   |          641    |         476     |          554      |         622       |       524         |      396          |
|CRHR1  |             22  |            18   |           32    |          40     |           24      |          43       |        36         |       26          |

## <div id='2'>2. 数据质量检测</div>

### 2.0 分组

该数据共分为 2 组，4 个对照组 Ctrl 与 4 个实验组 PPHLN1。

```R
group <- data.frame(
    group = rep(
        x = c("Ctrl", "PPHLN1"), # 2 组
        each = 4 # 各 4 个
    ),
    row.names = colnames(counts) # 行名即样本名
)
```

### 2.1 PCA

PCA （主成分分析）是一种数据降维方法。
1. 首先，所有的数据以其每个基因为基底，张成一个高维空间，而每个样本在这个高维空间中具有一个唯一的坐标 (x1, x2, ... , xn)：
    1. 计算高维空间中各个方向上样本间的方差，取方差最大的方向上的向量作为其主成分向量，记作 PC_1；
    2. 显然，在 PC_1 上，样本的分布最为分散，不同样本间的区分度也就最高；
2. 在垂直于 PC_1 的各个方向上，取方差最大的方向上的向量，记作 PC_2：
   1. PC_1 与 PC_2 共同张成一个正交的二维平面；
   2. 显然，在这一平面上，各个样本间的分布最为分散，不同样本间的区分度也就最高；
3. 在由 PC_1 与 PC_2 张成的二维平面上考察不同分组的样本的分布模式，理想的情况下：
   1. 处于同一分组的样本之间的位置最为接近；
   2. 处于不同分组的样本之间的位置相距最远。



``` R
pcaBefore <- PCA(
    X = t(counts),
    graph = FALSE # 计算时不画图
)
picPcaBefore <- fviz_pca_ind(
    X = pcaBefore,
    col.ind = group$group, # 规定各组的颜色
    repel = TRUE, # 文字标签之间不重叠
    addEllipses = TRUE # 加入椭圆形，每组不多于 3 个样本时无法加入椭圆
                       # （五点确定一个椭圆，因此至少需要五个自由度，而对于不精确的拟合，可以只用 4 个点）
)
print(picPcaBefore)
```

<img src="/img/RNAseq_数据预处理的常规流程/picPcaBefore.jpeg">

### 2.2 相关性热图
计算各个样本之间的相关性，以热图的形式展示，理想情况下：
1. 同属一组的样本之间相关性最高；
2. 属于不同组的样本之间相关性最低。

```R
picHeatmapBefore <- pheatmap(
    mat = counts %>% cor(), # 计算相关性矩阵
    annotation_row = group, # 行色块标签，表示分组
    annotation_col = group, # 列色块标签，表示分组
    legend = TRUE, # 显示色条图注
    annotation_legend = FALSE, # 不显示分组图注
    cutree_rows = 2, # 根据树状图分为两行
    cutree_cols = 2, # 根据树状图分为两列
    color = colorRampPalette(c("blue", "white", "red"))(100) # 规定颜色为 (蓝, 白, 红)，过渡为 100
)
plot.new() # 新建画布
print(picHeatmapBefore) # 画图
```

<img src="/img/RNAseq_数据预处理的常规流程/picHeatmapBefore.jpeg">

### 2.3 箱线图
通过箱线图表示各个样本的分布，以考察各个样本是否属于同分布；

若样本间的分布不同，不利于后续的统计学分析，需要进行标准化；

由于 Raw counts 符合负二项分布，因此需要取 log2(counts+1) 使其接近正态分布。

```R
picBoxplotBefore <- (counts+1) %>% log2() %>% # 取 log2(counts+1)
    melt() %>% # 将宽表形式的表达矩阵转换为长表形式，便于画图
    ggplot + # 画图开始
    geom_boxplot( # 画箱线图
        mapping = aes(
            x = variable, # x 轴为样本
            y = value # y 轴为基因表达量
        )
    ) +
    labs(x = "Sample", y = "Expression") + # 指定 x 轴与 y 轴的标题
    theme_bw() + # 方框空白主题
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1) # x 轴标签旋转 45 度并右对齐
    )
plot.new() # 新建画布
print(picBoxplotBefore) # 画图
```

<img src="/img/RNAseq_数据预处理的常规流程/picBoxplotBefore.jpeg">

## <div id='3'>3. 过滤低表达基因</div>

观察到箱线图中的中位数线位于 0，说明存在大量低表达甚至不表达的基因；

考虑到 $\frac{0}{0}$ 无意义，会影响后续的差异倍数的计算；

因此需要去除低表达的基因；

此处使用 edgeR 包中的 filterByExpr() 函数。

```R
# 过滤低表达基因
keep <- filterByExpr(
    counts,
    group = group$group
)

# 取出过滤后仍保留的行
countsRLE <- counts[keep, ]

# 制图
picBoxplotRLE <- (countsRLE+1) %>% log2() %>% # log2(counts+1)
    melt() %>%
    ggplot +
    geom_boxplot(
        mapping = aes(
            x = variable,
            y = value
        )
    ) +
    labs(x = "Sample", y = "Expression") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# 画图
print(picBoxplotRLE)
```

<img src="/img/RNAseq_数据预处理的常规流程/picBoxplotRLE.jpeg">

经过过滤低表达基因后，可以看到，箱线图中的中位数线位于箱子的正中间，说明现在的表达量分布接近正态分布。

## <div id='4'>4. 数据标准化</div>

观察经过过滤后的箱线图，可以发现，虽然大体上的分布趋势相似，但是箱子的高度仍然有细微的差别；

我们可以使用 limma 包中的 normalizeBetweenArrays() 函数进行优化

```R
countsNorm <- normalizeBetweenArrays(countsRLE)
picBoxplotAfter <- (countsNorm+1) %>% log2() %>%
    melt() %>%
    ggplot +
    geom_boxplot(
        mapping = aes(
            x = Var2,
            y = value
        )
    ) +
    labs(x = "Sample", y = "Expression") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
print(picBoxplotAfter)
```

<img src="/img/RNAseq_数据预处理的常规流程/picBoxplotAfter.jpeg">

可以看到，经过标准化后的数据分布完全一致，有利于我们下一步的分析

## <div id='5'>5. 去批次</div>

此次实验中的样本全部来自一次测序结果，因此不存在批次效应；

若实验中的样本是多次测序汇总的结果，则不同的测序批次之间存在批次效应，需要去除；

去除批次效应可以使用 sva 包中的 ComBat() 函数以及 limma 包中的 removeBatchEffect() 函数；

因为后续的差异分析主要基于 limma 包，因此此处也主要介绍 limma 包中的方法。

```R
# 将分组信息转化为实验设计矩阵
design <- model.matrix(~0+group$group)

# 提供批次信息
batch <- data.frame(
    batch = c(...) # 此处填入批次信息
    row.names = col.names(countsNorm) # 行名来自样本名
)

# 去批次
countsRBE <- removeBatchEffect(
    counts,
    batch = batch$batch, # 去除批次效应
    design = design # 保留生物学效应
)
```

去批次后，可以再使用 PCA 方法验证去批次的效果。

## <div id='6'>6. 处理后质量检测</div>

同样使用 PCA、相关性热图、箱线图三种方法；

代码不赘述

<img src="/img/RNAseq_数据预处理的常规流程/picPcaAfter.jpeg">

<img src="/img/RNAseq_数据预处理的常规流程/picHeatmapAfter.jpeg">

<img src="/img/RNAseq_数据预处理的常规流程/picBoxplotAfter.jpeg">