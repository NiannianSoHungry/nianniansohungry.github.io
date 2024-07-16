---
title: Limma 差异分析的一般流程
tags:
  - 生物信息学
  - bulk-seq
  - limma
categories:
  - 生物信息学
date: 2024-07-15 18:23:34
mathjax: true
cover: /img/Limma-差异分析的一般流程/picVolcano.jpeg
---

# Limma 差异分析的一般流程

## 目录

0. <a href="#0">环境配置</a>
1. <a href="#1">实验设计与数据标准化</a>
2. <a href="#1">limma 的标准流程</a>
3. <a href="#2">火山图的绘制</a>
4. <a href="#3">表达热图的绘制</a>

## 0. 环境配置

## 0.1 加载数据包

```R
library(tidyverse)
library(limma)
library(ggplo2)
library(pheatmap)
```

## 0.2 数据预览

在<a href="/post/RNAseq-数据预处理的常规流程">上一篇文章</a>中，我们已经将数据整理为了标准化后的表达矩阵

|       | HDMB03_shCTRL_1 | HDMB03_shCTRL_2 | HDMB03_shCTRL_3 | HDMB03_shCTRL_4 | HDMB03_shPPHLN1_1 | HDMB03_shPPHLN1_2 | HDMB03_shPPHLN1_3 | HDMB03_shPPHLN1_4 |
|-------|-----------------|-----------------|-----------------|-----------------|-------------------|-------------------|-------------------|-------------------|
|ARL17A |             84  |            90   |           88    |          70     |           98      |         103       |       102         |       78          |
|BAZ2B  |            694  |           730   |          674    |         734     |          672      |         930       |       778         |      752          |
|BOLA2  |              2  |             0   |            2    |           6     |            0      |           0       |         0         |        0          |
|CCNYL6 |              1  |             2   |            0    |           1     |            0      |           0       |         1         |        3          |
|CKMT1B |            490  |           726   |          641    |         476     |          554      |         622       |       524         |      396          |
|CRHR1  |             22  |            18   |           32    |          40     |           24      |          43       |        36         |       26          |

## 1. 实验设计与数据标准化

### 1.1 实验设计

所谓实验设计就是实验组与对照组的设置，通常 limma 处理的是单因素两组或多组的实验设计

```R
# 设置分组向量
group <- data.frame(
    group = rep(
        x = c("Ctrl", "PPHLN1"), # 2 组
        each = 4 # 各 4 个
    ),
    row.names = colnames(counts) # 行名即样本名
)

# 将分组向量转化为 dummy variables
design <- model.matrix(~0+group$group)

# 修改 design 表的行名与列名，与样本相对应
colnames(design) <- unique(group$group)
rownames(design) <- colnames(countsNorm)
```

### 1.2 数据标准化

之前已经提过，counts 数据符合负二项分布，而良好的差异分析的输入则要求符合正态分布，因此这里有三种方法进行标准化

#### 1.2.1 Voom

Voom 是 limma 自带的缩放方法，其本质是 $\log_{2}(CPM+0.5)$；

其中 CPM 的计算方法如下式, $N_{i}$ 代表第 i 个基因的表达量：

$$
CPM_{i} = \frac{N_{i}}{\sum{N_{i}}} * 10^6
$$

这相当于将某个基因的表达量转化为其在所有基因中的百分位数，实现将长拖尾数据转化为正态分布。

```R
countsNormVoom <- v0oom(
  counts = countsNorm,
  design = design
)
```

#### 1.2.2 Log2(counts+1)

Log2(counts+1) 就其形式而言与 voom 的公式如出一辙；

由于我们上一步已经进行过 normalizeBetweenArrays() ，因此可以认为 $\sum{N_i}$ 对于所有样本都是相等的，即测序深度归一化；

因此直接做 Log2(counts+1) 的结果与 voom 的结果仅相差一个比例倍数，其分布应当是完全一致的，计算差异倍数的结果也会是一致的。

```R
countsNormLog <- log2(countsNorm+1)
```

#### 1.2.3 TPM 或 FPKM

limma 理论上可以接受任何标准化方法的输入，因此通常不被 DESeq2 或 edgeR 所接受的 TPM 或 FPKM 数据会交给 limma 处理；

然而 limma 的官方文件也指出，与 DESeq2 和 edgeR 相同，其在处理 counts 数据时效果是最好的；

因此一般有 raw counts 数据时，不会转化为 TPM 或 FPKM；

因而此处只给出 TPM 与 FPKM 的公式，而不给出具体代码。

其中 $L_i$ 代表第 i 个基因转录本的长度，不特指时取最长转录本，若特指某个转录本则使用该转录本的长度

$$
TPM = \frac{\frac{N_i}{L_i} \times 10^6}{\sum{\frac{N_i}{L_i}}}
$$

$$
FPKM = \frac{N_i \times 10^6}{L_i \times \sum{N_i}}
$$

观察公式可知，TPM 计算时，求和符号中同时存在 $N_i$ 与 $L_i$， 因此同时对测序深度与基因长度进行了归一化；

而 FPKM 在计算时，求和符号中仅有 $N_i$，$L_i$ 在求和符号之外，因而其仅对测序深度进行了归一化，而未对基因长度进行归一化。

#### 1.2.4 标准化方法的选择

首先，在有 raw counts 数据的情况下，不选择 TPM 或 FPKM 标准化；

然后，在前期进行了 normalizeBetweenArrays() 的前提下， voom 和 log2(counts+1) 是等价的；

但是考虑到前期处理时可能不规范，没有进行 normalizeBetweenArrays()，此时选择 voom 相较之下更为保险；

之后的分析中也会主要使用 countsNormVoom。

## 2. Limma 的标准流程

我们现在已经准备好了 countsNormVoom 与 design 两个数据，可以进行 limma 的标准流程。

```R
# 建立对比矩阵
# 我们将对比设置为 PPHLN1-Ctrl，即 PPHLN1 为实验组， Ctrl 为对照组
# 最后得到的差异倍数也是 log(PPHLN1/Ctrl)
contrastsMatrix <- makeContrasts(
    contrasts = "PPHLN1-Ctrl",
    levels = design
)

# 线性拟合，通过拟合的线性函数代表基因的平均表达量
fit <- lmFit(
    object = countsNormVoom,
    design = design
)

# 计算倍数，输入拟合后的结果与对比矩阵，计算两组间的差异倍数
fit <- contrasts.fit(
    fit = fit,
    contrasts = contrastsMatrix
)

# 通过经验贝叶斯方法对计算结果的标准差进行同比的缩放，以缓和统计分析的结果，使得 t 检验，F 检验等的结果更显著
fit <- eBayes(fit)

# 从拟合结果中取出差异分析表
DEG <- topTable(
  fit,
  number = Inf
)
DEG$gene <- rownames(DEG) # 提取出基因列

# 根据矫正 P 值与差异倍数将差异基因分为上调或下调
DEG$regulation <- case_when(
    DEG$adj.P.Val < 0.05 & DEG$logFC > 1 ~ "up",
    DEG$adj.P.Val < 0.05 & DEG$logFC < -1 ~ "down",
    .default = "ns"
)

# 查看最终得到的差异基因的数量
table(DEG$regulation)
```

## 3. 火山图的绘制

差异基因通常通过火山图进行呈现；

此处实验设置只有两组，因此主要介绍两组火山图的绘制；

之后在讲解单细胞分析的时候会遇到多组间的差异基因，介时会讲解多组火山图的绘制。

### 3.1 认识火山图

1. 坐标轴：
   1. 横轴：log2FC，即差异倍数，越往右越上调，越往左越下调；
   2. 纵轴：-log10(adj.P.Value)，即矫正 P 值，越往上统计学意义越显著
2. 分隔线：
   1. $x= \pm 1$，标记了差异倍数为 1  的位置，在这两条线外侧的基因差异显著；
   2. $y = -\log_{10}{0.05}$，标记了矫正 P 值为 0.05 的位置，在这条线之上的基因具有统计学意义
3. 颜色：
   1. 红色：矫正 P 值 < 0.05 且 差异倍数 > 1，代表实验组相对与对照组上调的基因
   2. 蓝色：矫正 P 值 < 0.05 且 差异倍数 < -1，代表实验组相对与对照组下调的基因
   3. 灰色：不符合上述条件，代表无明显差异的基因

### 3.2 绘制

```R
picVolcano <- ggplot(DEG) +
    # 散点图，标记每个基因的位置
    geom_point(
        mapping = aes(
            x = logFC,
            y = -log10(adj.P.Val),
            color = regulation
        )
    ) +
    # 指定颜色
    scale_color_manual(
        values = c("blue", "gray", "red")
    ) +
    # 纵分隔线
    geom_vline(
        xintercept = c(-1, 1),
        lty = 5, # 虚线
        color = "gray"
    ) +
    # 横分隔线
    geom_hline(
        yintercept = -log10(0.05),
        lty = 5, # 虚线
        color = "gray"
    ) +
    # 横纵轴标签与标题
    labs(
        x = expression(log[2]*"(Fold Change)"),
        y = expression(-log[10]*"(Adjusted P Value)"),
        title = "PPHLN1 vs. Ctrl"
    ) +
    # 方框空白主题
    theme_bw() +
    # 其他主题设置
    theme(
        plot.title = element_text(
            face = "bold", # 标题加粗
            hjust = 0.5 # 标题居中
        ),
        legend.position = "bottom" # 图例置于底部
    )
print(picVolcano) # 画图
```

<img src="/img/Limma-差异分析的一般流程/picVolcano.jpeg">

## 4. 表达热图的绘制

热图表现了差异基因在不同样本间表达量的差异，通常选取与自己研究主题相关的基因进行展示；

此处因为仅做演示，因此选取高变 2000 个基因进行展示

```R
# 按照差异倍数正序与倒序各取前 1000 个基因，总共 2000 个
topVar2k <- rbind(
    DEG[order(DEG$logFC, decreasing = TRUE)[1:1000], ],
    DEG[order(DEG$logFC, decreasing = FALSE)[1:1000], ]
)$gene

# 取出这 2000 个基因的表达矩阵
dat <- countsNormVoom[topVar2k, ]

# 绘制热图
picHeatmap <- pheatmap(
    mat = dat,
    scale = "row", # 按行标准化
    color = colorRampPalette(c("blue", "white", "red"))(100), # 指定颜色为（蓝，白，红），过渡为 100 级
    treeheight_row = 0, # 基因方向上的聚类树高度设置为 0，因为此处聚类无意义，因而隐藏起来，
                        # 若在你的实验中需要看不同基因之间相关性的高低，可以不设置隐藏
    cutree_cols = 2, # 将样本分为两组，目的在于检查实验组与对照组间是否有混杂
    show_rownames = FALSE # 不显示基因名，因为此处展示了大量基因，若仅展示少量基因时，应当显示基因名
)
plot.new() # 新建画布
print(picHeatmap) # 画图
```

<img src="/img/Limma-差异分析的一般流程/picHeatmap.jpeg">