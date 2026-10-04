# scMEDAL 模型讲解（组会版）

## 1. 我们想解决什么问题

单细胞 RNA-seq 可以看成一个“细胞 × 基因”的表达矩阵。同一种细胞之间的差异，既可能来自真实生物学，例如细胞类型、疾病状态，也可能来自 donor、样本批次或实验技术。传统 batch correction 主要尝试把 batch 信息去掉；问题是 donor 与疾病状态可能纠缠在一起，过度校正可能连有意义的生物学信号也一起删掉。

scMEDAL 的核心想法不是只问“怎样去掉 batch”，而是把问题拆成两个互补部分：

- 哪些信息相对独立于 batch？
- 哪些信息是 batch/donor-specific，而且可能仍包含生物学意义？

## 2. 两个独立训练的子网络

输入包括每个细胞的基因表达向量 `x` 和 one-hot batch label `z`。

### scMEDAL-FE：学习 batch-invariant 表示

FE 是 conventional autoencoder 加 adversarial batch classifier：

- Autoencoder 需要从 latent representation 重建原始表达，因此不能丢掉所有信息。
- Adversarial classifier 尝试根据 encoder features 预测 batch。
- FE 的优化反过来惩罚能够预测 batch 的特征。

直观理解：FE 一边要“记住这个细胞是什么”，一边要“尽量忘记它来自哪个 donor/batch”。

### scMEDAL-RE：学习 batch-specific 表示

RE 是 Bayesian autoencoder：

- batch label 会进入网络的多层结构；
- Bayesian weights/variational inference 用分布而不是单一点估计表示不确定性；
- batch classifier 正向鼓励 latent space 保留 batch-specific information；
- KL divergence 对权重分布进行正则化，降低过拟合。

直观理解：RE 专门回答“这个 donor/batch 对表达模式造成了什么系统性偏移”。

FE 和 RE 不是两个互相竞争的 correction models。论文实际把它们分开训练，之后才把两个 latent representations 拼接给下游 classifier。模型受 mixed-effects 思想启发，但不是传统线性混合模型的联合估计。

## 3. Counterfactual reconstruction

RE 可以固定一个细胞的 latent state，只替换 decoder 使用的 batch label：

> 如果模型把这个细胞放到另一个 donor/batch 条件下，它会重建出怎样的表达？

这可以比较同一细胞在不同 batch 条件下的模型预测变化。结果是模型生成的 as-if reconstruction，表示给定替换后的 donor/batch 标签时，模型预测的表达模式。

## 4. 我们在 AML 上做了什么

- 使用作者处理好的 AML 数据，而不是从 FASTQ 开始。
- 2,916 HVGs、19 donor/batches、21 个配置的 cell-type categories。
- 使用作者提供的五折划分。
- 训练 scMEDAL-FE、scMEDAL-RE、Harmony、Scanorama、scVI、scANVI。
- Input PCA 作为未经整合的参考。
- 正式神经网络训练最多 500 epochs，并使用 validation-based early stopping，patience 30。
- 使用五折测试集中的 Mono 和 Mono-like 细胞做 RE counterfactual 分析；每折把同一组源细胞放入 19 种目标 donor/batch 条件。
- 完成五折基因差异方向、代表性 donor 图和逐一排除 control donor 的敏感性检查。

## 5. 怎么比较才公平

统一使用相同的 cells、genes、metadata、folds、50-dimensional latent space 和评价代码。

- Batch correction methods：希望 batch clustering 较低，同时 cell-type structure 保留。
- scMEDAL-RE：希望 batch/donor-specific structure 清楚，因此较高 batch clustering 是设计目标。
- ASW 是论文主要指标；CH 和 1/DB 提供补充视角。
- UMAP 只用于定性观察。不同模型分别计算的 UMAP 坐标方向、旋转和距离不能直接逐点比较。

## 6. Baselines 的角色

- Harmony：对低维表示迭代校正，使不同 batches 更好混合。
- Scanorama：利用跨数据集的匹配关系进行整合。
- scVI：无监督 probabilistic deep generative model。
- scANVI：在 scVI 基础上加入 cell-type labels 的半监督模型。

这些方法主要提供 batch-invariant 表示；scMEDAL-RE 的定位是提供互补的 batch-aware 信息，而不只是替代它们。

## 7. 当前能说与不能说的结论

可以说：

- AML 上六种方法的正式五折训练和产物检查已经完成。
- FE latent 的解释目标是 batch mixing 与 biological preservation。
- RE latent 的解释目标是 donor/batch modeling。
- 完整 counterfactual 输出链路已经跑通。
- SAMSN1 的模型重构差异在五折中方向一致，并且五折都进入差异幅度前 30；结果用于展示目标 donor/batch 相关模式。

暂时不能说：

- “scMEDAL 在所有数据集上最好”：我们目前只正式复现了 AML。
- “RE 的 batch separation 高，所以 correction 更差”：RE 本来就不是 correction latent。
- “每个基因都达到统计显著”：fold 1 的 2,916 个基因检验中，FDR < 0.05 的数量为 0。
- “只看一张 UMAP 就证明模型优越”：最终需要五折指标、置信区间和下游任务共同支持。

## 8. 一分钟口头版本

我们研究的是单细胞数据里的 batch effect。donor 和疾病状态可能纠缠在一起，因此需要同时观察共享的细胞状态与 donor/batch 特异变化。scMEDAL 训练两个互补子网络：FE 用 adversarial autoencoder 学习尽量独立于 batch 的主要细胞状态；RE 用 Bayesian autoencoder 学习 donor/batch-specific variation。两者独立训练，后续可以组合用于预测。RE 固定一个细胞的 latent state、替换 batch label，生成模型在另一目标条件下的 as-if 表达重建。我们用作者处理好的 AML 数据完成了 FE、RE 和四种比较方法的正式五折运行，并用五折指标、UMAP、counterfactual 基因效应图及 donor 图解释模型输出。SAMSN1 在五折目标条件比较中方向一致，是一个适合展示的模型预测模式。
