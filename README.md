J/ψ SPlot 分析与 Acceptance 计算

本项目用于 J/ψ → μ⁺μ⁻ 双微分截面分析中的 prompt J/ψ 信号提取与 Acceptance 计算。

 1. makeSPlot_converged.py
- 对mass和ctau进行二维拟合，区分 prompt J/ψ、non-prompt J/ψ 和 background。
- 使用 SPlot 方法提取 prompt J/ψ 的 sWeights。
- 构建 prompt J/ψ 的 (p_T)-(|y|) 参考分布。

输出： fitOutput_final.root，包含hPromptPtY、投影直方图及拟合检查图。

2. Acceptance_reweighted.py

- 使用 prompt J/ψ 的 (p_T)-(|y|) 参考分布对 Particle Gun MC 进行重加权。
- 根据生成级 J/ψ 信息及末态 μ 子的运动学要求，计算未加权和重加权的 Acceptance。
- 比较重加权前后的运动学分布及 Acceptance。

输出：Acceptance_reweighted.root、Acceptance_reweighted.txt及相关检查图。

分析流程
数据 → mass和ctau二维拟合 → SPlot → MC 重加权 → Acceptance 计算
