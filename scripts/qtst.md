# QTST 脚本系列说明

本系列脚本用于基于量子过渡态理论（QTST）计算反应速率、量子隧穿校正等物理量，适用于以氢及其同位素迁移为主的反应体系。
目录结构

```code
scripts/
├── qtst_read_freq.py        # 读取VASP输出的频率，生成YAML格式数据
├── qtst_calc_rate.py        # 基于频率文件计算TST速率常数及量子校正
├── qtst_tunneling.py        # 计算量子隧穿校正因子
├── qtst_input.sh            # 生成输入模板
├── qtst_input.yaml          # 输入参数样例
├── qtst_tunning.yaml        # 隧穿校正输出样例
```

## 一、脚本功能与使用流程

1. 频率提取（`qtst_read_freq.py`）

    - 功能：自动识别当前目录下 stable、saddle 等文件夹，解析 VASP 的 OUTCAR 文件，提取实部/虚部频率，保存为 {name}_f.yaml 和 {name}_fi.yaml。

    - 用法：

    ```bash
    python qtst_read_freq.py
    ```

    - 输出：每个稳定态/过渡态目录下生成 xxx_f.yaml（实频）和 xxx_fi.yaml（虚频）。

2. 速率与量子校正计算（`qtst_calc_rate.py`）

    - 功能：读取频率数据，计算零点能（ZPE）、TST下的速率常数和量子校正因子。

    - 用法：

    ```bash
    python qtst_calc_rate.py --input qtst_input.yaml
    ```

    - 输入：`qtst_input.yaml`（温度等参数）、`_f.yaml/_fi.yaml`（频率文件，自动检索）。

    - 输出：`qtst_rates.yaml`，内容示例：

    ```YAML
    - stable: stable
      saddle: saddle
      zpe_correction (eV): 0.123
      k0 (Classical TST rate constant (THz)): 4.56
      quantum_corrections:
        200: 0.98
        300: 1.02
        ...
    ```

3. 量子隧穿校正（qtst_tunneling.py）

    - 功能：基于反应势垒、特征虚频等参数，数值积分计算量子隧穿校正因子（适用于氢/同位素）。

    - 用法：

    ```bash
    python qtst_tunneling.py --input qtst_input.yaml
    ```

    - 输入：qtst_input.yaml，内容示例：

    ```YAML
    deltaE0: 0.65     # 势垒高度 (eV)
    i_omega: 46.41    # 虚频 (THz)
    #down_limit: -100 # 可选，积分下限
    temperatures: [200, 300, 400, 500, 600]
    #output_file: qtst_tunning.yaml
    ```

    - 输出：qtst_tunning.yaml，包含各温度下的隧穿校正因子、积分误差等。

4. 输入模板生成（`qtst_input.sh`）

    - 功能：快速生成 qtst_input.yaml 模板，便于批量参数测试。

    ```bash
    bash qtst_input.sh
    ```

## 二、输入与输出说明

### 输入

- 频率数据：通过 VASP 计算产生的 OUTCAR 文件，须放置于名为 stable、saddle 等文件夹下。
- 参数文件：qtst_input.yaml，包含势垒、虚频、温度等参数。

### 输出

- `{name}_f.yaml` / `{name}_fi.yaml`：实部/虚部频率列表
- `qtst_rates.yaml`：TST速率、ZPE校正及量子修正结果
- `qtst_tunning.yaml`：量子隧穿校正因子及积分误差

## 三、注意事项

- 仅适用于过渡态主导虚频明确、单一反应坐标的体系，主要为轻原子（如氢）隧穿过程。
- 输入参数单位需严格对应（如能垒为 eV，虚频为 THz）。
- 当出现积分溢出警告时，可适当调整 down_limit 参数。
- 建议结合文献和实际体系特征，合理选择输入参数。

## 四、参考

[1] J. T. Fermann and S. Auerbach, J. Chem. Phys. 112, 6787 (2000).
