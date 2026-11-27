# motion-control-parameters

PMSM（永磁同步电机）控制参数整定工具。输入电机参数、控制器采样设置和原理图硬件参数，自动计算电流环、速度环、位置环、观测器参数和 ADC/PWM 硬件换算，并在离散域 (z 域) 校核。

## 安装与使用

应用在 `app/` 目录，是一个 PWA（可安装的网页应用），Windows、iPhone、Android 共用一套代码，首次打开后可离线使用。

| 平台 | 方法 |
|---|---|
| Windows（最简单） | 直接双击 `app/index.html` 用 Edge/Chrome 打开，完全离线可用 |
| Windows（安装为应用） | 用 Edge/Chrome 打开托管地址 → 地址栏右侧“安装应用”，之后从开始菜单启动，独立窗口 |
| iPhone / iPad | 用 Safari 打开托管地址 → 分享按钮 →“添加到主屏幕” |
| Android | 用 Chrome 打开托管地址 → 菜单 →“安装应用” |

托管地址：仓库启用 GitHub Pages 后为 `https://<用户名>.github.io/motion-control-parameters/`。启用方法：仓库 Settings → Pages → Source 选 **GitHub Actions**，之后每次 `app/` 合并到 `main` 时由 `.github/workflows/pages.yml` 自动发布。（私有仓库使用 Pages 需要 GitHub 付费计划。）

本地调试：`cd app && python3 -m http.server 8000`，浏览器打开 `http://localhost:8000`。

## 功能

- **电机**：SPMSM / IPMSM，Rs、Ld、Lq、磁链（ψf、ψa、Ke、Kt 互换）、极对数、惯量、摩擦、额定/极限值；坐标变换可选等幅值或等功率（瑞萨 absolute）。
- **控制器**：开关频率、单采单更/双采双更、采样到占空比生效的延时（可填小数拍，修正 z 变换精确建模）、PI 实现（瑞萨 motor_pi_ctrl 位置式或增量式）、PI 离散化方法（前向/后向欧拉、Tustin）、速度环与位置环分频、测速滤波。
- **电流环**：瑞萨固件极点配置、连续域 IMC、z 域直接对消；离散 Bode、裕度、z 平面极点、采样阶跃仿真。
- **速度环**：瑞萨固件极点配置、带宽法、对称最优；速度反馈可选无感 PLL（本固件）或编码器差分；多速率离散校核与负载扰动仿真。
- **自动整定与鲁棒性**：按目标相位 / 增益裕度自动求最大带宽（可按标称参数或最差参数组合），带宽–裕度权衡曲线；电感、电阻、惯量、转矩常数失配时的最差裕度与带宽范围。
- **位置环**：P + 速度/加速度前馈，跟随误差。
- **观测器**：PLL、PI 型扩展反电势观测器、滑模观测器、负载转矩观测器。
- **导入与对比**：解析固件 `r_motor_*_cfg.h`（支持注释、条件编译、宏表达式）自动填入参数，并设为对比基准，各页以虚线和“基准 / 变化”列显示差异。
- **硬件换算**：每条采样链路提供多种电路形式（附基本单元电路图，R1/R2… 编号与输入一一对应，本板另附原理图位号对照）；ADC 分辨率、相电流采样增益/偏置/量程/建立时间、母线电压分压、硬件过流阈值、VOT 温度、PWM 计数/死区误差/低边采样占空比上限；与固件宏对照。
- **使用体验**：图表悬停读数、输入错误标红、手机端“参数 / 结果”切换、命名方案保存 / 载入 / 导出 JSON。
- **导出**：C 头文件、MATLAB/Simulink、瑞萨 RX24T `*_cfg.h` 宏、JSON。

## 测试

```
npm test                 # 计算核心的自动化测试（Node 18+，无依赖）
npm run reference        # 重新生成 python-control 参考值（需要 pip install control）
npm run snapshot         # 结果有意改变后，更新回归快照
```

测试覆盖：电流环裕度与闭环极点对照 python-control 的独立计算；小数拍延时下 z 域模型与时域仿真逐拍一致；瑞萨固件增益公式；硬件换算与容差分析；cfg 头文件解析；导出代码均为合法 C 字面量；关键结果回归快照。GitHub Actions 在每次推送时运行测试，Pages 发布前也会先跑测试。

## 目录

- `app/core.js` 计算核心（纯函数，无页面依赖，浏览器与 Node 测试共用）
- `app/index.html` 界面与绘图
- `app/manifest.webmanifest`、`app/sw.js`、`app/icons/` PWA 安装与离线缓存
- `app/vendor/mathjax/` 本地公式渲染库（Apache-2.0）
- `tests/` 自动化测试与参考数据；`tools/` 参考值 / 快照生成脚本、在线预览版构建脚本
