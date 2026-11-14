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
- **控制器**：开关频率、单采单更/双采双更、计算延时拍数、PI 离散化方法（前向/后向欧拉、Tustin）、速度环与位置环分频、测速滤波。
- **电流环**：瑞萨固件极点配置、连续域 IMC、z 域直接对消；离散 Bode、裕度、z 平面极点、采样阶跃仿真。
- **速度环**：瑞萨固件极点配置、带宽法、对称最优；多速率离散校核与负载扰动仿真。
- **位置环**：P + 速度/加速度前馈，跟随误差。
- **观测器**：PLL、PI 型扩展反电势观测器、滑模观测器、负载转矩观测器。
- **硬件换算**：每条采样链路提供多种电路形式（附基本单元电路图，R1/R2… 编号与输入一一对应，本板另附原理图位号对照）；ADC 分辨率、相电流采样增益/偏置/量程/建立时间、母线电压分压、硬件过流阈值、VOT 温度、PWM 计数/死区误差/低边采样占空比上限；与固件宏对照。
- **导出**：C 头文件、MATLAB/Simulink、瑞萨 RX24T `*_cfg.h` 宏、JSON。

## 目录

- `app/index.html` 应用本体（单文件）
- `app/manifest.webmanifest`、`app/sw.js`、`app/icons/` PWA 安装与离线缓存
- `app/vendor/mathjax/` 本地公式渲染库（Apache-2.0）
- `tools/build_artifact.py` 生成 claude.ai 在线预览版
