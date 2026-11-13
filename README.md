# P25 Phase 1 核心（IMBE 声码器 + 空中接口）

用在 Baofeng DM-32UV 的开源固件 [OpenDMR](https://github.com/va3bpf/OpenDMR) 里（它的子模块 `p25/`）。

APCO P25 Phase 1 常规（直通 / 中继）语音的两部分，自己实现，C 代码不依赖电台，在电脑上就能测试：

- **IMBE 声码器**（全速率 7200x4400：每 20 ms 一帧，88 bit 参数 + 56 bit 纠错 = 144 bit）
- **IMBE 编码器**（语音 → 88 bit 参数）
- **空中接口**：帧同步、NID（BCH）、HDU / LDU1 / LDU2 / TDU / TDULC 的收和发，Golay / Hamming / Reed-Solomon 纠错
  （Reed-Solomon 支持“错误 + 删除”：收不到的符号当删除），链路控制、加密同步

电台那边的部分（HR_C7000 调制解调器的原始收发、码表、界面）在 OpenDMR 里。

**现状**：在 OpenDMR 里两台 DM-32UV 之间直通收发、和 MMDVM 热点互通（见 OpenDMR 的
[docs/p25.md](https://github.com/va3bpf/OpenDMR/blob/develop/docs/p25.md)，包括已知问题）。计划和核对过的协议参数见
[docs/plan.md](docs/plan.md)。

## 目录

| 路径 | 内容 |
|---|---|
| `c/imbe/` | IMBE 声码器（定点 C，解码和编码） |
| `c/p25/` | 空中接口：`p25_frame.c` 收发帧，`p25_fec.c` 纠错码 |
| `c/test/` | 电脑上的测试和命令行工具，`c/build.sh` 编译（WinLibs MinGW-w64 GCC） |
| `py/` | Python：IMBE 表的生成、和参考实现的比较、C4FM 解调 |
| `docs/plan.md` | 计划、参考资料和许可证、已经核对过的协议参数 |
| `ref/` | 参考源码（mbelib、dsd），**不进 git**，`fetch_refs.sh` 拉取 |
| `bin/`、`out/` | 编译和测试的输出，**不进 git** |

## 许可证

见 [LICENSE](LICENSE)：非商业用途可以自由使用、修改、分发，商业用途要作者书面许可。移植自 mbelib 和 DSD 的部分同时适用
它们的 ISC 许可证。IMBE 是 DVSI 的技术：本许可证不授予任何专利权利。
