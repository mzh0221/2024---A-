# -*- coding: utf-8 -*-
"""
problem5.py —— 问题5：最大龙头速度求解（模块12，算法5，§10.2）
================================================
题目：舞龙队沿问题4设定的路径行进，龙头速度保持不变，确定龙头的最大行进速度，
     使舞龙队各把手速度均不超过 2 m/s。

缩放法原理（§10.1，公式28）：
  速度递推式(11)关于 dθ_i/dt 齐次线性：若龙头速度由 v0 变为 λv0，则所有把手速度
  同步变为原来的 λ 倍（位置仅由几何确定、与速度无关）。故以 v0=1 m/s 仿真问题4
  路径，记录全队最大速度 v_max^0，则龙头允许的最大行进速度为
      v_max · v_max^0 ≤ 2   ⟹   v_max = 2 / v_max^0    (28)

算法5（§10.2）：
  1. 在问题4的调头路径上，设 v0=1 m/s
  2. 计算调头全程所有把手速度
  3. 记录全队最大速度 v_max^0
  4. 龙头最大速度 v_max = 2 / v_max^0

注意：问题4 的 result4.xlsx 按 1 s 采样，会漏掉调头过程中的速度峰值（峰值出现
在非整数时刻）。因此问题5 需对调头过程（0≤t≤L）做精细扫描（0.05 s 步长）以
准确定位最大速度。最大速度通常出现在调头圆弧段（曲率变化剧烈处）。

参考结果：最大龙头速度约 1.246267 m/s（黄文杰等论文结果）。
"""

import time as _timer

import common as bd
import problem4 as p4
import excel_writer as ew


def solve_problem5():
    print("\n========== 问题5：最大龙头速度 ==========")
    p4_pitch = 1.7
    geom = p4.turn_path_geometry(p4_pitch, R=4.5)
    L = geom['L']                                  # 调头曲线长度

    # 精细扫描调头过程 0~L+1，步长 0.05 s（用解析速度，无需数值微分）
    t0 = _timer.time()
    v_max_0 = 0.0                                  # v0=1 时全队最大速度
    t_at_max = 0.0
    h_at_max = 0
    t = 0.0
    while t <= L + 1.0:
        s_head = bd.V0 * t
        xy = p4.compute_handles_turning(s_head, geom)
        v = p4.velocities_turning(s_head, xy, geom)
        v_max_t = v.max()
        if v_max_t > v_max_0:
            v_max_0 = v_max_t
            t_at_max = t
            h_at_max = int(v.argmax())
        t += 0.05
    print("  精细扫描完成（0~%.1fs, 步长0.05s），用时 %.1fs" % (L + 1.0, _timer.time() - t0))
    print("  峰值位置: t=%.2fs, 把手%d" % (t_at_max, h_at_max))

    v_max = 2.0 / v_max_0                          # 公式(28)
    print("  v0=1m/s 时全程最大把手速度 v_max^0=%.6f m/s" % v_max_0)
    print("  最大龙头速度 v_max=2/v_max^0=%.6f m/s (参考 1.246267 m/s)" % v_max)

    ew.write_result5(v_max, v_max_0, t_at_max, h_at_max)
    return v_max


if __name__ == "__main__":
    solve_problem5()
