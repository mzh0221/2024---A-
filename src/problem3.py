# -*- coding: utf-8 -*-
"""
problem3.py —— 问题3：最小螺距求解（模块10，算法4，§8.2）
================================================
题目：调头空间为以螺线中心为圆心、直径 9 m（半径 R=4.5 m）的圆形区域。
     确定最小螺距 p_min，使龙头前把手能沿相应螺线盘入到调头空间边界，
     且盘入全程（到达边界前）不发生碰撞。

算法4（§8.2，二分+碰撞检验）：
  1. 螺距搜索区间 [p_low, p_high]=[0.30, 0.60]，R=4.5，精度 ε_p=1e-5
  2. While p_high - p_low > ε_p:
       p_mid = (p_low+p_high)/2
       龙头到达边界极角 θ_in = 2πR/p_mid（取 0≤θ_in<2π 的等价角，起点 θ_start=32π）
       仿真龙头从 θ_start 盘入至 θ_in，逐时刻检查是否碰撞
       if 碰撞: p_low=p_mid（螺距太小）  else: p_high=p_mid（可继续缩小）
  3. p_min = (p_low+p_high)/2

说明：θ_in 的多圈处理用 θ_b=2πR/p 直接作为“到达边界时的极角”（螺线 r=aθ
     到达 r=R 对应 θ_b=2πR/p），检验全程 [0, t_b]（t_b 为到达边界时刻）是否碰撞。
     由于碰撞可能间歇性发生（不同板凳对在不同时刻碰撞后分离），需在整个区间
     以 dt 步长扫描，任一时刻碰撞即不可行。

参考结果：最小螺距约 0.45032 m（公开论文结果，量级 0.45 m）。
"""

import time as _timer
import numpy as np

import common as bd
import excel_writer as ew


def feasible_pitch(p, R=4.5, dt=0.5):
    """判定螺距 p 是否可行：龙头从 θ_start 盘入到调头空间边界（r=R）全程无碰撞。
    边界极角 θ_b = 2πR/p；到达边界时刻 t_b = (s(θ_start)-s(θ_b))/v0。
    在 [0, t_b] 区间以 dt 步长扫描，任一时刻碰撞即不可行。"""
    p = max(p, 1e-6)
    theta_b = 2.0 * np.pi * R / p
    if theta_b >= bd.THETA_START:
        return True                      # 螺距过大，起点已在边界内，必然可行
    s_start = bd.arc_length(bd.THETA_START, p)
    s_b = bd.arc_length(theta_b, p)
    t_b = (s_start - s_b) / bd.V0
    t = 0.0
    while t <= t_b + 1e-9:
        tt = min(t, t_b)
        if bd.collision_at(tt, p):
            return False
        t += dt
    return True


def solve_problem3():
    """二分法求最小螺距（可行=无碰撞）。
    feasible=True（可行）→ 尝试更小（p_hi=p_mid）；False（碰撞）→ 需更大（p_lo=p_mid）。"""
    print("\n========== 问题3：求最小螺距 ==========")
    R = 4.5
    p_lo, p_hi = 0.30, 0.60          # 下界碰撞（不可行），上界可行
    eps_p = 1e-5
    t0 = _timer.time()
    while p_hi - p_lo > eps_p:
        p_mid = (p_lo + p_hi) / 2.0
        if feasible_pitch(p_mid, R):
            p_hi = p_mid             # 可行，尝试更小
        else:
            p_lo = p_mid             # 碰撞，需更大
    p_min = (p_lo + p_hi) / 2.0
    print("  最小螺距 p_min≈%.5f m (参考 0.45032 m)，用时 %.1fs"
          % (p_min, _timer.time() - t0))

    # 边界极角与圈数
    theta_b = 2 * np.pi * R / p_min
    print("  到达边界极角 θ_b≈%.4f rad (%.2f 圈)" % (theta_b, theta_b / (2 * np.pi)))

    ew.write_result3(p_min, theta_b)
    return p_min


if __name__ == "__main__":
    solve_problem3()
