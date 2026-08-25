# -*- coding: utf-8 -*-
"""
problem2.py —— 问题2：盘入终止时刻求解（模块9，算法3，§7.1）
================================================
题目：舞龙队沿问题1的螺线继续盘入，确定盘入终止时刻（板凳间不能再不碰撞地
     继续盘入的时间），给出此时舞龙队的位置和速度，存入 result2.xlsx。

算法3（§7.1，二分法）：
  1. 下界 t_low=300（问题1已确认 300 s 内无碰撞），上界 t_high=500（保证已碰撞）
  2. While t_high - t_low > ε(=1e-3):
       t_mid = (t_low+t_high)/2
       计算 t_mid 时刻所有板凳顶点（算法1 + 公式16-19）
       执行碰撞检测（§11.2，含端点接触/共线判定）
       if 碰撞: t_high=t_mid  else: t_low=t_mid
  3. T_collision = t_low，输出该时刻位置速度至 result2.xlsx

参考结果：首次碰撞时间约 412.47 s（公开论文结果）。
碰撞检测假设（§5.3）：首次碰撞必涉及龙头与相邻外圈龙身，只检查龙头板凳与
不相邻外圈板凳，详见 common.head_collision。
"""

import time as _timer

import common as bd
import excel_writer as ew


def find_collision_time(p=0.55, t_low=300.0, t_high=None, eps=1e-3):
    """二分法求盘入终止时刻（碰撞临界时刻）。
    上界自适应：取龙头到达螺线中心总时间的 95%，确保不会因龙头过中心导致
    碰撞检测失效（th0 钳位到 0 时碰撞检测返回 False）。
    返回 T_collision（≈412.47 s）。"""
    import numpy as np
    if t_high is None:
        # 龙头从 θ_start 盘入到中心的总弧长时间
        s_total = bd.arc_length(bd.THETA_START, p)
        t_center = s_total / bd.V0                # 龙头到达中心时刻
        t_high = min(500.0, t_center * 0.97)      # 不超过中心，留 3% 余量
        # 确保 t_high 时确实碰撞
        while t_high > t_low + 1.0 and not bd.collision_at(t_high, p):
            t_high = (t_low + t_high) / 2.0
    t0 = _timer.time()
    while t_high - t_low > eps:
        t_mid = (t_low + t_high) / 2.0
        if bd.collision_at(t_mid, p):
            t_high = t_mid          # 已碰撞，终止时刻更早
        else:
            t_low = t_mid           # 未碰撞，终止时刻更晚
    t_coll = t_low
    print("  终止时刻 T≈%.4f s (参考 412.47 s)，用时 %.1fs"
          % (t_coll, _timer.time() - t0))
    return t_coll


def solve_problem2():
    print("\n========== 问题2：求盘入终止时刻 ==========")
    p = 0.55
    t_coll = find_collision_time(p)

    # 输出该时刻位置速度
    xy, v = bd.state_at(t_coll, p)
    print("  龙头坐标=(%.4f, %.4f)  速度=%.6f" % (xy[0, 0], xy[0, 1], v[0]))

    # 关键把手表格（龙头、第1/51/101/151/201节龙身前把手、龙尾后把手）
    key_h = [0, 1, 51, 101, 151, 201, 223]
    key_name = ["龙头", "第1节龙身", "第51节龙身", "第101节龙身",
                "第151节龙身", "第201节龙身", "龙尾（后）"]
    print("  --- 终止时刻位置速度 ---")
    print("  \t横坐标x\t纵坐标y\t速度")
    for name, h in zip(key_name, key_h):
        print("  %s\t%.6f\t%.6f\t%.6f"
              % (name, xy[h, 0], xy[h, 1], v[h]))

    ew.write_result2("result2.xlsx", "result2.xlsx", xy, v)
    return t_coll


if __name__ == "__main__":
    solve_problem2()
