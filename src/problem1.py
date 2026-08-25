# -*- coding: utf-8 -*-
"""
problem1.py —— 问题1：300 秒仿真（模块8，算法2，§6.1）
================================================
题目：舞龙队沿螺距 p=0.55 m 的等距螺线顺时针盘入，龙头前把手速度恒 1 m/s，
     初始龙头位于螺线第16圈 A 点（θ_start=32π）。给出 0~300 s 每秒整个
     舞龙队的位置和速度，存入 result1.xlsx。

算法2（§6.1）：
  1. 参数初始化：p=0.55, v0=1, N=223, T=300, θ_start=32π
  2. For t=0..300（步长 1 s）:
       a. 由公式(5)弧长反演求龙头极角 θ_head(t)        → common.head_theta
       b. 由算法1递推全部 224 个把手极角 θ[0:224]        → common.position_recurse
       c. 由公式(2)转直角坐标                          → common.handles_xy
       d. 由公式(11)(12)求各把手速度                    → common.velocity_recurse
       e. 保存本时刻结果
  3. 按模板输出 result1.xlsx（保留6位小数）

并在控制台输出 0/60/120/180/240/300 s 关键时刻：
  龙头前把手、龙头后第 1/51/101/151/201 节龙身前把手、龙尾后把手 的位置与速度。
"""

import time as _timer
import numpy as np

import common as bd
import excel_writer as ew


def solve_problem1():
    print("\n========== 问题1：300 s 仿真 ==========")
    p = 0.55                         # 螺距
    T = 300                          # 仿真时长 (s)
    times = list(range(0, T + 1))    # 0..300，共 301 个时刻

    pos_by_t, vel_by_t = [], []
    t0 = _timer.time()
    for t in times:
        xy, v = bd.state_at(t, p)    # 单时刻位置+速度
        pos_by_t.append(xy)
        vel_by_t.append(v)
    print("  仿真完成，用时 %.1fs" % (_timer.time() - t0))

    # ---- 校验 ----
    print("  t=0   龙头速度=%.6f (应=1)" % vel_by_t[0][0])
    print("  t=300 龙头速度=%.6f (应=1)" % vel_by_t[300][0])
    print("  t=300 龙头坐标=(%.4f, %.4f)" % (pos_by_t[300][0, 0], pos_by_t[300][0, 1]))

    # ---- 关键时刻表格（论文用，表1/表2 格式） ----
    key_t = [0, 60, 120, 180, 240, 300]
    # 把手编号：龙头(0)、第1/51/101/151/201节龙身前把手(1/51/101/151/201)、龙尾后把手(223)
    key_h = [0, 1, 51, 101, 151, 201, 223]
    key_name = ["龙头", "第1节龙身", "第51节龙身", "第101节龙身",
                "第151节龙身", "第201节龙身", "龙尾（后）"]

    print("\n  --- 位置表（关键时刻） ---")
    print("  " + "\t".join("%ds" % t for t in key_t))
    for name, h in zip(key_name, key_h):
        xs = ["%.6f" % pos_by_t[t][h, 0] for t in key_t]
        ys = ["%.6f" % pos_by_t[t][h, 1] for t in key_t]
        print("  %sx\t" % name + "\t".join(xs))
        print("  %sy\t" % name + "\t".join(ys))

    print("\n  --- 速度表（关键时刻） ---")
    print("  " + "\t".join("%ds" % t for t in key_t))
    for name, h in zip(key_name, key_h):
        vs = ["%.6f" % vel_by_t[t][h] for t in key_t]
        print("  %s\t" % name + "\t".join(vs))

    # ---- 输出 result1.xlsx ----
    ew.write_result1(times, pos_by_t, vel_by_t)
    return pos_by_t, vel_by_t


if __name__ == "__main__":
    solve_problem1()
