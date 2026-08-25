# -*- coding: utf-8 -*-
"""
problem4.py —— 问题4：调头路径建模（模块11，§9）
================================================
题目：盘入螺距 p=1.7 m，盘出螺线与盘入螺线关于中心对称；舞龙队在问题3设定的
     调头空间（半径 R=4.5 m 的圆）内完成调头；调头路径由两段相切圆弧连成 S 形，
     前段半径是后段 2 倍，与盘入/盘出螺线均相切。问：能否调整圆弧使调头曲线变短？
     并以调头开始时刻为零时刻，给出 −100 s~100 s 每秒的位置与速度，存 result4.xlsx。

调头路径四段（《建模思路》§9.4）：
  ① 盘入螺线（p=1.7，θ 由起点递减至 θ₁，交边界于 A）
  ② 第一段圆弧（O₁,r₁，A→C）
  ③ 第二段圆弧（O₂,r₂，C→A′，r₁=2r₂）
  ④ 盘出螺线（与盘入螺线中心对称，r(θ)=(p/2π)(θ+π)）

关键公式（§9.2，切向连续修正版）：
  (22) θ₁ = 2πR/p                       入口极角
  (23) n₁ = (−sinθ₁−θ₁cosθ₁, cosθ₁−θ₁sinθ₁)/√(1+θ₁²)   A 处内法向
  (24) O₁ = A + r₁·n₁；O₂ = A′ − r₂·n₁   圆心（O₁在A内侧，O₂在A′内侧，构成S形）
  (25) r₂ = R·√(1+θ₁²)/(3θ₁)，r₁=2r₂    由外切条件 |O₁−O₂|=r₁+r₂ 推导（见下）
       推导：O₁−O₂=2A+3r₂·n₁，|O₁−O₂|²=4R²−12r₂·R·θ₁/√(1+θ₁²)+9r₂²=9r₂²
             ⟹ 4R²=12r₂·R·θ₁/√(1+θ₁²) ⟹ r₂=R·√(1+θ₁²)/(3θ₁)
  (26) L = r₁α₁ + r₂α₂                   调头曲线长度

结论（§9.3）：在题设约束下两圆弧半径由相切条件唯一确定，调头曲线长度为定值，
不存在调整圆弧的自由度，故不能使调头曲线变短。

把手机制：龙头沿路径行进（速度恒 1 m/s），后续把手由刚体链接约束
|Pᵢ−P_{i+1}|=lᵢ 求解，把手 i+1 在把手 i 后方（弧长更小）。路径分段后用各段
解析方程（螺线余弦定理/圆弧距离）求解，避免数值搜索跳根。
"""

import time as _timer
import numpy as np
from scipy.optimize import brentq, fsolve

import common as bd
import excel_writer as ew


# ========================= 弧长反演工具 =========================
def _arc_length_inv(s_target, p, lo=1e-6, hi=None):
    """弧长反演：求 θ 使 arc_length(θ,p)=s_target。s 关于 θ 单调递增。"""
    if hi is None:
        hi = bd.THETA_START
    if s_target <= bd.arc_length(lo, p):
        return lo
    if s_target >= bd.arc_length(hi, p):
        return hi
    return brentq(lambda th: bd.arc_length(th, p) - s_target, lo, hi,
                  xtol=1e-13, rtol=1e-13)


def _arc_s_from_phi(phi, phi_start, r, direction):
    """从弧起点 phi_start 沿 direction 方向到 phi 的弧长（非负）。
    direction>0: CCW（角度递增）；direction<0: CW（角度递减，可能跨 2π）。"""
    delta = phi - phi_start
    if direction > 0:
        while delta < -1e-12:
            delta += 2 * np.pi
        while delta > 2 * np.pi + 1e-12:
            delta -= 2 * np.pi
        return delta * r
    else:
        while delta > 1e-12:
            delta -= 2 * np.pi
        while delta < -2 * np.pi - 1e-12:
            delta += 2 * np.pi
        return -delta * r


# ========================= 调头路径几何（公式22-26，切向连续修正版） =========================
def turn_path_geometry(p, R=4.5):
    """计算调头路径全部几何参数。
    S 形双圆弧调头路径：弧1(CW) 从 A 到 C，弧2(CCW) 从 C 到 A′，两圆外切于 C。
    三处切向连续：A点(盘入螺线→弧1)、C点(弧1→弧2)、A'点(弧2→盘出螺线)。

    圆心配置（修正）：O₁=A+r₁·n₁（A内侧），O₂=A′−r₂·n₁（A′内侧）。
    外切条件 |O₁−O₂|=r₁+r₂ 代入 r₁=2r₂ 得 r₂=R·√(1+θ₁²)/(3θ₁)≈1.50m。
    两弧转角相等（α₁=α₂≈3.03 rad），总转角0（S形，两端切向相同）。
    """
    a = p / (2.0 * np.pi)
    theta1 = 2.0 * np.pi * R / p                         # (22) 入口极角
    A = np.array([R * np.cos(theta1), R * np.sin(theta1)])
    Ap = -A                                              # 出口（中心对称）
    sq = np.sqrt(1.0 + theta1 ** 2)
    n1 = np.array([-(np.sin(theta1) + theta1 * np.cos(theta1)),   # (23) 内法向
                   np.cos(theta1) - theta1 * np.sin(theta1)]) / sq
    # 外切条件 |O₁−O₂|=r₁+r₂，O₁=A+r₁·n₁(内侧)，O₂=A′−r₂·n₁(内侧)
    # 推导：O₁−O₂=2A+3r₂·n₁，|O₁−O₂|²=4R²−12r₂·R·θ₁/sq+9r₂²=9r₂² → r₂=R·sq/(3θ₁)
    r2 = R * sq / (3.0 * theta1)
    r1 = 2.0 * r2
    O1 = A + r1 * n1                                     # 弧1 圆心（A 内侧）
    O2 = Ap - r2 * n1                                    # 弧2 圆心（A′ 内侧）
    # 切点 C（在 O₁−O₂ 连线上，外切：|O₁C|=r₁，|O₂C|=r₂）
    C = O1 + (r1 / (r1 + r2)) * (O2 - O1)
    # 各点相对圆心的极角
    ang_A = np.arctan2(A[1] - O1[1], A[0] - O1[0])
    ang_C1 = np.arctan2(C[1] - O1[1], C[0] - O1[0])
    ang_C2 = np.arctan2(C[1] - O2[1], C[0] - O2[0])
    ang_Ap = np.arctan2(Ap[1] - O2[1], Ap[0] - O2[0])
    # 弧1: CW(dir1=−1)，A→C；弧2: CCW(dir2=+1)，C→A′
    # CW+CCW（反向）+ 外切 → C处切向连续；A处切向=盘入方向，A′处=盘出方向
    dir1 = -1.0
    alpha1 = (ang_A - ang_C1) % (2.0 * np.pi)            # CW: 角度递减
    dir2 = 1.0
    alpha2 = (ang_Ap - ang_C2) % (2.0 * np.pi)           # CCW: 角度递增
    L1 = r1 * alpha1                                     # 弧1 长度
    L2 = r2 * alpha2                                     # 弧2 长度
    L = L1 + L2                                          # (26) 调头曲线总长
    return dict(p=p, R=R, a=a, theta1=theta1, A=A, Ap=Ap, n1=n1,
                r1=r1, r2=r2, O1=O1, O2=O2, C=C,
                alpha1=alpha1, alpha2=alpha2, L1=L1, L2=L2, L=L,
                dir1=dir1, dir2=dir2,
                ang_A_O1=ang_A, ang_C_O1=ang_C1,
                ang_C_O2=ang_C2, ang_Ap_O2=ang_Ap)


# ========================= 路径上位置（四段） =========================
def path_position(s, geom):
    """调头路径上弧长 s 处的 (x,y) 坐标。
    s=0 为入口 A（调头开始）；s>0 为行进方向（弧→盘出）；s<0 为盘入螺线方向。
    弧长定义：盘入段 s<0（θ>θ₁），弧1 段 0≤s<L1，弧2 段 L1≤s<L1+L2，盘出段 s≥L。
    """
    p = geom['p']; a = geom['a']; theta1 = geom['theta1']
    L1 = geom['L1']; L2 = geom['L2']; L = L1 + L2
    s_theta1 = bd.arc_length(theta1, p)
    if s < 0:
        # 盘入螺线段：θ > θ₁（入口外侧）
        s_target = s_theta1 - s            # s<0 → s_target>s(θ₁)
        theta = _arc_length_inv(s_target, p, theta1, theta1 + 50 * np.pi)
        return np.array([a * theta * np.cos(theta), a * theta * np.sin(theta)])
    elif s < L1:
        # 第一段圆弧
        phi = geom['ang_A_O1'] + (s / geom['r1']) * geom['dir1']
        return geom['O1'] + geom['r1'] * np.array([np.cos(phi), np.sin(phi)])
    elif s < L:
        # 第二段圆弧
        phi = geom['ang_C_O2'] + ((s - L1) / geom['r2']) * geom['dir2']
        return geom['O2'] + geom['r2'] * np.array([np.cos(phi), np.sin(phi)])
    else:
        # 盘出螺线段：Q(ψ) = −spiral_xy(ψ)，ψ随 s 递增（盘出方向=向外，半径增大）
        # 弧长关系：arc_length(ψ) = s(θ₁) + (s − L)，故 ψ > θ₁ 且随 s 增大
        s_from_Ap = s - L
        s_target = s_theta1 + s_from_Ap    # ψ递增（向外）
        psi = _arc_length_inv(s_target, p, theta1, theta1 + 50 * np.pi)
        return -np.array([a * psi * np.cos(psi), a * psi * np.sin(psi)])


def _path_tangent(s, geom):
    """调头路径上弧长 s 处的单位切向量（沿行进方向）。"""
    p = geom['p']; a = p / (2.0 * np.pi)
    L1 = geom['L1']; L2 = geom['L2']; L = L1 + L2
    theta1 = geom['theta1']
    s_theta1 = bd.arc_length(theta1, p)
    if s < 0:
        # 盘入螺线：s 增加时 θ 递减，切向 = −dP/dθ
        s_target = s_theta1 - s
        theta = _arc_length_inv(s_target, p, theta1, theta1 + 50 * np.pi)
        dx = -(a * (np.cos(theta) - theta * np.sin(theta)))
        dy = -(a * (np.sin(theta) + theta * np.cos(theta)))
        n = np.sqrt(dx * dx + dy * dy)
        return np.array([dx, dy]) / n
    elif s < L1:
        # 弧1：CCW 切向 = R_90CCW(rad)
        phi = geom['ang_A_O1'] + (s / geom['r1']) * geom['dir1']
        rad = np.array([np.cos(phi), np.sin(phi)])
        return np.array([-rad[1], rad[0]]) if geom['dir1'] > 0 else np.array([rad[1], -rad[0]])
    elif s < L:
        # 弧2：CW 切向 = R_90CW(rad)
        phi = geom['ang_C_O2'] + ((s - L1) / geom['r2']) * geom['dir2']
        rad = np.array([np.cos(phi), np.sin(phi)])
        return np.array([-rad[1], rad[0]]) if geom['dir2'] > 0 else np.array([rad[1], -rad[0]])
    else:
        # 盘出螺线：Q(ψ)=−a·ψ·(cos,sin)，ψ随 s 递增（向外），切向 = dQ/dψ 方向 = −dP/dθ 方向
        s_target = s_theta1 + (s - L)        # ψ递增（向外）
        psi = _arc_length_inv(s_target, p, theta1, theta1 + 50 * np.pi)
        # dP/dθ 分量（P=a·θ·(cosθ,sinθ)）
        dx = a * (np.cos(psi) - psi * np.sin(psi))
        dy = a * (np.sin(psi) + psi * np.cos(psi))
        n = np.sqrt(dx * dx + dy * dy)
        return np.array([-dx, -dy]) / n        # 切向 = −dP/dθ 方向


# ========================= 各段后把手求解器 =========================
def _solve_next_on_spiral_out(psi_i, p, l):
    """盘出螺线上：已知 ψ_i，求后把手 ψ_{i+1}（ψ_{i+1}<ψ_i，盘出方向后把手在内侧）。
    方程(9)同形式，在 (ψ_i−2π, ψ_i) 内求根。"""
    a = p / (2.0 * np.pi)

    def eq(psi):
        return l * l - a * a * (psi_i * psi_i + psi * psi
                                - 2.0 * psi_i * psi * np.cos(psi - psi_i))
    dpsi_est = l / (a * np.sqrt(1.0 + psi_i * psi_i))
    lo = psi_i - 2.0 * np.pi + 1e-9
    hi = psi_i - 1e-9
    prev_d, prev_g = 1e-9, eq(psi_i - 1e-9)
    for fac in [0.3, 1.0, 1.5, 2.0, 3.0, 6.0, 12.0]:
        d = fac * dpsi_est
        if d <= prev_d:
            d = prev_d + 1e-6
        g = eq(psi_i - d)
        if prev_g > 0 and g <= 0:
            return brentq(eq, psi_i - d, psi_i - prev_d, xtol=1e-13, rtol=1e-13)
        prev_d, prev_g = d, g
    glo, ghi = eq(lo), eq(hi)
    if glo * ghi < 0:
        return brentq(eq, lo, hi, xtol=1e-13, rtol=1e-13)
    sol = fsolve(eq, psi_i - dpsi_est)[0]
    if lo < sol < hi and abs(eq(sol)) < 1e-10:
        return sol
    raise ValueError("盘出螺线无有效根: psi_i=%r p=%r l=%r" % (psi_i, p, l))


def _solve_next_on_arc(P_i, O, r, phi_start, phi_end, direction, l, phi_hint=None):
    """圆弧段：已知把手 i 位置 P_i，在圆心 O、半径 r、起止角 phi_start/phi_end、
    方向 direction 的圆弧上找把手 i+1 使 |P_i−P_{i+1}|=l。
    phi_hint 提供估计角度，避免长弧上找到错误根。返回 φ 或 None。"""
    def pos(phi):
        return O + r * np.array([np.cos(phi), np.sin(phi)])

    def dist_sq(phi):
        d = P_i - pos(phi)
        return np.dot(d, d) - l * l

    def arc_s(phi):
        delta = phi - phi_start
        if direction > 0:
            while delta < 0:
                delta += 2 * np.pi
            while delta > 2 * np.pi + 1e-10:
                delta -= 2 * np.pi
            return delta * r
        else:
            while delta > 1e-10:
                delta -= 2 * np.pi
            while delta < -2 * np.pi - 1e-10:
                delta += 2 * np.pi
            return -delta * r

    # 实际弧角度跨度与弧长
    if direction > 0:
        dphi_arc = (phi_end - phi_start) % (2 * np.pi)
        if dphi_arc < 1e-10:
            dphi_arc = 2 * np.pi
        lo, hi = phi_start - 0.01, phi_start + dphi_arc + 0.01
    else:
        dphi_arc = (phi_start - phi_end) % (2 * np.pi)
        if dphi_arc < 1e-10:
            dphi_arc = 2 * np.pi
        lo, hi = phi_start - dphi_arc - 0.01, phi_start + 0.01
    L_arc = dphi_arc * r

    # 多点探测找所有变号区间
    span = abs(hi - lo)
    n_probe = max(200, int(np.ceil(span / (dphi_arc / 50 + 1e-10))) + 1)
    ts = np.linspace(lo, hi, n_probe)
    ds = [dist_sq(t) for t in ts]
    roots = []
    for k in range(n_probe - 1):
        if ds[k] * ds[k + 1] <= 0:
            root = brentq(dist_sq, ts[k], ts[k + 1], xtol=1e-13, rtol=1e-13)
            s_root = arc_s(root)
            if -0.01 <= s_root <= L_arc + 0.01:
                roots.append(root)
    for phi_endpt in [lo, hi]:
        if abs(dist_sq(phi_endpt)) < 0.05:
            s_root = arc_s(phi_endpt)
            if -0.01 <= s_root <= L_arc + 0.01:
                roots.append(phi_endpt)
    if not roots:
        return None
    # 去重
    unique = [roots[0]]
    for rr in roots[1:]:
        if min(abs(rr - ur) for ur in unique) > 1e-8:
            unique.append(rr)
    if phi_hint is not None:
        def ang_dist(a, b):
            d = abs(a - b) % (2 * np.pi)
            return min(d, 2 * np.pi - d)
        return min(unique, key=lambda rr: ang_dist(rr, phi_hint))
    return min(unique, key=lambda rr: arc_s(rr))


def _solve_next_on_spiral_from_point(P_i, p, theta_lo, theta_hi, l):
    """已知把手 i 位置 P_i（可在任意段），在盘入螺线 [theta_lo,theta_hi] 上
    找把手 i+1 使 |P_i−spiral_xy(θ)|=l。返回 θ 或 None。"""
    a = p / (2.0 * np.pi)

    def dist_sq(theta):
        x = a * theta * np.cos(theta)
        y = a * theta * np.sin(theta)
        return (P_i[0] - x) ** 2 + (P_i[1] - y) ** 2 - l * l
    n_probe = 100
    ts = np.linspace(theta_lo, theta_hi, n_probe)
    ds = [dist_sq(t) for t in ts]
    for k in range(n_probe - 1):
        if ds[k] * ds[k + 1] <= 0:
            return brentq(dist_sq, ts[k], ts[k + 1], xtol=1e-13, rtol=1e-13)
    return None


# ========================= 把手 i+1 弧长求解（分段解析） =========================
def _solve_next_handle(s_i, geom, l, P_i=None):
    """已知把手 i 弧长 s_i 与位置 P_i，求把手 i+1 的弧长 s_{i+1}（<s_i，后方）。
    根据把手 i 所在段使用对应解析方程，避免数值搜索跳根。"""
    p = geom['p']; a = p / (2.0 * np.pi)
    L1 = geom['L1']; L2 = geom['L2']; L = L1 + L2
    theta1 = geom['theta1']
    s_theta1 = bd.arc_length(theta1, p)
    s_est = s_i - l                                  # 估计 s_{i+1}
    P_pos = P_i if P_i is not None else path_position(s_i, geom)

    # --- 把手 i 在盘入螺线段 ---
    if s_i < 0:
        theta_i = _arc_length_inv(s_theta1 - s_i, p, theta1, theta1 + 50 * np.pi)
        if s_est < 0:                                # 后把手也在盘入螺线段
            theta_next = bd.solve_theta_next(theta_i, p, l)
            return -(bd.arc_length(theta_next, p) - s_theta1)
        # 后把手可能跨到圆弧1
        phi = _solve_next_on_arc(P_pos, geom['O1'], geom['r1'],
                                 geom['ang_A_O1'], geom['ang_C_O1'], geom['dir1'], l)
        if phi is not None:
            return _arc_s_from_phi(phi, geom['ang_A_O1'], geom['r1'], geom['dir1'])
        return _search_fallback(P_pos, s_est, geom, l)

    # --- 把手 i 在圆弧1 ---
    elif s_i < L1:
        O1 = geom['O1']; r1 = geom['r1']
        s_est_arc1 = max(s_est, 0)
        phi_hint1 = geom['ang_A_O1'] + (s_est_arc1 / r1) * geom['dir1']
        phi = _solve_next_on_arc(P_pos, O1, r1, geom['ang_A_O1'], geom['ang_C_O1'],
                                 geom['dir1'], l, phi_hint=phi_hint1)
        if phi is not None:
            s_cand = _arc_s_from_phi(phi, geom['ang_A_O1'], r1, geom['dir1'])
            if 0 <= s_cand < s_i:
                return s_cand
        theta_next = _solve_next_on_spiral_from_point(P_pos, p, theta1, theta1 + 2 * np.pi, l)
        if theta_next is not None:
            return -(bd.arc_length(theta_next, p) - s_theta1)
        return _search_fallback(P_pos, s_est, geom, l)

    # --- 把手 i 在圆弧2 ---
    elif s_i < L:
        O2 = geom['O2']; r2 = geom['r2']
        s_est_arc2 = max(s_est, L1)
        phi_hint2 = geom['ang_C_O2'] + ((s_est_arc2 - L1) / r2) * geom['dir2']
        phi = _solve_next_on_arc(P_pos, O2, r2, geom['ang_C_O2'], geom['ang_Ap_O2'],
                                 geom['dir2'], l, phi_hint=phi_hint2)
        if phi is not None:
            s_cand = L1 + _arc_s_from_phi(phi, geom['ang_C_O2'], r2, geom['dir2'])
            if L1 <= s_cand < s_i:
                return s_cand
        # 后把手可能在圆弧1
        s_est_arc1 = max(s_est, 0)
        phi_hint1 = geom['ang_A_O1'] + (s_est_arc1 / geom['r1']) * geom['dir1']
        phi1 = _solve_next_on_arc(P_pos, geom['O1'], geom['r1'], geom['ang_A_O1'],
                                  geom['ang_C_O1'], geom['dir1'], l, phi_hint=phi_hint1)
        if phi1 is not None:
            s_cand = _arc_s_from_phi(phi1, geom['ang_A_O1'], geom['r1'], geom['dir1'])
            if 0 <= s_cand < L1:
                return s_cand
        # 后把手可能在盘入螺线段
        theta_next = _solve_next_on_spiral_from_point(P_pos, p, theta1, theta1 + 0.5, l)
        if theta_next is None:
            theta_next = _solve_next_on_spiral_from_point(P_pos, p, theta1, theta1 + 2 * np.pi, l)
        if theta_next is not None:
            return -(bd.arc_length(theta_next, p) - s_theta1)
        return _search_fallback(P_pos, s_est, geom, l)

    # --- 把手 i 在盘出螺线段（ψ 随 s 递增，向外）---
    else:
        psi_i = _arc_length_inv(s_theta1 + (s_i - L), p, theta1, theta1 + 50 * np.pi)
        if s_est >= L:                                # 后把手也在盘出螺线段
            psi_next = _solve_next_on_spiral_out(psi_i, p, l)   # ψ_{i+1}<ψ_i（后方）
            return L + (bd.arc_length(psi_next, p) - s_theta1)
        # 后把手可能在圆弧2
        s_est_arc2 = max(s_est, L1)
        phi_hint2 = geom['ang_C_O2'] + ((s_est_arc2 - L1) / geom['r2']) * geom['dir2']
        phi2 = _solve_next_on_arc(P_pos, geom['O2'], geom['r2'], geom['ang_C_O2'],
                                  geom['ang_Ap_O2'], geom['dir2'], l, phi_hint=phi_hint2)
        if phi2 is not None:
            s_cand = L1 + _arc_s_from_phi(phi2, geom['ang_C_O2'], geom['r2'], geom['dir2'])
            if L1 <= s_cand < L:
                return s_cand
        # 后把手可能在圆弧1
        s_est_arc1 = max(s_est, 0)
        phi_hint1 = geom['ang_A_O1'] + (s_est_arc1 / geom['r1']) * geom['dir1']
        phi1 = _solve_next_on_arc(P_pos, geom['O1'], geom['r1'], geom['ang_A_O1'],
                                  geom['ang_C_O1'], geom['dir1'], l, phi_hint=phi_hint1)
        if phi1 is not None:
            s_cand = _arc_s_from_phi(phi1, geom['ang_A_O1'], geom['r1'], geom['dir1'])
            if 0 <= s_cand < L1:
                return s_cand
        theta_next = _solve_next_on_spiral_from_point(P_pos, p, theta1, theta1 + 0.5, l)
        if theta_next is None:
            theta_next = _solve_next_on_spiral_from_point(P_pos, p, theta1, theta1 + 2 * np.pi, l)
        if theta_next is not None:
            return -(bd.arc_length(theta_next, p) - s_theta1)
        return _search_fallback(P_pos, s_est, geom, l)


def _search_fallback(P_i, s_est, geom, l):
    """兜底：在路径上数值搜索 s 使 |P_i−P(s)|=l。"""
    def dist_sq(s):
        d = P_i - path_position(s, geom)
        return np.dot(d, d) - l * l
    lo = s_est - 2.0 * l
    hi = s_est + 0.5 * l
    n_probe = 30
    ss = np.linspace(lo, hi, n_probe)
    ds = [dist_sq(s) for s in ss]
    for k in range(n_probe - 1):
        if ds[k] * ds[k + 1] <= 0:
            return brentq(dist_sq, ss[k], ss[k + 1], xtol=1e-13, rtol=1e-13)
    best_k = int(np.argmin(np.abs(ds)))
    return ss[best_k]


# ========================= 全部把手位置与速度（调头路径） =========================
def compute_handles_turning(s_head, geom, lengths=None):
    """给定龙头弧长 s_head（从入口 A 起），递推全部 224 个把手位置 (224,2)。
    龙头在 s_head，后续把手依次在后方（弧长递减）。"""
    if lengths is None:
        lengths = bd.L_ARR
    xy = np.empty((bd.N_HANDLE, 2))
    xy[0] = path_position(s_head, geom)
    s_cur = s_head
    for i in range(bd.N_BENCH):
        l = lengths[i]
        s_next = _solve_next_handle(s_cur, geom, l, xy[i])
        xy[i + 1] = path_position(s_next, geom)
        s_cur = s_next
    return xy


def velocities_turning(s_head, xy, geom, lengths=None, v0=bd.V0):
    """解析法计算调头路径上各把手速度。
    刚体链接约束 |P_i−P_{i+1}|=l（常数），求导得
      (P_i−P_{i+1})·v_{i+1} = (P_i−P_{i+1})·v_i
    且 v_{i+1} 沿路径切向 T_{i+1}，故
      v_{i+1} = [(P_i−P_{i+1})·v_i] / [(P_i−P_{i+1})·T_{i+1}] · T_{i+1}
    龙头速度 = v0·T_0。返回各把手速度大小数组 (N_HANDLE,)。
    （解析法避免数值微分在切线衔接处的不稳定。）"""
    if lengths is None:
        lengths = bd.L_ARR
    n = len(xy)
    speeds = np.empty(n)
    T0 = _path_tangent(s_head, geom)
    v_vec = v0 * T0
    speeds[0] = v0
    s_cur = s_head
    for i in range(n - 1):
        l = lengths[i]
        s_next = _solve_next_handle(s_cur, geom, l, xy[i])
        T_next = _path_tangent(s_next, geom)
        link = xy[i] - xy[i + 1]
        proj_v = np.dot(link, v_vec)
        proj_T = np.dot(link, T_next)
        if abs(proj_T) < 1e-15:                 # 退化（链接垂直切向）
            speeds[i + 1] = 0.0
            v_vec = np.zeros(2)
        else:
            alpha = proj_v / proj_T
            v_vec = alpha * T_next
            speeds[i + 1] = abs(alpha)
        s_cur = s_next
    return speeds


# ========================= 问题4 主流程 =========================
def solve_problem4():
    """问题4：调头路径建模，输出 −100 s~100 s 每秒位置速度至 result4.xlsx。"""
    print("\n========== 问题4：调头路径建模 ==========")
    p4 = 1.7
    geom = turn_path_geometry(p4, R=4.5)
    print("  螺距 p=%.1fm, 调头空间 R=4.5m" % p4)
    print("  圆弧半径 r1=%.4fm, r2=%.4fm (r1=2r2)" % (geom['r1'], geom['r2']))
    print("  圆心角 α1=%.6frad, α2=%.6frad (S形两弧转角相等)" % (geom['alpha1'], geom['alpha2']))
    print("  调头曲线长度 L=%.4fm, 调头时间=%.4fs" % (geom['L'], geom['L']))
    print("  结论：半径由相切条件唯一确定，调头曲线长度为定值，不能变短。")

    times = list(range(-100, 101))               # −100..100，共 201 个时刻
    t0 = _timer.time()
    pos_by_t, vel_by_t = [], []
    for t in times:
        s_head = bd.V0 * t                        # 龙头弧长（调头开始为零时刻）
        xy = compute_handles_turning(s_head, geom)
        v = velocities_turning(s_head, xy, geom)
        pos_by_t.append(xy)
        vel_by_t.append(v)
    print("  位置+速度计算完成，用时 %.1fs" % (_timer.time() - t0))

    # 校验
    print("  t=0   龙头速度=%.6f (应=1)" % vel_by_t[100][0])
    print("  t=0   龙头坐标=(%.4f, %.4f)" % (pos_by_t[100][0, 0], pos_by_t[100][0, 1]))

    # 关键时刻表格（−100,−50,0,50,100 s）
    key_t_idx = [0, 50, 100, 150, 200]
    key_h = [0, 1, 51, 101, 151, 201, 223]
    print("\n  --- 速度表（关键时刻） ---")
    print("  " + "\t".join("%ds" % times[i] for i in key_t_idx))
    for h in key_h:
        vs = ["%.6f" % vel_by_t[i][h] for i in key_t_idx]
        print("  把手%3d\t" % h + "\t".join(vs))

    ew.write_result4(times, pos_by_t, vel_by_t)
    return pos_by_t, vel_by_t, geom


if __name__ == "__main__":
    solve_problem4()
