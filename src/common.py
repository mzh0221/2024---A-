# -*- coding: utf-8 -*-
"""
板凳龙（2024 国赛 A 题）核心数学模型模块 common.py
================================================
对应《代码需求表》模块 1~7，依据《建模思路（修订版）》公式 (1)~(19) 实现：

  模块1  solve_theta_next()      位置递推方程(9)求根（§11.1）
  模块2  cross/on_segment/        线段相交判定（§11.2），≤ 判相交 + 共线判接触
         segments_intersect
  模块3  head_theta(t)           弧长积分公式(5)反演龙头极角（§2.3，二分法）
  模块4  polar_to_xy()           公式(2)极坐标转直角坐标（§2.1）
  模块5  position_recurse()      算法1：逐节递推全部 224 个把手极角（§3.2）
  模块6  velocity_recurse()      速度递推公式(11)(12)（§4.2）
  模块7  vertex_coords()         公式(16)-(19) 板凳四顶点坐标（§5.1）

【方向说明（重要）】
  建模思路文档第 045 段写“盘入时 θᵢ₊₁ < θᵢ”，但据题目附件 result1.xlsx 校验：
  t=0 时龙头前把手 θ₀ = 32π ≈ 100.531 rad（半径 8.8 m），
  而龙头后把手（第1节龙身前把手）θ₁ ≈ 100.89 rad > θ₀。
  物理上盘入时龙头在最内圈（极角最小），龙身/龙尾逐节向外，θ 递增。
  故本实现采用校验过的正确方向：θ_{i+1} ∈ (θ_i, θ_i + 2π)，后把手在外圈。
  （建模思路文档此处的方向描述有误，特此修正。）
"""

import numpy as np
from scipy.optimize import brentq, fsolve


# ========================= 全局参数（《建模思路》§1.2 符号表） =========================
N_BENCH = 223          # 板凳总数：1 龙头 + 221 龙身 + 1 龙尾
N_HANDLE = 224         # 把手总数 = 板凳数 + 1（把手编号 0..223）
L_HEAD = 2.86          # 龙头前后把手距离 (m)：341cm - 2×27.5cm = 286cm
L_BODY = 1.65          # 龙身/龙尾前后把手距离 (m)：220cm - 2×27.5cm = 165cm
WIDTH = 0.30           # 板凳全宽 (m)
HALF_W = WIDTH / 2     # 半宽 w = 0.15 m（顶点公式中的 w）
D_HANDLE = 0.275      # 把手孔中心到最近板头距离 (m) = 27.5 cm
V0 = 1.0              # 龙头前把手基准行进速度 (m/s)
THETA_START = 32 * np.pi   # 龙头初始极角（螺线第16圈 A 点），θ_start = 32π


def bench_lengths():
    """每节板凳的把手距离数组（下标 = 板凳号 0..222）。
    l[0] = 龙头 = 2.86 m，其余 = 龙身/龙尾 = 1.65 m。
    板凳 i 连接把手 i（前把手）与把手 i+1（后把手）。"""
    L = np.full(N_BENCH, L_BODY)
    L[0] = L_HEAD
    return L


# 默认把手距离数组（模块内常用）
L_ARR = bench_lengths()


# ========================= 模块3：弧长积分与龙头位置（公式4、5） =========================
def arc_length(theta, p):
    """阿基米德螺线弧长积分（公式4）。
    螺线 r = (p/2π)·θ，记 a = p/(2π)，弧长
        s(θ) = (a/2)·[ θ·√(1+θ²) + ln(θ + √(1+θ²)) ]
    s 关于 θ 单调递增（θ>0），故可由 s 反解 θ。
    """
    a = p / (2.0 * np.pi)
    sq = np.sqrt(1.0 + theta * theta)
    return (a / 2.0) * (theta * sq + np.log(theta + sq))


def head_theta(t, p, v0=V0, theta_start=THETA_START):
    """由弧长守恒反演 t 时刻龙头前把手极角 θ₀(t)（公式5）。
    盘入时龙头向内运动，弧长守恒：
        s(θ_start) - s(θ₀(t)) = v0·t   →   s(θ₀(t)) = s(θ_start) - v0·t
    因 s 单调增、θ₀(t) 单调减，解唯一，用 brentq 二分求解（误差<1e-12）。
    若 t 超过螺线总弧长时间（理论已盘到中心），钳位到极小 θ 避免越界。
    """
    s_start = arc_length(theta_start, p)
    s_target = s_start - v0 * t
    lo = 1e-6
    hi = theta_start
    if s_target <= arc_length(lo, p):
        return lo                 # 龙头已逼近中心，钳位
    return brentq(lambda th: arc_length(th, p) - s_target, lo, hi,
                  xtol=1e-13, rtol=1e-13)


# ========================= 模块5：逐节把手极角递推（算法1，公式9） =========================
def _eq_theta_next(th, theta_i, p, l):
    """方程(9)残差：l² - a²·[θ_i² + th² - 2·θ_i·th·cos(th - θ_i)] = 0。
    其中 a = p/(2π)。"""
    a = p / (2.0 * np.pi)
    return l * l - a * a * (theta_i * theta_i + th * th
                            - 2.0 * theta_i * th * np.cos(th - theta_i))


def solve_theta_next(theta_i, p, l):
    """已知前把手极角 θ_i，求后把手极角 θ_{i+1}（公式9求根，模块1）。
    正确方向：盘入时后把手在外圈，θ_{i+1} ∈ (θ_i, θ_i + 2π)。
    方程(9)在该区间通常有两个根：同圈近邻根（物理正确，Δθ≈l/(a·√(1+θ²))）
    与隔圈根（Δθ≈2π）。两端残差同号无法直接括号，故用弧长估计定位“第一个根”。
    """
    a = p / (2.0 * np.pi)
    # 弧长近似估计 Δθ（同圈近邻根的近似位置）
    dth_est = l / (a * np.sqrt(1.0 + theta_i * theta_i))

    def eq(th):
        return _eq_theta_next(th, theta_i, p, l)

    # 在估计根附近取若干采样点，找首个 +→- 变号区间（即第一个根）
    probes = [0.3, 1.0, 1.5, 2.0, 3.0, 6.0, 12.0]
    prev_d, prev_g = 1e-9, eq(theta_i + 1e-9)              # g(θ_i+)≈l²>0
    for fac in probes:
        d = fac * dth_est
        if d <= prev_d:
            d = prev_d + 1e-6
        g = eq(theta_i + d)
        if prev_g > 0 and g <= 0:                           # 首个变号区间
            return brentq(eq, theta_i + prev_d, theta_i + d,
                          xtol=1e-13, rtol=1e-13)
        prev_d, prev_g = d, g
    # 兜底：大区间扫描
    lo = theta_i + 1e-9
    hi = theta_i + 2.0 * np.pi - 1e-9
    glo, ghi = eq(lo), eq(hi)
    if glo * ghi < 0:
        return brentq(eq, lo, hi, xtol=1e-13, rtol=1e-13)
    sol = fsolve(eq, theta_i + dth_est)[0]
    if theta_i < sol < theta_i + 2 * np.pi and abs(eq(sol)) < 1e-10:
        return sol
    raise ValueError("无有效根: theta_i=%r p=%r l=%r" % (theta_i, p, l))


def position_recurse(theta_head, p, lengths=L_ARR):
    """算法1：给定龙头极角，递推全部 224 个把手极角（模块5）。
    返回 theta[0..223]：theta[0]=龙头前把手（最小），theta[223]=龙尾后把手（最大）。
    """
    thetas = np.empty(N_HANDLE)
    thetas[0] = theta_head
    for i in range(N_BENCH):
        thetas[i + 1] = solve_theta_next(thetas[i], p, lengths[i])
    return thetas


# ========================= 模块4：极坐标↔直角坐标（公式2） =========================
def polar_to_xy(theta, p):
    """公式(2)：x = a·θ·cosθ, y = a·θ·sinθ，其中 a = p/(2π)。
    输入标量或数组 θ，返回 (x, y)。"""
    a = p / (2.0 * np.pi)
    return a * theta * np.cos(theta), a * theta * np.sin(theta)


def handles_xy(thetas, p):
    """全部把手的直角坐标，返回 (N_HANDLE, 2) 数组。"""
    x, y = polar_to_xy(thetas, p)
    return np.column_stack([x, y])


# ========================= 模块6：速度解析递推（公式11、12） =========================
def velocity_recurse(thetas, p, v0=V0):
    """速度解析递推（公式11、12，模块6）。
      龙头: dθ₀/dt = -v0 / [a·√(1+θ₀²)]   （盘入时 θ 递减，取负）
      递推: dθ_{i+1}/dt = f(θ_i, θ_{i+1})·dθ_i/dt   （公式11）
      速度大小: |v_i| = a·√(1+θ_i²)·|dθ_i/dt|   （公式12）
    返回各把手速度大小数组 (N_HANDLE,)。
    注：公式(11)关于 dθ_i/dt 齐次线性，是问题5速度缩放法成立的基础。
    """
    a = p / (2.0 * np.pi)
    n = len(thetas)
    dthetadt = np.empty(n)
    dthetadt[0] = -v0 / (a * np.sqrt(1.0 + thetas[0] ** 2))
    for i in range(n - 1):
        ti, tip1 = thetas[i], thetas[i + 1]
        dt = tip1 - ti
        c, s = np.cos(dt), np.sin(dt)
        num = ti - tip1 * c - ti * tip1 * s
        den = ti * c - tip1 - ti * tip1 * s
        dthetadt[i + 1] = (num / den) * dthetadt[i]
    speeds = a * np.sqrt(1.0 + thetas ** 2) * np.abs(dthetadt)
    return speeds


# ========================= 模块7：板凳顶点坐标（公式16-19） =========================
def vertex_coords(Pi, Pip1, l, d=D_HANDLE, w=HALF_W):
    """计算第 i 节板凳四个顶点 (Q,R,S,T)（公式16-19，模块7）。
    Pi=(xi,yi) 前把手, Pip1=(xip1,yip1) 后把手, l=把手距离, d=孔心距板头, w=半宽。
    返回 4×2 数组，顺序：Q(后右)、R(后左)、S(前左)、T(前右)。
    （这里“后/前”指板凳后端/前端，后端含后把手 P_{i+1} 一侧。）
    """
    xi, yi = Pi
    xip1, yip1 = Pip1
    dx, dy = xip1 - xi, yip1 - yi       # 轴线方向（未归一化）
    k_back = (l + d) / l                # 后端缩放系数 (l+d)/l
    k_front = d / l                     # 前端缩放系数 d/l
    ox = w * dy / l                     # 右法向偏移 x 分量
    oy = -w * dx / l                    # 右法向偏移 y 分量
    # 后端中心 = P_i + ((l+d)/l)·(P_{i+1}-P_i)
    bx = xi + k_back * dx
    by = yi + k_back * dy
    # 前端中心 = P_i - (d/l)·(P_{i+1}-P_i)
    fx = xi - k_front * dx
    fy = yi - k_front * dy
    Q = np.array([bx + ox, by + oy])   # 后右
    R = np.array([bx - ox, by - oy])   # 后左
    S = np.array([fx - ox, fy - oy])   # 前左
    T = np.array([fx + ox, fy + oy])   # 前右
    return np.array([Q, R, S, T])


def all_bench_vertices(xy, lengths=L_ARR):
    """计算全部 223 节板凳的顶点，返回 list of (4,2) 数组。"""
    return [vertex_coords(xy[i], xy[i + 1], lengths[i]) for i in range(N_BENCH)]


# ========================= 模块2：碰撞检测（线段相交法，§11.2） =========================
def cross(o, a, b):
    """向量 OA × OB 的 z 分量。"""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def on_segment(o, a, b, eps=1e-12):
    """点 o 是否在线段 ab 上（含端点），eps=1e-12。"""
    return (min(a[0], b[0]) - eps <= o[0] <= max(a[0], b[0]) + eps and
            min(a[1], b[1]) - eps <= o[1] <= max(a[1], b[1]) + eps)


def segments_intersect(A, B, C, D, eps=1e-12):
    """线段 AB 与 CD 是否相交（含端点接触与共线），碰撞用 ≤ 判定（§11.2，模块2）。
    严格相交用叉积变号判定；临界接触（端点在线段上、共线）单独判为相交。"""
    d1, d2 = cross(A, B, C), cross(A, B, D)
    d3, d4 = cross(C, D, A), cross(C, D, B)
    if ((d1 > eps and d2 < -eps) or (d1 < -eps and d2 > eps)) and \
       ((d3 > eps and d4 < -eps) or (d3 < -eps and d4 > eps)):
        return True
    if abs(d1) <= eps and on_segment(C, A, B, eps):
        return True
    if abs(d2) <= eps and on_segment(D, A, B, eps):
        return True
    if abs(d3) <= eps and on_segment(A, C, D, eps):
        return True
    if abs(d4) <= eps and on_segment(B, C, D, eps):
        return True
    return False


def benches_collide(verts_a, verts_b, eps=1e-12):
    """两节板凳是否碰撞：4×4 边组合（共16对），任一边相交即碰撞。"""
    ea = [(verts_a[i], verts_a[(i + 1) % 4]) for i in range(4)]
    eb = [(verts_b[i], verts_b[(i + 1) % 4]) for i in range(4)]
    for A, B in ea:
        for C, D in eb:
            if segments_intersect(A, B, C, D, eps):
                return True
    return False


def head_collision(thetas, p, lengths=L_ARR, eps=1e-12):
    """单时刻碰撞检查：龙头板凳（板凳0，把手0-1）是否与不相邻外圈板凳碰撞（§5.3）。
    盘入时龙头在内圈（θ最小），碰撞发生在龙头与相邻外圈（θ≈θ_0+2π）板凳之间。
    搜索窗口取 θ∈[θ_0+2π-1.5, θ_1+2π+1.5]，覆盖相邻外圈所有可能碰撞的板凳。
    跳过相邻板凳1（与龙头铰接，不判碰撞）。返回是否碰撞。
    """
    xy = handles_xy(thetas, p)
    verts = all_bench_vertices(xy, lengths)
    th0, th1 = thetas[0], thetas[1]
    lo_t, hi_t = th0 + 2 * np.pi - 1.5, th1 + 2 * np.pi + 1.5
    head = verts[0]
    for j in range(2, N_BENCH):           # 跳过相邻板凳1
        th_j, th_jp1 = thetas[j], thetas[j + 1]
        if th_jp1 < lo_t or th_j > hi_t:
            continue
        if benches_collide(head, verts[j], eps):
            return True
    return False


# ========================= 通用：单时刻状态（位置+速度） =========================
def state_at(t, p, v0=V0):
    """返回 t 时刻全部把手坐标 (224,2) 与速度 (224,)。问题1~3 通用。"""
    th_head = head_theta(t, p, v0)
    ths = position_recurse(th_head, p)
    xy = handles_xy(ths, p)
    v = velocity_recurse(ths, p, v0)
    return xy, v


def collision_at(t, p):
    """t 时刻是否发生碰撞（龙头 vs 相邻外圈）。问题2、3 通用。"""
    th_head = head_theta(t, p)
    ths = position_recurse(th_head, p)
    return head_collision(ths, p)


# ========================= 自检 =========================
if __name__ == "__main__":
    print("=== 自检：龙头速度应为 1 m/s，把手极角单调递增 ===")
    p = 0.55
    th_head = head_theta(0.0, p)
    ths = position_recurse(th_head, p)
    sp = velocity_recurse(ths, p)
    xy = handles_xy(ths, p)
    print("t=0: 龙头θ=%.6f (%.4fπ), 期望32π=%.1fπ" % (th_head, th_head / np.pi,
                                                    THETA_START / np.pi))
    print("龙头速度=%.6f m/s (应=1)" % sp[0])
    print("θ单调递增:", bool(np.all(np.diff(ths) > 0)))
    print("把手0坐标=(%.6f, %.6f) (期望约 (8.8, 0))" % (xy[0, 0], xy[0, 1]))
    print("把手1坐标=(%.6f, %.6f) (参考 (8.363824, 2.826544))" % (xy[1, 0], xy[1, 1]))
