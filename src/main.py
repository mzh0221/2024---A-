# -*- coding: utf-8 -*-
"""
main.py —— 板凳龙问题求解主程序（问题1~5）
================================================
依序求解 2024 国赛 A 题“板凳龙闹元宵”的全部 5 个问题：

  问题1：300 s 仿真，输出 result1.xlsx      → problem1.solve_problem1()
  问题2：二分法求盘入终止时刻，输出 result2.xlsx → problem2.solve_problem2()
  问题3：二分+碰撞检验求最小螺距，输出 result3   → problem3.solve_problem3()
  问题4：调头路径建模，输出 result4.xlsx      → problem4.solve_problem4()
  问题5：缩放法求最大龙头速度，输出 result5     → problem5.solve_problem5()

运行：在 answer/ 目录下执行  python main.py
依赖：numpy、scipy、openpyxl、python-docx（仅文档读取，求解不需要）
附件：result1/2/4.xlsx 模板需位于 ../2024---A--main/A题 copy/附件/（用于输出格式）
"""

import time as _timer


def main():
    print("=" * 60)
    print("2024 国赛 A 题“板凳龙”求解（问题1~5）")
    print("=" * 60)
    t_total = _timer.time()

    import problem1
    import problem2
    import problem3
    import problem4
    import problem5

    problem1.solve_problem1()
    problem2.solve_problem2()
    problem3.solve_problem3()
    problem4.solve_problem4()
    problem5.solve_problem5()

    print("\n" + "=" * 60)
    print("全部问题求解完成，总用时 %.1fs" % (_timer.time() - t_total))
    print("结果文件 result1.xlsx ~ result5.xlsx 已保存到 answer/ 目录")
    print("=" * 60)


if __name__ == "__main__":
    main()
