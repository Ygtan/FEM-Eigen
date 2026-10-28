// 解析解断言测试公共辅助：test_2d.cpp / test_3d.cpp 共用。
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

inline int g_pass = 0;
inline int g_fail = 0;

void check(const std::string& name, double got, double exact, double tol = 1e-8) {
    const double denom = std::max(1.0, std::abs(exact));
    const double rel = std::abs(got - exact) / denom;
    if (rel <= tol) {
        ++g_pass;
        std::printf("  [PASS] %-44s got=%14.8g  exact=%14.8g  rel=%.2e\n",
                    name.c_str(), got, exact, rel);
    } else {
        ++g_fail;
        std::printf("  [FAIL] %-44s got=%14.8g  exact=%14.8g  rel=%.2e (tol=%.0e)\n",
                    name.c_str(), got, exact, rel, tol);
    }
}
