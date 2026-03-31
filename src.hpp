#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
#define IGNORE_MATRIX 0
// #define IGNORE_MATRIX 1

#if !IGNORE_MATRIX

class matrix {
private:

    // m行n列的矩阵，用动态二维数组存储，每个元素是分数类实例
    int m, n;
    fraction **data;

public:

    // 默认构造函数
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    // 构造函数，构建 m_*n_ 的矩阵，矩阵元素设为0。
    matrix(int m_, int n_) {
        m = m_;
        n = n_;
        if (m <= 0 || n <= 0) {
            data = nullptr;
            return;
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = fraction(0);
            }
        }
    }

    // 拷贝构造函数，构建与 obj 完全相同的矩阵。
    matrix(const matrix &obj) {
        m = obj.m;
        n = obj.n;
        if (m <= 0 || n <= 0 || obj.data == nullptr) {
            data = nullptr;
            return;
        }
        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
    }

    // 移动拷贝构造函数。
    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    // 析构函数。
    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    // 重载赋值号。
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;

        // 释放原有内存
        if (data != nullptr) {
            for (int i = 0; i < m; i++) {
                delete[] data[i];
            }
            delete[] data;
        }

        m = obj.m;
        n = obj.n;
        if (m <= 0 || n <= 0 || obj.data == nullptr) {
            data = nullptr;
            return *this;
        }

        data = new fraction*[m];
        for (int i = 0; i < m; i++) {
            data[i] = new fraction[n];
            for (int j = 0; j < n; j++) {
                data[i][j] = obj.data[i][j];
            }
        }
        return *this;
    }

    // 重载括号，返回矩阵的第i行(1-based)、第j列(0-based)的元素的引用。如果 i、j 不合法，抛出 matrix_error 错误。
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) {
            throw matrix_error();
        }
        return data[i-1][j];
    }

    // 重载乘号，返回矩阵乘法 lhs * rhs 的结果。如果 lhs 的列数与 rhs 的行数不相等，抛出 matrix_error 错误。
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m || lhs.data == nullptr || rhs.data == nullptr) {
            throw matrix_error();
        }
        matrix result(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; i++) {
            for (int j = 0; j < rhs.n; j++) {
                result.data[i][j] = fraction(0);
                for (int k = 0; k < lhs.n; k++) {
                    result.data[i][j] = result.data[i][j] + lhs.data[i][k] * rhs.data[k][j];
                }
            }
        }
        return result;
    }

    // 返回矩阵的转置。若矩阵为空，抛出 matrix_error 错误。
    matrix transposition() {
        if (data == nullptr || m <= 0 || n <= 0) {
            throw matrix_error();
        }
        matrix result(n, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // 返回矩阵的行列式。建议用高斯消元实现。若矩阵不是方阵或为空，抛出 matrix_error 错误。
    fraction determination() {
        if (data == nullptr || m != n || m <= 0) {
            throw matrix_error();
        }

        // 创建矩阵副本用于高斯消元
        matrix temp(*this);
        fraction det(1);

        // 高斯消元
        for (int i = 0; i < m; i++) {
            // 寻找主元
            int pivot = i;
            for (int j = i + 1; j < m; j++) {
                if (temp.data[j][i] == fraction(0)) continue;
                if (temp.data[pivot][i] == fraction(0) ||
                    !(temp.data[j][i] == fraction(0))) {
                    pivot = j;
                }
            }

            // 如果主元为0，行列式为0
            if (temp.data[pivot][i] == fraction(0)) {
                return fraction(0);
            }

            // 交换行
            if (pivot != i) {
                for (int j = 0; j < m; j++) {
                    fraction tmp = temp.data[i][j];
                    temp.data[i][j] = temp.data[pivot][j];
                    temp.data[pivot][j] = tmp;
                }
                det = det * fraction(-1);
            }

            // 消元
            for (int j = i + 1; j < m; j++) {
                if (!(temp.data[j][i] == fraction(0))) {
                    fraction factor = temp.data[j][i] / temp.data[i][i];
                    for (int k = i; k < m; k++) {
                        temp.data[j][k] = temp.data[j][k] - factor * temp.data[i][k];
                    }
                }
            }
        }

        // 计算对角线元素的乘积
        for (int i = 0; i < m; i++) {
            det = det * temp.data[i][i];
        }

        return det;
    }

    // 辅助函数：获取去除第 row 行和 col 列后的子矩阵 (1-based indices)
    matrix submatrix(int row, int col) const {
        if (data == nullptr || m <= 1 || n <= 1 || row < 1 || row > m || col < 1 || col > n) {
            throw matrix_error();
        }
        matrix result(m - 1, n - 1);
        int ri = 0;
        for (int i = 0; i < m; i++) {
            if (i == row - 1) continue;
            int rj = 0;
            for (int j = 0; j < n; j++) {
                if (j == col - 1) continue;
                result.data[ri][rj] = data[i][j];
                rj++;
            }
            ri++;
        }
        return result;
    }

    // 辅助函数：获取去除第 row1, row2 行和 col1, col2 列后的子矩阵 (1-based indices)
    matrix submatrix2(int row1, int row2, int col1, int col2) const {
        if (data == nullptr || m <= 2 || n <= 2) {
            throw matrix_error();
        }
        matrix result(m - 2, n - 2);
        int ri = 0;
        for (int i = 0; i < m; i++) {
            if (i == row1 - 1 || i == row2 - 1) continue;
            int rj = 0;
            for (int j = 0; j < n; j++) {
                if (j == col1 - 1 || j == col2 - 1) continue;
                result.data[ri][rj] = data[i][j];
                rj++;
            }
            ri++;
        }
        return result;
    }

    // 辅助函数：用向量替换第col列，然后去除第row行和第delcol列 (1-based indices)
    matrix replace_col_and_remove(int col, fraction* vec, int vecsize, int row, int delcol) const {
        if (data == nullptr || col < 1 || col > n || row < 1 || row > m || delcol < 1 || delcol > n) {
            throw matrix_error();
        }

        // 先替换列
        matrix temp(*this);
        for (int i = 0; i < m && i < vecsize; i++) {
            temp.data[i][col - 1] = vec[i];
        }

        // 然后去除行和列
        return temp.submatrix(row, delcol);
    }
};

#endif

class resistive_network {
private:

    // 节点数量 和 接线数量
    int interface_size, connection_size;

    // 邻接矩阵A，电导矩阵C，Laplace矩阵(A^tCA) (具体定义见 reading_materials.pdf)
    matrix adjacency, conduction, laplace;

public:

    // 设置电阻网络。节点数量为interface_size_，接线数量为connection_size_。
    //       对于 1<=i<=connection_size_，从节点from[i-1]到节点to[i-1]有接线，对应电阻为resistance[i-1]。
    //       保证接线使得电阻网络联通，from[i-1] < to[i-1]，resitance[i-1] > 0，均合法。
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;

        // 构建邻接矩阵 A (m x n)
        adjacency = matrix(connection_size, interface_size);
        for (int i = 1; i <= connection_size; i++) {
            int from_node = from[i - 1];
            int to_node = to[i - 1];
            adjacency(i, from_node - 1) = fraction(1);
            adjacency(i, to_node - 1) = fraction(-1);
        }

        // 构建电导矩阵 C (m x m)
        conduction = matrix(connection_size, connection_size);
        for (int i = 1; i <= connection_size; i++) {
            conduction(i, i - 1) = fraction(1) / resistance[i - 1];
        }

        // 构建 Laplace 矩阵 L = A^T * C * A
        matrix A_T = adjacency.transposition();
        matrix temp = A_T * conduction;
        laplace = temp * adjacency;
    }

    ~resistive_network() = default;

    // 返回节点 interface_id1 和 interface_id2 (1-based)之间的等效电阻。
    //       保证 interface_id1 <= interface_id2 均合法。
    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) {
            return fraction(0);
        }

        // R_ij = det(M_ij) / det(M_i)
        // M_i: 去除第i行、第i列
        // M_ij: 去除第i行、第i列、第j行、第j列

        matrix M_i = laplace.submatrix(interface_id1, interface_id1);
        fraction det_M_i = M_i.determination();

        matrix M_ij = laplace.submatrix2(interface_id1, interface_id2, interface_id1, interface_id2);
        fraction det_M_ij = M_ij.determination();

        return det_M_ij / det_M_i;
    }

    // 在给定节点电流I的前提下，返回节点id(1-based)的电压。认为节点interface_size(1-based)的电压为0。
    //       对于 1<=i<=interface_size，节点i(1-based)对应电流为 current[i-1]。
    //       保证 current 使得电阻网络有解，id < interface_size 合法。
    fraction get_voltage(int id, fraction current[]) {
        // u_i = det(M_i^n) / det(M_n)
        // M_n: 去除第n行、第n列
        // M_i^n: 将A^TCA的第i列替换成I，然后去除第n行、第n列

        matrix M_n = laplace.submatrix(interface_size, interface_size);
        fraction det_M_n = M_n.determination();

        matrix M_i_n = laplace.replace_col_and_remove(id, current, interface_size, interface_size, interface_size);
        fraction det_M_i_n = M_i_n.determination();

        return det_M_i_n / det_M_n;
    }

    // 在给定节点电压U的前提下，返回电阻网络的功率。
    //       对于 1<=i<=interface_size，节点i (1-based) 对应电压为 voltage[i-1]。
    //       保证 voltage 合法。
    fraction get_power(fraction voltage[]) {
        // P = Σ(u_wi^2 / r_i)
        // u_w = A * U^T

        // 构建电压向量矩阵 (n x 1)
        matrix U(interface_size, 1);
        for (int i = 1; i <= interface_size; i++) {
            U(i, 0) = voltage[i - 1];
        }

        // 计算 u_w = A * U
        matrix u_w = adjacency * U;

        fraction power(0);
        for (int i = 1; i <= connection_size; i++) {
            fraction u_wi = u_w(i, 0);
            fraction r_i = fraction(1) / conduction(i, i - 1);
            power = power + (u_wi * u_wi) / r_i;
        }

        return power;
    }
};


#endif //SRC_HPP
