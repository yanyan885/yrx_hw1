/*
 * @Author: 鄢瑞贤 1493754709@qq.com
 * @Date: 2025-05-9 19:31:03
 * @LastEditors: 鄢瑞贤 2324935200@qq.com
 * @LastEditTime: 2025-05-10 14:00:00
 * @FilePath: \yrx_hw1\src\algebra.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "algebra.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> // 需要 malloc, free

#define EPSILON 1e-9 // 用于浮点数比较的小值


Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows) {
        printf("Error: Matrix a's columns must match Matrix b's rows.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}


// --- 辅助函数 ---

// 创建子矩阵：移除指定行和列
Matrix create_submatrix(Matrix a, int skip_row, int skip_col) {
    if (a.rows <= 0 || a.cols <= 0) {
         return create_matrix(0,0); // 不能从空矩阵创建子矩阵
    }
    if (a.rows == 1 || a.cols == 1) {
        // 从 1xN 或 Mx1 矩阵移除一行一列会得到空矩阵或无效矩阵
        // 根据行列式/伴随矩阵定义，1x1 矩阵的代数余子式涉及 0x0 矩阵，其行列式为 1
        return create_matrix(0, 0);
    }
    Matrix sub = create_matrix(a.rows - 1, a.cols - 1);
    int current_row = 0;
    for (int i = 0; i < a.rows; i++) {
        if (i == skip_row) continue;
        int current_col = 0;
        for (int j = 0; j < a.cols; j++) {
            if (j == skip_col) continue;
            // 检查边界，防止访问越界 (虽然理论上不应发生)
            if (current_row < sub.rows && current_col < sub.cols) {
                 sub.data[current_row][current_col] = a.data[i][j];
            }
            current_col++;
        }
        current_row++;
    }
    return sub;
}

// 复制矩阵
Matrix copy_matrix(Matrix a) {
    Matrix copy = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            copy.data[i][j] = a.data[i][j];
        }
    }
    return copy;
}

// 交换矩阵的两行
void swap_rows(Matrix *m, int r1, int r2) {   
        if (!m || r1 < 0 || r2 < 0 || r1 >= m->rows || r2 >= m->rows || r1 == r2) return;
        for (int j = 0; j < m->cols; j++) {
            double tmp = m->data[r1][j];
            m->data[r1][j] = m->data[r2][j];
            m->data[r2][j] = tmp;
        }
}



double det_matrix(Matrix a) // 矩阵行列式
{
    if (a.rows != a.cols) {
        printf("Error: Determinant can only be calculated for square matrices.\n");
        return 0;
    }
    int n = a.rows;

    // 基本情况
    if (n == 0) return 1.0; // 0x0 矩阵行列式定义为 1
    if (n == 1) return a.data[0][0];
    // n=2 的情况可以作为优化，但递归也能处理
    // if (n == 2) {
    //     return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    // }

    // 递归步骤：沿第一行展开 (i=0)
    double det = 0.0;
    for (int j = 0; j < n; j++) {
        Matrix submatrix = create_submatrix(a, 0, j);
        double sub_det = det_matrix(submatrix); // 递归计算子行列式
        det += pow(-1.0, 0 + j) * a.data[0][j] * sub_det;
        
    }
    return det;
}

Matrix inv_matrix(Matrix a)// 矩阵求逆
{
    if (a.rows != a.cols) {
        printf("Error: Inverse can only be calculated for square matrices.\n");
        return create_matrix(0, 0);
    }
    int n = a.rows;
    if (n == 0) {
         printf("Error: Cannot invert an empty matrix.\n");
         return create_matrix(0,0);
    }

    double det = det_matrix(a);

    // 检查行列式是否接近零
    if (fabs(det) < EPSILON) {
        printf("Error: Matrix is singular (determinant is approximately zero) and cannot be inverted.\n");
        return create_matrix(0, 0);
    }

    // 计算伴随矩阵 A*
    Matrix adjugate = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // 计算代数余子式 C_ji = (-1)^(j+i) * M_ji
            Matrix submatrix = create_submatrix(a, j, i); // 注意是 A_ji
            double cofactor = pow(-1.0, i + j) * det_matrix(submatrix);
            adjugate.data[i][j] = cofactor; // 伴随矩阵元素 a*_ij = C_ji
            
        }
    }

    // 逆矩阵 = (1/行列式) * 伴随矩阵
    Matrix inverse = scale_matrix(adjugate, 1.0 / det);


    return inverse;
}

int rank_matrix(Matrix a)
{
    int rows = a.rows;
    int cols = a.cols;
    if (rows == 0 || cols == 0) return 0;

    // 创建一个副本进行高斯消元，不修改原矩阵
    Matrix temp_matrix = copy_matrix(a);

    int rank = 0;
    int pivot_row = 0; // 当前处理的主元行

    // 遍历所有列（或直到所有行都成为主元行）
    for (int j = 0; j < cols && pivot_row < rows; j++) {
        // 1. 选主元：在当前列(j)的 pivot_row 及下方行中，找到绝对值最大的元素
        int max_row = pivot_row;
        for (int k = pivot_row + 1; k < rows; k++) {
            if (fabs(temp_matrix.data[k][j]) > fabs(temp_matrix.data[max_row][j])) {
                max_row = k;
            }
        }

        // 2. 检查主元是否过小（视为零）
        if (fabs(temp_matrix.data[max_row][j]) < EPSILON) {
            // 当前列在 pivot_row 及下方没有非零元素，跳到下一列
            continue;
        }

        // 3. 交换行：将找到的主元行换到 pivot_row
        swap_rows(&temp_matrix, pivot_row, max_row);

        // 4. 消元：将 pivot_row 下方所有行的第 j 列元素变为 0
        // (可选：可以将主元行的主元变为1，但仅为求秩则非必需)
        // double pivot_value = temp_matrix.data[pivot_row][j]; // 主元值
        for (int i = pivot_row + 1; i < rows; i++) {
            double factor = temp_matrix.data[i][j] / temp_matrix.data[pivot_row][j];
            // 对该行从第 j 列开始的所有元素进行操作
            temp_matrix.data[i][j] = 0.0; // 直接设为0，避免精度问题
            for (int k = j + 1; k < cols; k++) {
                temp_matrix.data[i][k] -= factor * temp_matrix.data[pivot_row][k];
            }
        }

        // 5. 移动到下一主元行，秩增加
        pivot_row++;
        rank++;
    }

    return rank; // rank 即为主元的数量
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols) {
        printf("Error: Trace can only be calculated for square matrices.\n");
        return 0;
    }
    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}