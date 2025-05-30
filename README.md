# 矩阵运算程序设计文档

## 功能模块实现

### 行列式计算（`det_matrix`）
- **方法**：采用递归展开法
- **特殊情况处理**：
  - 对于 0x0 矩阵，行列式值为 1
  - 1x1 矩阵的行列式就是其唯一元素值
- **计算流程**：
  1. 验证矩阵为方阵
  2. 对第 0 行的每个元素 `a.data[0][j]`：
     - 生成移除第 0 行和第 j 列的子矩阵 `submatrix`
     - 递归调用 `det_matrix(submatrix)` 计算子行列式 `sub_det`
     - 通过公式累加：`det += (-1)^(0 + j) * a.data[0][j] * sub_det`
- **公式**：
  $$\text{det}(A) = \sum_{j = 0}^{n - 1} (-1)^{0 + j} \cdot a_{0, j} \cdot \text{det}(A_{0, j})$$

### 逆矩阵计算（`inv_matrix`）
- **方法**：使用伴随矩阵法
- **计算流程**：
  1. 计算行列式 `det`
  2. 若 `|det|` 近似为 0，矩阵奇异不可逆，返回 0x0 矩阵
  3. 否则遍历矩阵每个位置 `(i, j)`：
     - 计算代数余子式：$C_{ji} = (-1)^{j + i} \cdot \text{det}(A_{ji})$（`A_{ji}` 是移除第 j 行和第 i 列的子矩阵）
     - 赋值到伴随矩阵 `adjugate.data[i][j]`
  4. 将伴随矩阵乘以 $1.0 / det$ 得到逆矩阵
- **公式**：
  $$A^{-1} = \frac{1}{\text{det}} \cdot \text{adjugate}(A)$$

### 矩阵秩计算（`rank_matrix`）
- **方法**：列主元高斯消元法
- **计算流程**：
  1. 对每列 `j`：
     - 从 `pivot_row` 开始向下找绝对值最大元素所在行 `max_row`
     - 若该列 `pivot_row` 下方无主元（最大值近似为 0）则跳过
     - 交换 `max_row` 与 `pivot_row` 行
     - 对 `pivot_row` 下方的每行 `i`：
       - 计算因子 `factor = a[i][j] / a[pivot_row][j]`
       - 执行行操作：$R_i = R_i - factor \cdot R_{pivot_row}$
       - （为避免精度误差，可直接将操作后的元素设为 0.0 并对后续列执行减法操作）
  2. 每成功消去一列，秩 `rank` 加 1

## 技术难点与解决方案

### 静态数组行交换问题
- **问题描述**：
  - 使用静态二维数组 `double data[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE]` 存储矩阵时
  - 代码 `m->data[r1] = m->data[r2];` 会报错"表达式必须是可修改的左值"
  - 原因：静态数组的行地址不可直接赋值
- **解决方案**：
  - 不使用 `swap_rows` 交换指针
  - 改为逐元素交换行数据
### 还有个死int因为c89的问题天天阴我
### 差点忘了，还有那个输出，那个流程也是先清空这个再清空那个的，燃尽了
![加](picture/jia.png)

![减](picture/jian.png)
![叉乘](picture/cha.png)
![行列式](picture/hang.png)
![迹](picture/ji.png)
![数乘](picture/shu.png)
![秩序](picture/zhi.png)
![转置](picture/zhuan.png)
![求逆](picture/ni.png)