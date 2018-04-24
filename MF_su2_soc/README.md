> First trial with jupyter notebook.

# Mean field solution on spin-$\frac{1}{2}$ bond-bond interaction Hamiltonian

Outline:
- HS decomposition; Mean-field approximation;
- Define the order parameter
- Programming it.

## Mean field approach in lattice space
The restricted mean-field staggered spin flux order parameter $N$ is defined as
$$ \vec N_{ij}=2×N(-1)^{i_x+i_y}(\delta_{j-i,\hat x}-\delta_{j-i,\hat y} )\hat z \tag 1$$
### 1. Order parameter against interaction strength
![](./fig/OP_gtcur.png)

### 2. Order parameter against doping strength
![](./fig/Dope_gtcur.png)


## Useful advices:

- [你为什么使用 jupyter ，进行分析，而不是用 python 脚本或仅仅利用 excel ？](https://www.zhihu.com/question/37490497)
- [Git基本操作](http://www.runoob.com/git/git-basic-operations.html)
