> First trial with jupyter notebook.

# Mean field solution on spin-1/2 bond-bond interaction Hamiltonian

Outline:
- HS decomposition; Mean-field approximation;
- Define the order parameter
- Programming it.

## Mean field approach in lattice space
First, the triplet current order in the mean field approximation is

![](https://latex.codecogs.com/gif.latex?\vec%10N_{ij}=i\langle%10c_{i\sigma}^{\dagger}(\frac{\vec\sigma}{2})_{\sigma\sigma'}c_{j\sigma'}-H.c.\rangle)

The restricted mean-field staggered spin flux order parameter N is defined as

![](https://latex.codecogs.com/gif.latex?\vec%10N_{ij}=2×N(-1)^{i_x+i_y}(\delta_{j-i,\hat%10x}-\delta_{j-i,\hat%10y})\hat%10z\tag1)

This represents a flux phase. This phase breaks the parity, translation by one site and conjugation symmetry.
### 1. Order parameter against interaction strength
![](./fig/OP_gtcur.png)

### 2. Order parameter against doping strength
![](./fig/Dope_gtcur.png)

---
The order parameter M is defined in the way,

![](https://latex.codecogs.com/gif.latex?M_{ij}=\langle%10c_{i\sigma}^{\dagger}c_{j\sigma}+H.c.\rangle)

and for the mean field order parameter, I chose

![](https://latex.codecogs.com/gif.latex?M_{ij}=M(\delta_{j-i,\hat%10x}+\delta_{j-i,\hat%10y})\tag2)

but it doesn't work well in the simulation. So I digged some old papers which discussed the similar model.

In Ref. 3, the SU(n) Hubbard-Heisenberg model can be restored to our SU(2) model when n=2, and the ground state is constructed out by the valence bonds. The normal site-centered charge density wave becomes the bond-centered density wave. However, the nearest neighber state is formed under the condition that n is approching infinity. Those processes that create non-nearest-neighbor bonds have amplitudes which are of order 1/n. The singlet bond operators also correspond to a bond-centered spin-density wave(the Heisenberg exchange term),

![](https://latex.codecogs.com/gif.latex?\langle%10\mid%10c_{i\sigma}^{\dagger}c_{j\sigma}\mid^2\rangle\propto\langle%10S_i\cdot%10S_j\rangle)

maybe the order parameter in mean-field calculation should be defined on each bond,

![](https://latex.codecogs.com/gif.latex?M_{ij}=(M_x\delta_{j-i,\hat%10x}+M_y\delta_{j-i,\hat%10y})\tag3)

and the new unit cell consists of the four sites at the corners of a square in the origin lattice.
## Useful advices:

- [你为什么使用 jupyter ，进行分析，而不是用 python 脚本或仅仅利用 excel ？](https://www.zhihu.com/question/37490497)
- [Git基本操作](http://www.runoob.com/git/git-basic-operations.html)
- [Emacs: Python最好的编译器？](https://segmentfault.com/a/1190000004165173)
