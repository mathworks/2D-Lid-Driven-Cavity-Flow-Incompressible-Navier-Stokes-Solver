# 2D Lid Driven Cavity Flow
[![View 2D-Lid-Driven-Cavity-Flow---Incompressible-Navier-Stokes-Sol on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74483-2d-lid-driven-cavity-flow-incompressible-navier-stokes-sol)

This repo provides a MATLAB example code for the lid-driven cavity flow where incompressible 
Navier Stokes equation is numerically solved using a simple 2nd order finite difference scheme on a staggered grid system.

![sample](./gif/animation_sample1e2and1e4.gif)
(Left: Re = 100, Right: Re = 10,000)

The arrow denotes the velocity field, and the contour denotes its magnitude.

## Part 1


- Click [here](./docs_part1/vanilaCavityFlow_EN.md) for detailed documentation in English.
- 日本語のドキュメントは[こちら](./docs_part1/vanilaCavityFlow_JP.md) から

The numerical scheme is kept primitive; the explicit treatment of viscous term (the solution diverges at low Reynolds number), and the time integration is Euler.

まずは単純な手法でキャビティ流れのシミュレーションを実施します。

## Part 2

- Click [here](./docs_part2/vanilaCavityFlowImplicit_EN.md) for detailed documentation in English.
- 日本語のドキュメントは[こちら](./docs_part2/vanilaCavityFlowImplicit_JP.md) から

The implicit treatments for viscous terms are implemented at low Reynolds number, namely the Crank-Nicolson method. For better stability for non-linear terms, Adams-Bashforth, and 3 steps-Runge-Kutta is also implemented. 

拡散項に対して陰解法を実装しました。対流項へアダムス・バッシュフォースを使用したもの、3段階のルンゲクッタ法の２つの時間発展を実装しています。


## Next to come

The plan is to allow arbitrary boundary conditions for more fun simulations.




# Environment

- MATLAB R2019b
- Signal Processing Toolbox if you use dct in solving Poisson eqn.

# ToDo

1. Implement implicit treatment of viscous terms
2. Implement crank-Nicolson for the non-linear terms
3. Allow obstacles within the domain
4. Allow inflow from the wall
5. Make it to 3D

--

Copyright (c) 2020, The MathWorks, Inc.
