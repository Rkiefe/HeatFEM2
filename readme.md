# HeatFEM
This is a finite element implementation of the 2D heat equation, with Neumann boundary conditions.

## Model
The heat equation is as follows,

$\rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k\nabla T) = \frac{d q_V}{dt}$

With Neumann boundary conditions: 

$k\nabla T \cdot \vec{n} = h(T_{ext}-T)$


![Example](https://github.com/user-attachments/assets/199e0c59-01db-4181-b984-45a5d7fb5712)
