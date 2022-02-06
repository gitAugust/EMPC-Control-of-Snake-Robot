# EMPC-Control-of-Snale-Robot
Reproduction of paper ['Economic model predictive control for snake robot locomotion'](https://ieeexplore.ieee.org/abstract/document/9029627)

## Run Dependencies
this code is relay on the [casadi-windows-matlabR2016a-v3.5.5](https://github.com/casadi/casadi/releases) matlab toolbox, to run this code the toolbox should be installed to your matlab of place in the same path of the code as an .zip file.

## The state of the snake robot is setted as following:
Nl is the number of the links of snake robot.

| $\Phi(t)$ | $\Theta(t)$ | $\P_x(t)$ | $\P_y(t)$ | $\dot{Phi(t)}$ | $\dot{Phi(t)}$ | $\dot{\Theta(t)}$ | $\dot{\P_x(t)}$ | $\dot{\P_y(t)}$ |
| :-----: |:--:| :--: | :--: | :--: | :--------: | :---: | :---: | :--: |
| 1：Nl-1 | Nl | Nl+1 | Nl+2 | Nl+3 | Nl+3：2Nl+1 | 2Nl+2 | 2Nl+3 |2Nl+4 |


## *Reference*  
[1]Nonhoff M, Köhler P N, Kohl A M, et al. Economic model predictive control for snake robot locomotion[C] 2019 IEEE 58th Conference on Decision and Control. IEEE, 2019: 8329-8334.  
[2]Liljebäck P, Pettersen K Y, Stavdahl Ø, et al. Snake robots: modelling, mechatronics, and control[M]. Springer Science & Business Media, 2012.
