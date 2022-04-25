# Optimal-Power-Flow
Solve optimal power flow problem by Gurobi.

<img src="https://github.com/yijunchengg/Optimal-Power-Flow/blob/main/power%20flow.png">


The constraints (1b) and (1c) are for each bus. The power loss is calculated for each brunch as (1c).
The objective as (1a) is the sum of the cost function of generator _Pg_ and substation _Ps_, and the power loss.

## Data
The branch data and load data of 4-bus system is as follows,

```
branch = pd.DataFrame(columns=('bus1', 'bus2', 'resistance', 'reactance'),
                      data=[
                      [1,2,0.0922,0.047],
                      [2,3,0.493,0.2511],
                      [3,4,0.366,0.1864],
                      ])

load = pd.DataFrame(columns=('activepower', 'reactivepower'),
                    data=[
                    [0,0],
                    [0.1,0.06],
                    [0.09,0.04],
                    [0.12,0.08],
                    ])
```
Bus 1 is the slack bus, and there is a generator at bus 4. They can provide active power and reactive power.

## Calculate the admittance matrix Y

## Create the optimization model
For each bus, there are four variables, which are active power _P_ and reactive power _Q_, voltage magnitude |_V_| and phase angle <img src="https://render.githubusercontent.com/render/math?math=\delta">.

We need to introduce new varibales to support constraints containing more general multilinear terms. Further information: [How do I model multilinear terms in Gurobi?](https://support.gurobi.com/hc/en-us/articles/360049744691).


## Results
Here are the optimal solution.
```
substation[0] 0.03801858663939228
substation[1] 0.1942591006904868
active power of generator 0.3
reactive power of generator 1.61247068928283e-05
voltage magnitude of load[0] 1.0
voltage magnitude of load[1] 0.9875049256162309
voltage magnitude of load[2] 0.990141853291291
voltage magnitude of load[3] 1.0374120696892812
phase angle of load[0] 0.0
phase angle of load[1] 0.01632876407190631
phase angle of load[2] 0.10003770115765243
phase angle of load[3] 0.1612398131926849
running time: 188.73901224136353
active generation 0.33801858663939227
reactive generation 0.19427522539737963
active load 0.31
reactive load 0.18
```
