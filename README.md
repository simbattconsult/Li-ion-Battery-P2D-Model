# Open-source DFN P2D model for Li-ion cell in MATLAB/Simulink



Introduction



This repository contains the simulation code for Doyle-Fuller-Newman (DFN or P2D) model of Li-ion battery cell associated with the following research work (see the "reference" folder):



Ref 1

Guo, M., Jin, X. and White, R.E., 2017. "Nonlinear State-Variable Method (NSVM) for Li-Ion Batteries: Finite-Element Method and Control Mode". Journal of The Electrochemical Society, 164(11), p.E3200.



Ref 2

Guo, M., Jin, X. and White, R.E., 2017. "Nonlinear state-variable method for solving physics-based li-ion cell model with high-frequency inputs". Journal of The Electrochemical Society, 164(11), p.E3001.



Ref 3

Guo, M., Jin, X. and White, R.E., 2017. "An adaptive reduced-order-modeling approach for simulating real-time performances of li-ion battery systems". Journal of The Electrochemical Society, 164(14), p.A3602.



wherein the open-source code solves full-order P2D model PDEs in discrete state space.







MATLAB version



In the main directory, the script "main\_sim\_p2d\_matlab.m" runs the original MATLAB code that generates the plots in Figure 5 (a) and (b) of Ref 1. The parameters and constant matrix operators are loaded from "\\workspace\\p2d\_workspace\_LIB4.mat".


<img width="481" height="655" alt="image" src="https://github.com/user-attachments/assets/eb7310d5-69e8-48f3-9d0d-011f65259338" />

(Figure 5 (a) & (b) in paper "Nonlinear State-Variable Method (NSVM) for Li-Ion Batteries: Finite-Element Method and Control Mode". Journal of The Electrochemical Society, 164(11), p.E3200.)




Simulink version



An equivalent Simulink version of this P2D model is also included here and is implemented by the script "main\_sim\_p2d\_simulink.m". In this Simulink model, original MATLAB code are slightly modified and embedded into User-Defined-Function (UDF) blocks. The Simulink model uses explicit solution method described by Eq(36) through Eq(43) in Ref 2 with current control mode only. The constant voltage simulation is managed by an external feedback cycle controller connected to the P2D model block. The control gain signal "dVexdI" is evaluated by Eq(31) of Ref 3 plus an external resistance term "Rex". The Simulink code runs under fixed-step discrete-time simulation settings. Not just reporducing Figure 5 (a) and (b) of Ref 1, the Simulink model can also plot electrolyte concentration/potential and reaction current vs x-dim location at selected time.



The Simulink version of P2D model has the following benefits:

1\) More consistent with battery control strategy used in realistic dynamic systems

2\) Easy integration into large-scale simulation systems like MIL, HIL, SIL

3\) Convenient C/C++ code building for micro-controller or BMS execution

4\) Allowing single precision and integer datatype to save memory for hardware implementation


Contact info
Email : simbattconsult0821@outlook.com
