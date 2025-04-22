# Networked-Sensing-PCRB-Analysis

# What is this?

This project contains scripts to reproduce experiments from the paper
"[Joint Transmission and Compression Optimization for Networked Sensing with Limited-Capacity Fronthaul Links](https://ieeexplore.ieee.org/document/10948152)"
by 
[Weifeng Zhu](mailto://eee-wf.zhu@polyu.edu.hk)
,
[Shuowen Zhang](mailto://shuowen.zhang@polyu.edu.hk)
, 
and [Liang Liu](mailto://liangeie.liu@polyu.edu.hk).
Published in IEEE Transactions on Wireless Communications.
See also the related [preprint](https://arxiv.org/abs/2408.03174).

If this project help your researches, it is appreciated that this work is added in the reference of your papers.

# The Problem of Interest

Briefly, the Posterior Cramer-Rao Bound Analysis and Optimization (transmit and quantization covariance matrix) in the networked sensing system with fronthaul capacity limit are provided. 

# Overview

The included scripts 
- are generally written in matlab and require cvx.


# Description of Files

## [PCRB_MCJBQ_opt.m](PCRB_MCJBQ_opt.m) 

Main function to optimize the two kinds of covariance matrices.
