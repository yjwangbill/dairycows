稀疏众数可加模型

*****************************************************************************************************************
* This toolbox is designed to work with Matlab 2018b   *
*********************************************************
------------------------------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides an efficient way to learn the SpMAM. 
------------------------------------------------------------------------------------------------------------------------------------------------
Specifications for Using it

One demo file 'Demo_SpMAM.m' is proposed to illustrate the principle of the method with dynamic displays     
-----------------------------------------------------------------------------------------------------------------------------------------------


To evaluate the model(s) under Example A (Chi-squared noise), run Demo_SpMAM.m;

Note: we here use a hard-threshold to screen out the informative variable for simplicity. An effective way to select a stability-inducing threshold is proposed in [Sun et al 2013]. 

Some parameters:
n                      the sample size
p                      the input dimension 
lambda                 the regularization parameter                
norm      = 1 or 2     1 means SpMAM with \ell_1-norm; 2 means SpMAM with \ell_2-norm 
para.kerOpt            the modal kernel
-----------------------------------------------------------------------------------------------------------------------------------------------
Results: 

ASE: 0.2662
----------------------------------------------------------------------------------------------------------------------------------------------

Ref:
Sun W, Wang J, and Fang Y. Consistent selection of tuning parameters via variable selection stability. Journal of Machine Learning Research, 14(9):3419–3440, 2013.