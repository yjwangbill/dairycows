稀疏众数可加模型

*****************************************************************************************************************
* This toolbox is designed to work with Matlab 2018b   *
*********************************************************
------------------------------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides an efficient way to learn the modal regression. 
------------------------------------------------------------------------------------------------------------------------------------------------
Specifications for Using it

One demo file 'demo_simulation.m' is proposed to illustrate the principle of the method with dynamic displays     
-----------------------------------------------------------------------------------------------------------------------------------------------
Experimental settings
   
Example A:                              Simulate_data.m                 Example 1--Function 1                               
Example B:                              Simulate_data.m                 Example 2--Function 2 
Gaussian Noise:    	               Simulate_data.m                 Y  = true_function(X)+r*normrnd(0,1,n,1);
Student   Noise:                         Simulate_data.m                 Y  = true_function(X)+r*nctrnd(3,0,n,1);
Chi-square Noise:               	       Simulate_data.m                 Y  = true_function(X)+r*random('chi',1,n,1);
------------------------------------------------------------------------------------------------------------------------------------------------

To evaluate the model(s) with q=1 in the Sect.4, run Demo_SpMAM1.m;
To evaluate the model(s) with q=2 in the Sect.4, run Demo_SpMAM2.m;
-----------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------
Results: 

Size=8   TP=7.96   FP=0.04
C=24  U=1  O=0     ASE=0.19756
----------------------------------------------------------------------------------------------------------------------------------------------