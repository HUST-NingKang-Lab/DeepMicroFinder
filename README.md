# DeepMicroFinder
DeepMicroFinder is a deep learning framework, which integrated the neural network and transfer learning and could effectively reduce the regional effects for microbial-based cross-regional diagnosis of T2D.

<img src="https://github.com/HUST-NingKang-Lab/DeepMicroFinder/blob/main/figure1.png" style="zoom:150%;" />

``
Samples of each cohort were randomly divided into the training subset and the testing subset, and four models were constructed for assessment: (1) Independent disease neural network (DNN) model: ab initio training the DNN model on the training subset and testing on the testing subset of SGMP cohort, respectively. (2) Regional DNN model: ab initio training the DNN model using the training subset of GGMP cohort and testing it on the testing subset of SGMP cohort. (3) Regional+ DNN model: ab initio training the DNN model using the training subset of GGMP cohort as well as the training subset of SGMP cohort, then testing it on the testing subset of SGMP cohort. (4) Transfer DNN model: ab initio training the DNN model using the training subset of GGMP cohort, followed by applying transfer learning to a certain proportion (from 20% to 80%) of samples from SGMP cohort to generate the transfer DNN model, and then testing the transfer DNN model on the testing subset of SGMP cohort. The three boxes on the right represent the evaluation and applications of DeepMciroFinder, including cross-regional diagnosis of T2D and biomarker discovery.
``

## Get and use
To learn how to install the model and how to use it, click [here](https://github.com/HUST-NingKang-Lab/EXPERT)
