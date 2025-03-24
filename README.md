# DeepMicroFinder
DeepMicroFinder is a deep learning framework, which integrated the neural network and transfer learning and could effectively reduce the regional effects for microbial-based cross-regional diagnosis of T2D.

<img src="https://github.com/HUST-NingKang-Lab/DeepMicroFinder/blob/main/figure1.png" style="zoom:150%;" />

``
Samples of each cohort were randomly divided into the training subset and the testing subset, and four models were constructed for assessment: (1) Independent disease neural network (DNN) model: ab initio training the DNN model on the training subset and testing on the testing subset of SGMP cohort, respectively. (2) Regional DNN model: ab initio training the DNN model using the training subset of GGMP cohort and testing it on the testing subset of SGMP cohort. (3) Regional+ DNN model: ab initio training the DNN model using the training subset of GGMP cohort as well as the training subset of SGMP cohort, then testing it on the testing subset of SGMP cohort. (4) Transfer DNN model: ab initio training the DNN model using the training subset of GGMP cohort, followed by applying transfer learning to a certain proportion (from 20% to 80%) of samples from SGMP cohort to generate the transfer DNN model, and then testing the transfer DNN model on the testing subset of SGMP cohort. The three boxes on the right represent the evaluation and applications of DeepMciroFinder, including cross-regional diagnosis of T2D and biomarker discovery.
``

## Get and use
To learn how to install the model and how to use it, click [here](https://github.com/HUST-NingKang-Lab/EXPERT)

## Example
The example data for DeepMicroFinder:

Species abundance tables([reference format](https://github.com/HUST-NingKang-Lab/EXPERT)): shandong_train.tsv shandong_test.tsv       
Disease models: GGMP_Disease_Model.h5     

#### Transfer learning
- Transfer the knowledge of shandong to the guandong DNN model for better performance in disease diagnosis on shandong. You'll see running log and training process in the printed message.
```
expert transfer -i shandong_trainCM.h5 -l shandong_train_labels.h5 -t ontology.pkl -m GGMP_Disease_Model.h5  -o SGMP_Disease_Model 
```
#### Search
- Search the test set of shandong against the transferred DNN model.
```
expert search -i shandong_testCM.h5 -m SGMP_Disease_Model.h5 -o Search_shandong
```
#### Evaluation
- Evaluate the performance of the Transferred DNN model. You'll obtain a performance report.
```
expert evaluate -i Search_shandong -l shandong_test_labels.h5 -o Evaluation
```
