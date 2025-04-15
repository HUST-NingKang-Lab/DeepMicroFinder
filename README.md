# DeepMicroFinder
DeepMicroFinder is a deep learning framework, which integrated the neural network and transfer learning and could effectively reduce the regional effects for microbial-based cross-regional diagnosis of T2D.

- /fig1 # The figure for illustrating the framework of DeepMicroFinder
- /fig2 # The data for conducting the GBTM grouping trajectory
- /model # The DNN models of GGMP and SGMP, and the ontology file
- supplementary_figures_codes # General R scripts used in this study

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

## Maintainer

|   Name    | Email                 | Organization                                                 |
| :-------: | --------------------- | ------------------------------------------------------------ |
| Nan Wang | wangnan123@hust.edu.cn | Phd student, School of Life Science and Technology, Huazhong University of Science & Technology|
| Kang Ning | ningkang@hust.edu.cn  | Professor, School of Life Science and Technology, Huazhong University of Science & Technology |
