# DeepMicroFinder
DeepMicroFinder is a deep learning framework, which integrated the neural network and transfer learning and could effectively reduce the regional effects for microbial-based cross-regional diagnosis of T2D.

<img src="https://github.com/HUST-NingKang-Lab/DeepMicroFinder/blob/main/figure1.png" style="zoom:150%;" />

``
Samples of each cohort were randomly divided into the training subset and the testing subset, and four models were constructed for assessment: (1) Independent disease neural network (DNN) model: ab initio training the DNN model on the training subset and testing on the testing subset of SGMP cohort, respectively. (2) Regional DNN model: ab initio training the DNN model using the training subset of GGMP cohort and testing it on the testing subset of SGMP cohort. (3) Regional+ DNN model: ab initio training the DNN model using the training subset of GGMP cohort as well as the training subset of SGMP cohort, then testing it on the testing subset of SGMP cohort. (4) Transfer DNN model: ab initio training the DNN model using the training subset of GGMP cohort, followed by applying transfer learning to a certain proportion (from 20% to 80%) of samples from SGMP cohort to generate the transfer DNN model, and then testing the transfer DNN model on the testing subset of SGMP cohort. The three boxes on the right represent the evaluation and applications of DeepMciroFinder, including cross-regional diagnosis of T2D and biomarker discovery.
``

## Get and use
To learn how to install the model and how to use it, click [here](https://github.com/HUST-NingKang-Lab/EXPERT)

## Example
Here we choose Guangdong(Region A) and Shandong(Region B) as the source region and target region, and we obtained genus-level species abundance tables for these two regions, 
individuals with T2D are the positive samples and controls are the negitive samples.

#### Prepared files
Genus-level species abundance tables([reference format](https://github.com/HUST-NingKang-Lab/EXPERT)): region1_train.tsv  region1_test.tsv  region2_train.tsv  region2_test.tsv       
Biome profiles([reference format](https://github.com/HUST-NingKang-Lab/EXPERT)): biome.tsv      
Mapper profiles([reference format](https://github.com/HUST-NingKang-Lab/EXPERT)): region1_train_mapper.tsv  region1_test_mapper.tsv  region2_train_mapper.tsv region2_test_mapper.tsv 

#### Ontology construct
- Construct a biome ontology representing stages of T2D. You'll see constructed ontology like a tree in the printed message.
([Reference](https://github.com/HUST-NingKang-Lab/EXPERT))
#### Source mapping
- Map microbial community samples to the biome ontology to obtain hierarchical labels. You'll see counts of the samples on each biome ontology layer in the printed message.
([Reference](https://github.com/HUST-NingKang-Lab/EXPERT))
#### Data convert
- Convert input abundance data to model-acceptable hdf file. The EXPERT model only accepts standardized abundance data. Here we standardize the abundance data using convert mode.
([Reference](https://github.com/HUST-NingKang-Lab/EXPERT))

#### Ab initio training 
- Train the disease neural network model from scratch. Here we will use ontology.pkl and hdf files.
```
expert train -i region1_trainCM.h5 -l region1_train_labels.h5 -t ontology.pkl -o region1_DNN
```
#### Transfer learning
- Transfer the knowledge of region B to the DNN model of region A for better performance in disease diagnosis on region B. You'll see running log and training process in the printed message.
```
expert transfer -i region2_trainCM.h5 -l region2_train_labels.h5 -t ontology.pkl -m  region1_DNN -o Transfer_DNN
```
#### Search
- Search the test set of region B against the transferred DNN model.
```
expert search -i region2_testCM.h5 -m Transfer_DNN -o Search_Transfer_DNN
```
#### Evaluation
- Evaluate the performance of the Transferred DNN model. You'll obtain a performance report.
```
expert evaluate -i Search_Transfer_DNN -l region2_test_labels.h5 -o Evaluation
```

## Maintainers
| Name | Email | Organization |
| ---- | ----- | ------------ |
|Nan Wang|[wangnan123@hust.edu.cn](mailto:wangnan123@hust.edu.cn)|PhD Student, School of Life Science and Technology, Huazhong University of Science & Technology|
|Kang Ning  | [ningkang@hust.edu.cn](mailto:ningkang@hust.edu.cn)       | Professor, School of Life Science and Technology, Huazhong University of Science & Technology |
