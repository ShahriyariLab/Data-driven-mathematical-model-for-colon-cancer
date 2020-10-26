# Data-driven-mathematical-model-for-colon-cancer

Every colon cancer has its own unique characteristics, and therefore may respond differently to identical treatments. Here, we develop a data driven mathematical model for the interaction network of key components of immune microenvironment in colon cancer. We estimate the relative abundance of each immune cell from gene expression profiles of tumors, and group patients based on their immune patterns. Then we compare the tumor sensitivity and progression in each of these groups of patients, and observe differences in the patterns of tumor growth between the groups.

This repository contains the following scripts:

qspmodel.py: python classes and methods for analysis of 
          Mathematical model of immune response in cancer

QSPClusterAnalysis.py: python script for analysis of  parameter sensitivity and dynamics of 
                    Mathematical model of immune response in colon cancer

ParseQSPData.m: MatLab script to parse the results obtained by running 
			  QSPClusterAnalysis.py for ploting by PlotQSPData.m
		
PlotQSPData.m: MatLab script for plotting the results obtained by running 
			 QSPClusterAnalysis.py and ParseQSPData.m
			 
If using this or related code please cite

```
 Kirshtein, A.; Akbarinejad, S.; Hao, W.; Le, T.; Aronow, R.A.; Shahriyari, L. 
  Data driven mathematical model of colon cancer progression. (Manuscript submitted for publication).
```
