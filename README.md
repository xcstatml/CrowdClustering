Main program:
BDP_cluster_reliable.m: The Opt-KG policy for crowdsourced clustering with reliable workers (Algorithm 1 in a batch mode)
BDP_cluster_unreliable.m: Bayesian decision process with unreliable workers based on the Opt-KG policy (Algorithm 2 in a batch mode)

Test program:
test_BDP_cluster_reliable.m
test_BDP_cluster_unreliable.m

Misc:
We use grPartition (J. Hespanha. An efficient MATLAB Algorithm for Graph Partitioning.
http://www.ece.ucsb.edu/~hespanha/techrep.html, 2004.) when implementing our algorithms.

Remark:
In main programs, the matrix of similarity parameters is required for generating similarity labels in simulation studies. In real applications, similarity parameters are unknown and the similarity labels of the selected pairs of items are required. 


Please cite the following paper if the package was used:

Xiaozhou Wang, Xi Chen, Qihang Lin and Weidong Liu. Bayesian Decision Process for Budget-efficient Crowdsourced Clustering. In Proceedings of the Twenty-Ninth International Joint Conference on
 Artificial Intelligence, pages 2044-2050, 2020.

@inproceedings{ijcai2020-283,
  
  title     = {Bayesian Decision Process for Budget-efficient Crowdsourced Clustering},
  
  author    = {Wang, Xiaozhou and Chen, Xi and Lin, Qihang and Liu, Weidong},
  
  booktitle = {Proceedings of the Twenty-Ninth International Joint Conference on
 Artificial Intelligence, {IJCAI-20}},
  
  publisher = {International Joint Conferences on Artificial Intelligence Organization},
  
  editor    = {Christian Bessiere},	
  
  pages     = {2044--2050},
  
  year      = {2020},
  
  month     = {7},
  
  note      = {Main track}
  
  doi       = {10.24963/ijcai.2020/283},
  
  url       = {https://doi.org/10.24963/ijcai.2020/283},

}