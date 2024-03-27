# Research project overview
This project is a part of the published research work "The Impact of Physicality on Network Structure", done in the emerging field of Physical Networks, which aims to understand the properties of three-dimensional networked systems (such as biological neural networks).

In this paper, a linear physical network model is introduced, along with the collision-avoiding representation graph, called a "metagraph" representation. The concept of "generalized meta-graph" is introduced as well, which can be applied to empirical data, by encoding information on Euclidean distances between neighboring edges.

In my contribution to this project, I applied this new representation to efficiently solve a computationally challenging task: finding out how many unique neuron-to-neuron physical collisions are obtained if the thickness of the neurons is increased for 20 different additive factors. The results of this analysis have shown that biological neural networks are composed of highly confined objects, which have many neighbors in their local physical neighborhood, which cannot be said for mitochondrial, vascular, and plant root networks.
I further developed the representation to analyze not the entire dataset, but each individual neuron and its physical confinement.

If you want to read more about this research, please check out the links below:

Publication link: [https://www.nature.com/articles/s41567-023-02267-1](https://www.nature.com/articles/s41567-023-02267-1)

arXiv link: [https://arxiv.org/abs/2211.13265](https://arxiv.org/abs/2211.13265)

# Technical project overview
This project is divided into four different parts:

[**1. data**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/1.%20data) - I perform data processing of physical network skeletons, where I perform skeleton healing and remove redundancies in the physical structure. There is a more detailed set of notebooks where I query and download neurons and connect them via their synapses to form a single network.

[**2. basic_measures**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/2.%20basic_measures) - I consider physical aspects of my data, so I analyze their fractal dimension and space-filling properties (physical density).

[**3. abstract_network**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/3.%20abstract_network)- With a graph structure to work with, I perform standard network analysis (degree distribution and motifs).

[**4. link_confinement_measure**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/4.%20link_confinement_measure) - Finally, I apply the labeled point cloud algorithm to find the number of physically close neurons (objects) as their radial thickness is increased, both for individual neurons and all neurons in the dataset at once.

# Data
[**neuron_regions_information**](https://github.com/lukablagoje/the-impact-of-physicality-on-network-structure/tree/main/neuron_regions_information) - Folder for storing information about the neurons in specific regions, which are in CSV format.

[**neuron_regions_points**](https://github.com/lukablagoje/the-impact-of-physicality-on-network-structure/tree/main/neuron_regions_points) - Folder for storing point clouds of neuron data, which are in CSV format.

[**querying results**](https://github.com/lukablagoje/the-impact-of-physicality-on-network-structure/tree/main/querying_results) - Folder for storing labeled kd-tree results, which are in pickle format (the number in the filename represents querying radius)

To access the original dataset, you will need to use the Python library provided by the Janelia project, which you can install from https://pypi.org/project/neuprint-python/

If you want to understand this library more, you can read through their documentation https://connectome-neuprint.github.io/neuprint-python/docs/index.html
