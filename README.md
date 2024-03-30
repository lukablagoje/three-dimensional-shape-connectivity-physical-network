# Three-dimensional shape and connectivity of physical networks
Data describing the three-dimensional structure of physical networks is increasingly available, leading to a surge of interest in network science to explore the relationship between the shape and connectivity of physical networks. We contribute to this effort by standardizing and analyzing 15 data sets from different domains. Treating junction points as nodes and connections between them as links, we divide the networks into three categories: lattice-like networks, trees, and linked trees. We find that the degree distribution of physical networks is bounded, with most nodes having degree one or three. Characterizing the physical properties of links, we show that links have an elongated shape and tend to follow a nearly straight trajectory, while a small fraction of links follow a winding path. These typical node and link properties must be reflected by physical network models. We also measure how confined a link is in space by comparing its trajectory to a randomized null-model, showing that links that are central in the abstract network tend to be physically confined by their neighbors. The fact that the shape and connectivity of the physical networks are intertwined highlights that their three-dimensional layout must be taken into account to understand the evolution and function of physical networks.

arXiv link: [https://arxiv.org/abs/2211.13265](https://arxiv.org/abs/2403.19333))

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
