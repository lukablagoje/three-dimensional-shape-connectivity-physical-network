# Three-dimensional shape and connectivity of physical networks
Data describing the three-dimensional structure of physical networks is increasingly available, leading to a surge of interest in network science to explore the relationship between the shape and connectivity of physical networks. We contribute to this effort by standardizing and analyzing 15 data sets from different domains. Treating junction points as nodes and connections between them as links, we divide the networks into three categories: lattice-like networks, trees, and linked trees. We find that the degree distribution of physical networks is bounded, with most nodes having degrees one or three. Characterizing the physical properties of links, we show that links have an elongated shape and tend to follow a nearly straight trajectory, while a small fraction of links follow a winding path. These typical node and link properties must be reflected by physical network models. We also measure how confined a link is in space by comparing its trajectory to a randomized null-model, showing that links that are central in the abstract network tend to be physically confined by their neighbors. The fact that the shape and connectivity of the physical networks are intertwined highlights that their three-dimensional layout must be taken into account to understand the evolution and function of physical networks. Below is a figure from the paper explaining how link confinement measure is computed:

![image](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/assets/52599010/616c67f1-e5c9-4b89-b926-db7740d10f80)


arXiv link: [https://arxiv.org/abs/2211.13265](https://arxiv.org/abs/2403.19333)

# Technical project overview
This project is divided into four different parts:

[**1. data**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/1.%20data) - I perform data processing of physical network skeletons, where I perform skeleton healing and merging of parallel segments. There are two folders, one for the fruit fly datasets, and one for other skeletons.

[**2. basic_measures**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/2.%20basic_measures) - I consider physical aspects of my data, so I analyze their fractal dimension and space-filling (physical density).

[**3. abstract_network**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/3.%20abstract_network) - With a graph structure to work with, I perform standard network analysis (degree distribution and motifs).

[**4. link_confinement_measure**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/4.%20link_confinement_measure) - I apply collision detection algorithms to compute a link-level measure of physical confinement.

# Data
All processed datasets can be accessed here: [**1. data/ 3. final data**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-network/tree/main/1.%20data/3.%20final_data)

If you are interested in the raw data (before processing), please take a look at the supplemental information of the article.
