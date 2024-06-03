# Three-dimensional shape and connectivity of physical networks
Data describing the three-dimensional structure of physical networks is increasingly available, leading to a surge of interest in network science to explore the relationship between the shape and connectivity of physical networks. We contribute to this effort by standardizing and analyzing 15 data sets from different domains. Treating junction points as nodes and connections between them as links, we divide the networks into three categories: lattice-like networks, trees, and linked trees. We find that the degree distribution of physical networks is bounded, with most nodes having degrees one or three. Characterizing the physical properties of links, we show that links have an elongated shape and tend to follow a nearly straight trajectory, while a small fraction of links follow a winding path. These typical node and link properties must be reflected by physical network models. We also measure how confined a link is in space by comparing its trajectory to a randomized null-model, showing that links that are central in the abstract network tend to be physically confined by their neighbors. The fact that the shape and connectivity of the physical networks are intertwined highlights that their three-dimensional layout must be taken into account to understand the evolution and function of physical networks. Below is a figure from the paper explaining how link confinement measure is computed:

![image](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/assets/52599010/616c67f1-e5c9-4b89-b926-db7740d10f80)


arXiv link: [https://arxiv.org/abs/2211.13265](https://arxiv.org/abs/2403.19333)

# Technical project overview
This project is divided into four different parts:

[**1. data**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/1.%20data) - I perform data processing of physical network skeletons, where I perform skeleton healing and merging of parallel segments. There are two folders, one for the fruit fly datasets, and one for other skeletons.

[**2. basic_measures**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/2.%20basic_measures) - I consider physical aspects of my data, so I analyze their fractal dimension and space-filling (physical density).

[**3. abstract_network**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/3.%20abstract_network)- With a graph structure to work with, I perform standard network analysis (degree distribution and motifs).

[**4. link_confinement_measure**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-networks/tree/main/4.%20link_confinement_measure) - I apply collision detection algorithms to compute a link-level measure of physical confinement.

# Data
All processed datasets can be accessed here: [**1. data/ 3. final data**](https://github.com/lukablagoje/three-dimensional-shape-connectivity-physical-network/tree/main/1.%20data/3.%20final_data)

\begin{table}[h]
\begin{center}
\begin{tabular}{| c | c | c | c |c|c|c| c| c|} 
 \hline
Data set  & \textbf{$N_{\T{link}}$} & \textbf{$N_{\T{node}}$} & \textbf{$N_{\T{seg}}$} & \textbf{$l_{\T{seg}}$} & \textbf{$\rho_{\T{seg}}$} & \textbf{$a_{\T{seg}}$} & Category & Source\\ [0.5ex] 
\hline
h\_neuron     & 631   & 632   & 38340 & $1.16 \pm 0.07$ & $0.11 \pm 0.00$ & $0.1 \pm 0.01$ & tree & \cite{koch2016big}\\ 
\hline
r\_neuron       & 185   & 186   & 4536 &  $2.21 \pm 1.98$ &  $0.46 \pm 0.19$ & $0.22 \pm 0.27$ & tree & \cite{bartho2007cortical}\\
\hline
m\_neuron    & 154   & 155   & 16814 & $0.36 \pm 0.28$ & $0.23 \pm 0.15$ & $0.62 \pm 0.3$ & tree & \cite{coskren2015functional}\\
\hline
z\_neuron & 307   & 308   & 2867 &  $1.51 \pm 1.24$ &  $0.34 \pm 0.12$ & $0.22 \pm 0.19$ & tree & \cite{torvund2017cone}\\
\hline
anthill          & 15240 & 15241 & 29387 & $16.4 \pm 8.38$ & $9.61 \pm 3.17$ & $0.59 \pm 0.3$ & tree & \cite{anthill}\\
\hline
root\_1           & 975   & 976   & 5621  &  $25.66 \pm 0.61$  & $4.69 \pm 5.11$ & $0.18 \pm 0.2$ & tree & \cite{ohashi2019reconstruction}\\
\hline
root\_2           & 410   & 411   & 2132   &  $38.49 \pm 12.56$ & $7.2 \pm 7.47$ & $0.18 \pm 0.2$ & tree & \cite{ohashi2019reconstruction} \\
\hline
fruit\_fly\_1      & 100388 & 97588 & 535611 & $48.5 \pm 49.2$ & $19.8 \pm 17.38$ & $0.41 \pm 0.19$ & linked tree & \cite{clements2020neuprint}\\
\hline
fruit\_fly\_2      & 32121 & 31408 & 181068 & $48.99 \pm 46.06$ & $20.0 \pm 15.86$ & $0.41 \pm 0.2$ & linked tree & \cite{clements2020neuprint}\\
\hline
fruit\_fly\_3      & 49599 & 49233 & 121318 & $39.6 \pm 41.97$ & $18.38 \pm 20.2$ & $0.49 \pm 0.27$ & linked tree & \cite{clements2020neuprint}\\
\hline
fruit\_fly\_4      & 34987 & 32749 & 138488 & $41.57 \pm 41.57$ & $16.97 \pm 17.85$ & $0.41 \pm 0.19$ & linked tree & \cite{clements2020neuprint}\\
\hline
vascular\_1       & 2359 & 1558 & 17935 &  $4.69 \pm 1.62$ & $3.0 \pm 1.22$ & $0.62 \pm 0.33$ & lattice & \cite{gagnon2015quantifying} \\
\hline
vascular\_2       & 1300 & 862 & 16078 &  $3.91 \pm 1.11$ & $3.17 \pm 1.41$ & $0.83 \pm 0.45$ & lattice & \cite{gagnon2015quantifying}\\
\hline
vascular\_3       & 1181 & 789 & 12487 &  $5.1 \pm 1.95$ & $2.96 \pm 0.97$ & $0.61 \pm 0.31$ & lattice & \cite{gagnon2015quantifying}\\
\hline
mitochon    & 73 & 59 & 847 & $0.09 \pm 0.06$ & $0.11 \pm 0.02$ & $1.17 \pm 0.67$ & lattice & \cite{viana2020mitochondrial}\\
\hline
\end{tabular}
\caption{\textbf{Data sets summary:} For each data set, we provide the total number of physical nodes $N_\T{node}$ (i.e., junction and terminal points in the skeletonized representation), physical links $N_\T{link}$ and skeleton segments $N_\T{seg}$. We also provide the segment statistics - segment length $l_{\T{seg}}$, segment radius $\rho_{\T{seg}}$, and segment aspect ratio $a_{\T{seg}}$, along with the abstract network categories and sources. For quantities with $\pm$, we used the median and interquartile range (difference between the 75th and 25th percentile) to quantify their variation.}
\label{table_datasets_SI}
\end{center}
\end{table}
