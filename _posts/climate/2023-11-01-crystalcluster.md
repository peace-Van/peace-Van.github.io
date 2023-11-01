---
layout: post
title: A novel clustering method - Crystal clustering
author: Van
category: climate
---

Code is on [GitHub](github.com/peace-Van/crystal-cluster)   
*Note: This clustering method, proposed by myself, is used in my climate classification.*   
   
## Phenomenon
Imagine how a grain of sugar dissolves in water. Before it dissolves, the sucrose molecules form hydrogen bond (also van der Waals force but much weaker, only mention hydrogen bond in the following discussion) in between, so they condense as a macroscopic solid bulk.   
   
![p1](/assets/crystalcluster/sucrose01.png)   
The sucrose molecule (credit: [UCLA chem](https://www.chem.ucla.edu/~harding/IGOC/S/sucrose.html))   
   
![p2](/assets/crystalcluster/sucrose-crystal.png)   
The crystal structure of sucrose (credit: [Ji-Hun An et al.](https://www.mdpi.com/2073-4352/7/10/284)). The dashed lines show the hydrogen bonds.

When met with water, some of the hydrogen bonds break for water molecules cut in. Say there's not enough water, and not all the bonds are broken. In this case, at the macroscopic level, the sugar is partially dissolved and there's no ongoing net dissolution in the system, while at the microscopic level, the sucrose molecules in the aquatic phase occasionally collides with the ones on the edge of the solid phase, forming hydrogen bond and precipitates. Meanwhile, on the other side of the solid water molecules drag out the solid-phase-locked sucrose molecules. Microscopic dissolution and precipitation are happening simultaneously at the same speed, so the macroscopic observation is no net dissolution or precipitation.
 
> If you can't imagine that, [this video](https://www.youtube.com/watch?v=BavKFkSzNFE) (credit: [Berean Builders](https://www.youtube.com/@Bereanbuilders/about)) may help.  

As the partial dissolution system stablizes (aka reaches the **thermodynamically stable state**), it naturally forms several clusters of solid blocks. If we view the sucrose molecules as the data observations to cluster, now we get the result with outliers dissolved away (or as singleton clusters).   

## Modelling
But how can we model this process by computer algorithm? We need to dive deeper into the motive of this spontaneous process - what drives the bond to form/break?   

> For the sake of simplicity we will make certain assumptions, e.g. there's no kinetic barrier, no pressure-volume work in the process.   

It's the **Gibbs Free Energy**. The definition is

$$G=H-TS$$

where $H$ is the enthalpy (bond energy), $T$ is the thermodynamics temperature (always positive), $S$ is the entropy. A process is spontaneous in thermodynamics only if it **results in a decrease in the system's Gibbs free energy**.

We keep the system at a pre-set temperature, so $T$ is a constant here. Whenever a bond breaks, the system absorbs energy ($\Delta H$ is positive), and the degrees of freedom increases ($\Delta S$ is positive), and vice versa. Whether this is thermodynamically possible depends on the condition - if $\Delta G = \Delta H - T\Delta S < 0$, it is. It's easy to know that **higher temperature favors bond breaking, leading to more clusters**. 

Next let's define $H$ and $S$ in terms of the data clustering system. We have to define the system itself first.   
   
> This is a blog post, so I prefer to make it as simple as possible, using minimal math representations.
   
 - **Data obervations** to be clustered: a **$n$ by $m$ matrix $X$**, rows are the observations and columns are feature dimensions. The observations can have different weights, $w_i$ for the $i$-th observation. This can be plotted as a graph $G$ with $n$ nodes. At this point let's assume the graph is fully connected, and the edge weight is the distance in the feature space between the connected nodes.   

 > Distance measures are not limited. Here we assume Euclidean distance.   
   
 - **Clusters: connected components of the graph $G$**. In this way **only the minimum spanning edges (edges on the minimum spanning tree) matter**, because non minimum spanning edges do not affect the connected components of the graph. From a viewpoint of the physical world, if a molecule in the aquatic phase choose to form a bond with a molecule in the solid phase, it would always prioritize the nearest one, resulting in a minimum spanning forest on the graph representation. So the fully connected graph reduces to **a minimum spanning tree of $G$, denoted as $G_{MST}$**.
 - **State**: The system evolves with bonds forming and breaking during the process before the stable state is reached. During the process, get a snapshot of the system, and we will have a subgraph of $G_{MST}$, with the same set of nodes and a subset of edges.
 - **Action**: A bond forming or breaking is called an action. An action that results in the decrease of the Gibbs free energy is called a **feasible action**.   

![p3](/assets/crystalcluster/adjacency.png)   
A state of the data clustering system   

### Enthalpy
For two data observations $X_i, X_j$ with weight $w_i, w_j$ and the distance between them is $d_{i, j}$, how do we define the bond energy between them? Here I would employ the concept in electrostatics:

$$\Delta H_{i,j}^{form} = - \frac{w_i w_j}{d_{i, j}^2}$$

This is analogous to the energy stored in a system of two point charges. The further two points are, the less potential they form a bond in between ($\Delta H$ is less negative). This value **is not** dependent on the state of the system.  

>Under common scenarios we do not add up the weights of data observations in the same connected component when calculating enthalpy. If we do so, it would produce a strong gravitational effect that draws a large number of data observations into one cluster, which, is a more faithful reproduction of the crystallization process (crystal growth after nucleation), but is not desired under the data clustering setting. What we want to simulate here is the system kept at the nucleation stage.

Enthalpy is summable with respect to bonds. For a certain state of the system, add up the $\Delta H^{form}$ of all existing bonds and we get the enthalpy of the state. That is
$$H = -\sum_{i, j \in e} \frac{w_i w_j}{d_{i, j}^2}$$

### Entropy
For a state with $k$ clusters and each cluster $i$ has data observations whose weights sum up to $c_i$, the entropy is calculated as
$$S = \frac{1}{N} \sum_{i=1}^k c_i \log \frac{c_i}{N}$$
where $N = c_1 + c_2 + ... + c_k = \sum_{j=1}^n w_j$. This is from the definition of information entropy.  
   
It's not easy to calculate the entropy change of an action ($\Delta S_{i,j}$) without using some approximation method. We have to calculate both the entropy of the state before the action and the entropy of the state after the action, and then subtract. Unlike enthalpy change, this value **is** dependent on the state of the system.

## Algorithm
We are going to solve the optimization problem on $G_{MST}$. Remove some edges of $G_{MST}$, so that
$$f(G_{MST}, T) = H - TS$$
is minimized. $T$ is a given positive constant (hyperparamter). $H$ and $S$ are defined above.  
   
The code provides four algorithms. Here I would only introduce the `greedy-backtrack` algorithm, which should be used in most cases. This algorithm simulates the precipitation-solubility equilibrium.  

![p3](/assets/crystalcluster/algorithm.PNG)

The description is self-explanatory and I would rather not dive into the details here. One thing to note is that the algorithm does not guarantee to converge nor to find a global optimum, but in most cases it converges fast and the result is good enough to use. It may not be applicable to large datasets because it's $O(n^2)$ in both time and space complexity.

## Test
On the [Iris flower dataset](https://archive.ics.uci.edu/dataset/53/iris), the clustering results are shown as follows. It's clear that number of clusters increases with rise of temperature. 
![pp1](/assets/crystalcluster/exp1.png)
![pp2](/assets/crystalcluster/exp2.png)
![pp3](/assets/crystalcluster/exp3.png)
![pp4](/assets/crystalcluster/exp4.png)
![pp5](/assets/crystalcluster/exp5.png)
![pp6](/assets/crystalcluster/exp6.png)

## Further Information
This method actually belongs to a group of clustering algorithms of *optimising different criteria – divisive strategy over MST* named by [Marek Gagolewski et al.](https://arxiv.org/pdf/2303.05679.pdf). In the article they extensively studied the clustering algorithms based on MST and compared with traditional algorithms like $k$-means, BIRCH, linkage, PAM, etc. They concluded that *MSTs are suitable for representing dense clusters of arbitrary shapes, and are relatively robust to outliers in the outer parts of the feature space*.   

This method is out of their scope of study because it does not have a fixed number of clusters $k$ as input. Rather, we use a *relaxation strength* parameter $T$. This can have its own pros and cons. $k$ is more intuitive to humans, while the good things are 
- the scale is continuous, unlike $k$ has to be an integer
- easier to deal with outliers, no need to investigate the number of outliers for a good guess of $k$
   
It may need more trial and error to get the best $T$ for the desired outcome.   

> The most similar method discussed by [Marek Gagolewski et al.](https://arxiv.org/pdf/2303.05679.pdf) is [ITM by Andreas C. M¨uller et al.](https://www.nowozin.net/sebastian/papers/mueller2012itclustering.pdf). Their method requires the number of clusters $k$ as input and targets at maximizing the mutual information between input data and output labels. It approximates on entropy calculation and thus may run much faster.

This method may also have more advanced applications, such as a multi-stage clustering process which simulates the fine fractionation operation in chemical production: at a higher temperature, the first aggregated cluster is removed from the system, and then gradually cool down to remove the subsequent clusters step by step at different levels of temperature. The gravitational effect of agglomerated weight is yet to be further explored and may have its potential applications.
