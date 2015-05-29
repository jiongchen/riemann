#surf-param-deform#

trying to get familiar with libigl, and may implement some classical algorithms for surface parameterization and deformation techniques gradually.

## Log ##

[2015-4-15] Implement LSCM/DCP, see [[Mullen08]](https://hal.inria.fr/inria-00334477/document).

[2014-4-18] Implement calculating Green Coordinates in 2D for points inside cage, and slight modification of GC for exterior and boundary points will be included in future, as well as Green Coordinates in 3D. 

[2015-4-20] Finish Green Coordinates in 3D and up till now [[Lipman08]](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.150.2116&rep=rep1&type=pdf) is partially reproducted. The computation for GC  of exterior points of cage may be  included in future, which allows for controlling deformation using a partial cage.

[2015-5-22] Partially implement [[Funck06]](https://isgwww.cs.uni-magdeburg.de/visual/files/publications/Archive/Funck_2006_SIGGRAPH.pdf). More modeling metaphors may be included in future.

[2015-5-27] Implement twisting and bending effects based on vector field.