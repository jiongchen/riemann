#surf-param-deform#

trying to get familiar with libigl, and may implement some classical algorithms for surface parameterization and deformation techniques gradually.

## Log ##

[2015-4-15] Implement LSCM/DCP, see [[Mullen08]](https://hal.inria.fr/inria-00334477/document).

[2014-4-18] Implement calculating Green Coordinates in 2D for points inside cage, and slight modification of GC for exterior and boundary points will be included in future, as well as Green Coordinates in 3D. 

[2015-4-20] Finish Green Coordinates in 3D and up till now [[Lipman08]](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.150.2116&rep=rep1&type=pdf) is partially reproducted. The computation for GC  of exterior points of cage may be  included in future, which allows for controlling deformation using a partial cage.

[2015-5-22] Partially implement [[Funck06]](https://isgwww.cs.uni-magdeburg.de/visual/files/publications/Archive/Funck_2006_SIGGRAPH.pdf). More modeling metaphors may be included in future.

[2015-5-27] Implement twisting and bending effects based on vector field.

[2015-6-13] Implement frame-driven deformation proposed by [[Panozzo14]](http://igl.ethz.ch/projects/frame-fields/frame-fields.pdf). The other part about frame-field aligned quadrangulation may included later, which consistis of processes to paramterize and remesh the deformed surface for obtaining a uniform quadrilateral mesh. 

[2015-6-25] Implement deformation transfer for triangle meshes [[Sumner04]](http://people.csail.mit.edu/sumner/research/deftransfer/Sumner2004DTF.pdf), which uses gradient based deformation method.

[2015-7-1] Implement another method for solving triangle correspondence guided by harmonic fields in deformation transfer, which is proposed by [[Zayer05]](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.226.1150&rep=rep1&type=pdf).

[2015-8-27] Implement construction of standing wave on surface proposed by [[Zhang10]](http://www.cad.zju.edu.cn/home/hj/10/Huang10WaveQuad.pdf).

[2015-12-1] Implement a 3D bounded distortion solver according to [[Kovalsky15]](http://www.wisdom.weizmann.ac.il/~ylipman/2015_LargeScaleBD.pdf). 