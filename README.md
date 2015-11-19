#LinearFW

This is the code to reproduce all the experiments in our 
[NIPS 2015 paper](http://arxiv.org/pdf/1511.05932v1):
```
On the Global Linear Convergence of Frank-Wolfe Optimization Variants
Simon Lacoste-Julien and Martin Jaggi
NIPS 2015
```
which covers the global linear convergence rate of Frank-Wolfe
optimization variants for problems described as in Eq. (1) in the paper.
It contains the implementation of Frank-Wolfe, 
away-steps Frank-Wolfe and pairwise Frank-Wolfe on two applications:

1. l1-constrained least-square regression (lasso);
2. a QP on the flow polytope coming from a video co-localization application.

The code runs in Matlab (was tested in Matlab 2014 on Linux, Windows and Mac).
But for the first two folders below, it also runs in Octave easily by
removing the line initializing the random seed.

There are three folders:
* `FW_lasso` contains the Lasso experiment to produce the top figure in 
   Figure 2. Launch `run_FW.m` in the folder to produce the plot (takes a few seconds).
* `triangle_FW_experiment` contains the empirical tigthness of the linear rate
   constant experiment, from Appendix E (Figure 5). Launch `run_triangles.m`
   in the folder to produce the plots (takes about 30 seconds).
* `FW_video_colocalization` containts the video co-localization QP experiment
   to produce the bottom figure in Figure 2. Launch `run_FW.m` in the folder
   to produce the plot (takes less than 1 minute). If you get the 
   `Undefined function 'solver_video_mex'` error, you need to mex the following
   file to get the correct LMO; go to the `solvers` subfolder, and then
   run `mex solver_video_mex.cpp` in it.

##Credits

The video co-localization code was written by 
[Armand Joulin](http://ai.stanford.edu/~ajoulin/) and 
[Kevin Tang](http://ai.stanford.edu/~kdtang/).
We obtained it [here](http://ai.stanford.edu/~ajoulin/code/FW.zip)
and modified the Frank-Wolfe code by adding a hashing function
to make the active set maintenance more efficient, as well as
added the pairwise FW variant. Their
video co-localization approach is described in the paper:
```
Efficient Image and Video Co-localization with Frank-Wolfe Algorithm
Armand Joulin, Kevin Tang and Li Fei-Fei
ECCV 2014
```

##Extending the code

You can easily re-use the code in the `FW_video_colocalization` folder
to adapt it to other QPs with different domains. For this, you mainly need
to modify two things:

1. You need to implement your own Linear Minimization Oracle (LMO) on your 
   domain. You pass it as the `fun_optim` argument to the FW functions
   (`FW` for standard FW; `AFW` for away-steps FW and `PFW` for pairwise
    FW). This function takes a vector as argument and returns an atom
    that minimizes the dot product with this vector over the domain.
2. If your atoms are something different than just 0-1 vectors, then
   you also need to modify the `hashing` internal function for `AFW`
   and `PFW` so that you can properly encode the atoms in your domain
   in a unique string (or number) (this is used for the efficient 
   maintenance of the active set). The current hashing function only
   supports 0-1 vectors.

##Disclaimer

This code is not meant to be the most efficient. For example, an implementation
handling the sparsity of the atoms could be faster for the video
co-localization application. Our goal was to simply demonstrate how the FW
variants work on practical applications.

##Citation

Please use the following BibTeX entry to cite this software in your work:

    @inproceedings{LacosteJulien2015linearFW,
      author    = {Simon Lacoste-Julien and Martin Jaggi},
      title     = {On the Global Linear Convergence of {F}rank-{W}olfe Optimization Variants},
      booktitle = {Advances in Neural Information Processing Systems (NIPS)},
      year      = {2015},
    }

And if you use the video co-localization LMO, you also need to cite:

    @inproceedings{JouTangFeiECCV14,
      title     = {Efficient Image and Video Co-localization with {F}rank-{W}olfe Algorithm},
      author    = {Armand Joulin and Kevin Tang and Li Fei-Fei},
      booktitle = {European Conference on Computer Vision (ECCV)},
      year      = {2014},
    }


##Authors

* [Simon Lacoste-Julien](http://www.di.ens.fr/~slacoste/)
* [Martin Jaggi](http://people.inf.ethz.ch/jaggim/)
