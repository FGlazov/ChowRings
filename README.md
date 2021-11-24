# ChowRings

This is a small Julia Library that computes the Chow Ring and augmented Chow Rings of matroids. This library is based on the OSCAR.jl project. In particular, many of its inputs expect matroids as BigObject from the Polymake.jl subsystem, and it represents the chow rings with the type MPolyQuoRing from Oscar.jl. This library can also compute the decompositions in Theorem 1.1 and Theorem 1.2 in the paper 

@misc{braden2020semismall,
      title={A semi-small decomposition of the Chow ring of a matroid}, 
      author={Tom Braden and June Huh and Jacob P. Matherne and Nicholas Proudfoot and Botong Wang},
      year={2020},
      eprint={2002.03341},
      archivePrefix={arXiv},
      primaryClass={math.AG}
}

- which can be found here https://arxiv.org/abs/2002.03341. The notation used inside the code is inspired by the notation used in that paper.

# Installation

You must have Julia installed, along with Oscar.jl. Until there is a stable release, currently the best way to install the library is to clone this repository, and then activate it locally on your computer. I.e.

```
> git clone https://github.com/FGlazov/ChowRings.git
> cd ChowRings
> julia
( Press ] once in Julia to go into package manager)
> pkg> activate .
( Press backspace to exit package manager and to return to the julia shell)
```
