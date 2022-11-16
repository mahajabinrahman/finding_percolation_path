# finding_percolation_path
When coffee is strained through a series of filters, it has successfully percolated. Like coffee, other things also percolate (water, diseases, blood, polymers, networks of neurons, etc). Read this great Scientific American article for a fun intro: https://www.scientificamerican.com/article/the-mathematics-of-how-connections-become-global/. 

Physicists care about percolation because it typically marks a phase transition, where the nature of microscopic interactions fundamentally change to allow long range correlations. 

"main.py" uses a .txt file of indices as its input, converts those indices to 3D coordinates. From the coordinates, the burning algorithm is used to identify separate networks. The code will go through each of the networks and identify whether any of them has spanned from wall to wall. 

[percolation_path.pdf](https://github.com/mahajabinrahman/finding_percolation_path/files/10025680/percolation_path.pdf)
