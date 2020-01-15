# Exercises and lectures on NEURON and NetPyNE (Salvador and Bill) for Lascon 2020 

## Software sources

# Contents of this repo:

<!-- github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet -->

## [Talk pdfs](https://drive.google.com/drive/folders/1U8UJvg4OIHsrzWqTatu4ByoV7uwNrY_n)

- Day1Talk-NeuralCircuits.pdf `General introduction to computational neuroscience + examples of neural circuits`

- Day2Talk-HopBrody.pdf `Building a Hopfield-Brody network in NEURON` (this talk contains additional slides that we did not get to in the time period)

- Day3Talk-netpyne.pdf `Introduction to netpyne`

- Day4Talk-MPI.pdf `HPC, MPI and ParallelContext`

## Code

**test0.py** `simple test of mpi` **mpiexec -n 4 --oversubscribe nrniv -python -mpi test0.py**

**HopBrod.py**

```
presented on Day 2,4 -- Hopfield-Brody network using INVL ArtificialCell or HH cell (with useHHCell=True)

uses IClamp for HH cell rate randomization

mpirun() to run under MPI (using HHcell only)
```

**invlfire.mod**

`NMODL source file for INVL ArtificialCell used by HopBrod.py`

**HopBrodNetPyNE.py**

```
shortened and simplified version of netpyne/doc/source/code/tut7.py
uses HH with NetStims for HH cell rate randomization
```

## Additional Exercises

- `Accessing and making github repos`

    get a github account from http://github.com

    set up a repo for your LASCON project

    access github repos for [NetPyNE](https://github.com/Neurosim-lab/netpyne) and [NEURON](https://github.com/neuronsimulator/nrn)
    

