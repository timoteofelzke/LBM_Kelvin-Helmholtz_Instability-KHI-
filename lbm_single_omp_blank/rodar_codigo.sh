#!/bin/bash

# Limpa a compilação anterior
make clean

# Compila o código utilizando múltiplos núcleos
make -j

# Executa o programa gerado
./grad-lbm

# executa o paraview
paraview
