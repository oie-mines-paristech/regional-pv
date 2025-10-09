# regional-pv

Code to convert weather into regional photovoltaic estimates.

## Background

The code was originally developed as part of the [PhD research](https://kobra.uni-kassel.de/items/cbf2d921-2685-4143-b915-ed3033da902f)
of Yves-Marie Saint Drenan at the University of Kassel, Germany.

Following his move to Mines Paris PSL, Yves-Marie continued to enhance the
codebase, applying it in several European research initiatives, including the
[European Climate Energy Mixes (ECEM)](https://climate.copernicus.eu/european-climate-energy-mixes) in 2018
and [Clim2Power](https://jpi-climate.eu/project/clim2power/) in 2020.

Between 2022-2025, Rodrigo Amaro e Silva (ULisbon & Mines Paris – PSL) joined
the development team, contributing to the integration of the code into the [Pan-European Climate Database (PECD)](https://climate.copernicus.eu/powering-europe-through-climate-uncertainty),
funded by Copernicus Climate Change Services (C3S).

## Workflow for users

After cloning the repository, define it as your current directory.

For best experience create a new conda environment (e.g. `regional-pv`) with Python 3.12:

```
conda create -n regional-pv -c conda-forge python=3.12
conda activate regional-pv
conda env create -f environment.yml
pip install .
```

Then, in a Python script, import the package and run `spv_workflow`:

```python
import regional_pv

out = regional_pv.spv_workflow(...)
```

## Publications

TBD.

## License

```
Copyright 2018-2025, ARMINES.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```
