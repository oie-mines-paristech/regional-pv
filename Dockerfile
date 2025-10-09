FROM continuumio/miniconda3

WORKDIR /src/regional-pv

COPY environment.yml /src/regional-pv/

RUN conda install -c conda-forge gcc python=3.12 \
    && conda env update -n base -f environment.yml

COPY . /src/regional-pv

RUN pip install --no-deps -e .
