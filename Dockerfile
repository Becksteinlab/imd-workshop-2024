FROM ghcr.io/becksteinlab/streaming-md-docker:main

# COPY --from=tianon/gosu /gosu /usr/local/bin/

# Create user and directory with the correct permissions
# ARG NB_USER=conda
ARG NB_USER=conda
ENV USER ${NB_USER}
ARG NB_UID=1000
ENV HOME /home/${NB_USER}
RUN useradd --shell /bin/bash -u ${NB_UID} -G lucky -o -c "" -m conda

COPY . .
RUN chown -R ${NB_UID} ${HOME}

# Modify the entrypoint script
COPY entrypoint /opt/docker/bin/entrypoint
RUN chmod +x /opt/docker/bin/entrypoint

# RUN source /opt/conda/etc/profile.d/conda.sh && \
#     conda create -n workshop python notebook jupyterlab jupyterhub imdclient MDAnalysis
RUN source /opt/conda/etc/profile.d/conda.sh && \
        conda env create -n workshop --file=env.yaml

# # Install JupyterLab as root
# RUN source /opt/conda/etc/profile.d/conda.sh && \
#     conda activate env && \
#     python3 -m pip install --no-cache-dir notebook jupyterlab jupyterhub imdclient MDAnalysis



CMD [ "/bin/bash" ]