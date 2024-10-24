FROM ghcr.io/becksteinlab/streaming-md-docker:main

# COPY --from=tianon/gosu /gosu /usr/local/bin/

# Create user and directory with the correct permissions
# ARG NB_USER=conda
ARG NB_USER=conda
ENV USER ${NB_USER}
ARG NB_UID=1000
ENV HOME /home/${NB_USER}
RUN useradd --shell /bin/bash -u ${NB_UID} -G lucky -o -c "" -m conda

# ENV USER ${NB_USER}
# ENV NB_UID ${NB_UID}
# ENV HOME /home/${NB_USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN chmod +x ${HOME}/entrypoint

# # Ensure correct ownership of the home directory
# RUN chown -R ${NB_UID}:${NB_UID} ${HOME}

# USER ${NB_USER}
# COPY . ${HOME}


# USER root
# RUN chown -R ${NB_UID} ${HOME}

# Install JupyterLab as root
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate env && \
    python3 -m pip install --no-cache-dir notebook jupyterlab

ENTRYPOINT [ "/opt/conda/bin/tini", "--", "/home/conda/entrypoint" ]
CMD [ "/bin/bash" ]

