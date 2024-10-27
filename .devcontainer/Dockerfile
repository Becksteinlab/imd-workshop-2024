FROM ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

# Modify the entrypoint script so that it activates the workshop environment
COPY create_conda_user /opt/docker/bin/create_conda_user
RUN chmod +x /opt/docker/bin/create_conda_user
RUN /opt/docker/bin/create_conda_user

RUN source /opt/conda/etc/profile.d/conda.sh && \
        conda env create -n workshop --file=env.yaml

# Codespaces does not use these!
# ENTRYPOINT [ "/opt/conda/bin/tini", "--", "/opt/docker/bin/workshop_entrypoint" ]
# CMD [ "/bin/bash" ]

USER conda 