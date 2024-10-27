FROM ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

# Modify the entrypoint script so that it activates the workshop environment
COPY entrypoint /opt/docker/bin/entrypoint
RUN chmod +x /opt/docker/bin/entrypoint

RUN source /opt/conda/etc/profile.d/conda.sh && \
        conda env create -n workshop --file=env.yaml

# Use inherited entrypoint

CMD [ "/bin/bash" ]
