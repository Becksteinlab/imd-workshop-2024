FROM ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

# Modify the entrypoint script so that it activates the workshop environment
COPY workshop_entrypoint /opt/docker/bin/workshop_entrypoint
RUN chmod +x /opt/docker/bin/workshop_entrypoint

RUN source /opt/conda/etc/profile.d/conda.sh && \
        conda env create -n workshop --file=env.yaml

ENTRYPOINT [ "/opt/conda/bin/tini", "--", "/opt/docker/bin/workshop_entrypoint" ]

CMD [ "/bin/bash" ]
