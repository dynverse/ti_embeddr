FROM dynverse/dynwrapr:v0.1.0

RUN R -e 'devtools::install_github("dynverse/scaterlegacy")'
RUN R -e 'devtools::install_github("dynverse/embeddr")'

COPY definition.yml run.R example.sh /code/

ENTRYPOINT ["/code/run.R"]
