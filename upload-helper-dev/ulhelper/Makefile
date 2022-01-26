# makefile for running tests in upload-helper
#
TOP_DIR = ${CURDIR}

PROJECT_DIR  =${TOP_DIR}/projects
ULHELPER_SRC =${TOP_DIR}/ulhelper
WEBSERVER    =${TOP_DIR}/webserver


default: help

# NOTE: this code taken from https://gist.github.com/rcmachado/af3db315e31383502660
help: ## This Makefile can be used to run the upload-helper unit tests
	$(info Available targets:)
	@awk '/^[a-zA-Z\-\_0-9]+:/ {                                   \
          nb = sub( /^## /, "", helpMsg );                             \
          if(nb == 0) {                                                \
            helpMsg = $$0;                                             \
            nb = sub( /^[^:]*:.* ## /, "", helpMsg );                  \
          }                                                            \
          if (nb)                                                      \
            printf "\033[1;31m%-" width "s\033[0m %s\n", $$1, helpMsg; \
        }                                                              \
        { helpMsg = $$0 }'                                             \
        width=$$(grep -o '^[a-zA-Z_0-9]\+:' $(MAKEFILE_LIST) | wc -L)  \
	$(MAKEFILE_LIST)
test:  ## Run all unit tests
	pytest
cov: ## produce test coverage
	pytest --cov=./
cov-verbose: ## produce a verbose test coverage
	pytest -v --cov-report term-missing --cov=./
