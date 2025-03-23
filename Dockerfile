ARG REMOTE_SRC=1
ARG GIT_BRANCH=develop
ARG BUILDPLATFORM=linux/amd64
FROM --platform=${BUILDPLATFORM} sagemath/sagemath:10.5 AS base
RUN sudo apt-get update
RUN sudo dpkg --configure -a
RUN sudo apt-get -f -y install git

FROM base AS use-git-1
RUN git clone https://github.com/fredstro/hilbertmodgroup.git
WORKDIR "/home/sage/hilbertmodgroup"
RUN git config pull.rebase false && git checkout ${GIT_BRANCH}

FROM base AS use-git-0
ARG GIT_BRANCH=''
COPY --chown=sage . hilbertmodgroup
WORKDIR "/home/sage/hilbertmodgroup"

FROM use-git-${REMOTE_SRC} AS final
RUN sudo apt-get -y install make
COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
CMD ["run"]