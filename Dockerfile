FROM sagemath/sagemath:9.4
RUN sudo apt-get update
RUN sudo apt-get -y install git
RUN git clone https://github.com/fredstro/hilbertmodgroup.git
WORKDIR "hilbertmodgroup"
RUN git config pull.rebase false && git checkout -b develop
RUN sage setup.py install
COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
CMD ["run"]