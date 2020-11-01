FROM ubuntu:20.04

# [Option] Install zsh
ARG INSTALL_ZSH="true"
# [Option] Upgrade OS packages to their latest versions
ARG UPGRADE_PACKAGES="true"

# Install needed packages and setup non-root user. Use a separate RUN statement to add your own dependencies.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

COPY library-scripts/common-debian.sh /tmp/library-scripts/
RUN apt-get update \
    && /bin/bash /tmp/library-scripts/common-debian.sh "${INSTALL_ZSH}" "${USERNAME}" "${USER_UID}" "${USER_GID}" "${UPGRADE_PACKAGES}" \
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts/

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y gcc-10 g++-10 cmake ninja-build ccache wget python3 python3-pip git libblas-dev liblapack-dev libeigen3-dev software-properties-common liblapacke-dev

COPY library-scripts/install-dune.sh /tmp/library-scripts/
RUN apt-get update \
    && /bin/bash /tmp/library-scripts/install-dune.sh \
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts/
RUN echo "#define HAVE_UG 1" >> /usr/local/include/dune/grid/config.h && echo "#define UG_USE_NEW_DIMENSION_DEFINES" >> /usr/local/include/dune/grid/config.h

RUN SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.zsh_history" \
    && echo $SNIPPET >> "/root/.zshrc" \
    && mkdir /commandhistory \
    && touch /commandhistory/.zsh_history \
    && chown -R $USERNAME /commandhistory \
    && echo $SNIPPET >> "/home/$USERNAME/.zshrc"
RUN mkdir /ccache && chown -R $USERNAME /ccache
ENV CCACHE_DIR=/ccache

# Setup default python tools in a venv via pipx to avoid conflicts
ENV PIPX_HOME=/usr/local/py-utils \
    PIPX_BIN_DIR=/usr/local/py-utils/bin
ENV PATH=${PATH}:${PIPX_BIN_DIR}
COPY library-scripts/python-debian.sh /tmp/library-scripts/
RUN apt-get update && DEBIAN_FRONTEND=noninteractive sudo apt-get install -y python3-venv libatlas-base-dev fenics clang-format clang-tidy
RUN bash /tmp/library-scripts/python-debian.sh "none" "/usr/local" "${PIPX_HOME}" "${USERNAME}" "false" \ 
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts/

# conan default setup
RUN conan profile new default --detect
RUN conan profile update settings.compiler=gcc default
RUN conan profile update settings.compiler.libcxx=libstdc++11 default
RUN conan profile update settings.build_type=Release default
RUN mkdir -p /home/$USERNAME/.conan/data && sudo chown -R $USERNAME /home/$USERNAME/.conan

COPY library-scripts/install-kaskade-deps.sh /tmp/library-scripts/
RUN bash /tmp/library-scripts/install-kaskade-deps.sh

COPY Kaskade7.4.tar.gz /tmp/deps/Kaskade7.4.tar.gz
COPY library-scripts/install-kaskade.sh /tmp/library-scripts/
RUN bash /tmp/library-scripts/install-kaskade.sh \ 
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* /tmp/library-scripts/

USER $USERNAME
RUN conan profile new default --detect
RUN conan profile update settings.compiler=gcc default
RUN conan profile update settings.compiler.libcxx=libstdc++11 default
RUN conan profile update settings.build_type=Release default
