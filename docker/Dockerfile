FROM ubuntu:20.04

# Éviter les interactions pendant l'installation
ENV DEBIAN_FRONTEND=noninteractive

# Installer les dépendances nécessaires
RUN apt-get update -y && \
    apt-get install -y \
    gcc \
    libblas-dev \
    liblapacke-dev \
    make \
    vim

# Définir le répertoire de travail
WORKDIR /app

# Copier les fichiers du projet
COPY . .

# Créer le dossier bin s'il n'existe pas
RUN mkdir -p bin

# Utiliser la configuration Docker pour la compilation
ENV HOSTNAME=docker

# Compiler et exécuter les tests
CMD make all && make run