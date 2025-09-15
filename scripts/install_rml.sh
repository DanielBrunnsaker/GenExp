# Make directory /opt/tools
RMLMAPPER_JAR="https://github.com/RMLio/rmlmapper-java/releases/download/v7.3.3/rmlmapper-7.3.3-r374-all.jar"
mkdir -p /opt/tools && cd /opt/tools
curl -LO $RMLMAPPER_JAR