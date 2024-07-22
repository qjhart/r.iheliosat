FROM osgeo/grass-gis:releasebranch_8_3-debian as grass

RUN mkdir /usr/local/grass83/raster && cd /usr/local/grass83/raster && \
    git clone --branch=grass83 https://github.com/qjhart/r.iheliosat.git && \
    cd r.iheliosat && make

WORKDIR /grassdb

#RUN rm -rf /usr/local/grass83/raster
