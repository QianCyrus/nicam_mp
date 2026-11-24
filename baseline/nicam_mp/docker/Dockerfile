FROM nvcr.io/nvidia/nvhpc:26.3-devel-cuda_multi-ubuntu24.04

ENV NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
ENV PATH=${NVHPC_ROOT}/compilers/bin:${PATH}
ENV LD_LIBRARY_PATH=${NVHPC_ROOT}/compilers/lib:${LD_LIBRARY_PATH}

# Build dependencies
RUN apt-get update && apt-get install -y \
    wget make cmake m4 curl git libxml2-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# ------- HDF5 (release asset - URL is fine) -------
RUN curl -L -o hdf5-1.14.5.tar.gz \
    https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz \
    && tar xzf hdf5-1.14.5.tar.gz && cd hdf5-1.14.5 \
    && CC=nvc CXX=nvc++ FC=nvfortran ./configure \
       --prefix=/usr/local --enable-fortran --enable-shared \
    && make -j$(nproc) && make install \
    && cd / && rm -rf hdf5-*

# ------- NetCDF-C (use archive URL) -------
RUN curl -L -o netcdf-c-4.9.2.tar.gz \
    https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz \
    && tar xzf netcdf-c-4.9.2.tar.gz && cd netcdf-c-4.9.2 \
    && CC=nvc CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --disable-dap \
    && make -j$(nproc) && make install \
    && cd / && rm -rf netcdf-c-*

# ------- NetCDF-Fortran (use archive URL) -------
RUN curl -L -o netcdf-fortran-4.6.1.tar.gz \
    https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz \
    && tar xzf netcdf-fortran-4.6.1.tar.gz && cd netcdf-fortran-4.6.1 \
    && CC=nvc FC=nvfortran \
       CFLAGS="-fPIC" FCFLAGS="-fPIC" \
       CPPFLAGS="-I/usr/local/include" \
       LDFLAGS="-L/usr/local/lib" \
       LD_LIBRARY_PATH="/usr/local/lib:${NVHPC_ROOT}/compilers/lib" \
       ./configure --prefix=/usr/local \
    && make -j$(nproc) && make install \
    && cd / && rm -rf netcdf-fortran-*

# Update library cache
RUN ldconfig

# Dev tools (optional: fortls for VSCode)
RUN apt-get update && apt-get install -y python3-pip \
    && pip install --break-system-packages fortls fprettify \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
