FROM ubuntu:18.04
# FROM python:2.7-stretch
# ADD . /temp
COPY . /tmp
# RUN mkdir new
WORKDIR /tmp
RUN cp test.csv imusim-py3/imusim/maths
ENV LC_ALL=C.UTF-8
ENV TERM=xterm
# RUN dpkg --add-architecture i386
# RUN apt-get update && apt-get upgrade -y && apt-get install -y \
#         python3.7 python3-pip libsm6 libx11-6 libx11-dev \
#         libxcb1-dev libxrender1 libxext6 libfontconfig1 \
#         libgl1-mesa-glx:i386 libgl1-mesa-glx git cython3 libxkbcommon-x11-0 \
#         xvfb x11-apps xauth xorg openbox fluxbox

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
        python3.7 python3-pip libsm6 libx11-6 libx11-dev \
        libxcb1-dev libxrender1 libxext6 libfontconfig1 \
        libgl1-mesa-glx git cython3 nano

# RUN pip install --upgrade pip

# ENV LD_LIBRARY_PATH=/usr/local/lib/python3.7/site-packages/:$LD_LIBRARY_PATH
# ENV PATH=/usr/include/:$PATH
# RUN export OPENLGL_INCLUDE_DIR=/usr/include/GL
# ENV PYTHONPATH=/usr/local/lib/python3.7/site-packages/:$PYTHONPATH
# ENV PATH=/usr/local/lib/python3.7/site-packages/vtk:$PATH
# ENV PATH=/usr/local/VTKBuild/bin:$PATH
# ENV PYTHONPATH=/usr/local/VTKBuild/Wrapping/Python/vtk:/usr/local/VTKBuild/Wrapping/Python:$PYTHONPATH
# RUN pip install --trusted-host pypi.python.org -r requirements_imu.txt
# RUN TMPDIR=new
# RUN apt-get update  && apt-get install -y apt-utils && \
#     apt-get install -y python-numpy python-scipy python-matplotlib cython\
#     python-dev mayavi2 ipython python-setuptools python-simpy python-pyparsing
# RUN apt-get update && apt-get install -y python-pip
# RUN pip install -U -t /usr/local/lib/python3.7/site-packages/ -r old_req.txt -b /usr/local/lib/python3.7/site-packages/

# ENV DEBIAN_FRONTEND=noninteractive

# ENV ETS_TOOLKIT=null
# ENV DISPLAY=:0
# ENV XSOCK=/tmp/.X11-unix
# ENV XAUTH=/tmp/.docker.xauth
# ENV QT_X11_NO_MITSHM=1
# ENV QT_DEBUG_PLUGINS=1
# ENV QT_QPA_PLATFORM=offscreen
# ENTRYPOINT ["tini", "-g", "--", "xvfb-run"]
# RUN apt-get update && apt-get upgrade -y && apt-get\
#      install -y libgl1-mesa-glx x11-apps xvfb xserver xauth xorg openbox
# RUN apt-get update && apt-get upgrade -y \
#     && apt-get install -y libgl1-mesa-glx paraview

RUN pip3 install -U -r reqs.txt -b /tmp

# ENV VTK_OPENGL_HAS_OSMESA=ON

RUN pip3 install -U mayavi -b /tmp

# RUN python -c "import vtk"
# RUN python -c "print(hasattr(vtk, 'vtkOSOpenGLRenderWindow'))"
# RUN pip install -U -r py2_req.txt -b /new
# RUN apt-get update && apt-get install -y python-numpy
# RUN python -c "import numpy"
# RUN python -c "import vtk"
# RUN pip install -U -t /usr/local/lib/python3.7/site-packages/ mayavi -b /usr/local/lib/python3.7/site-packages/
# RUN python -c "import mayavi"
# RUN pip install -r old_req.txt
# RUN pip install --upgrade -t /usr/local/lib/python3.7/site-packages/ -r requirements_imu.txt -b /usr/local/lib/python3.7/site-packages/
# RUN pip install python-vtk
# RUN python -c "import vtk"
# RUN pip install mayavi
# RUN python -c "import vtk"
# RUN git clone git://github.com/enthought/mayavi
# WORKDIR mayavi
# RUN python setup.py install
# RUN apt-get update && \
#     apt-get install -y \
#     apt-utils
# RUN apt-get update && \
#     apt-get install -y \
#     mayavi2
# RUN pip wheel -r requirements_imu.txt
# RUN pip install mayavi
# RUN pip wheel mayavi
# ENV passenv = VTK_PATH LD_LIBRARY_PATH
# RUN pip install mayavi
# RUN pip install --trusted-host pypi.python.org -r requirements_seg.txt
# COPY . /seglearn
# COPY requirements.txt /temp
# RUN mkdir seg
# RUN cd seg
# RUN git clone https://github.com/dmbee/seglearn /seg
# RUN git clone git://github.com/martinling/imusim
# RUN git clone git://github.com/Shmuma/imusim/tree/py3 /imu
# RUN Xvfb :1 -screen 0 800x600x16
# RUN apt-get update && apt-get upgrade -y && apt-get install -y git cython3

# RUN git clone git://github.com/Shmuma/imusim -b py3
# WORKDIR imusim

WORKDIR imusim-py3
RUN cython3 -3 -a imusim/maths/*.pyx

# RUN cython3 -2 -a imusim/maths/*.pyx
# ENV DISPLAY=:1
# EXPOSE 5901

# RUN cython language_level=Py3
# RUN cython -a imusim/maths/*.pyx
# RUN apt-get update && \
#     apt-get install -y \
#     libvtk6-dev
# RUN apt-get update && \
#     apt-get install -y \
#     mayavi2
# RUN pip install --trusted-host pypi.python.org -r requirements_imu.txt
# RUN git clone https://github.com/martinling/imusim /imu

RUN python3 setup.py install
WORKDIR imusim/maths
RUN mkdir mount
# RUN jupyter notebook --allow-root --ip='0.0.0.0' --port=8888