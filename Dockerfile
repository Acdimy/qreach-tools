# syntax=docker/dockerfile:1

# Comments are provided throughout this file to help you get started.
# If you need more help, visit the Dockerfile reference guide at
# https://docs.docker.com/go/dockerfile-reference/

# Want to help us make this template better? Share your feedback here: https://forms.gle/ybq9Krt8jtBL3iCk7

ARG PYTHON_VERSION=3.9.17
FROM python:${PYTHON_VERSION}-slim as base

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update
RUN apt-get install -y make
RUN apt-get install -y gcc g++


# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
# ARG UID=10001
# RUN adduser \
#     --disabled-password \
#     --gecos "" \
#     --home "/nonexistent" \
#     --shell "/sbin/nologin" \
#     --no-create-home \
#     --uid "${UID}" \
#     appuser

# Download dependencies as a separate step to take advantage of Docker's caching.
# Leverage a cache mount to /root/.cache/pip to speed up subsequent builds.
# Leverage a bind mount to requirements.txt to avoid having to copy them into
# into this layer.
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

# Switch to the non-privileged user to run the application.
# USER appuser

# Copy the source code into the container.
COPY . /app

RUN export PYTHON_INCLUDE=`python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"`
RUN export BOOST_PATH=$(pwd)/boost_1_81_0/
# RUN export BOOST_PATH=/app/boost_1_81_0/

WORKDIR /app/python_pkg
RUN invoke build-qreach
RUN invoke build-pybind11

# Expose the port that the application listens on.
EXPOSE 8000

# Run the application.
CMD python3 /app/python_pkg/test_new.py qrw 2 1
# CMD python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"
# WORKDIR /app
# CMD ["bash"]