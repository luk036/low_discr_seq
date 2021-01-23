# low_discr_seq

[![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/luk036/low_discr_seq)
[![Build Status](https://travis-ci.com/luk036/low_discr_seq.svg?branch=main)](https://travis-ci.com/luk036/low_discr_seq)
[![Documentation Status](https://readthedocs.org/projects/low_discr_seq/badge/?version=latest)](https://low_discr_seq.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c7e8c69a3335427aa2f08e3e2d455552)](https://app.codacy.com/app/luk036/low_discr_seq?utm_source=github.com&utm_medium=referral&utm_content=luk036/low_discr_seq&utm_campaign=badger)
[![CodeFactor](https://www.codefactor.io/repository/github/luk036/low_discr_seq/badge)](https://www.codefactor.io/repository/github/luk036/low_discr_seq)
[![codecov](https://codecov.io/gh/luk036/low_discr_seq/branch/main/graph/badge.svg?token=Kxl10DrV6g)](https://codecov.io/gh/luk036/low_discr_seq)

Low Discrepancy Sequence C++ Code

To run in gitpod.io:

    ./envconfig.sh  # first time when gitpod image is built

To build with Ninja:

    mkdir build && cd build
    cmake -GNinja ..
    ninja all

To run CTest:

    ninja test
