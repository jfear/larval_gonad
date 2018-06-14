#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

#requirements = [i.strip() for i in open('requirements.txt').readlines()]

setup(
    name='larval_gonad',
    version='0.0.1',
    description="Local library for the larval gonad project",
    author="Justin M Fear",
    author_email='justin.m.fear@gmail.com',
    url='https://github.com/jfear/larval_gonad',
    packages=['larval_gonad'],
#    install_requires=requirements,
    license="MIT license",
)
