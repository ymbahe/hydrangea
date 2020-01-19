# -----------------------------------------------------------------------
# This file is part of the hydrangea tools package
# Copyright (C) 2019 Yannick Bahe (bahe@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hydrangea",
    version="0.1",
    author="Yannick Bahe",
    author_email="bahe@strw.leidenuniv.nl",
    description="Tools for working with the Hydrangea/C-EAGLE simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://www.leidenuniv.nl",
    include_package_data=True,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: LGPL-3.0-or-later",
        "Operating System :: OS independent"
    ],
    python_requires='>=3.6'
)

